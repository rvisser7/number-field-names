#!/usr/bin/env python3
"""Name LMFDB fields directly: for each field, test every family and keep the best name.

This is the backward direction (see recognise.py). Instead of generating fields from
parameters up to a bound, it goes through the LMFDB fields of the given degrees and, for each,
decides which families it belongs to. For cyclotomic, real cyclotomic, multiquadratic, pure,
Kummer, imprimitive quartic, nested radical, nested biquadratic and nested quartic (octics),
and nested pure (quadratic over pure) fields the answer is complete, with no bound. Cubic
fields with no other name get a Cardano name when it is small enough (--cardano-bound)
(nested radicals with too many candidates are skipped, see recognise.py). Composita come from
named subfields, with two or more generators.

The names are merged into data/degree<n>.txt every --batch fields, keeping the better name
when a field already has one (see name_key in nf_names.py), so an interrupted run keeps
what it has done, and runs can be split with --shard or disc ranges and run in parallel.

Fields come from the LMFDB database (`from lmfdb import db`, or lmfdb-lite), or from a file
with one `label|coeffs` line per field, e.g. `3.1.108.1|[-2,0,0,1]`, optionally followed by
`|ramps` (the ramified primes, e.g. `[2,3]`). Degree 1 and 2 fields are left out (they are
named from their label, see generate_degree2.py).

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    sage -python scripts/name_lmfdb_fields.py --degrees 3-6 --max-disc 10^6
    sage -python scripts/name_lmfdb_fields.py --degrees 6 --galois 6T3 --families pure,kummer
    sage -python scripts/name_lmfdb_fields.py --degrees 8 --shard 0/4    # and 1/4, 2/4, 3/4
    sage -python scripts/name_lmfdb_fields.py --input fields.txt --print
"""

import argparse
import json
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nf_names import FAMILIES, data_path, family_rank, merge_data, read_data
from pari_fields import pari  # noqa: F401 (checks cypari2 is there)
from cypari2 import PariError
from progress import Progress, fmt_time
import recognise
from recognise import FAMILY_TESTS, best_name, from_label


def parse_degrees(s):
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return sorted(set(d for d in out if d >= 3))


def parse_int(s):
    """An integer, also written as 10^6 or 1e6."""
    if s is None:
        return None
    if "^" in s:
        b, e = s.split("^")
        return int(b) ** int(e)
    return int(float(s)) if "e" in s.lower() else int(s)


def parse_list(s):
    s = s.strip()
    return [int(c) for c in (json.loads(s) if s.startswith("[") else s.split(","))]


def db_rows(degree, args):
    """(total, iterator of (label, coeffs, ramps)) for the LMFDB fields to name."""
    try:
        from lmfdb import db
    except ImportError:
        from lmf import db
    query = {"degree": degree}
    disc = {}
    if args.min_disc is not None:
        disc["$gte"] = args.min_disc
    if args.max_disc is not None:
        disc["$lte"] = args.max_disc
    if disc:
        query["disc_abs"] = disc
    if args.galois:
        query["galois_label"] = args.galois
    total = db.nf_fields.count(query)
    # ramps (the ramified primes) saves factoring the discriminant; computed if not there
    has_ramps = "ramps" in getattr(db.nf_fields, "search_cols", ["ramps"])
    cols = ["label", "coeffs"] + (["ramps"] if has_ramps else [])
    rows = db.nf_fields.search(query, cols)
    return total, ((r["label"], [int(c) for c in r["coeffs"]],
                    [int(p) for p in r["ramps"]] if r.get("ramps") is not None else None)
                   for r in rows)


def file_rows(path, degree, args):
    rows = []
    with open(path) as f:
        for line in f:
            cells = line.strip().split("|")
            if len(cells) < 2 or not cells[0][:1].isdigit():
                continue  # blank line or header
            label = cells[0]
            n, _, disc_abs, _ = (int(x) for x in label.split("."))
            if n != degree:
                continue
            if args.min_disc is not None and disc_abs < args.min_disc:
                continue
            if args.max_disc is not None and disc_abs > args.max_disc:
                continue
            rows.append((label, parse_list(cells[1]), parse_list(cells[2]) if len(cells) > 2 else None))
    return len(rows), iter(rows)


def family(name):
    rank = family_rank(name)
    return FAMILIES[rank][0] if rank < len(FAMILIES) else "other"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", required=True, help="e.g. 3-8 or 4,6")
    parser.add_argument("--input", help="read fields from this file instead of the LMFDB")
    parser.add_argument("--source", choices=["db"], default="db",
                        help="where to get fields when there is no --input (only db: the API "
                             "is too slow to go through whole degrees)")
    parser.add_argument("--min-disc", type=parse_int, help="smallest |disc|, e.g. 10^6")
    parser.add_argument("--max-disc", type=parse_int, help="largest |disc|, e.g. 10^8")
    parser.add_argument("--galois", help="only fields with this LMFDB Galois group label, e.g. 6T3")
    parser.add_argument("--families", default=",".join(FAMILY_TESTS).replace(" ", "_"),
                        help="which tests to run (default all): " + ", ".join(
                            f.replace(" ", "_") for f in FAMILY_TESTS))
    parser.add_argument("--shard", default="0/1", help="i/k: only every k-th field, from the i-th")
    parser.add_argument("--only-unnamed", action="store_true",
                        help="skip fields that already have a name in the data files")
    parser.add_argument("--batch", type=int, default=1000,
                        help="merge into the data files after this many fields (default 1000)")
    parser.add_argument("--print", action="store_true",
                        help="print label|name for each named field instead of merging")
    parser.add_argument("--replace", action="store_true",
                        help="overwrite names already in the data files")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    parser.add_argument("--cardano-bound", type=int, default=recognise.MAX_CARDANO_SIZE,
                        help="largest number in a Cardano cubic name (default %(default)s; 0 for "
                             "no bound), see cardano.py")
    args = parser.parse_args()
    recognise.MAX_CARDANO_SIZE = args.cardano_bound or None

    families = [f.replace("_", " ") for f in args.families.split(",")]
    unknown = set(families) - set(FAMILY_TESTS)
    if unknown:
        sys.exit(f"unknown families: {', '.join(sorted(unknown))}")
    shard, nshards = (int(x) for x in args.shard.split("/"))
    progress = Progress(every=args.progress_every)

    for n in parse_degrees(args.degrees):
        t0 = time.time()
        total, rows = file_rows(args.input, n, args) if args.input else db_rows(n, args)
        named_before = {}
        if args.only_unnamed and os.path.exists(data_path(n)):
            named_before = dict(read_data(data_path(n)))
        found, batch = Counter(), {}
        merged = Counter()
        checked = errors = 0
        progress.start(f"degree {n}, naming LMFDB fields")

        def flush():
            if batch and not args.print:
                merged.update(dict(zip(["added", "replaced", "kept"],
                                       merge_data(data_path(n), batch, args.replace))))
            batch.clear()

        for i, (label, coeffs, ramps) in enumerate(rows):
            if i % nshards != shard or label in named_before:
                continue
            checked += 1
            try:
                name = best_name(from_label(label, coeffs, ramps), families)
            except PariError as e:
                errors += 1
                print(f"PARI error on {label}: {e}", file=sys.stderr)
                name = None
            if name:
                found[family(name)] += 1
                batch[label] = name
                if args.print:
                    print(f"{label}|{name}", flush=True)
            if checked % args.batch == 0:
                flush()
            progress(i + 1, total, f", {sum(found.values())} named")
        flush()
        print(f"degree {n}: {checked} fields checked, {sum(found.values())} named"
              + (f" ({', '.join(f'{f}: {c}' for f, c in found.most_common())})" if found else "")
              + (f"; {merged['added']} added, {merged['replaced']} replaced, "
                 f"{merged['kept']} kept existing name" if not args.print else "")
              + (f"; {errors} PARI errors" if errors else "")
              + f", in {fmt_time(time.time() - t0)}", file=sys.stderr if args.print else sys.stdout)


if __name__ == "__main__":
    main()
