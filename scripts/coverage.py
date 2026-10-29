#!/usr/bin/env python3
"""How many fields have a name, per degree and per family.

Reads data/degree<n>.txt. With --db, also counts the LMFDB fields of each degree (through the
lmfdb package or lmfdb-lite) to show what fraction is named.

With --by-galois, goes through the LMFDB fields of the given degrees instead and counts the
named and unnamed fields per Galois group, most unnamed first, with the smallest unnamed
fields of each group as examples: this shows which families would name the most fields.
Going through a whole degree reads every field's label, so use --max-disc for a quicker look.

Usage:
    python3 scripts/coverage.py
    python3 scripts/coverage.py --db
    python3 scripts/coverage.py --by-galois --degrees 3-8 --max-disc 10^8
"""

import argparse
import glob
import heapq
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nf_names import DATA_DIR, FAMILIES, data_path, family_rank, label_key, read_data


def parse_degrees(s):
    """'3-6,9' -> [3, 4, 5, 6, 9]"""
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return out


def parse_int(s):
    """An integer, also written as 10^6 or 1e6."""
    if "^" in s:
        b, e = s.split("^")
        return int(b) ** int(e)
    return int(float(s)) if "e" in s.lower() else int(s)


def galois_counts(rows, named, examples=3):
    """From (label, galois_label) rows, {group: (fields, named, smallest unnamed labels)}."""
    total, have, unnamed = Counter(), Counter(), defaultdict(list)
    for label, group in rows:
        group = group or "?"
        total[group] += 1
        if label in named:
            have[group] += 1
        else:  # keep the `examples` smallest, in LMFDB order
            heap = unnamed[group]
            item = tuple(-x for x in label_key(label)), label
            if len(heap) < examples:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)
    return {g: (total[g], have[g], [lab for _, lab in sorted(unnamed[g], reverse=True)])
            for g in total}


def group_names(db, groups):
    """{LMFDB transitive group label: its name, e.g. S3}, or {} if gps_transitive is not there."""
    try:
        return {r["label"]: r["name"] for r in db.gps_transitive.search(
            {"label": {"$in": sorted(groups)}}, ["label", "name"])}
    except Exception:
        return {}


def print_galois_report(db, degree, data_dir, max_disc=None, top=20, examples=3):
    path = data_path(degree, data_dir)
    named = {label for label, _ in read_data(path)} if os.path.exists(path) else set()
    query = {"degree": degree}
    if max_disc is not None:
        query["disc_abs"] = {"$lte": max_disc}
    rows = ((r["label"], r["galois_label"])
            for r in db.nf_fields.search(query, ["label", "galois_label"]))
    counts = galois_counts(rows, named, examples)
    fields = sum(c[0] for c in counts.values())
    have = sum(c[1] for c in counts.values())
    print(f"degree {degree}" + (f", |disc| <= {max_disc}" if max_disc is not None else "")
          + f": {fields} fields, {have} named ({100 * have / fields if fields else 0:.2f}%)")
    order = sorted(counts, key=lambda g: (-(counts[g][0] - counts[g][1]), g))
    names = group_names(db, order[:top])
    print(f"  {'group':<14} {'fields':>9} {'named':>9} {'named %':>8}  smallest unnamed")
    for g in order[:top]:
        total, n, ex = counts[g]
        group = g + (f" {names[g]}" if g in names else "")
        print(f"  {group:<14} {total:>9} {n:>9} {100 * n / total:>7.2f}%  {', '.join(ex)}")
    if len(order) > top:
        rest = [counts[g] for g in order[top:]]
        print(f"  ({len(rest)} more groups, {sum(c[0] - c[1] for c in rest)} unnamed fields)")
    print()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", action="store_true", help="compare with the LMFDB field counts")
    parser.add_argument("--by-galois", action="store_true",
                        help="count named and unnamed LMFDB fields per Galois group (needs --degrees)")
    parser.add_argument("--degrees", help="for --by-galois, e.g. 3-8 or 4,6")
    parser.add_argument("--max-disc", type=parse_int, help="for --by-galois, largest |disc|, e.g. 10^8")
    parser.add_argument("--top", type=int, default=20, help="for --by-galois, groups to show (default 20)")
    parser.add_argument("--examples", type=int, default=3,
                        help="for --by-galois, unnamed fields to show per group (default 3)")
    parser.add_argument("-d", "--data-dir", default=DATA_DIR)
    args = parser.parse_args()
    if args.by_galois and not args.degrees:
        parser.error("--by-galois needs --degrees")

    if args.db or args.by_galois:
        try:
            from lmfdb import db
        except ImportError:
            from lmf import db

    if args.by_galois:
        for n in parse_degrees(args.degrees):
            print_galois_report(db, n, args.data_dir, args.max_disc, args.top, args.examples)
        return

    families = [name for name, _ in FAMILIES] + ["other"]
    total = Counter()
    print(f"{'degree':>6} {'named':>9}" + (f" {'in LMFDB':>10} {'named %':>8}" if args.db else "")
          + "  by family")
    paths = glob.glob(os.path.join(args.data_dir, "degree*.txt"))
    for n, path in sorted((int(re.search(r"degree(\d+)", p)[1]), p) for p in paths):
        by_family = Counter(families[family_rank(name)] for _, name in read_data(path))
        named = sum(by_family.values())
        total.update(by_family)
        line = f"{n:>6} {named:>9}"
        if args.db:
            count = db.nf_fields.count({"degree": n})
            line += f" {count:>10} {100 * named / count:>7.2f}%"
        print(line + "  " + ", ".join(f"{f}: {c}" for f, c in by_family.most_common()))
    print(f"\n{'all':>6} {sum(total.values()):>9}  "
          + ", ".join(f"{f}: {c}" for f, c in total.most_common()))


if __name__ == "__main__":
    main()
