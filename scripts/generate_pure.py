#!/usr/bin/env python3
"""Add pure fields Q(\\sqrt[n]{m}) to data/degree<n>.txt.

For each degree n, runs over 1 <= |m| <= bound, finds the distinct fields (see pure.py),
looks up their labels in the LMFDB, and merges the names into data/degree<n>.txt.
Fields already named in the data file keep their name, unless it is from a family that
pure fields win over (see FAMILIES in nf_names.py) or --replace is given.
Fields that are not in the LMFDB are listed at the end.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_pure.py --degrees 3-10 --bound 100            # LMFDB API
    python3 scripts/generate_pure.py --degrees 3,5 --bound 1000 --source db  # lmfdb / lmfdb-lite
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lmfdb_lookup import api_source, db_source, label_fields
from nf_names import data_path, merge_data
from progress import Progress
from pure import is_isomorphic, pure_fields


def parse_degrees(s):
    """'3-6,9' -> [3, 4, 5, 6, 9]"""
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", required=True, help="e.g. 3-10 or 3,5,7")
    parser.add_argument("--bound", type=int, required=True, help="largest |m|")
    parser.add_argument("--source", choices=["api", "db"], default="api")
    parser.add_argument("--replace", action="store_true",
                        help="overwrite names already in the data files")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    args = parser.parse_args()
    source = api_source if args.source == "api" else db_source()

    progress = Progress(every=args.progress_every)
    missing = []
    for n in parse_degrees(args.degrees):
        if n < 3:
            sys.exit("pure fields start at degree 3 (degree 2 is generate_degree2.py)")
        progress.start(f"degree {n}, computing fields")
        fields = pure_fields(n, args.bound,
                             lambda i, total, k: progress(i, total, f", {k} fields so far"))
        progress.start(f"degree {n}, LMFDB lookups")
        found, not_found = label_fields(fields, source, is_isomorphic, progress)
        missing += not_found
        added, replaced, kept = merge_data(data_path(n), found, args.replace)
        print(f"degree {n}: {len(found)} pure fields in the LMFDB; "
              f"{added} added, {replaced} replaced, {kept} kept existing name")

    if missing:
        print(f"\n{len(missing)} pure fields not found in the LMFDB:")
        for K in missing:
            print(f"  {K.name}  disc {K.disc}  polredabs {list(K.coeffs)}")


if __name__ == "__main__":
    main()
