#!/usr/bin/env python3
"""Add pure fields Q(\\sqrt[n]{m}) to data/degree<n>.txt.

For each degree n, runs over 1 <= |m| <= bound, finds the distinct fields (see pure.py),
looks up their labels in the LMFDB, and merges the names into data/degree<n>.txt.
Fields already named in the data file keep their name unless --replace is given.
Fields that are not in the LMFDB are listed at the end.

Needs cypari2 (pip install cypari2).

Usage:
    python3 scripts/generate_pure.py --degrees 3-10 --bound 100            # LMFDB API
    python3 scripts/generate_pure.py --degrees 3,5 --bound 1000 --source db  # lmfdb / lmfdb-lite
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lmfdb_lookup import api_source, db_source, find_label
from nf_names import data_path, read_data, write_data
from pure import is_isomorphic, pure_fields


def parse_degrees(s):
    """'3-6,9' -> [3, 4, 5, 6, 9]"""
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return out


def merge(path, new_rows, replace):
    """Merge {label: name} into the data file at path; return (added, replaced, kept)."""
    rows = dict(read_data(path)) if os.path.exists(path) else {}
    added = replaced = kept = 0
    for label, name in new_rows.items():
        if label not in rows:
            added += 1
        elif rows[label] == name:
            continue
        elif replace:
            replaced += 1
        else:
            kept += 1
            continue
        rows[label] = name
    write_data(path, rows.items())
    return added, replaced, kept


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", required=True, help="e.g. 3-10 or 3,5,7")
    parser.add_argument("--bound", type=int, required=True, help="largest |m|")
    parser.add_argument("--source", choices=["api", "db"], default="api")
    parser.add_argument("--replace", action="store_true",
                        help="overwrite names already in the data files")
    args = parser.parse_args()
    source = api_source if args.source == "api" else db_source()

    missing = []
    for n in parse_degrees(args.degrees):
        if n < 3:
            sys.exit("pure fields start at degree 3 (degree 2 is generate_degree2.py)")
        found = {}
        for K in pure_fields(n, args.bound):
            label = find_label(K.coeffs, K.r1, source(n, abs(K.disc)), is_isomorphic)
            if label is None:
                missing.append(K)
            else:
                found[label] = K.name
        added, replaced, kept = merge(data_path(n), found, args.replace)
        print(f"degree {n}: {len(found)} pure fields in the LMFDB; "
              f"{added} added, {replaced} replaced, {kept} kept existing name")

    if missing:
        print(f"\n{len(missing)} pure fields not found in the LMFDB:")
        for K in missing:
            print(f"  {K.name}  disc {K.disc}  polredabs {list(K.coeffs)}")


if __name__ == "__main__":
    main()
