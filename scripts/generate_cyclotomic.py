#!/usr/bin/env python3
"""Add cyclotomic fields Q(\\zeta_n) and real cyclotomic fields Q(\\zeta_n)^+ to the data files.

Computes every such field of degree 3 to --max-degree (see cyclotomic.py), looks up their
labels in the LMFDB, and merges the names into data/degree<d>.txt. Cyclotomic names win
over pure names (e.g. 4.0.256.1 becomes \\Q(\\zeta_{8}) rather than \\Q(\\sqrt[4]{-1})),
following the LMFDB's field_pretty; other existing names are kept unless --replace is given.
Fields that are not in the LMFDB are listed at the end.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_cyclotomic.py --max-degree 47 --source db
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cyclotomic import cyclotomic_fields, real_cyclotomic_fields
from lmfdb_lookup import api_source, db_source, label_fields
from nf_names import data_path, merge_data
from pari_fields import is_isomorphic
from progress import Progress


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-degree", type=int, default=47,
                        help="largest field degree (default 47, the largest in the LMFDB)")
    parser.add_argument("--source", choices=["api", "db"], default="api")
    parser.add_argument("--replace", action="store_true",
                        help="overwrite names already in the data files")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    args = parser.parse_args()
    source = api_source if args.source == "api" else db_source()

    progress = Progress(every=args.progress_every)
    found, missing = {}, []
    for what, compute in [("cyclotomic", cyclotomic_fields),
                          ("real cyclotomic", real_cyclotomic_fields)]:
        progress.start(f"{what}, computing fields")
        fields = compute(args.max_degree, progress)
        progress.start(f"{what}, LMFDB lookups")
        f, m = label_fields(fields, source, is_isomorphic, progress)
        print(f"{what}: {len(fields)} fields of degree 3 to {args.max_degree}, {len(f)} in the LMFDB")
        found.update(f)
        missing += m

    by_degree = defaultdict(dict)
    for label, name in found.items():
        by_degree[int(label.split(".")[0])][label] = name
    for d in sorted(by_degree):
        added, replaced, kept = merge_data(data_path(d), by_degree[d], args.replace)
        print(f"degree {d}: {len(by_degree[d])} fields; "
              f"{added} added, {replaced} replaced, {kept} kept existing name")

    if missing:
        print(f"\n{len(missing)} fields not found in the LMFDB:")
        for K in missing:
            print(f"  {K.name}  degree {K.degree}  disc {K.disc}")


if __name__ == "__main__":
    main()
