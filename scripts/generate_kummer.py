#!/usr/bin/env python3
"""Add Kummer fields Q(\\zeta_n, \\sqrt[n]{a}) to the data files.

For each n, runs over 1 <= |a| <= bound (see radical.py), looks the fields up in the LMFDB
and merges the names into data/degree<d>.txt. Names of families that win over these (see
FAMILIES in nf_names.py) are kept unless --replace is given.

The full degree is phi(n) * n: 6, 8, 20, 12, 42, 32 for n = 3, 4, 5, 6, 7, 8. Degrees 20
and up are much slower.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_kummer.py --n 3,4,6 --bound 100 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline import add_common_args, lookup_and_merge, new_progress
from radical import kummer_fields


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", default="3,4,6", help="which n, e.g. 3,4,6 (the default)")
    parser.add_argument("--bound", type=int, required=True, help="largest |a|")
    add_common_args(parser)
    args = parser.parse_args()
    ns = [int(n) for n in args.n.split(",")]
    if any(n < 3 for n in ns):
        sys.exit("n must be at least 3")

    progress = new_progress(args)
    progress.start("Kummer fields, computing")
    fields = kummer_fields(ns, args.bound, progress)
    lookup_and_merge(fields, args, progress, "Kummer fields")


if __name__ == "__main__":
    main()
