#!/usr/bin/env python3
"""Add nested radicals Q(\\sqrt[n]{a + b\\sqrt{d}}), n >= 3, to the data files.

For each n, runs over |a|, |b|, |d| <= bound (see nested_radical.py), looks the fields up in
the LMFDB and merges the names into data/degree<2n>.txt. Names of families that win over
these (see FAMILIES in nf_names.py) are kept unless --replace is given.

There are about 2.5 bound^3 triples (a, b, d) for each n, e.g. 20000 for bound 20.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_nested_radical.py --n 3 --bound 20 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nested_radical import nested_radical_fields
from pipeline import add_common_args, lookup_and_merge, new_progress


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", default="3,4", help="which n, e.g. 3,4 (the default)")
    parser.add_argument("--bound", type=int, required=True, help="largest |a|, |b| and |d|")
    add_common_args(parser)
    args = parser.parse_args()
    ns = [int(n) for n in args.n.split(",")]
    if any(n < 3 for n in ns):
        sys.exit("n must be at least 3 (n = 2 is generate_imprimitive_quartic.py)")

    progress = new_progress(args)
    fields = []
    for n in ns:
        progress.start(f"nested radicals, n = {n}, computing (by d)")
        fields += nested_radical_fields(n, args.bound, progress)
    lookup_and_merge(fields, args, progress, "nested radicals")


if __name__ == "__main__":
    main()
