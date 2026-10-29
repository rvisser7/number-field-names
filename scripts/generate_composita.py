#!/usr/bin/env python3
"""Add radical composita Q(\\sqrt[m]{a}, \\sqrt[n]{b}) (m, n coprime) to the data files.

For each pair m < n, runs over |a|, |b| <= bound (see radical.py), looks the fields up in the
LMFDB and merges the names into data/degree<mn>.txt. Names of families that win over these
(see FAMILIES in nf_names.py, e.g. \\Q(\\sqrt[6]{2}) for Q(\\sqrt{2}, \\sqrt[3]{2})) are kept
unless --replace is given.

The number of fields is about (2 * bound)^2 per pair, and degrees 20 and up are much slower,
so use a smaller bound for bigger m * n.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_composita.py --pairs 2x3 --bound 50 --source db
    python3 scripts/generate_composita.py --pairs 2x5,3x4,2x7 --bound 15 --source db
"""

import argparse
import os
import sys
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline import add_common_args, lookup_and_merge, new_progress
from radical import composita_fields


def parse_pairs(s):
    pairs = []
    for part in s.split(","):
        m, n = sorted(int(x) for x in part.split("x"))
        if m < 2 or gcd(m, n) != 1:
            sys.exit(f"{part}: need coprime degrees m, n >= 2")
        pairs.append((m, n))
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pairs", default="2x3", help="degree pairs, e.g. 2x3,2x5 (default 2x3)")
    parser.add_argument("--bound", type=int, required=True, help="largest |a| and |b|")
    add_common_args(parser)
    args = parser.parse_args()

    progress = new_progress(args)
    progress.start("radical composita, computing")
    fields = composita_fields(parse_pairs(args.pairs), args.bound, progress)
    lookup_and_merge(fields, args, progress, "radical composita")


if __name__ == "__main__":
    main()
