#!/usr/bin/env python3
"""Add cubic fields Q(\\sqrt[3]{u + e\\sqrt{d}} + \\sqrt[3]{u - e\\sqrt{d}}) (Cardano's formula) to data/degree3.txt.

Runs over the polynomials x^3 - 3kx - 2u with 0 < u <= bound and 0 < |k| <= --k-bound (see
cardano.py), keeps the names whose numbers are all at most --max-size, looks the fields up
in the LMFDB and merges the names into data/degree3.txt. Names of every other family win
over these (e.g. \\Q(\\sqrt[3]{2}), \\Q(\\zeta_{7})^+), unless --replace is given.
name_lmfdb_fields.py finds these names too, with no bound on u and k.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_cardano.py --bound 100 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cardano import cardano_fields, cardano_size
from nf_names import parse_cardano
from pipeline import add_common_args, lookup_and_merge, new_progress


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bound", type=int, required=True, help="largest u")
    parser.add_argument("--k-bound", type=int, help="largest |k| (default: --bound)")
    parser.add_argument("--max-size", type=int, default=100,
                        help="largest number in a name (default 100; 0 for no bound)")
    add_common_args(parser)
    args = parser.parse_args()

    progress = new_progress(args)
    progress.start("cardano cubics, computing (by u)")
    fields = cardano_fields(args.bound, args.k_bound, progress)
    if args.max_size:
        fields = [K for K in fields if cardano_size(*parse_cardano(K.name)) <= args.max_size]
    lookup_and_merge(fields, args, progress, "cardano cubics")


if __name__ == "__main__":
    main()
