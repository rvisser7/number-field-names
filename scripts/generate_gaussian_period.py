#!/usr/bin/env python3
"""Add abelian fields named by a Gaussian period, Q(\\zeta_f + \\zeta_f^{h_2} + ... + \\zeta_f^{h_k}), to the data files.

For each degree n, lists every field Q(\\zeta_f)^H of degree n, f its conductor, whose Gaussian
period sum_{h in H} \\zeta_f^h generates it and has |H| <= --max-terms terms (see
gaussian_period.py; there are finitely many, as phi(f) <= |H| n), looks them up in the LMFDB
and merges the names into data/degree<n>.txt. Names of other families win over these (e.g.
\\Q(\\zeta_{5}, \\sqrt{2})), unless --replace is given. name_lmfdb_fields.py finds these names too.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_gaussian_period.py --degrees 3-24 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gaussian_period import MAX_TERMS, period_fields
from pipeline import add_common_args, lookup_and_merge, new_progress


def parse_degrees(s):
    """'3-5,8' -> [3, 4, 5, 8]"""
    out = set()
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out.update(range(int(lo), int(hi or lo) + 1))
    return sorted(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", required=True, help="e.g. 3-24 or 4,6 (each at least 3)")
    parser.add_argument("--max-terms", type=int, default=MAX_TERMS,
                        help="most terms in a name (default %(default)s)")
    add_common_args(parser)
    args = parser.parse_args()
    degrees = parse_degrees(args.degrees)
    if degrees[0] < 3:
        parser.error("degrees must be at least 3")

    progress = new_progress(args)
    progress.start("gaussian periods, computing (by degree)")
    fields = period_fields(degrees, args.max_terms, progress)
    lookup_and_merge(fields, args, progress, "gaussian periods")


if __name__ == "__main__":
    main()
