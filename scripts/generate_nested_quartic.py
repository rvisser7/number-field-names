#!/usr/bin/env python3
"""Add octics over an imprimitive quartic, Q(\\sqrt{a + b\\sqrt{D} + c\\sqrt{A + B\\sqrt{D}} + d\\sqrt{D}\\sqrt{A + B\\sqrt{D}}}),
to data/degree8.txt.

Runs over |a|, |b|, |c|, |d| <= bound and the imprimitive quartic fields Q(\\sqrt{A + B\\sqrt{D}})
with |A|, |B|, |D| <= --quartic-bound (see nested_quartic.py), looks the fields up in the
LMFDB and merges the names into data/degree8.txt. Names of every other family win over these
(see FAMILIES in nf_names.py), so only octics with no other name get one, unless --replace is
given. name_lmfdb_fields.py finds these names too, with no bound.

There are about (bound + 1)(2 bound + 1)^3 quadruples for each quartic field, e.g. 500 for
bound 2, and each takes about 1 ms.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_nested_quartic.py --bound 2 --quartic-bound 5 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nested_quartic import nested_quartic_fields
from pipeline import add_common_args, lookup_and_merge, new_progress


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bound", type=int, required=True, help="largest |a|, |b|, |c| and |d|")
    parser.add_argument("--quartic-bound", type=int, default=5,
                        help="largest |A|, |B| and |D| of the quartic field (default 5)")
    add_common_args(parser)
    args = parser.parse_args()

    progress = new_progress(args)
    progress.start("nested quartic, computing (by quartic field)")
    fields = nested_quartic_fields(args.bound, args.quartic_bound, progress)
    lookup_and_merge(fields, args, progress, "nested quartic fields")


if __name__ == "__main__":
    main()
