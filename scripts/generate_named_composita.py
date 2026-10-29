#!/usr/bin/env python3
"""Add composita of two named fields, e.g. Q(\\zeta_{5}, \\sqrt{2}), to the data files.

See named_composita.py for which fields are combined and how they are named. --bound limits
the quadratic d and pure a, --quartic-bound the A, B, D of the \\sqrt{A + B\\sqrt{D}} blocks,
and --max-degree the degree of the compositum. Names of other families win over these (see
FAMILIES in nf_names.py), so this is best run after the other generate_*.py scripts, but
the order does not change the result.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_named_composita.py --max-degree 16 --bound 20 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from named_composita import named_composita_fields
from pipeline import add_common_args, lookup_and_merge, new_progress


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-degree", type=int, default=16)
    parser.add_argument("--bound", type=int, required=True, help="largest |d| and |a|")
    parser.add_argument("--quartic-bound", type=int, default=5,
                        help="largest |A|, |B|, |D| (default 5)")
    add_common_args(parser)
    args = parser.parse_args()

    progress = new_progress(args)
    progress.start("composita, computing")
    fields = named_composita_fields(args.max_degree, args.bound, args.quartic_bound, progress)
    lookup_and_merge(fields, args, progress, "composita of named fields")


if __name__ == "__main__":
    main()
