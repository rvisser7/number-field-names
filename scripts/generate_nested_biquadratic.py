#!/usr/bin/env python3
"""Add octic fields Q(\\sqrt{A + B\\sqrt{d_1} + C\\sqrt{d_2} + D\\sqrt{d_3}}) to data/degree8.txt.

Runs over |A|, |B|, |C|, |D| <= bound and the biquadratic fields Q(\\sqrt{d_1}, \\sqrt{d_2})
with |d_1|, |d_2| <= --basis-bound (see nested_biquadratic.py), looks the fields up in the
LMFDB and merges the names into data/degree8.txt. Names of every other family win over these
(see FAMILIES in nf_names.py), so only octics with no other name get one, unless --replace
is given. name_lmfdb_fields.py finds these names too, with no bound.

There are about 2 (bound + 1)^2 (2 bound + 1)^2 quadruples for each biquadratic field, e.g.
about 1600 for bound 3, and each takes about 1 ms.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_nested_biquadratic.py --bound 3 --basis-bound 10 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nested_biquadratic import nested_biquadratic_fields
from pipeline import add_common_args, lookup_and_merge, new_progress


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bound", type=int, required=True, help="largest |A|, |B|, |C| and |D|")
    parser.add_argument("--basis-bound", type=int, default=10,
                        help="largest |d_1| and |d_2| (default 10)")
    add_common_args(parser)
    args = parser.parse_args()

    progress = new_progress(args)
    progress.start("nested biquadratic, computing (by biquadratic field)")
    fields = nested_biquadratic_fields(args.bound, args.basis_bound, progress)
    lookup_and_merge(fields, args, progress, "nested biquadratic fields")


if __name__ == "__main__":
    main()
