#!/usr/bin/env python3
"""Write data/degree2.txt: pretty names Q(sqrt(D)) for quadratic fields.

Quadratic names can be read off the label, so no field data is needed.

Usage:
    python3 scripts/generate_degree2.py --bound 1000        # all |D| <= 1000
    python3 scripts/generate_degree2.py --labels labels.txt # one LMFDB label per line
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nf_names import (data_path, is_fundamental_discriminant, quadratic_label,
                      quadratic_name, write_data)


def labels_up_to(bound):
    return [quadratic_label(D) for D in range(-bound, bound + 1)
            if is_fundamental_discriminant(D)]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--bound", type=int, help="include every |D| up to this bound")
    src.add_argument("--labels", help="file with one quadratic field label per line")
    parser.add_argument("-o", "--output", default=data_path(2))
    args = parser.parse_args()

    if args.bound is not None:
        labels = labels_up_to(args.bound)
    else:
        with open(args.labels) as f:
            labels = [ln.strip() for ln in f if ln.strip()]

    write_data(args.output, [(label, quadratic_name(label)) for label in labels])
    print(f"wrote {len(labels)} rows to {args.output}")


if __name__ == "__main__":
    main()
