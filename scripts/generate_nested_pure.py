#!/usr/bin/env python3
"""Add quadratic extensions of pure fields Q(\\sqrt{A_0 + A_1\\sqrt[n]{m_1} + ...}) to the data files.

For each odd prime n, runs over |A_k| <= bound and the pure fields Q(\\sqrt[n]{m}) with
m <= --m-bound (see nested_pure.py), looks the fields up in the LMFDB and merges the names
into data/degree<2n>.txt. Names of every other family win over these (see FAMILIES in
nf_names.py), so only fields with no other name get one, unless --replace is given.
name_lmfdb_fields.py finds these names too, with no bound.

There are about (2 bound + 1)^n coefficient vectors for each m, e.g. 343 for n = 3 and
bound 3, and each takes about 1 ms.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_nested_pure.py --n 3 --bound 5 --m-bound 100 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nested_pure import nested_pure_fields
from pipeline import add_common_args, lookup_and_merge, new_progress


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", default="3", help="which odd primes n, e.g. 3,5 (default 3)")
    parser.add_argument("--bound", type=int, required=True, help="largest |A_k|")
    parser.add_argument("--m-bound", type=int, default=50,
                        help="largest radicand m of the pure field (default 50)")
    add_common_args(parser)
    args = parser.parse_args()
    ns = [int(n) for n in args.n.split(",")]
    if any(n < 3 or any(n % p == 0 for p in range(2, n)) for n in ns):
        sys.exit("n must be an odd prime (n = 2 is generate_imprimitive_quartic.py)")

    progress = new_progress(args)
    fields = []
    for n in ns:
        progress.start(f"nested pure, n = {n}, computing (by m)")
        fields += nested_pure_fields(n, args.bound, args.m_bound, progress)
    lookup_and_merge(fields, args, progress, "nested pure fields")


if __name__ == "__main__":
    main()
