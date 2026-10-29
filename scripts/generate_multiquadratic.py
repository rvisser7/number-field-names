#!/usr/bin/env python3
"""Add multiquadratic fields Q(\\sqrt{d_1}, ..., \\sqrt{d_k}) to the data files.

For each degree 2^k, runs over the fields whose greedy basis (see multiquadratic.py) has
every |d_i| <= bound, looks up their labels in the LMFDB, and merges the names into
data/degree<2^k>.txt. Names of families that win over these (e.g. \\Q(\\zeta_{8}) for
4.0.256.1, see FAMILIES in nf_names.py) are kept unless --replace is given.

The discriminant of a multiquadratic field is known without computing the field, so we
only build a polynomial for fields that have an LMFDB field of the same degree,
discriminant and signature. The rest are counted as not in the LMFDB.

The number of fields grows quickly with the rank, so use a smaller bound for higher degrees,
e.g. run degree 4 with bound 200 but degree 32 with bound 15.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_multiquadratic.py --degrees 4 --bound 200 --source db
    python3 scripts/generate_multiquadratic.py --degrees 8,16 --bound 40 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lmfdb_lookup import api_source, db_source, find_label
from multiquadratic import disc_and_r1, greedy_bases, multiquadratic_coeffs, multiquadratic_name
from nf_names import data_path, merge_data
from pari_fields import is_isomorphic
from progress import Progress


def parse_degrees(s):
    degrees = [int(d) for d in s.split(",")]
    for d in degrees:
        if d not in (4, 8, 16, 32):
            sys.exit(f"degree {d}: multiquadratic fields here have degree 4, 8, 16 or 32")
    return degrees


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", required=True, help="e.g. 4 or 4,8,16")
    parser.add_argument("--bound", type=int, required=True, help="largest |d_i|")
    parser.add_argument("--source", choices=["api", "db"], default="api")
    parser.add_argument("--replace", action="store_true",
                        help="overwrite names already in the data files")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    args = parser.parse_args()
    source = api_source if args.source == "api" else db_source()

    progress = Progress(every=args.progress_every)
    unmatched = []
    for n in parse_degrees(args.degrees):
        bases = list(greedy_bases(n.bit_length() - 1, args.bound))
        progress.start(f"degree {n}, LMFDB lookups")
        found, absent = {}, []
        for i, basis in enumerate(bases, start=1):
            disc, r1 = disc_and_r1(basis)
            candidates = [(label, c) for label, c in source(n, abs(disc))
                          if label.split(".")[1] == str(r1)]
            if not candidates:
                absent.append(basis)
            else:
                label = find_label(multiquadratic_coeffs(basis), r1, candidates, is_isomorphic)
                if label is None:
                    unmatched.append(basis)
                else:
                    found[label] = multiquadratic_name(basis)
            progress(i, len(bases), f", {len(found)} found")
        added, replaced, kept = merge_data(data_path(n), found, args.replace)
        print(f"degree {n}: {len(bases)} fields, {len(found)} in the LMFDB, {len(absent)} with no "
              f"LMFDB field of that discriminant; {added} added, {replaced} replaced, "
              f"{kept} kept existing name")
        if absent:
            print("  e.g. " + ", ".join(multiquadratic_name(b) for b in absent[:5]))

    if unmatched:
        # An LMFDB field has the right discriminant and signature but is not this field.
        print(f"\n{len(unmatched)} fields not found in the LMFDB:")
        for basis in unmatched:
            print(f"  {multiquadratic_name(basis)}  disc {disc_and_r1(basis)[0]}")


if __name__ == "__main__":
    main()
