#!/usr/bin/env python3
"""Add imprimitive quartic fields Q(\\sqrt{A + B\\sqrt{D}}) to data/degree4.txt.

Runs over |A|, |B|, |D| <= bound, finds the distinct fields and the simplest (A, B, D) for
each (see imprimitive_quartic.py), looks up their labels in the LMFDB, and merges the names
into data/degree4.txt. Existing names of families that win over these (see FAMILIES in
nf_names.py, e.g. \\Q(\\zeta_{5}) for 4.0.125.1) are kept unless --replace is given.
Fields that are not in the LMFDB are listed at the end.

Needs cypari2: pip install cypari2, or run with sage -python.

Usage:
    python3 scripts/generate_imprimitive_quartic.py --bound 50 --source db
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imprimitive_quartic import imprimitive_quartic_fields
from lmfdb_lookup import api_source, db_source, label_fields
from nf_names import data_path, merge_data
from pari_fields import is_isomorphic
from progress import Progress


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bound", type=int, required=True, help="largest |A|, |B| and |D|")
    parser.add_argument("--source", choices=["api", "db"], default="api")
    parser.add_argument("--replace", action="store_true",
                        help="overwrite names already in the data files")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    args = parser.parse_args()
    source = api_source if args.source == "api" else db_source()

    progress = Progress(every=args.progress_every)
    progress.start("computing fields (by D)")
    fields = imprimitive_quartic_fields(args.bound, progress)
    progress.start("LMFDB lookups")
    found, missing = label_fields(fields, source, is_isomorphic, progress)
    added, replaced, kept = merge_data(data_path(4), found, args.replace)
    print(f"{len(fields)} imprimitive quartic fields, {len(found)} in the LMFDB; "
          f"{added} added, {replaced} replaced, {kept} kept existing name")

    if missing:
        print(f"\n{len(missing)} fields not found in the LMFDB:")
        for K in missing:
            print(f"  {K.name}  disc {K.disc}")


if __name__ == "__main__":
    main()
