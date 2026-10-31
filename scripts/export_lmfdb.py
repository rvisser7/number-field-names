#!/usr/bin/env python3
"""Export the data files with full LaTeX names, for uploading to the LMFDB.

The data files store names in a compact form (Q(S{1+2i})); this writes them out with the
full macros and spacing (\\Q(\\sqrt{1 + 2i})), in the same psycodict upload format
(label|pretty_name, text|text, a blank line, then the rows in LMFDB order). The data
files are not changed.

Usage:
    # one file per degree, export/degree<n>.txt
    python3 scripts/export_lmfdb.py -o export
    # also every degree in one file
    python3 scripts/export_lmfdb.py -o export --combined export/pretty_names.txt
    # only some degrees
    python3 scripts/export_lmfdb.py -o export --degrees 2-8
"""

import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nf_names import DATA_DIR, data_path, read_data, write_data

FILE_RE = re.compile(r"^degree(\d+)\.txt$")


def parse_degrees(s):
    """'2-5,8' -> {2, 3, 4, 5, 8}"""
    out = set()
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out.update(range(int(lo), int(hi or lo) + 1))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--data-dir", default=DATA_DIR)
    parser.add_argument("-o", "--output-dir", help="write <dir>/degree<n>.txt for each data file")
    parser.add_argument("--combined", help="also write every degree into this one file")
    parser.add_argument("--degrees", help="only these degrees, e.g. 2-8 or 3,4,6 (default: all)")
    args = parser.parse_args()
    if not args.output_dir and not args.combined:
        parser.error("give --output-dir and/or --combined")
    if args.output_dir and os.path.abspath(args.output_dir) == os.path.abspath(args.data_dir):
        parser.error("the output directory must not be the data directory")

    degrees = parse_degrees(args.degrees) if args.degrees else None
    files = []
    for path in glob.glob(os.path.join(args.data_dir, "degree*.txt")):
        m = FILE_RE.match(os.path.basename(path))
        if m and (degrees is None or int(m[1]) in degrees):
            files.append((int(m[1]), path))
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
    every = []
    for n, path in sorted(files):
        rows = read_data(path)
        every += rows
        if args.output_dir:
            write_data(data_path(n, args.output_dir), rows, compact=False)
            print(f"degree {n}: {len(rows)} names -> {data_path(n, args.output_dir)}")
    if args.combined:
        if os.path.dirname(args.combined):
            os.makedirs(os.path.dirname(args.combined), exist_ok=True)
        write_data(args.combined, every, compact=False)
        print(f"{len(every)} names -> {args.combined}")


if __name__ == "__main__":
    main()
