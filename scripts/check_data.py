#!/usr/bin/env python3
"""Check that every data/degree<n>.txt is well formed.

For each file: the header and types lines are right, every row has two columns,
labels are valid, of degree n, sorted and unique, names are non-empty, have balanced
brackets and no \\( \\) delimiters. Exits nonzero if anything is wrong.

Usage:
    python3 scripts/check_data.py              # check <repo>/data
    python3 scripts/check_data.py -d /path/to/data
"""

import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nf_names import DATA_DIR, label_key, read_data

LABEL_RE = re.compile(r"^(\d+)\.(\d+)\.(\d+)\.(\d+)$")
FILE_RE = re.compile(r"^degree(\d+)\.txt$")


def brackets_balanced(s):
    """True if the (), [] and {} in s are properly nested."""
    pairs = {")": "(", "]": "[", "}": "{"}
    stack = []
    for c in s:
        if c in "([{":
            stack.append(c)
        elif c in pairs:
            if not stack or stack.pop() != pairs[c]:
                return False
    return not stack


def check_file(path, degree):
    """Return a list of error strings for one data file."""
    try:
        rows = read_data(path)
    except ValueError as e:
        return [str(e)]
    errors = []
    seen = set()
    prev = None
    for label, name in rows:
        m = LABEL_RE.match(label)
        if not m:
            errors.append(f"{path}: bad label {label!r}")
            continue
        if int(m.group(1)) != degree:
            errors.append(f"{path}: {label} is not of degree {degree}")
        if label in seen:
            errors.append(f"{path}: duplicate label {label}")
        seen.add(label)
        if prev is not None and label_key(label) < label_key(prev):
            errors.append(f"{path}: {label} is out of order (after {prev})")
        prev = label
        if not name:
            errors.append(f"{path}: {label} has an empty name")
        elif not brackets_balanced(name):
            errors.append(f"{path}: {label} has unbalanced brackets: {name}")
        elif r"\(" in name or r"\)" in name:
            errors.append(f"{path}: {label} name should not contain \\( \\): {name}")
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--data-dir", default=DATA_DIR)
    args = parser.parse_args()

    errors = []
    nfiles = 0
    for path in sorted(glob.glob(os.path.join(args.data_dir, "degree*.txt"))):
        m = FILE_RE.match(os.path.basename(path))
        if not m:
            continue
        nfiles += 1
        errors += check_file(path, int(m.group(1)))
    for e in errors:
        print(e)
    print(f"checked {nfiles} files, {len(errors)} errors")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
