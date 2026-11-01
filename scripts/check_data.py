#!/usr/bin/env python3
"""Check that every data/degree<n>.txt is well formed.

For each file: the header and types lines are right, every row has two columns,
labels are valid, of degree n, sorted and unique, names are non-empty, stored in the
compact form (no backslashes or spaces; run compact_data.py to convert older files), have
balanced brackets and no \\( \\) delimiters. Exits nonzero if anything is wrong.

With --fields it also checks that each name describes the field of its label: the name is
read as LaTeX (independently of the scripts that wrote it, see name_field.py) and the field it
gives must have the degree, signature and |discriminant| of the label; with --source it must
also be isomorphic to the LMFDB field of the label. This needs cypari2 (sage -python), and
--drop-wrong removes the rows whose name gives another field.

Usage:
    python3 scripts/check_data.py              # check <repo>/data
    python3 scripts/check_data.py -d /path/to/data
    sage -python scripts/check_data.py --fields --degrees 3-8               # offline
    sage -python scripts/check_data.py --fields --source db --shard 0/4     # against the LMFDB
"""

import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nf_names import DATA_DIR, compress_name, is_compact, label_key, merge_lock, read_data, write_data

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


def _stored_form(name):
    """How write_data would store name (None if it cannot)."""
    try:
        return compress_name(name)
    except ValueError:
        return None


def check_file(path, degree):
    """Return a list of error strings for one data file."""
    try:
        rows = read_data(path, expand=False)
    except ValueError as e:
        return [str(e)]
    errors = []
    seen = set()
    prev = None
    not_compact = 0
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
        elif not is_compact(name) and _stored_form(name) != name:
            if not not_compact:
                first = f"{label}|{name}"
            not_compact += 1
    if not_compact:  # one line, as an old file has every name in the full form
        errors.append(f"{path}: {not_compact} names are not in the compact form, e.g. {first} "
                      f"(run scripts/compact_data.py to convert them)")
    return errors


def parse_degrees(s):
    """'3-5,8' -> {3, 4, 5, 8}"""
    out = set()
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out.update(range(int(lo), int(hi or lo) + 1))
    return out


def check_fields(path, args, progress):
    """Check that the names in one data file describe the fields of their labels; return
    (number checked, [(label, name, why)] for the wrong ones, [(label, name, why)] for those
    that could not be checked: unreadable names, labels not in the LMFDB)."""
    from name_field import check_name  # needs cypari2
    rows = read_data(path)
    shard, nshards = (int(x) for x in args.shard.split("/"))
    rows = rows[shard::nshards]
    source = None
    if args.source:
        from lmfdb_lookup import api_source, db_source
        source = api_source if args.source == "api" else db_source()
    cache = {}
    wrong, unchecked = [], []
    progress.start(os.path.basename(path))
    for k, (label, name) in enumerate(rows, start=1):
        lmfdb = None
        if source:
            n, _, D, _ = (int(x) for x in label.split("."))
            if (n, D) not in cache:
                cache = {(n, D): dict(source(n, D))}  # rows are in |disc| order
            lmfdb = cache[(n, D)].get(label)
            if lmfdb is None:
                unchecked.append((label, name, "the label is not in the LMFDB"))
                continue
        why = check_name(label, name, lmfdb)
        if why is not None:
            (unchecked if why.startswith("cannot read") else wrong).append((label, name, why))
        progress(k, len(rows), f", {len(wrong)} wrong")
    return len(rows), wrong, unchecked


def drop_rows(path, labels):
    with merge_lock(path):
        rows = [(label, name) for label, name in read_data(path) if label not in labels]
        write_data(path, rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--data-dir", default=DATA_DIR)
    parser.add_argument("--fields", action="store_true",
                        help="also check that each name gives the field of its label (needs cypari2)")
    parser.add_argument("--source", choices=["api", "db"],
                        help="with --fields, also compare with the LMFDB polynomial of each label")
    parser.add_argument("--degrees", help="with --fields, only these degrees, e.g. 3-8")
    parser.add_argument("--shard", default="0/1", help="with --fields, i/k: only every k-th row, from the i-th")
    parser.add_argument("--drop-wrong", action="store_true",
                        help="with --fields, remove the rows whose name gives another field "
                             "(not those that could not be checked)")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS")
    args = parser.parse_args()
    if (args.source or args.degrees or args.drop_wrong) and not args.fields:
        parser.error("--source, --degrees and --drop-wrong go with --fields")

    errors = []
    nfiles = 0
    paths = []
    for path in sorted(glob.glob(os.path.join(args.data_dir, "degree*.txt"))):
        m = FILE_RE.match(os.path.basename(path))
        if not m:
            continue
        nfiles += 1
        errors += check_file(path, int(m.group(1)))
        paths.append((int(m.group(1)), path))
    for e in errors:
        print(e)
    print(f"checked {nfiles} files, {len(errors)} errors")

    bad = 0
    if args.fields:
        from progress import Progress
        progress = Progress(every=args.progress_every)
        degrees = parse_degrees(args.degrees) if args.degrees else None
        for n, path in sorted(paths):
            if degrees is not None and n not in degrees:
                continue
            checked, wrong, unchecked = check_fields(path, args, progress)
            for label, name, why in wrong + unchecked:
                print(f"{path}: {label}|{name}: {why}")
            print(f"degree {n}: checked {checked} names, {len(wrong)} wrong, {len(unchecked)} could not be checked"
                  + ("" if args.source else " (degree, signature and |disc| only)"))
            bad += len(wrong) + len(unchecked)
            if args.drop_wrong and wrong:
                drop_rows(path, {label for label, _, _ in wrong})
                print(f"  removed the {len(wrong)} wrong rows from {path}")
    sys.exit(1 if errors or bad else 0)


if __name__ == "__main__":
    main()
