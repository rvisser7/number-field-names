#!/usr/bin/env python3
"""Rewrite the data files with their names in the compact form.

The data files store names compactly (Q, S and Z for \\Q, \\sqrt and \\zeta, and no spaces,
e.g. Q(S{1+2i}) for \\Q(\\sqrt{1 + 2i}); see nf_names.compress_name). This converts data files
written before that, in place. Any merge also converts a file it writes to, so this is only
needed once, for files no script has written to since. Running it again changes nothing.

Usage:
    python3 scripts/compact_data.py              # convert <repo>/data/degree*.txt
    python3 scripts/compact_data.py -d /path/to/data
"""

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nf_names import DATA_DIR, compress_name, is_compact, merge_lock, read_data, write_data


def compact_file(path):
    """Convert one data file, and return (names compressed, names, names kept as they are)."""
    with merge_lock(path):
        raw = read_data(path, expand=False)
        rows = read_data(path)
        changed = sum(1 for (_, old), (_, name) in zip(raw, rows) if compress_name(name) != old)
        kept = sum(1 for _, name in rows if not is_compact(compress_name(name)))
        if changed:
            write_data(path, rows)
    return changed, len(rows), kept


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--data-dir", default=DATA_DIR)
    args = parser.parse_args()

    paths = sorted(glob.glob(os.path.join(args.data_dir, "degree*.txt")))
    total_changed = total = 0
    for path in paths:
        before = os.path.getsize(path)
        changed, n, kept = compact_file(path)
        total_changed += changed
        total += n
        note = f", {kept} kept as they are (not of the compact form)" if kept else ""
        print(f"{os.path.basename(path)}: {changed} of {n} names compressed{note}, "
              f"{before} -> {os.path.getsize(path)} bytes")
    print(f"{len(paths)} files, {total_changed} of {total} names compressed")


if __name__ == "__main__":
    main()
