#!/usr/bin/env python3
"""How many fields have a name, per degree and per family.

Reads data/degree<n>.txt. With --db, also counts the LMFDB fields of each degree (through the
lmfdb package or lmfdb-lite) to show what fraction is named.

Usage:
    python3 scripts/coverage.py
    python3 scripts/coverage.py --db
"""

import argparse
import glob
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nf_names import DATA_DIR, FAMILIES, family_rank, read_data


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", action="store_true", help="compare with the LMFDB field counts")
    parser.add_argument("-d", "--data-dir", default=DATA_DIR)
    args = parser.parse_args()

    if args.db:
        try:
            from lmfdb import db
        except ImportError:
            from lmf import db

    families = [name for name, _ in FAMILIES] + ["other"]
    total = Counter()
    print(f"{'degree':>6} {'named':>9}" + (f" {'in LMFDB':>10} {'named %':>8}" if args.db else "")
          + "  by family")
    paths = glob.glob(os.path.join(args.data_dir, "degree*.txt"))
    for n, path in sorted((int(re.search(r"degree(\d+)", p)[1]), p) for p in paths):
        by_family = Counter(families[family_rank(name)] for _, name in read_data(path))
        named = sum(by_family.values())
        total.update(by_family)
        line = f"{n:>6} {named:>9}"
        if args.db:
            count = db.nf_fields.count({"degree": n})
            line += f" {count:>10} {100 * named / count:>7.2f}%"
        print(line + "  " + ", ".join(f"{f}: {c}" for f, c in by_family.most_common()))
    print(f"\n{'all':>6} {sum(total.values()):>9}  "
          + ", ".join(f"{f}: {c}" for f, c in total.most_common()))


if __name__ == "__main__":
    main()
