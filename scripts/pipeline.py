"""The shared tail of the generate_*.py scripts: look Fields up in the LMFDB and merge
their names into the data files."""

from collections import defaultdict

from lmfdb_lookup import api_source, db_source, label_fields
from nf_names import data_path, merge_data
from pari_fields import is_isomorphic
from progress import Progress


def add_common_args(parser):
    parser.add_argument("--source", choices=["api", "db"], default="api")
    parser.add_argument("--replace", action="store_true",
                        help="overwrite names already in the data files")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")


def lookup_and_merge(fields, args, progress, what):
    """Look up fields (a list of Fields) in the LMFDB, merge the names found into the data
    files, print a summary, and list the fields that are not in the LMFDB."""
    source = api_source if args.source == "api" else db_source()
    progress.start(f"{what}, LMFDB lookups")
    found, missing = label_fields(fields, source, is_isomorphic, progress)
    print(f"{what}: {len(fields)} fields, {len(found)} in the LMFDB")

    by_degree = defaultdict(dict)
    for label, name in found.items():
        by_degree[int(label.split(".")[0])][label] = name
    for d in sorted(by_degree):
        added, replaced, kept = merge_data(data_path(d), by_degree[d], args.replace)
        print(f"  degree {d}: {len(by_degree[d])} fields; "
              f"{added} added, {replaced} replaced, {kept} kept existing name")

    if missing:
        print(f"\n{len(missing)} fields not found in the LMFDB:")
        for K in missing:
            print(f"  {K.name}  degree {K.degree}  disc {K.disc}")


def new_progress(args):
    return Progress(every=args.progress_every)
