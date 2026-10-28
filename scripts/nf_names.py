"""Shared helpers for the pretty name data files.

A data file data/degree<n>.txt looks like

    label|pretty_name
    text|text

    2.0.3.1|\\Q(\\sqrt{-3})
    2.0.4.1|\\Q(\\sqrt{-1})
"""

import os
import re

HEADER = "label|pretty_name"
TYPES = "text|text"

# Resolve data/ from this file, not from the cwd.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "data")


def data_path(degree, data_dir=DATA_DIR):
    return os.path.join(data_dir, f"degree{degree}.txt")


def label_key(label):
    """Sort key for an LMFDB number field label n.r.D.i, matching the LMFDB's default
    order: degree, |disc|, disc sign, iso number. The sign is (-1)^s with s = (n - r)/2
    complex places, so e.g. 2.0.8.1 comes before 2.2.8.1. Fields of the same degree,
    |disc|, sign and iso number but different signature (e.g. 4.0.D.i and 4.4.D.i)
    have no LMFDB order; we break the tie by r.
    """
    n, r, D, i = (int(x) for x in label.split("."))
    sign = -1 if (n - r) // 2 % 2 else 1
    return (n, D, sign, i, r)


def read_data(path):
    """Return a list of (label, pretty_name) rows from a data file."""
    with open(path) as f:
        lines = f.read().split("\n")
    if lines[:3] != [HEADER, TYPES, ""]:
        raise ValueError(f"{path}: expected header {HEADER!r}, {TYPES!r} and a blank line")
    rows = []
    for n, ln in enumerate(lines[3:], start=4):
        if ln == "":
            continue
        cells = ln.split("|")
        if len(cells) != 2:
            raise ValueError(f"{path}:{n}: expected 2 columns, got {len(cells)}: {ln!r}")
        rows.append((cells[0], cells[1]))
    return rows


def write_data(path, rows):
    """Write (label, pretty_name) rows to a data file, sorted by label."""
    rows = sorted(rows, key=lambda row: label_key(row[0]))
    with open(path, "w", newline="\n") as f:
        f.write(HEADER + "\n" + TYPES + "\n\n")
        for label, name in rows:
            f.write(f"{label}|{name}\n")


# When a field has names from several families, the earlier family wins. This follows the
# order of the cases in the LMFDB's field_pretty, e.g. 4.0.256.1 is \Q(\zeta_{8}), not
# \Q(\sqrt[4]{-1}), and 2.0.3.1 is \Q(\sqrt{-3}), not \Q(\zeta_{3}).
# Names that match none of these rank last, so add each new family to this list.
FAMILIES = [
    ("rational", re.compile(r"^\\Q$")),
    ("quadratic", re.compile(r"^\\Q\(\\sqrt\{-?\d+\}\)$")),
    ("cyclotomic", re.compile(r"^\\Q\(\\zeta_\{\d+\}\)$")),
    ("real cyclotomic", re.compile(r"^\\Q\(\\zeta_\{\d+\}\)\^\+$")),
    ("pure", re.compile(r"^\\Q\(\\sqrt\[\d+\]\{-?\d+\}\)$")),
    ("imprimitive quartic", re.compile(r"^\\Q\(\\sqrt\{-?\d+ \+ \d*(i|\\sqrt\{-?\d+\})\}\)$")),
]


def family_rank(name):
    """Position of the name's family in FAMILIES (lower wins); unknown names rank last."""
    for i, (_, rx) in enumerate(FAMILIES):
        if rx.match(name):
            return i
    return len(FAMILIES)


def merge_data(path, new_rows, replace=False):
    """Merge {label: name} into the data file at path, and return (added, replaced, kept).

    A label already in the file keeps its name, unless the new name is from a family that
    wins over it (see FAMILIES) or replace is True.
    """
    rows = dict(read_data(path)) if os.path.exists(path) else {}
    added = replaced = kept = 0
    for label, name in new_rows.items():
        old = rows.get(label)
        if old is None:
            added += 1
        elif old == name:
            continue
        elif replace or family_rank(name) < family_rank(old):
            replaced += 1
        else:
            kept += 1
            continue
        rows[label] = name
    write_data(path, rows.items())
    return added, replaced, kept


def is_squarefree(n):
    n = abs(n)
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return False
        if n % p == 0:
            n //= p
        p += 1
    return True


def is_fundamental_discriminant(D):
    if D in (0, 1):
        return False
    if D % 4 == 1:
        return is_squarefree(D)
    if D % 4 == 0:
        return (D // 4) % 4 in (2, 3) and is_squarefree(D // 4)
    return False


def quadratic_label(D):
    """LMFDB label of Q(sqrt(D)) for a fundamental discriminant D."""
    return f"2.{2 if D > 0 else 0}.{abs(D)}.1"


def quadratic_name(label):
    """Pretty name of a quadratic field from its label, matching LMFDB's field_pretty.

    Note 2.0.4.1 is \\Q(\\sqrt{-1}), not \\Q(i).
    """
    n, r, disc, _ = label.split(".")
    if n != "2" or r not in ("0", "2"):
        raise ValueError(f"not a quadratic field label: {label}")
    D = int(disc) if r == "2" else -int(disc)
    if not is_fundamental_discriminant(D):
        raise ValueError(f"not a fundamental discriminant: {label}")
    return r"\Q(\sqrt{%d})" % (D if D % 4 else D // 4)
