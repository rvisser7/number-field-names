"""Shared helpers for the pretty name data files.

A data file data/degree<n>.txt looks like

    label|pretty_name
    text|text

    2.0.3.1|\\Q(\\sqrt{-3})
    2.0.4.1|\\Q(\\sqrt{-1})
"""

import fcntl
import math
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
    """Write (label, pretty_name) rows to a data file, sorted by label.

    The rows go to a temporary file that then replaces the data file, so a reader never
    sees a half-written file.
    """
    rows = sorted(rows, key=lambda row: label_key(row[0]))
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", newline="\n") as f:
        f.write(HEADER + "\n" + TYPES + "\n\n")
        for label, name in rows:
            f.write(f"{label}|{name}\n")
    os.replace(tmp, path)


# When a field has names from several families, the earlier family wins. This follows the
# order of the cases in the LMFDB's field_pretty, e.g. 4.0.256.1 is \Q(\zeta_{8}), not
# \Q(\sqrt[4]{-1}), and 2.0.3.1 is \Q(\sqrt{-3}), not \Q(\zeta_{3}).
# Names that match none of these rank last, so add each new family to this list.
class _Checked:
    """A regex whose match must also pass a check on its groups."""

    def __init__(self, pattern, check):
        self.rx, self.check = re.compile(pattern), check

    def match(self, s):
        m = self.rx.match(s)
        return m if m and self.check(m) else None


FAMILIES = [
    ("rational", re.compile(r"^\\Q$")),
    ("quadratic", re.compile(r"^\\Q\(\\sqrt\{-?\d+\}\)$")),
    ("cyclotomic", re.compile(r"^\\Q\(\\zeta_\{\d+\}\)$")),
    ("real cyclotomic", re.compile(r"^\\Q\(\\zeta_\{\d+\}\)\^\+$")),
    ("multiquadratic", re.compile(r"^\\Q\((i|\\sqrt\{-?\d+\})(, \\sqrt\{-?\d+\})+\)$")),
    ("pure", re.compile(r"^\\Q\(\\sqrt\[\d+\]\{-?\d+\}\)$")),
    ("imprimitive quartic", re.compile(r"^\\Q\(\\sqrt\{-?\d+ \+ \d*(i|\\sqrt\{-?\d+\})\}\)$")),
    # \zeta_{n} with \sqrt[n], or \zeta_{n/2} with \sqrt[n] for n = 2 mod 4, as radical.py writes them
    ("kummer", _Checked(r"^\\Q\((i, \\sqrt\[4\]|\\zeta_\{(\d+)\}, \\sqrt\[(\d+)\])\{-?\d+\}\)$",
                        lambda m: not m[2] or int(m[3]) in (int(m[2]), 2 * int(m[2]) if int(m[2]) % 2 else 0))),
    # radicals of coprime degrees (\sqrt{a} and i have degree 2)
    ("radical compositum",
     _Checked(r"^\\Q\((i|\\sqrt\{-?\d+\}|\\sqrt\[(\d+)\]\{-?\d+\}), \\sqrt\[(\d+)\]\{-?\d+\}\)$",
              lambda m: math.gcd(int(m[2] or 2), int(m[3])) == 1)),
    # Anything else of the form \Q(g_1, g_2), as written by named_composita.py
    ("compositum", re.compile(r"^\\Q\(.+, .+\)$")),
    ("nested radical", re.compile(r"^\\Q\(\\sqrt\[\d+\]\{-?\d+ \+ \d*(i|\\sqrt\{-?\d+\})\}\)$")),
    # \Q(\sqrt{A + B\sqrt{d_1} + C\sqrt{d_2} + D\sqrt{d_3}}), as nested_biquadratic.py writes them
    ("nested biquadratic",
     _Checked(r"^\\Q\(\\sqrt\{(-?\d+|-?\d*(?:i|\\sqrt\{-?\d+\}))( [+-] \d*(?:i|\\sqrt\{-?\d+\}))+\}\)$",
              lambda m: parse_nested_biquadratic(m[0]) is not None)),
]


def family_rank(name):
    """Position of the name's family in FAMILIES (lower wins); unknown names rank last."""
    for i, (_, rx) in enumerate(FAMILIES):
        if rx.match(name):
            return i
    return len(FAMILIES)


_ROOT = r"\\sqrt\{(-?\d+) \+ (\d*)(i|\\sqrt\{(-?\d+)\})\}"  # \sqrt{A + B\sqrt{D}}
_NESTED_RE = re.compile(r"^\\sqrt\[(\d+)\]\{(-?\d+) \+ (\d*)(i|\\sqrt\{(-?\d+)\})\}$")
_QUARTIC_RE = re.compile("^" + _ROOT + "$")


def _abd(m, first):
    """(a, b, d) from a match of _QUARTIC_RE or _NESTED_RE, whose a is group `first`."""
    a, b = int(m[first]), int(m[first + 1] or 1)
    return a, b, (-1 if m[first + 2] == "i" else int(m[first + 3]))


def _nested_key(a, b, d):
    # as in imprimitive_quartic.py and nested_radical.py
    return (max(abs(a), abs(b)), abs(b), abs(a), abs(d), a < 0, d < 0)


def _totient(n):
    out, p, m = n, 2, n
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            out -= out // p
        p += 1
    return out - out // m if m > 1 else out


def generator_key(g):
    """Order of the generators of a \\Q(g_1, g_2) name, as in named_composita.py: cyclotomic
    (by degree, i first), then quadratic, pure and imprimitive quartic, each simplest first."""
    if g == "i":
        return (0, 2, False, 4)
    m = re.match(r"^\\zeta_\{(\d+)\}$", g)
    if m:
        return (0, _totient(int(m[1])), True, int(m[1]))
    m = re.match(r"^\\sqrt\{(-?\d+)\}$", g)
    if m:
        return (1, abs(int(m[1])), int(m[1]) < 0)
    m = re.match(r"^\\sqrt\[(\d+)\]\{(-?\d+)\}$", g)
    if m:
        return (2, int(m[1]), abs(int(m[2])), int(m[2]) < 0)
    m = _QUARTIC_RE.match(g)
    if m:
        return (3,) + _nested_key(*_abd(m, 1))
    return (4, g)


# Families whose names are built from radicals. Between these, the name with the smallest
# radicands wins first, and only then the family order: Q(\zeta_{3}, \sqrt[3]{2}) rather than
# its pure name Q(\sqrt[6]{-108}), but Q(\sqrt[6]{-3}) rather than Q(\zeta_{3}, \sqrt[3]{3}).
RADICAL_FAMILIES = {"pure", "imprimitive quartic", "kummer", "radical compositum", "compositum",
                    "nested radical"}


def radicand_size(name):
    """The largest |integer| in a name, leaving out root indices and \\zeta subscripts."""
    s = re.sub(r"\\sqrt\[\d+\]|\\zeta_\{\d+\}", "", name)
    return max((abs(int(k)) for k in re.findall(r"-?\d+", s)), default=0)


def name_key(name):
    """Sort key for names of one field: the smaller key is the better name.

    Names of different families are ordered as in FAMILIES, except that between the
    RADICAL_FAMILIES the smallest radicand_size comes first. Within a family the name with
    the simplest parameters wins, by the rule the generate_*.py scripts use to choose (e.g.
    smallest |m|, positive first, for \\Q(\\sqrt[n]{m})), so that the name kept does not depend
    on which script found it first. Names of no known family are never replaced.
    """
    rank = family_rank(name)
    family = FAMILIES[rank][0] if rank < len(FAMILIES) else None
    if family is None:
        return (rank, 0, rank, ())
    tier = (min(i for i, (f, _) in enumerate(FAMILIES) if f in RADICAL_FAMILIES), radicand_size(name)) \
        if family in RADICAL_FAMILIES else (rank, 0)
    return tier + (rank, _family_key(name, family))


def nested_biquadratic_key(A, B, C, D, d1, d2):
    # as in nested_biquadratic.py: the smallest coefficients, then the smallest basis
    return (max(abs(A), abs(B), abs(C), abs(D)), abs(D), abs(C), abs(B), abs(A),
            (abs(d1), -d1), (abs(d2), -d2), B < 0, C < 0, A < 0, D < 0)


def parse_nested_biquadratic(name):
    """(A, B, C, D, d_1, d_2) from a name written by nested_biquadratic.py, or None."""
    from multiquadratic import order_key, squarefree_part  # (multiquadratic imports this module)
    m = re.match(r"^\\Q\(\\sqrt\{(.*)\}\)$", name)
    if not m:
        return None
    terms = re.findall(r"(^-?| [+-] )(\d*)(i|\\sqrt\{(-?\d+)\})?", m[1])
    if "".join(sign + c + (sym or "") for sign, c, sym, _ in terms) != m[1]:
        return None
    A, coeff = 0, {}
    for sign, c, sym, d in terms:
        k = (-1 if "-" in sign else 1) * int(c or 1)
        if not sym:
            if A or coeff or not c:
                return None  # the rational term comes first, once
            A = k
            continue
        d = -1 if sym == "i" else int(d)
        if d in coeff or d in (0, 1) or squarefree_part(d) != d or c in ("0", "1") or c.startswith("0"):
            return None
        coeff[d] = k
    if len(coeff) < 2 or len(coeff) > 3:
        return None
    x, y = list(coeff)[:2]
    ds = sorted({x, y, squarefree_part(x * y)}, key=order_key)
    if len(ds) != 3 or set(coeff) - set(ds) or list(coeff) != [d for d in ds if d in coeff]:
        return None  # not in the order sqrt(d_1), sqrt(d_2), sqrt(d_3)
    return (A,) + tuple(coeff.get(d, 0) for d in ds) + (ds[0], ds[1])


def _family_key(name, family):
    inner = name[3:-1] if name.startswith("\\Q(") and name.endswith(")") else name
    if family == "nested biquadratic":
        return nested_biquadratic_key(*parse_nested_biquadratic(name))
    if family == "imprimitive quartic":
        return _nested_key(*_abd(_QUARTIC_RE.match(inner), 1))
    if family == "nested radical":
        m = _NESTED_RE.match(inner)
        return (int(m[1]),) + _nested_key(*_abd(m, 2))
    if family == "compositum":  # fewer generators first
        gens = inner.split(", ")
        return (len(gens),) + tuple(generator_key(g) for g in gens)
    return tuple((abs(k), k < 0) for k in map(int, re.findall(r"-?\d+", name)))


def merge_data(path, new_rows, replace=False):
    """Merge {label: name} into the data file at path, and return (added, replaced, kept).

    A label already in the file keeps its name, unless the new name is better (see
    name_key: a family that wins over it, or a simpler name of the same family) or
    replace is True. So without replace, the result does not depend on the order of merges.

    The read-merge-write holds a lock on the data directory, so several scripts can run
    at once without losing each other's rows.
    """
    with open(os.path.join(os.path.dirname(os.path.abspath(path)), ".merge.lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)  # released when the file is closed
        return _merge_data(path, new_rows, replace)


def _merge_data(path, new_rows, replace):
    rows = dict(read_data(path)) if os.path.exists(path) else {}
    added = replaced = kept = 0
    for label, name in new_rows.items():
        old = rows.get(label)
        if old is None:
            added += 1
        elif old == name:
            continue
        elif replace or name_key(name) < name_key(old):
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
