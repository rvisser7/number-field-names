"""Pure fields Q(\\sqrt[n]{m}), computed with PARI (via cypari2).

For a degree n and a bound B we run over every m with 1 <= |m| <= B, keep those
with x^n - m irreducible, and identify the field by its polredabs polynomial.
Different m can give the same field (e.g. Q(\\sqrt[3]{4}) = Q(\\sqrt[3]{2}), and
Q(\\sqrt[3]{-2}) = Q(\\sqrt[3]{2})); each field is named by the m of smallest |m|,
preferring m > 0, so the name is the same however the field was reached.
"""

from pari_fields import is_isomorphic, make_field, pari, polredabs_coeffs  # noqa: F401


def pure_name(n, m):
    return r"\Q(\sqrt[%d]{%d})" % (n, m)


def pure_fields(n, bound, progress=None):
    """Return the distinct pure fields Q(\\sqrt[n]{m}) with |m| <= bound, as Fields.

    If given, progress(done, total, nfields) is called after each m.
    """
    fields = {}
    ms = [m for m in sorted(range(-bound, bound + 1), key=lambda m: (abs(m), -m)) if m != 0]
    for i, m in enumerate(ms, start=1):
        if progress:
            progress(i, len(ms), len(fields))
        f = pari(f"x^{n} - ({m})")
        if not pari.polisirreducible(f):
            continue
        coeffs = polredabs_coeffs(f)
        if coeffs in fields:
            continue  # already reached with a smaller |m|
        fields[coeffs] = make_field(pure_name(n, m), coeffs)
    return list(fields.values())
