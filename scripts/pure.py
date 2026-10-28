"""Pure fields Q(\\sqrt[n]{m}), computed with PARI (via cypari2).

For a degree n and a bound B we run over every m with 1 <= |m| <= B, keep those
with x^n - m irreducible, and identify the field by its polredabs polynomial.
Different m can give the same field (e.g. Q(\\sqrt[3]{4}) = Q(\\sqrt[3]{2}), and
Q(\\sqrt[3]{-2}) = Q(\\sqrt[3]{2})); each field is named by the m of smallest |m|,
preferring m > 0, so the name is the same however the field was reached.
"""

from dataclasses import dataclass

try:
    import cypari2
except ImportError:
    raise SystemExit("cypari2 is needed: pip install cypari2, or run this script with sage -python")

pari = cypari2.Pari()


@dataclass(frozen=True)
class PureField:
    n: int
    m: int
    coeffs: tuple    # polredabs polynomial, constant term first (as in the LMFDB)
    disc: int        # field discriminant
    r1: int          # number of real places

    @property
    def name(self):
        return pure_name(self.n, self.m)


def pure_name(n, m):
    return r"\Q(\sqrt[%d]{%d})" % (n, m)


def pure_fields(n, bound, progress=None):
    """Return the distinct pure fields Q(\\sqrt[n]{m}) with |m| <= bound, as PureFields.

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
        g = pari.polredabs(f)
        coeffs = tuple(int(c) for c in pari.Vecrev(g))
        if coeffs in fields:
            continue  # already reached with a smaller |m|
        fields[coeffs] = PureField(n, m, coeffs, int(pari.nfdisc(g)), int(pari.polsturm(g)))
    return list(fields.values())


def is_isomorphic(coeffs1, coeffs2):
    """True if the polynomials (constant term first) define isomorphic fields."""
    return pari.nfisisom(pari.Polrev(list(coeffs1)), pari.Polrev(list(coeffs2))).type() != "t_INT"
