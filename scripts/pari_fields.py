"""PARI (via cypari2) and a small record for a field we have computed."""

from dataclasses import dataclass

try:
    import cypari2
except ImportError:
    raise SystemExit("cypari2 is needed: pip install cypari2, or run this script with sage -python")

pari = cypari2.Pari()


@dataclass(frozen=True)
class Field:
    name: str        # pretty name, e.g. \Q(\zeta_{7})
    coeffs: tuple    # polredabs polynomial, constant term first (as in the LMFDB)
    disc: int        # field discriminant
    r1: int          # number of real places

    @property
    def degree(self):
        return len(self.coeffs) - 1


def polredabs_coeffs(f):
    """Coefficients (constant term first) of the polredabs polynomial of the PARI polynomial f."""
    return tuple(int(c) for c in pari.Vecrev(pari.polredabs(f)))


def make_field(name, coeffs):
    """Field record for the polredabs polynomial with the given coefficients."""
    g = pari.Polrev(list(coeffs))
    return Field(name, tuple(coeffs), int(pari.nfdisc(g)), int(pari.polsturm(g)))


def is_isomorphic(coeffs1, coeffs2):
    """True if the polynomials (constant term first) define isomorphic fields."""
    return pari.nfisisom(pari.Polrev(list(coeffs1)), pari.Polrev(list(coeffs2))).type() != "t_INT"
