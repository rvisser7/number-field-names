"""Octic fields Q(\\sqrt{A + B\\sqrt{d_1} + C\\sqrt{d_2} + D\\sqrt{d_3}}), computed with PARI (via cypari2).

These are the quadratic extensions K = L(\\sqrt{beta}) of a biquadratic field
L = Q(\\sqrt{d_1}, \\sqrt{d_2}), the analogue of the imprimitive quartics Q(\\sqrt{A + B\\sqrt{D}})
(imprimitive_quartic.py) one level up. Every octic field with a biquadratic subfield is one,
e.g. the quaternion field Q(\\sqrt{(2 + \\sqrt{2})(3 + \\sqrt{3})}) = Q(\\sqrt{6 + 3\\sqrt{2} + 2\\sqrt{3} + \\sqrt{6}}).

Conventions:

  * d_1, d_2 is the greedy basis of L (as in multiquadratic.py, so |d_1| <= |d_2|), and d_3 is
    the squarefree part of d_1 d_2; \\sqrt{-1} is written i. So \\sqrt{d_1}\\sqrt{d_2} = g\\sqrt{d_3}
    with d_1 d_2 = g^2 d_3, and the coefficient D is that of \\sqrt{d_3}.
  * beta = A + B\\sqrt{d_1} + C\\sqrt{d_2} + D\\sqrt{d_3} must not lie in a quadratic subfield of
    L, so at least two of B, C, D are nonzero (otherwise Q(\\sqrt{beta}) has degree 4).
  * Of the (A, B, C, D) giving the same field (beta times squares of L, and the conjugates of
    beta, which flip the signs of two of B, C, D), we keep the simplest: smallest
    max(|A|, |B|, |C|, |D|), then |D|, |C|, |B|, |A|, then the smallest basis, then B >= 0,
    C >= 0, A >= 0, D >= 0 (see simplicity).

Many of these fields also have a name of another family (multiquadratic, cyclotomic,
composita such as Q(i, \\sqrt{1 + \\sqrt{2}}), ...). This family comes last in FAMILIES, so
those names win, and only octics with no other name are named this way.
"""

import itertools

from multiquadratic import greedy_bases, sqrt_symbol, squarefree_part
from nf_names import nested_biquadratic_key
from pari_fields import make_field, pari, polredabs_coeffs


def nested_biquadratic_name(A, B, C, D, d1, d2):
    d3 = squarefree_part(d1 * d2)
    s = ""
    for c, d in [(A, None), (B, d1), (C, d2), (D, d3)]:
        if c == 0:
            continue
        body = str(abs(c)) if d is None else ("" if abs(c) == 1 else str(abs(c))) + sqrt_symbol(d)
        s += ("-" if c < 0 else "") + body if not s else (" - " if c < 0 else " + ") + body
    return r"\Q(\sqrt{%s})" % s


def simplicity(A, B, C, D, d1, d2):
    return nested_biquadratic_key(A, B, C, D, d1, d2)


def conjugates(A, B, C, D):
    """beta and its conjugates over Q: sqrt(d_1) -> -sqrt(d_1) and/or sqrt(d_2) -> -sqrt(d_2)."""
    return [(A, B, C, D), (A, -B, C, -D), (A, B, -C, -D), (A, -B, -C, D)]


def is_valid(A, B, C, D):
    """Whether beta lies in no quadratic subfield of L."""
    return sum(1 for c in (B, C, D) if c) >= 2


def nested_biquadratic_poly(A, B, C, D, d1, d2):
    """The minimal polynomial (if irreducible) of \\sqrt{beta}, of degree 8."""
    d3 = squarefree_part(d1 * d2)
    g = pari.sqrtint(d1 * d2 // d3)
    f = pari(f"x^2 - ({A}) - ({B})*s - ({C})*t - ({D})*s*t/{g}")
    f = pari.polresultant(f, pari(f"s^2 - ({d1})"), pari("s"))
    return pari.polresultant(f, pari(f"t^2 - ({d2})"), pari("t"))


def nested_biquadratic_fields(bound, basis_bound, progress=None):
    """The fields Q(\\sqrt{A + B\\sqrt{d_1} + C\\sqrt{d_2} + D\\sqrt{d_3}}) of degree 8 with
    |A|, |B|, |C|, |D| <= bound and the greedy basis d_1, d_2 of L with |d_i| <= basis_bound,
    each named by its simplest parameters (see above), as Fields.

    If given, progress(done, total, extra) is called after each basis.
    """
    bases = list(greedy_bases(2, basis_bound))
    reps = {}  # polredabs coeffs -> best (A, B, C, D, d1, d2)
    r = range(-bound, bound + 1)
    for i, (d1, d2) in enumerate(bases, start=1):
        for A, B, C, D in itertools.product(r, range(bound + 1), range(bound + 1), r):
            if not is_valid(A, B, C, D) or _has_square_factor(A, B, C, D):
                continue
            f = nested_biquadratic_poly(A, B, C, D, d1, d2)
            if not pari.polisirreducible(f):
                continue
            coeffs = polredabs_coeffs(f)
            p = (A, B, C, D, d1, d2)
            if coeffs not in reps or simplicity(*p) < simplicity(*reps[coeffs]):
                reps[coeffs] = p
        if progress:
            progress(i, len(bases), f", {len(reps)} fields so far")
    return [make_field(nested_biquadratic_name(*p), coeffs) for coeffs, p in reps.items()]


def _has_square_factor(*cs):
    """Whether some p^2 divides every coefficient (dividing by it gives the same field)."""
    g = 0
    for c in cs:
        g = int(pari.gcd(g, c))
    return g > 1 and any(e >= 2 for e in pari.factor(g)[1])
