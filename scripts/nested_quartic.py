"""Octic fields Q(\\sqrt{a + b\\sqrt{D} + c\\sqrt{A + B\\sqrt{D}} + d\\sqrt{D}\\sqrt{A + B\\sqrt{D}}}), computed
with PARI (via cypari2).

These are the quadratic extensions K = F(\\sqrt{beta}) of an imprimitive quartic field
F = Q(\\sqrt{A + B\\sqrt{D}}) (see imprimitive_quartic.py): every octic field with a C4 or D4
quartic subfield is one. E.g. \\Q(\\sqrt{1 + \\sqrt{1 + \\sqrt{2}}}).

Conventions:

  * (A, B, D) is F's own imprimitive quartic name, and beta is written in the basis 1, \\sqrt{D},
    gamma, \\sqrt{D} gamma of F, with gamma = \\sqrt{A + B\\sqrt{D}}; \\sqrt{-1} is written i.
  * beta does not lie in Q(\\sqrt{D}) (so c or d is nonzero, and Q(\\sqrt{beta}) = K), and has at
    least two terms (\\sqrt{c gamma} gives a nested radical Q(\\sqrt[4]{c^2 (A + B\\sqrt{D})})).
  * Of the (a, b, c, d) giving the same field (beta times squares of F, and the conjugates of
    beta under the automorphisms of F, e.g. gamma -> -gamma) we keep the simplest: smallest
    max(|a|, |b|, |c|, |d|), then |d|, |c|, |b|, |a|, then the simplest quartic, then c >= 0,
    d >= 0, b >= 0, a >= 0 (see nf_names.nested_quartic_key).

Many of these fields also have a name of another family (e.g. composita such as
Q(i, \\sqrt{1 + \\sqrt{2}}), or Q(\\sqrt{A + B\\sqrt{d_1} + C\\sqrt{d_2} + D\\sqrt{d_3}}) when K also has a
biquadratic subfield). This family comes last in FAMILIES, so those names win, and only
octics with no other name are named this way.
"""

import itertools

from imprimitive_quartic import imprimitive_quartic_fields, quartic_name
from nf_names import nested_quartic_key, parse_quartic
from pari_fields import make_field, pari, polredabs_coeffs


def nested_quartic_name(a, b, c, d, A, B, D):
    s_sym = "i" if D == -1 else r"\sqrt{%d}" % D
    g_sym = quartic_name(A, B, D)[3:-1]
    out = ""
    for k, syms in [(a, ""), (b, s_sym), (c, g_sym), (d, s_sym + g_sym)]:
        if k == 0:
            continue
        body = str(abs(k)) if not syms else ("" if abs(k) == 1 else str(abs(k))) + syms
        out += ("-" if k < 0 else "") + body if not out else (" - " if k < 0 else " + ") + body
    return r"\Q(\sqrt{%s})" % out


def simplicity(a, b, c, d, A, B, D):
    return nested_quartic_key(a, b, c, d, A, B, D)


def is_valid(a, b, c, d):
    return (c or d) and sum(1 for k in (a, b, c, d) if k) >= 2


def nested_quartic_poly(a, b, c, d, A, B, D):
    """The minimal polynomial (if irreducible) of \\sqrt{beta}, of degree 8."""
    f = pari(f"x^2 - ({a}) - ({b})*s - ({c})*g - ({d})*s*g")
    f = pari.polresultant(f, pari(f"g^2 - ({A}) - ({B})*s"), pari("g"))
    f = pari.polresultant(f, pari(f"s^2 - ({D})"), pari("s"))
    return f / pari.pollead(f)


def nested_quartic_fields(bound, quartic_bound, progress=None):
    """The fields Q(\\sqrt{a + b\\sqrt{D} + c\\sqrt{A + B\\sqrt{D}} + d\\sqrt{D}\\sqrt{A + B\\sqrt{D}}}) of
    degree 8 with |a|, |b|, |c|, |d| <= bound, over the imprimitive quartic fields of
    imprimitive_quartic_fields(quartic_bound), each named by its simplest parameters, as Fields.

    If given, progress(done, total, extra) is called after each quartic field.
    """
    quartics = [parse_quartic(F.name) for F in imprimitive_quartic_fields(quartic_bound)]
    reps = {}  # polredabs coeffs -> best (a, b, c, d, A, B, D)
    r = range(-bound, bound + 1)
    for i, ABD in enumerate(quartics, start=1):
        for a, b, c, d in itertools.product(r, r, range(bound + 1), r):
            if (c == 0 and d <= 0) or not is_valid(a, b, c, d) or _has_square_factor(a, b, c, d):
                continue  # gamma -> -gamma flips the signs of c and d
            f = nested_quartic_poly(a, b, c, d, *ABD)
            if not pari.polisirreducible(f):
                continue
            key = polredabs_coeffs(f)
            p = (a, b, c, d) + ABD
            if key not in reps or simplicity(*p) < simplicity(*reps[key]):
                reps[key] = p
        if progress:
            progress(i, len(quartics), f", {len(reps)} fields so far")
    return [make_field(nested_quartic_name(*p), key) for key, p in reps.items()]


def _has_square_factor(*cs):
    """Whether some p^2 divides every coefficient (dividing by it gives the same field)."""
    g = 0
    for c in cs:
        g = int(pari.gcd(g, c))
    return g > 1 and any(e >= 2 for e in pari.factor(g)[1])
