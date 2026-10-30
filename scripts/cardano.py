"""Cubic fields by Cardano's formula, Q(\\sqrt[3]{u + e\\sqrt{d}} + \\sqrt[3]{u - e\\sqrt{d}}), computed
with PARI (via cypari2).

Every cubic field K = Q(alpha) has a generator of trace 0 whose minimal polynomial is
x^3 - 3kx - 2u with integers k, u (scale any generator of trace 0 by 6), and then, by
Cardano's formula,

    alpha = \\sqrt[3]{u + \\sqrt{u^2 - k^3}} + \\sqrt[3]{u - \\sqrt{u^2 - k^3}},

with the cube roots chosen so that their product is k. We write u^2 - k^3 = e^2 d with d
squarefree, which K fixes: the quadratic resolvent Q(\\sqrt{-3 disc K}) is Q(\\sqrt{-3d}). So d = 1
for the pure cubics (which keep their pure name, see pure.py) and d = -3 for the cyclic
ones. For totally real K, d < 0 (the casus irreducibilis), e.g. the cube roots of 1 +- 2\\sqrt{-2}.
\\sqrt{-1} is written i.

Conventions:

  * u > 0 (alpha -> -alpha changes the sign of u), e > 0, and no p^3 divides both u and e
    (dividing by it scales alpha by 1/p);
  * of the (u, e) giving the same field we keep the one with the smallest max(u, e), then
    e, then u (see nf_names.cardano_key).

The numbers grow with the discriminant (disc(x^3 - 3kx - 2u) = -108 e^2 d), so
name_lmfdb_fields.py only keeps names up to a size (MAX_CARDANO_SIZE in recognise.py).
This family comes last in FAMILIES, so every other name wins (e.g. \\Q(\\zeta_{7})^+).
"""

from nf_names import cardano_key
from pari_fields import make_field, pari, polredabs_coeffs


def cardano_name(u, e, d):
    r = ("i" if d == -1 else r"\sqrt{%d}" % d)
    r = ("" if e == 1 else str(e)) + r
    return r"\Q(\sqrt[3]{%d + %s} + \sqrt[3]{%d - %s})" % (u, r, u, r)


def simplicity(u, e, d):
    return cardano_key(u, e, d)


def cardano_size(u, e, d):
    """The largest number in the name."""
    return max(u, e, abs(d))


def cardano_poly(u, k):
    return pari(f"x^3 - 3*({k})*x - 2*({u})")


def cardano_k(u, e, d):
    """The k with k^3 = u^2 - e^2 d."""
    v = u * u - e * e * d
    k = (1 if v > 0 else -1) * int(pari.sqrtnint(abs(v), 3))
    assert k ** 3 == v, (u, e, d)
    return k


def cardano_params(u, k):
    """(u, e, d) for the polynomial x^3 - 3kx - 2u, after removing cube factors, or None if
    it does not define a cubic field with d != 1."""
    v = u * u - k ** 3
    if u == 0 or v == 0:
        return None
    d = int(pari.core(v))
    if d == 1:
        return None  # pure (or reducible)
    e = int(pari.sqrtint(v // d))
    if u < 0:
        u = -u
    g = int(pari.gcd(u, e))
    for p in [int(p) for p in pari.factor(g)[0]] if g > 1 else []:
        while u % p ** 3 == 0 and e % p ** 3 == 0:
            u, e = u // p ** 3, e // p ** 3
    return u, e, d


def cardano_fields(bound, k_bound=None, progress=None):
    """The cubic fields x^3 - 3kx - 2u with 0 < u <= bound and 0 < |k| <= k_bound (default
    bound), each named by its simplest (u, e, d), as Fields.

    If given, progress(done, total, extra) is called after each u.
    """
    k_bound = bound if k_bound is None else k_bound
    reps = {}
    for u in range(1, bound + 1):
        for k in range(-k_bound, k_bound + 1):
            p = cardano_params(u, k) if k else None
            if not p:
                continue
            f = cardano_poly(u, k)
            if not pari.polisirreducible(f):
                continue
            key = polredabs_coeffs(f)
            if key not in reps or simplicity(*p) < simplicity(*reps[key]):
                reps[key] = p
        if progress:
            progress(u, bound, f", {len(reps)} fields so far")
    return [make_field(cardano_name(*p), key) for key, p in reps.items()]
