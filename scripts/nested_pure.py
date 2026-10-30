"""Quadratic extensions of pure fields, Q(\\sqrt{A_0 + A_1\\sqrt[n]{m_1} + ... + A_{n-1}\\sqrt[n]{m_{n-1}}}),
computed with PARI (via cypari2).

For an odd prime n and a pure field F = Q(\\sqrt[n]{m}), these are the fields K = F(\\sqrt{beta}),
beta in F, of degree 2n: every field of degree 2n with a pure subfield of degree n is one.
E.g. for n = 3 the sextics \\Q(\\sqrt{1 + \\sqrt[3]{2}}) and \\Q(\\sqrt{1 + \\sqrt[3]{2} + \\sqrt[3]{4}}).

Conventions:

  * m is the radicand of F's pure name (smallest |m|, m > 0; see pure.py), and the basis of F
    is 1, \\sqrt[n]{m_1}, ..., \\sqrt[n]{m_{n-1}}, where m_k is the n-th power free part of m^k (so
    \\sqrt[n]{m^k} = c_k \\sqrt[n]{m_k}, see nf_names.pure_basis), with real roots: for m = 2 that is
    1, \\sqrt[3]{2}, \\sqrt[3]{4}, and for m = 12 it is 1, \\sqrt[3]{12}, \\sqrt[3]{18}.
  * beta is not rational (then Q(\\sqrt{beta}) = K, as n is prime), and has at least two terms:
    \\sqrt{c\\sqrt[n]{m_k}} is pure, Q(\\sqrt[2n]{c^n m_k}), which generate_pure.py names.
  * Of the (A_0, ..., A_{n-1}) giving the same field (beta times squares of F) we keep the
    simplest: smallest max |A_k|, then |A_{n-1}|, ..., |A_1|, |A_0|, then A_1 >= 0, ...,
    A_{n-1} >= 0, A_0 >= 0 (see nf_names.nested_pure_key).

Many of these fields also have a name of another family (pure Q(\\sqrt[2n]{a}), Kummer,
radical composita such as Q(\\sqrt{2}, \\sqrt[3]{3}), nested radicals, ...). This family comes
last in FAMILIES, so those names win, and only fields with no other name are named this way.
"""

import itertools

from nf_names import nested_pure_key, pure_basis
from pari_fields import make_field, pari, polredabs_coeffs
from pure import pure_fields


def nested_pure_name(n, m, coeffs):
    s = ""
    for (_, mk), c in zip(pure_basis(n, m), coeffs):
        if c == 0:
            continue
        body = str(abs(c)) if mk == 1 else ("" if abs(c) == 1 else str(abs(c))) + r"\sqrt[%d]{%d}" % (n, mk)
        s += ("-" if c < 0 else "") + body if not s else (" - " if c < 0 else " + ") + body
    return r"\Q(\sqrt{%s})" % s


def simplicity(n, m, coeffs):
    return nested_pure_key(n, m, coeffs)


def is_valid(coeffs):
    """beta is not rational, and not a single term (see above)."""
    return any(coeffs[1:]) and sum(1 for c in coeffs if c) >= 2


def nested_pure_poly(n, m, coeffs):
    """The minimal polynomial (if irreducible) of \\sqrt{beta}, of degree 2n."""
    beta = sum(pari(c) * pari(f"y^{k}") / ck for k, ((ck, _), c) in enumerate(zip(pure_basis(n, m), coeffs)))
    f = pari.polresultant(pari("x^2") - beta, pari(f"y^{n} - {m}"), pari("y"))
    return f / pari.pollead(f)


def pure_radicands(n, bound):
    """The m > 0 of the pure names \\Q(\\sqrt[n]{m}) with m <= bound."""
    return sorted(int(F.name.split("{")[-1].rstrip("})")) for F in pure_fields(n, bound)
                  if not F.name.endswith("{-1})"))


def nested_pure_fields(n, bound, m_bound, progress=None):
    """The fields Q(\\sqrt{A_0 + A_1\\sqrt[n]{m_1} + ...}) of degree 2n with every |A_k| <= bound
    and pure name radicand m <= m_bound, each named by its simplest parameters, as Fields.

    If given, progress(done, total, extra) is called after each m.
    """
    ms = [m for m in pure_radicands(n, m_bound) if m > 0]
    reps = {}  # polredabs coeffs -> best (n, m, coeffs)
    r = range(-bound, bound + 1)
    for i, m in enumerate(ms, start=1):
        for coeffs in itertools.product(r, repeat=n):
            if not is_valid(coeffs) or _has_square_factor(coeffs):
                continue
            f = nested_pure_poly(n, m, coeffs)
            if not pari.polisirreducible(f):
                continue
            key = polredabs_coeffs(f)
            p = (n, m, coeffs)
            if key not in reps or simplicity(*p) < simplicity(*reps[key]):
                reps[key] = p
        if progress:
            progress(i, len(ms), f", {len(reps)} fields so far")
    return [make_field(nested_pure_name(*p), key) for key, p in reps.items()]


def _has_square_factor(cs):
    """Whether some p^2 divides every coefficient (dividing by it gives the same field)."""
    g = 0
    for c in cs:
        g = int(pari.gcd(g, c))
    return g > 1 and any(e >= 2 for e in pari.factor(g)[1])
