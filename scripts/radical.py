"""Kummer fields Q(\\zeta_n, \\sqrt[n]{a}) and radical composita Q(\\sqrt[m]{a}, \\sqrt[n]{b}).

Kummer fields: for each n and 1 <= |a| <= bound with x^n - a irreducible, K is the compositum
of Q(\\zeta_n) and Q(\\sqrt[n]{a}) (unique up to isomorphism, as Q(\\zeta_n) is Galois). Fields
are identified by their polredabs polynomial, and a field reached several ways is named by the
smallest n, then the smallest |a|, preferring a > 0. As for cyclotomic fields we write
\\zeta_{n/2} when n = 2 mod 4 (Q(\\zeta_6) = Q(\\zeta_3)), and i for \\zeta_4, so
e.g. \\Q(\\zeta_{3}, \\sqrt[3]{2}), \\Q(i, \\sqrt[4]{2}), \\Q(\\zeta_{3}, \\sqrt[6]{2}).

Radical composita: for coprime m < n, K is the compositum of a pure field Q(\\sqrt[m]{a}) and a
pure field Q(\\sqrt[n]{b}) (for m = 2, a quadratic field Q(\\sqrt{a})), which has degree mn and
is unique up to isomorphism since the degrees are coprime. Each factor is named as in pure.py
(smallest |a|, a > 0 first) and \\sqrt{-1} is written i, e.g. \\Q(\\sqrt{2}, \\sqrt[3]{3}).
"""

from math import gcd

from nf_names import is_squarefree
from pari_fields import make_field, pari, polredabs_coeffs
from pure import pure_fields


def zeta_symbol(n):
    if n == 4:
        return "i"
    return r"\zeta_{%d}" % (n // 2 if n % 4 == 2 else n)


def root_symbol(m, a):
    if m == 2:
        return "i" if a == -1 else r"\sqrt{%d}" % a
    return r"\sqrt[%d]{%d}" % (m, a)


def kummer_name(n, a):
    return r"\Q(%s, %s)" % (zeta_symbol(n), root_symbol(n, a))


def compositum_name(m, a, n, b):
    return r"\Q(%s, %s)" % (root_symbol(m, a), root_symbol(n, b))


def compositum(f, g):
    """polredabs coefficients of the compositum of the fields defined by f and g, assuming
    it is unique up to isomorphism (one of them Galois, or coprime degrees)."""
    return polredabs_coeffs(pari.polredbest(pari.polcompositum(f, g)[0]))


def kummer_fields(ns, bound, progress=None):
    """The distinct fields Q(\\zeta_n, \\sqrt[n]{a}) for n in ns and 1 <= |a| <= bound, as Fields.

    Fields where \\zeta_n is already in Q(\\sqrt[n]{a}) are left to the pure and cyclotomic
    names (e.g. n = 4, a = -1 gives Q(\\sqrt[4]{-1}) = Q(\\zeta_8)).
    """
    jobs = [(n, a) for n in sorted(ns)
            for a in sorted(range(-bound, bound + 1), key=lambda a: (abs(a), -a)) if a != 0]
    fields = {}
    for i, (n, a) in enumerate(jobs, start=1):
        if progress:
            progress(i, len(jobs), f", {len(fields)} fields so far")
        f = pari(f"x^{n} - ({a})")
        if not pari.polisirreducible(f):
            continue
        coeffs = compositum(pari.polcyclo(n), f)
        if len(coeffs) - 1 == n or coeffs in fields:
            continue  # zeta_n is already in Q(\sqrt[n]{a}), or reached with a smaller n or |a|
        fields[coeffs] = (n, a)
    return [make_field(kummer_name(n, a), c) for c, (n, a) in fields.items()]


def radical_factors(m, bound):
    """[(a, polynomial)] for the distinct Q(\\sqrt[m]{a}), |a| <= bound, named as in pure.py."""
    if m == 2:
        return [(d, pari(f"x^2 - ({d})")) for d in sorted(range(-bound, bound + 1),
                                                         key=lambda d: (abs(d), -d))
                if d not in (0, 1) and is_squarefree(d)]
    out = []
    for K in pure_fields(m, bound):
        a = int(K.name.split("{")[-1].rstrip("})"))
        if a != -1:  # Q(\sqrt[m]{-1}) is cyclotomic, and named so in composita
            out.append((a, pari.Polrev(list(K.coeffs))))
    return out


def composita_fields(degree_pairs, bound, progress=None):
    """The fields Q(\\sqrt[m]{a}, \\sqrt[n]{b}) for (m, n) in degree_pairs (coprime, m < n)
    and |a|, |b| <= bound, as Fields."""
    jobs = []
    for m, n in degree_pairs:
        assert m < n and gcd(m, n) == 1
        bs = radical_factors(n, bound)
        jobs += [(m, a, f, n, b, g) for a, f in radical_factors(m, bound) for b, g in bs]
    fields = {}
    for i, (m, a, f, n, b, g) in enumerate(jobs, start=1):
        if progress:
            progress(i, len(jobs), f", {len(fields)} fields so far")
        coeffs = compositum(f, g)
        # Factors are listed simplest first, so the first name found for a field is kept.
        fields.setdefault(coeffs, compositum_name(m, a, n, b))
    return [make_field(name, c) for c, name in fields.items()]
