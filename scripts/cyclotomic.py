"""Cyclotomic fields Q(\\zeta_n) and their maximal real subfields Q(\\zeta_n)^+.

Q(\\zeta_n) has degree phi(n), and Q(\\zeta_n)^+ = Q(\\zeta_n + \\zeta_n^{-1}) has degree
phi(n)/2. For odd n, Q(\\zeta_{2n}) = Q(\\zeta_n), so we only use n not 2 mod 4, as the
LMFDB does. Fields of degree 1 or 2 are left out: they are named as Q or as quadratic
fields (e.g. Q(\\zeta_3) is \\Q(\\sqrt{-3}) and Q(\\zeta_5)^+ is \\Q(\\sqrt{5})).
"""

from pari_fields import make_field, pari, polredabs_coeffs


def cyclotomic_name(n):
    return r"\Q(\zeta_{%d})" % n


def real_cyclotomic_name(n):
    return r"\Q(\zeta_{%d})^+" % n


def conductors(max_degree, real):
    """The n (not 2 mod 4) for which Q(\\zeta_n), or Q(\\zeta_n)^+ if real, has
    degree between 3 and max_degree."""
    # phi(n) >= sqrt(n/2), so phi(n) <= D implies n <= 2 D^2.
    max_phi = 2 * max_degree if real else max_degree
    out = []
    for n in range(3, 2 * max_phi ** 2 + 1):
        if n % 4 == 2:
            continue
        d = int(pari.eulerphi(n)) // (2 if real else 1)
        if 3 <= d <= max_degree:
            out.append(n)
    return out


def cyclotomic_fields(max_degree, progress=None):
    ns = conductors(max_degree, real=False)
    out = []
    for i, n in enumerate(ns, start=1):
        if progress:
            progress(i, len(ns))
        out.append(make_field(cyclotomic_name(n), polredabs_coeffs(pari.polcyclo(n))))
    return out


def real_cyclotomic_fields(max_degree, progress=None):
    ns = conductors(max_degree, real=True)
    out = []
    for i, n in enumerate(ns, start=1):
        if progress:
            progress(i, len(ns))
        f = pari(f"minpoly(Mod(x + 1/x, polcyclo({n})))")  # minimal polynomial of zeta_n + zeta_n^{-1}
        out.append(make_field(real_cyclotomic_name(n), polredabs_coeffs(f)))
    return out
