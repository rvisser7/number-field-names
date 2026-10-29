"""Recognise the families a given number field belongs to, and pick its best name.

This is the backward direction. The generate_*.py scripts run over small parameters and look
the fields they get up in the LMFDB; here we start from a field (e.g. an LMFDB field, from
its polynomial and ramified primes) and test each family. Every name is certified: it is
only returned once PARI has shown the field it describes is K. The best name is the one with
the smallest name_key (see nf_names.py), so it is the name the forward scripts would give,
without their bounds:

  * cyclotomic, real cyclotomic: K is compared with the Q(\\zeta_m), Q(\\zeta_m)^+ of its degree
    whose discriminant is that of K.
  * multiquadratic: K of degree 2^k is multiquadratic when it has 2^k - 1 quadratic subfields;
    it is named by the greedy basis, as in multiquadratic.py.
  * pure Q(\\sqrt[n]{a}): every prime dividing a ramifies in K, and for p not dividing n,
    v_p(disc K) = n - gcd(v_p(a), n). So a = +-prod p^e over the ramified p, 0 <= e < n, with
    gcd(e, n) fixed by the discriminant for p not dividing n: a finite list, which we sieve
    with primes q where K has a degree 1 prime (a must be an n-th power mod q) before testing
    the survivors with nfroots. The smallest |a| (a > 0 first) gives the name, as in pure.py.
  * Kummer Q(\\zeta_n, \\sqrt[n]{a}): the same, for the Galois fields containing \\zeta_n.
  * imprimitive quartic Q(\\sqrt{A + B\\sqrt{D}}): K = L(\\sqrt{beta}) with L = Q(\\sqrt{D}) its
    quadratic subfield. The A + B\\sqrt{D} that give K are the integral elements of
    beta L^{*2} (and of its conjugate). Each has ideal C^2 R, for a fixed squarefree R and C in
    a fixed ideal class, so running over the C of small norm in that class and over the units
    finds the simplest (A, B), as in imprimitive_quartic.py.
  * radical compositum, compositum: pairs of subfields of K with names (cyclotomic, quadratic,
    pure, imprimitive quartic) whose compositum is K, with the rules of radical.py and
    named_composita.py.

Not done yet: nested radicals Q(\\sqrt[n]{a + b\\sqrt{d}}) with n >= 3 (generate_nested_radical.py).
"""

import itertools
import math
from dataclasses import dataclass

from cyclotomic import conductors, cyclotomic_name, real_cyclotomic_name
from imprimitive_quartic import quartic_name
from multiquadratic import multiquadratic_name, order_key, span
from nf_names import name_key
from pari_fields import is_isomorphic, pari, polredabs_coeffs
from pure import pure_name
from radical import compositum_name, kummer_name, root_symbol, zeta_symbol

FAMILY_TESTS = ["cyclotomic", "real cyclotomic", "multiquadratic", "pure", "kummer",
                "imprimitive quartic", "composita"]
MAX_CANDIDATES = 10 ** 6  # radicands to sieve for one field before giving up


@dataclass
class NumberField:
    coeffs: tuple    # defining polynomial, constant term first
    r1: int
    disc_abs: int
    ramps: tuple     # ramified primes
    label: str = None

    @property
    def degree(self):
        return len(self.coeffs) - 1

    def poly(self, var="y"):
        # Number fields are in y, so that polynomials over them can be in x (which PARI
        # gives higher priority).
        return pari.Polrev(list(self.coeffs), var)


def from_label(label, coeffs, ramps=None):
    """The LMFDB field with this label and polynomial; ramps (its ramified primes) are
    computed if not given."""
    n, r1, disc_abs, _ = (int(x) for x in label.split("."))
    assert len(coeffs) == n + 1, label
    if ramps is None:
        ramps = [int(p) for p in pari.factor(disc_abs)[0]] if disc_abs > 1 else []
    return NumberField(tuple(int(c) for c in coeffs), r1, disc_abs, tuple(sorted(ramps)), label)


def from_coeffs(coeffs, primes=None):
    """The field of the polynomial coeffs. If given, primes contains its ramified primes."""
    f = pari.Polrev(list(coeffs))
    disc_abs = abs(int(pari.nfdisc(f)))
    if primes is None:
        ramps = [int(p) for p in pari.factor(disc_abs)[0]] if disc_abs > 1 else []
    else:
        ramps = [p for p in primes if disc_abs % p == 0]
    return NumberField(tuple(int(c) for c in coeffs), int(pari.polsturm(f)), disc_abs,
                       tuple(sorted(ramps)))


def valuation(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


# ---------------------------------------------------------------- small helpers

def _unramified_primes(K, count):
    """Up to count primes q, not dividing the discriminant of K's polynomial."""
    pdisc = abs(int(pari.poldisc(K.poly())))
    out = []
    for q in pari.primes(400):
        q = int(q)
        if pdisc % q:
            out.append(q)
            if len(out) == count:
                break
    return out


def _factor_degrees(K, q):
    fac = pari.factormod(K.poly(), q)
    return [int(pari.poldegree(g)) for g in fac[0]]


def _degree_one_primes(K, count=12):
    """Primes q not dividing disc(f) with f = K's polynomial having a root mod q, so K has a
    prime of degree 1 above q, whose residue field is F_q."""
    out = []
    for q in _unramified_primes(K, 200):
        if len(pari.polrootsmod(K.poly(), q)):
            out.append(q)
            if len(out) == count:
                break
    return out


def _has_root(K, g):
    """Whether the polynomial g (in x, rational coefficients) has a root in K."""
    return len(pari.nfroots(K.poly(), g)) > 0


def _sieve(K, n, options, signs):
    """The radicands a = s * prod p^e (s in signs, (p, es) in options, e in es) that are
    n-th powers modulo every q in _degree_one_primes(K), sorted by (|a|, a < 0)."""
    if math.prod(len(es) for _, es in options) * len(signs) > MAX_CANDIDATES:
        return None
    # a is an n-th power mod q iff a^((q-1)/g) = 1, g = gcd(n, q-1). Writing each p^((q-1)/g)
    # as h^L_p for a generator h of the g-th roots of unity, that is a linear condition on the
    # exponents e_p mod g.
    conds = []
    for q in _degree_one_primes(K):
        g = math.gcd(n, q - 1)
        if g == 1:
            continue
        h = pow(int(pari.znprimroot(q)), (q - 1) // g, q)
        dlog = {pow(h, j, q): j for j in range(g)}
        conds.append((g, [dlog[pow(p, (q - 1) // g, q)] for p, _ in options],
                      dlog[pow(q - 1, (q - 1) // g, q)]))
    out = []
    for s in signs:
        for es in itertools.product(*(es for _, es in options)):
            if all((sum(e * l for e, l in zip(es, ls)) + (lm if s < 0 else 0)) % g == 0
                   for g, ls, lm in conds):
                a = s * math.prod(p ** e for (p, _), e in zip(options, es))
                if a != 1:
                    out.append(a)
    return sorted(out, key=lambda a: (abs(a), a < 0))


# ---------------------------------------------------------------- the families

_CYCLO = {}


def _cyclotomic_list(n, real):
    """[(m, polredabs coeffs, |disc|)] for Q(\\zeta_m) (or Q(\\zeta_m)^+) of degree n."""
    key = (n, real)
    if key not in _CYCLO:
        out = []
        for m in conductors(n, real):
            if int(pari.eulerphi(m)) // (2 if real else 1) != n:
                continue
            f = pari(f"minpoly(Mod(x + 1/x, polcyclo({m})))") if real else pari.polcyclo(m)
            out.append((m, polredabs_coeffs(f), abs(int(pari.nfdisc(f)))))
        _CYCLO[key] = out
    return _CYCLO[key]


def cyclotomic_conductor(K, real=False):
    """The m with K = Q(\\zeta_m) (or Q(\\zeta_m)^+), m not 2 mod 4, or None. Degree >= 3."""
    n = K.degree
    if n < 3 or K.r1 != (n if real else 0):
        return None
    for m, coeffs, disc_abs in _cyclotomic_list(n, real):
        if disc_abs == K.disc_abs and (coeffs == K.coeffs or is_isomorphic(coeffs, K.coeffs)):
            return m
    return None


def multiquadratic_basis(K):
    """The greedy basis (see multiquadratic.py) if K is multiquadratic, else None."""
    n = K.degree
    k = n.bit_length() - 1
    if n < 4 or n != 2 ** k or K.r1 not in (0, n):
        return None
    # every unramified prime has residue degree 1 or 2
    if any(max(_factor_degrees(K, q)) > 2 for q in _unramified_primes(K, 8)):
        return None
    subs = pari.nfsubfields(K.poly(), 2)
    if len(subs) != n - 1:
        return None
    ds = sorted({int(pari.core(pari.poldisc(s[0]))) for s in subs}, key=order_key)
    basis = []
    for d in ds:
        if d not in span(basis):
            basis.append(d)
    return basis if len(basis) == k else None


def pure_radicand(K):
    """The a with K = Q(\\sqrt[n]{a}), smallest |a| and a > 0 first, or None. Degree >= 3."""
    n = K.degree
    if n < 3 or K.r1 != (1 if n % 2 else K.r1) or (n % 2 == 0 and K.r1 not in (0, 2)):
        return None
    options = []
    for p in K.ramps:
        if n % p == 0:
            options.append((p, range(n)))
            continue
        g = n - valuation(K.disc_abs, p)
        if g < 1 or n % g:
            return None
        options.append((p, [e for e in range(1, n) if math.gcd(e, n) == g]))
    signs = [1] if n % 2 else [1 if K.r1 == 2 else -1]
    for a in _sieve(K, n, options, signs) or []:
        g = pari(f"x^{n} - ({a})")
        if pari.polisirreducible(g) and _has_root(K, g):
            return a
    return None


def _is_galois(K):
    """Cheap necessary test: f mod q has factors of one degree, for some unramified q."""
    return all(len(set(_factor_degrees(K, q))) == 1 for q in _unramified_primes(K, 8))


def kummer_params(K):
    """(n, a) with K = Q(\\zeta_n, \\sqrt[n]{a}) of degree phi(n) n, smallest |a| and a > 0
    first, or None."""
    N = K.degree
    ns = [n for n in range(3, N + 1) if n * int(pari.eulerphi(n)) == N]
    if not ns or K.r1 != 0 or not _is_galois(K):
        return None
    n = ns[0]  # n phi(n) is injective
    cyc = pari.polcyclo(n)
    if not _has_root(K, cyc):
        return None
    options = [(p, range(n)) for p in K.ramps]
    for a in _sieve(K, n, options, [1, -1]) or []:
        g = pari(f"x^{n} - ({a})")
        if len(pari.nffactor(pari.polcyclo(n, "y"), g)[0]) == 1 and _has_root(K, g):
            return n, a
    return None


def _to_L(z, s, f):
    """(u, v) rational with z = u + v s, for z, s in K = Q[y]/f, when z is in Q(s)."""
    n = int(pari.poldegree(f))
    zv = pari.Vecrev(pari.lift(z), n)
    sv = pari.Vecrev(pari.lift(s), n)
    # z - v s is constant: pick the first index > 0 where s has a nonzero coordinate
    j = next(i for i in range(1, n) if sv[i] != 0)
    v = zv[j] / sv[j]
    return zv[0] - v * sv[0], v


def quartic_params(K):
    """(A, B, D) with K = Q(\\sqrt{A + B\\sqrt{D}}) simplest as in imprimitive_quartic.py and
    A != 0 (A = 0 is the pure quartic), or None if K is not a C4 or D4 quartic."""
    if K.degree != 4:
        return None
    fy = K.poly()
    subs = pari.nfsubfields(fy, 2)
    if len(subs) != 1:
        return None  # primitive, or biquadratic
    g, h = subs[0]
    disc = pari.poldisc(g)
    D = int(pari.core(disc))
    r = pari.sqrtint(disc / D)
    b = pari.polcoef(g, 1)
    s = pari.Mod((2 * h + b) / r, fy)  # sqrt(D) in K
    theta = pari.Mod(pari("y"), fy)
    for tau in pari.nfgaloisconj(fy):
        if tau == pari("y"):
            continue
        if pari.Mod(pari.subst(pari.lift(s), "y", tau), fy) == s:
            break
    else:
        return None
    gamma = theta - pari.Mod(tau, fy)  # tau(gamma) = -gamma, so gamma^2 is in L
    u, v = _to_L(gamma ** 2, s, fy)
    return _simplest_abd(u, v, D)


def _simplest_abd(u, v, D):
    """The simplest (A, B) (see imprimitive_quartic.simplicity) with A + B sqrt(D) in
    (u + v sqrt(D)) L^{*2} or its conjugate, L = Q(sqrt(D)), A != 0; returns (A, B, D)."""
    bnf = pari.bnfinit(pari("w^2") - D, 1)
    w = pari.Mod(pari("w"), pari("w^2") - D)
    beta = u + v * w
    fac = pari.idealfactor(bnf, beta)
    I, R = 1, 1
    for P, e in zip(fac[0], fac[1]):
        I = pari.idealmul(bnf, I, pari.idealpow(bnf, P, int(e) // 2))
        if int(e) % 2:
            R = pari.idealmul(bnf, R, P)
    norm_R = int(pari.idealnorm(bnf, R))
    if D > 0:
        eps = pari.Mod(pari.lift(pari("(b) -> b.fu[1]")(bnf)), pari("w^2") - D)
        unit_sq = None
    else:
        tu = pari("(b) -> b.tu")(bnf)
        unit_sq = [pari.Mod(pari.lift(tu[1]), pari("w^2") - D) ** (2 * j) for j in range(int(tu[0]))]
    sqrtD = math.sqrt(abs(D))

    def ab(z):
        """Integers (A, B) with z = (A + B sqrt D) / c^2 for the smallest c, B >= 0."""
        a0, b0 = pari.polcoef(pari.lift(z), 0), pari.polcoef(pari.lift(z), 1)
        den = int(pari.lcm(pari.denominator(a0), pari.denominator(b0)))
        A, B = int(a0 * den * den), int(b0 * den * den)
        g = math.gcd(A, B)
        for p in [int(p) for p in pari.factor(g)[0]] if g > 1 else []:
            while A % (p * p) == 0 and B % (p * p) == 0:
                A //= p * p
                B //= p * p
        return (A, abs(B)) if B else None

    def key(A, B):
        return (max(abs(A), abs(B)), abs(B), abs(A), A < 0)

    def candidates(max_norm):
        for C in itertools.chain.from_iterable(pari.ideallist(bnf, max_norm)):
            cls, nu = pari.bnfisprincipal(bnf, pari.idealdiv(bnf, C, I), 1)
            if any(int(c) for c in cls):
                continue
            z0 = beta * pari.nfbasistoalg(bnf, nu) ** 2
            if unit_sq is not None:
                zs = [z0 * t for t in unit_sq]
            else:
                zs = [z0 * eps ** (2 * k) for k in _unit_range(z0, eps, bound[0] * (1 + sqrtD))]
            for z in zs:
                res = ab(z)
                if res and res[0] != 0:
                    yield res

    bound = [math.inf]
    best = None
    # First a few small ideals for a good starting bound, then everything below it.
    for max_norm in (16, None):
        if max_norm is None:
            if best is None:
                return None
            M = key(*best)[0]
            max_norm = int(math.isqrt(int(M * M * (1 + abs(D)) // norm_R))) + 1
            bound[0] = M
        for A, B in candidates(max_norm):
            if best is None or key(A, B) < key(*best):
                best = (A, B)
    return best + (D,)


def _unit_range(z, eps, T):
    """The k with |sigma(z eps^{2k})| <= T at both real places (or all k near the balance
    point when T is infinite)."""
    x1, x2 = [abs(float(c)) for c in pari.conjvec(z)]
    e1 = abs(float(pari.conjvec(eps)[0]))
    le = 2 * math.log(e1)
    mid = (math.log(x2) - math.log(x1)) / (2 * le)  # balances the two conjugates
    if math.isinf(T):
        return range(math.floor(mid) - 1, math.ceil(mid) + 2)
    lo = (math.log(x2) - math.log(T)) / le
    hi = (math.log(T) - math.log(x1)) / le
    if le < 0:
        lo, hi = -hi, -lo
    return range(math.floor(lo) - 1, math.ceil(hi) + 2)


# ---------------------------------------------------------------- composita from subfields

_BLOCKS = {}


def _block(coeffs, primes):
    """How a subfield (polredabs coeffs) can appear in a \\Q(g_1, g_2) name: (kind, degree, data)
    with kind "quadratic" (data d), "cyclotomic" (m), "pure" (a) or "quartic" ((A, B, D)), or None."""
    if coeffs in _BLOCKS:
        return _BLOCKS[coeffs]
    d = len(coeffs) - 1
    out = None
    if d == 2:
        out = ("quadratic", 2, int(pari.core(pari.poldisc(pari.Polrev(list(coeffs))))))
    else:
        S = from_coeffs(coeffs, primes)
        m = cyclotomic_conductor(S)
        if m:
            out = ("cyclotomic", d, m)
        else:
            a = pure_radicand(S)
            if a:
                out = ("pure", d, a)
            elif d == 4:
                abd = quartic_params(S)
                if abd:
                    out = ("quartic", 4, abd)
    _BLOCKS[coeffs] = out
    return out


def _symbol(block):
    kind, d, data = block
    if kind == "quadratic":
        return zeta_symbol({-1: 4, -3: 3}[data]) if data in (-1, -3) else root_symbol(2, data)
    if kind == "cyclotomic":
        return zeta_symbol(data)
    if kind == "pure":
        return root_symbol(d, data)
    return quartic_name(*data)[3:-1]


def _poly(block):
    kind, d, data = block
    if kind == "quadratic":
        return pari(f"x^2 - ({data})")
    if kind == "cyclotomic":
        return pari.polcyclo(data)
    if kind == "pure":
        return pari(f"x^{d} - ({data})")
    A, B, D = data
    return pari(f"(x^2 - ({A}))^2 - ({B})^2 * ({D})")


def _named_kind(block):
    """The kind of a block in named_composita.py, where i and zeta_3 are cyclotomic."""
    if block[0] == "quadratic" and block[2] in (-1, -3):
        return "cyclotomic"
    return block[0]


KIND_ORDER = ["cyclotomic", "quadratic", "pure", "quartic"]


def composita_names(K):
    """The radical compositum and \\Q(g_1, g_2) names (see radical.py, named_composita.py) of K,
    from pairs of its subfields."""
    n = K.degree
    divisors = [d for d in range(2, n) if n % d == 0]
    if not divisors:
        return []
    subs = {d: [] for d in divisors}
    for d in divisors:
        for g, _ in pari.nfsubfields(K.poly(), d):
            b = _block(polredabs_coeffs(g), K.ramps)
            if b:
                subs[d].append(b)
    names = []
    for d1 in divisors:
        d2 = n // d1
        if d1 > d2:
            break
        for i, b1 in enumerate(subs[d1]):
            for b2 in subs[d2][i + 1:] if d1 == d2 else subs[d2]:
                names += _pair_names(b1, b2, n)
    return names


def _pair_names(b1, b2, n):
    out = []
    # radical compositum: Q(\sqrt[m]{a}, \sqrt[n]{b}), m < n coprime, quadratic or pure
    if (b1[0] in ("quadratic", "pure") and b2[0] == "pure" and b1[1] < b2[1]
            and math.gcd(b1[1], b2[1]) == 1):
        out.append(compositum_name(b1[1], b1[2], b2[1], b2[2]))
    # named compositum: the first generator cyclotomic or quadratic, kinds in order
    for c1, c2 in [(b1, b2), (b2, b1)]:
        k1, k2 = _named_kind(c1), _named_kind(c2)
        if k1 not in ("cyclotomic", "quadratic") or KIND_ORDER.index(k1) >= KIND_ORDER.index(k2):
            continue
        if c1[1] == c2[1] == 2 or (c1[1] == 2 and k2 == "pure" and c2[1] % 2):
            continue
        comp = pari.polcompositum(_poly(c1), _poly(c2))
        if int(pari.poldegree(comp[0])) == n:  # linearly disjoint, as c1 is Galois
            out.append(r"\Q(%s, %s)" % (_symbol(c1), _symbol(c2)))
    return out


# ---------------------------------------------------------------- all together

def all_names(K, families=FAMILY_TESTS):
    """Every name of K found by the tests in families (see FAMILY_TESTS)."""
    names = []
    if "cyclotomic" in families:
        m = cyclotomic_conductor(K)
        if m:
            names.append(cyclotomic_name(m))
    if "real cyclotomic" in families:
        m = cyclotomic_conductor(K, real=True)
        if m:
            names.append(real_cyclotomic_name(m))
    if "multiquadratic" in families:
        basis = multiquadratic_basis(K)
        if basis:
            names.append(multiquadratic_name(basis))
    if "pure" in families:
        a = pure_radicand(K)
        if a:
            names.append(pure_name(K.degree, a))
    if "kummer" in families:
        na = kummer_params(K)
        if na:
            names.append(kummer_name(*na))
    if "imprimitive quartic" in families:
        abd = quartic_params(K)
        if abd:
            names.append(quartic_name(*abd))
    if "composita" in families:
        names += composita_names(K)
    return list(dict.fromkeys(names))


def best_name(K, families=FAMILY_TESTS):
    """The best name of K (smallest name_key) among all_names, or None."""
    names = all_names(K, families)
    return min(names, key=name_key) if names else None
