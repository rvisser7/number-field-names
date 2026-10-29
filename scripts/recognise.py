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
    finds the simplest (A, B), as in imprimitive_quartic.py. Pure fields are left out (they
    keep their pure name).
  * nested radical Q(\\sqrt[n]{a + b\\sqrt{d}}), n >= 3: for each quadratic subfield L, the beta with
    K = L(\\sqrt[n]{beta}) lie in the S-units of L modulo n-th powers (S: the primes above the
    ramified primes and n, and enough primes for the class group). Being an n-th power at the
    degree 1 primes of K is a linear condition mod n on the exponents; the survivors are checked
    with nffactor and nfroots, and each beta simplified as for the quartics. Pure fields are left
    out, and fields with more than MAX_NESTED_CANDIDATES candidates are skipped.
  * radical compositum, compositum: pairs of subfields of K with names (cyclotomic, quadratic,
    pure, imprimitive quartic) whose compositum is K, with the rules of radical.py and
    named_composita.py. Also names with three or more generators, from a Galois named subfield
    (quadratic, cyclotomic, multiquadratic, Kummer) and another named subfield, e.g.
    Q(i, \\sqrt[3]{3}, \\sqrt[4]{2}); fewer generators win.
"""

import itertools
import math
import re
from dataclasses import dataclass, field

from cyclotomic import conductors, cyclotomic_name, real_cyclotomic_name
from imprimitive_quartic import quartic_name
from multiquadratic import multiquadratic_name, order_key, span, sqrt_symbol
from nested_radical import nested_name
from nf_names import generator_key, name_key, radicand_size
from pari_fields import is_isomorphic, pari, polredabs_coeffs
from pure import pure_name
from radical import compositum_name, kummer_name, root_symbol, zeta_symbol

FAMILY_TESTS = ["cyclotomic", "real cyclotomic", "multiquadratic", "pure", "kummer",
                "imprimitive quartic", "nested radical", "composita"]
MAX_CANDIDATES = 10 ** 6  # radicands to sieve for one field before giving up
# nested radicals: candidates for beta to test for one quadratic subfield before giving up. Only
# fields with many radicals get there (e.g. Q(\\zeta_5, \\sqrt{2}, \\sqrt{3})), and those have
# composita names.
MAX_NESTED_CANDIDATES = 20000  # e.g. 6^5 (about 1s), not 8^5


@dataclass
class NumberField:
    coeffs: tuple    # defining polynomial, constant term first
    r1: int
    disc_abs: int
    ramps: tuple     # ramified primes
    label: str = None
    _subfields: dict = field(default_factory=dict, repr=False, compare=False)

    @property
    def degree(self):
        return len(self.coeffs) - 1

    def subfields(self, d):
        """nfsubfields of degree d, computed once."""
        if d not in self._subfields:
            self._subfields[d] = pari.nfsubfields(self.poly(), d)
        return self._subfields[d]

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
    subs = K.subfields(2)
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
    """(A, B, D) with K = Q(\\sqrt{A + B\\sqrt{D}}) simplest as in imprimitive_quartic.py, or None
    if K is not a C4 or D4 quartic, or is pure (it keeps its pure name)."""
    if K.degree != 4 or pure_radicand(K):
        return None
    fy = K.poly()
    subs = K.subfields(2)
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


_BNF = {}


def _bnf(D):
    """bnfinit of L = Q(sqrt(D)), in the variable w."""
    if D not in _BNF:
        _BNF[D] = pari.bnfinit(pari("w^2") - D, 1)
    return _BNF[D]


def _simplest_abd(u, v, D, n=2):
    """The simplest (A, B) (see _abd_key) with A + B sqrt(D) in (u + v sqrt(D)) L^{*n} or its
    conjugate, L = Q(sqrt(D)), A != 0 and B != 0; returns (A, B, D), or None.

    Every such z has ideal C^n R, with R the n-th power free part of (u + v sqrt(D)) and C an
    integral ideal in a fixed class, and |N(z)| = |A^2 - B^2 D| <= M^2 (1 + |D|) when
    max(|A|, |B|) <= M. So once some (A, B) is known, running over the C of norm up to
    (M^2 (1 + |D|) / N(R))^(1/n) and over the units finds the simplest one.
    """
    bnf = _bnf(D)
    mod = pari("w^2") - D
    beta = u + v * pari.Mod(pari("w"), mod)
    fac = pari.idealfactor(bnf, beta)
    I, R = 1, 1
    for P, e in zip(fac[0], fac[1]):
        I = pari.idealmul(bnf, I, pari.idealpow(bnf, P, int(e) // n))
        if int(e) % n:
            R = pari.idealmul(bnf, R, pari.idealpow(bnf, P, int(e) % n))
    norm_R = int(pari.idealnorm(bnf, R))
    if D > 0:
        eps = pari.Mod(pari.lift(pari("(b) -> b.fu[1]")(bnf)), mod)
        unit_powers = None
    else:
        tu = pari("(b) -> b.tu")(bnf)
        unit_powers = [pari.Mod(pari.lift(tu[1]), mod) ** (n * j) for j in range(int(tu[0]))]
    sqrtD = math.sqrt(abs(D))

    def ab(z):
        """Integers (A, B) with z = (A + B sqrt D) / c^n for the smallest c, B >= 0."""
        a0, b0 = pari.polcoef(pari.lift(z), 0), pari.polcoef(pari.lift(z), 1)
        den = int(pari.lcm(pari.denominator(a0), pari.denominator(b0)))
        c = math.prod(int(p) ** -(-int(e) // n) for p, e in zip(*pari.factor(den))) if den > 1 else 1
        A, B = int(a0 * c ** n), int(b0 * c ** n)
        g = math.gcd(A, B)
        for p in [int(p) for p in pari.factor(g)[0]] if g > 1 else []:
            while A % p ** n == 0 and B % p ** n == 0:
                A //= p ** n
                B //= p ** n
        return (A, abs(B)) if B else None

    def candidates(max_norm, T):
        for C in itertools.chain.from_iterable(pari.ideallist(bnf, max_norm)):
            cls, nu = pari.bnfisprincipal(bnf, pari.idealdiv(bnf, C, I), 1)
            if any(int(c) for c in cls):
                continue
            z0 = beta * pari.nfbasistoalg(bnf, nu) ** n
            if unit_powers is not None:
                zs = [z0 * t for t in unit_powers]
            else:
                zs = [z0 * eps ** (n * k) for k in _unit_range(z0, eps, T, n)]
                if n % 2:  # -1 = (-1)^n is an n-th power too (whatever the sign of eps)
                    zs += [-z for z in zs]
            for z in zs:
                res = ab(z)
                if res and res[0] != 0:
                    yield res

    best = None
    # First a few small ideals for a good starting bound, then everything below it.
    for A, B in candidates(16, math.inf):
        if best is None or _abd_key(A, B) < _abd_key(*best):
            best = (A, B)
    if best is None:
        return None
    M = _abd_key(*best)[0]
    max_norm = int(((M * M * (1 + abs(D))) / norm_R) ** (1 / n)) + 1
    for A, B in candidates(max_norm, M * (1 + sqrtD)):
        if _abd_key(A, B) < _abd_key(*best):
            best = (A, B)
    return best + (D,)


def _abd_key(A, B):
    # as in imprimitive_quartic.simplicity and nested_radical.simplicity, for a fixed D
    return (max(abs(A), abs(B)), abs(B), abs(A), A < 0)


def _log_embeddings(z):
    """(log|u + v sqrt D|, log|u - v sqrt D|) for z = u + v sqrt(D) in a real quadratic field.

    The larger one has no cancellation, and the other comes from the exact norm, so this
    works however large u and v are (z often is a big power)."""
    lz = pari.lift(z)
    u, v = pari.polcoef(lz, 0), pari.polcoef(lz, 1)
    D = -pari.polcoef(z.mod(), 0)
    big = float(pari.log(pari.abs(u) + pari.abs(v) * pari.sqrt(D)))
    small = float(pari.log(pari.abs(u * u - v * v * D))) - big
    return (big, small) if u * v >= 0 else (small, big)


def _unit_range(z, eps, T, n=2):
    """The k with |sigma(z eps^{nk})| <= T at both real places (or the k near the balance
    point when T is infinite)."""
    l1, l2 = _log_embeddings(z)
    le = n * _log_embeddings(eps)[0]
    if math.isinf(T):
        mid = (l2 - l1) / (2 * le)  # balances the two conjugates
        return range(math.floor(mid) - 1, math.ceil(mid) + 2)
    lo = (l2 - math.log(T)) / le
    hi = (math.log(T) - l1) / le
    if le < 0:
        lo, hi = hi, lo
    return range(math.floor(lo) - 1, math.ceil(hi) + 2)


def _quadratic_subfields(K):
    """[(D, sqrt(D) in K as a Mod)] for the quadratic subfields Q(sqrt(D)) of K."""
    fy = K.poly()
    out = []
    for g, h in K.subfields(2):
        disc = pari.poldisc(g)
        D = int(pari.core(disc))
        r = pari.sqrtint(disc / D)
        out.append((D, pari.Mod((2 * h + pari.polcoef(g, 1)) / r, fy)))
    return out


def nested_params(K, max_size=None):
    """(n, a, b, d) with K = Q(\\sqrt[n]{a + b\\sqrt{d}}) of degree 2n, n >= 3, simplest as in
    nested_radical.py, or None. Pure fields are left out (they keep their pure name).

    For each quadratic subfield L = Q(sqrt(d)): if K = L(alpha) with alpha^n = beta in L, then
    v_P(beta) is divisible by n at every prime P of L not above a ramified prime of K or n.
    So, adding primes that generate the class group of L to the set S of those primes, beta
    is an S-unit times an n-th power. The S-units mod n-th powers are a finite group; we
    sieve it with primes of degree 1 of K (as for pure fields, now a linear system mod n),
    test the survivors with nfroots, and find the simplest a + b sqrt(d) for each.

    If given, only names with radicand_size below max_size are looked for (so d with
    |d| >= max_size are skipped).
    """
    if K.degree < 6 or K.degree % 2 or pure_radicand(K):
        return None
    n = K.degree // 2
    fy = K.poly()
    best = None
    for D, s in _quadratic_subfields(K):
        if max_size is not None and abs(D) >= max_size:
            continue
        for beta in _nested_betas(K, n, D, s):
            u, v = beta
            abd = _simplest_abd(u, v, D, n)
            if abd and (best is None or _nested_full_key(*abd) < _nested_full_key(*best)):
                best = abd
    return (n,) + best if best else None


def _nested_full_key(a, b, d):
    return (max(abs(a), abs(b)), abs(b), abs(a), abs(d), a < 0, d < 0)


def _nested_betas(K, n, D, s):
    """The (u, v) with beta = u + v sqrt(D) such that K = L(\\sqrt[n]{beta}), up to n-th powers."""
    bnf = _bnf(D)
    mod = pari("w^2") - D
    primes = set(K.ramps) | {int(p) for p in pari.factor(n)[0]}
    S = [P for p in sorted(primes) for P in pari.idealprimedec(bnf, p)]
    for g in pari("(b) -> b.gen")(bnf):  # so that S-units cover the class group
        for P in pari.idealfactor(bnf, g)[0]:
            if P not in S:
                S.append(P)
    tu = pari("(b) -> b.tu")(bnf)
    gens = [pari.Mod(pari.lift(tu[1]), mod)]
    if D > 0:
        gens.append(pari.Mod(pari.lift(pari("(b) -> b.fu[1]")(bnf)), mod))
    gens += [pari.Mod(pari.lift(pari.nfbasistoalg(bnf, x)), mod) for x in pari.bnfsunit(bnf, S)[0]]
    fy = K.poly()

    def in_K(z):  # u + v w  ->  u + v s
        return pari.polcoef(pari.lift(z), 0) + pari.polcoef(pari.lift(z), 1) * s

    gens_K = [pari.lift(in_K(g)) for g in gens]
    # Conditions: at a degree 1 prime P of K above q, beta = prod g_i^e_i must be an n-th
    # power in F_q, a linear condition on the e_i mod gcd(n, q - 1).
    rows = []
    pdisc = abs(int(pari.poldisc(fy)))
    for q in pari.primes(2000):
        q, g = int(q), math.gcd(n, int(q) - 1)
        if g == 1 or n % q == 0 or pdisc % q == 0:
            continue  # only primes with q = 1 mod a divisor of n give conditions
        h = pow(int(pari.znprimroot(q)), (q - 1) // g, q)
        dlog = {pow(h, j, q): j for j in range(g)}
        for r in pari.polrootsmod(fy, q):
            r = int(pari.lift(r))
            vals = []
            for gk in gens_K:
                den = int(pari.denominator(pari.content(gk)))
                x = int(pari.subst(gk * den, "y", r)) % q
                if den % q == 0 or x == 0:
                    vals = None  # some generator is not a unit at this prime
                    break
                x = x * pow(den, -1, q) % q
                vals.append(dlog[pow(x, (q - 1) // g, q)] * (n // g))
            if vals:
                rows.append(vals)
        if len(rows) >= 3 * len(gens) + 10:
            break
    if not rows:
        return []
    ker = pari.matkermod(pari.matrix(len(rows), len(gens), [x for row in rows for x in row]), n)
    cols = [[int(ker[i, j]) for i in range(len(gens))] for j in range(int(pari.matsize(ker)[1]))]
    if n ** len(cols) > MAX_NESTED_CANDIDATES:
        return []
    w = int(tu[0])  # gens[0] is a root of unity of order w, the others are free
    units = [k for k in range(1, n) if math.gcd(k, n) == 1]
    primes_n = [int(p) for p in pari.factor(n)[0]]
    out, seen = [], set()
    for cs in itertools.product(range(n), repeat=len(cols)):
        e = tuple(sum(c * col[i] for c, col in zip(cs, cols)) % n for i in range(len(gens)))
        if e in seen:
            continue
        # beta and beta^k (k prime to n) give the same field: test one of them
        orbit = {tuple(k * x % n for x in e) for k in units}
        seen |= orbit
        # x^n - beta is reducible when beta is a p-th power in L for a prime p | n, which
        # (up to the root of unity) is when every free exponent is divisible by p
        if any(all(x % p == 0 for x in e[1:]) and (w % p or e[0] % p == 0) for p in primes_n):
            continue
        beta = math.prod((gens[i] ** e[i] for i in range(len(gens))), start=pari.Mod(1, mod))
        # K = L(alpha) when x^n - beta is irreducible over L and has a root in K
        if (len(pari.nffactor(mod, pari("x^%d" % n) - beta)[0]) == 1
                and len(pari.nfroots(fy, pari("x^%d" % n) - in_K(beta)))):
            for k in units:
                lb = pari.lift(beta ** k)
                out.append((pari.polcoef(lb, 0), pari.polcoef(lb, 1)))
    return out


# ---------------------------------------------------------------- composita from subfields

_BLOCKS = {}


def _blocks(coeffs, primes):
    """How a subfield (polredabs coeffs) can appear in a \\Q(g_1, g_2) name: a list of
    (kind, degree, data) with kind "quadratic" (data d), "cyclotomic" (m), "pure" (a) or
    "quartic" ((A, B, D)). As in named_composita.py, each subfield is used as the first of
    cyclotomic, pure and imprimitive quartic that it is."""
    if coeffs in _BLOCKS:
        return _BLOCKS[coeffs]
    d = len(coeffs) - 1
    out = []
    if d == 2:
        out.append(("quadratic", 2, int(pari.core(pari.poldisc(pari.Polrev(list(coeffs)))))))
    else:
        S = from_coeffs(coeffs, primes)
        m = cyclotomic_conductor(S)
        if m:
            out.append(("cyclotomic", d, m))
        else:
            a = pure_radicand(S)
            abd = quartic_params(S) if d == 4 and not a else None
            if a:
                out.append(("pure", d, a))
            elif abd:
                out.append(("quartic", 4, abd))
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
        for g, _ in K.subfields(d):
            subs[d].append(_blocks(polredabs_coeffs(g), K.ramps))
    names = []
    for d1 in divisors:
        d2 = n // d1
        if d1 > d2:
            break
        for i, bs1 in enumerate(subs[d1]):
            for bs2 in subs[d2][i + 1:] if d1 == d2 else subs[d2]:
                for b1 in bs1:
                    for b2 in bs2:
                        names += _pair_names(b1, b2, n)
    return names + _long_composita_names(K, divisors)


_GENERATORS = {}


def _generators(g, primes):
    """For a subfield (polynomial g): (Galois, generators of its best name as a list, poly),
    when it is quadratic, cyclotomic, multiquadratic, Kummer, pure or an imprimitive quartic;
    else None."""
    coeffs = polredabs_coeffs(g)
    if coeffs in _GENERATORS:
        return _GENERATORS[coeffs]
    d = len(coeffs) - 1
    S = from_coeffs(coeffs, primes)
    galois = []  # names of S as a Galois field
    if d == 2:
        galois.append([_symbol(_blocks(coeffs, primes)[0])])
    m = cyclotomic_conductor(S)
    if m:
        galois.append([zeta_symbol(m)])
    basis = multiquadratic_basis(S)
    if basis:
        galois.append([sqrt_symbol(b) for b in basis])
    na = kummer_params(S)
    if na:
        galois.append([zeta_symbol(na[0]), root_symbol(*na)])
    if galois:
        out = (True, min(galois, key=lambda gs: name_key(r"\Q(%s)" % ", ".join(gs))))
    else:
        blocks = _blocks(coeffs, primes) if d > 2 else []
        out = (False, [_symbol(blocks[0])]) if blocks else None
    if out:
        out = out + (pari.Polrev(list(coeffs)),)
    _GENERATORS[coeffs] = out
    return out


def _long_composita_names(K, divisors):
    """\\Q(g_1, ..., g_k) names with k >= 3 generators: K = K_1 K_2 with K_1 a Galois subfield
    (quadratic, cyclotomic, multiquadratic or Kummer), K_2 any named subfield, and
    [K_1][K_2] = [K]; the generators of both best names, in the order of generator_key."""
    n = K.degree
    gens = {d: [x for x in (_generators(g, K.ramps) for g, _ in K.subfields(d)) if x]
            for d in divisors}
    names = []
    for d1 in divisors:
        for g1 in gens[d1]:
            if not g1[0]:
                continue
            for g2 in gens[n // d1]:
                all_gens = sorted(set(g1[1] + g2[1]), key=generator_key)
                if len(all_gens) < 3 or all(re.fullmatch(r"i|\\sqrt\{-?\d+\}", g) for g in all_gens):
                    continue  # two generators: composita_names; all quadratic: multiquadratic
                if int(pari.poldegree(pari.polcompositum(g1[2], g2[2])[0])) == n:
                    names.append(r"\Q(%s)" % ", ".join(all_gens))
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

def all_names(K, families=FAMILY_TESTS, best_only=False):
    """Every name of K found by the tests in families (see FAMILY_TESTS).

    With best_only, tests that cannot give a better name than one already found are
    skipped: nothing after a cyclotomic, real cyclotomic or multiquadratic name (those
    families win over all the radical ones), and nested radicals only below the smallest
    radicand_size found (they rank last among the radical families).
    """
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
    if best_only and names:
        return names
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
    if "nested radical" in families:
        max_size = min(map(radicand_size, names)) if best_only and names else None
        params = nested_params(K, max_size)
        if params:
            names.append(nested_name(*params))
    return list(dict.fromkeys(names))


def best_name(K, families=FAMILY_TESTS):
    """The best name of K (smallest name_key) among all_names, or None."""
    names = all_names(K, families, best_only=True)
    return min(names, key=name_key) if names else None
