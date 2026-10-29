"""Composita of two named fields, e.g. Q(\\zeta_{5}, \\sqrt{2}) or Q(i, \\sqrt{1 + \\sqrt{2}}).

The building blocks are fields we already name, each written by one generator:

  * cyclotomic Q(\\zeta_n)               -> \\zeta_{n}  (Galois)
  * quadratic Q(\\sqrt{d})               -> \\sqrt{d}, or i for d = -1  (Galois)
  * pure Q(\\sqrt[m]{a})                 -> \\sqrt[m]{a}
  * imprimitive quartic Q(\\sqrt{A + B\\sqrt{D}}) -> \\sqrt{A + B\\sqrt{D}}

and we take K = K_1 K_2 where K_1 is cyclotomic or quadratic. As K_1 is Galois the compositum
is unique up to isomorphism, and we keep it only when K_1 and K_2 are linearly disjoint
(degree [K_1 : Q][K_2 : Q]), naming it \\Q(g_1, g_2). Pairs that another family already
covers are skipped (two quadratics: multiquadratic; two cyclotomic: cyclotomic; quadratic and
odd degree pure: radical compositum). A field
reached several ways keeps the first name, with blocks ordered as listed above and then by
size, and where families overlap the other family's name wins (see FAMILIES in nf_names.py),
e.g. Q(\\zeta_{5}, \\sqrt{-3}) = Q(\\zeta_{15}) and Q(\\zeta_{3}, \\sqrt[3]{2}) is Kummer.
"""

from dataclasses import dataclass

from cyclotomic import conductors
from imprimitive_quartic import imprimitive_quartic_fields
from multiquadratic import sqrt_symbol
from nf_names import is_squarefree
from pari_fields import make_field, pari, polredabs_coeffs
from pure import pure_fields


@dataclass
class Block:
    kind: str        # "cyclotomic", "quadratic", "pure" or "quartic"
    symbol: str      # the generator, as LaTeX
    poly: object     # a PARI polynomial defining the field
    degree: int


def blocks(max_degree, bound, quartic_bound):
    """The building blocks, in the order used to pick names."""
    out = []
    for n in [4, 3] + conductors(max_degree // 2, real=False):  # i before zeta_3
        out.append(Block("cyclotomic", "i" if n == 4 else r"\zeta_{%d}" % n,
                         pari.polcyclo(n), int(pari.eulerphi(n))))
    for d in sorted(range(-bound, bound + 1), key=lambda d: (abs(d), -d)):
        if d not in (0, 1) and is_squarefree(d) and d not in (-1, -3):  # i, zeta_3 above
            out.append(Block("quadratic", sqrt_symbol(d), pari(f"x^2 - ({d})"), 2))
    for m in range(3, max_degree // 2 + 1):
        for K in pure_fields(m, bound):
            out.append(Block("pure", K.name[3:-1], pari.Polrev(list(K.coeffs)), m))
    for K in imprimitive_quartic_fields(quartic_bound):
        if r"\sqrt[4]" not in K.name:  # pure quartics are already above
            out.append(Block("quartic", K.name[3:-1], pari.Polrev(list(K.coeffs)), 4))
    # cyclotomic blocks by degree, so the smallest come first
    cyc = sorted((b for b in out if b.kind == "cyclotomic"), key=lambda b: b.degree)
    return cyc + [b for b in out if b.kind != "cyclotomic"]


def named_composita_fields(max_degree, bound, quartic_bound, progress=None):
    """Composita K_1 K_2 as described above, of degree at most max_degree, as Fields."""
    bs = blocks(max_degree, bound, quartic_bound)
    jobs = []
    for i, b1 in enumerate(bs):
        if b1.kind not in ("cyclotomic", "quadratic"):
            break
        for b2 in bs[i + 1:]:
            if b1.degree * b2.degree > max_degree:
                continue
            if b1.kind == b2.kind or b1.degree == b2.degree == 2:
                continue  # multiquadratic or cyclotomic
            if b1.degree == 2 and b2.kind == "pure" and b2.degree % 2:
                continue  # radical compositum
            jobs.append((b1, b2))
    fields = {}
    for i, (b1, b2) in enumerate(jobs, start=1):
        if progress:
            progress(i, len(jobs), f", {len(fields)} fields so far")
        comp = pari.polcompositum(b1.poly, b2.poly)
        f = comp[0]
        if len(comp) > 1 or pari.poldegree(f) != b1.degree * b2.degree:
            continue  # not linearly disjoint
        coeffs = polredabs_coeffs(pari.polredbest(f))
        fields.setdefault(coeffs, r"\Q(%s, %s)" % (b1.symbol, b2.symbol))
    return [make_field(name, c) for c, name in fields.items()]
