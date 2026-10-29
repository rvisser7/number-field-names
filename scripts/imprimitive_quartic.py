"""Imprimitive quartic fields Q(\\sqrt{A + B\\sqrt{D}}), computed with PARI (via cypari2).

Every quartic field K with a quadratic subfield Q(\\sqrt{D}) is Q(\\sqrt{A + B\\sqrt{D}}) for
integers A, B and squarefree D (clear denominators by multiplying by a rational square),
with minimal polynomial (x^2 - A)^2 - B^2 D. We run over |A|, |B|, |D| <= bound and identify
each field by its polredabs polynomial.

Fields with Galois group V4 have three quadratic subfields and are left out: they are
biquadratic, Q(\\sqrt{a}, \\sqrt{b}), which is a nicer name (the LMFDB's field_pretty does
the same). The others (C4 and D4) have just one quadratic subfield, so D is determined by K,
and we choose (A, B) to be as simple as possible:

  * B > 0 (B and -B give conjugate, so isomorphic, fields);
  * smallest max(|A|, |B|), then smallest |B|, then smallest |A|, then A > 0.

Pure fields are left out, as generate_pure.py names them: those with A = 0, where
Q(\\sqrt{B\\sqrt{D}}) = Q(\\sqrt[4]{B^2 D}), and those like Q(\\sqrt[4]{12}) = Q(\\sqrt{3 + 2\\sqrt{3}})
that only have A = 0 beyond the bound (found with recognise.pure_radicand).
"""

from math import gcd

from nf_names import is_squarefree
from pari_fields import make_field, pari, polredabs_coeffs


def quartic_name(A, B, D):
    if A == 0:
        return r"\Q(\sqrt[4]{%d})" % (B * B * D)
    sqrt_D = "i" if D == -1 else r"\sqrt{%d}" % D
    return r"\Q(\sqrt{%d + %s%s})" % (A, "" if B == 1 else B, sqrt_D)


def simplicity(ABD):
    A, B, D = ABD
    return (max(abs(A), abs(B)), abs(B), abs(A), A < 0)


def quartic_poly(A, B, D):
    return pari(f"(x^2 - ({A}))^2 - ({B})^2 * ({D})")


def imprimitive_quartic_fields(bound, progress=None):
    """Return the distinct non-V4 fields Q(\\sqrt{A + B\\sqrt{D}}) with |A|, |B|, |D| <= bound,
    each named by its simplest (A, B, D), as Fields.

    If given, progress(done, total, extra) is called after each D.
    """
    Ds = [D for D in range(-bound, bound + 1) if D not in (0, 1) and is_squarefree(D)]
    reps = {}  # polredabs coeffs -> list of (A, B, D)
    for i, D in enumerate(Ds, start=1):
        for A in range(-bound, bound + 1):
            for B in range(1, bound + 1):
                g = gcd(A, B)
                if any(g % (p * p) == 0 for p in range(2, int(g ** 0.5) + 1)):
                    continue  # (A/p^2, B/p^2) gives the same field and is simpler
                f = quartic_poly(A, B, D)
                if pari.polisirreducible(f):
                    reps.setdefault(polredabs_coeffs(f), []).append((A, B, D))
        if progress:
            progress(i, len(Ds), f", {len(reps)} fields so far")
    out = []
    for coeffs, abds in reps.items():
        if len(pari.nfsubfields(pari.Polrev(list(coeffs)), 2)) == 3:
            continue  # V4: biquadratic
        if not any(A == 0 for A, _, _ in abds) and not _is_pure(coeffs):
            out.append(make_field(quartic_name(*min(abds, key=simplicity)), coeffs))
    return out


def _is_pure(coeffs):
    from recognise import from_coeffs, pure_radicand  # recognise imports this module
    return pure_radicand(from_coeffs(coeffs)) is not None
