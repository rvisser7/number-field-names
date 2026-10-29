"""Nested radicals Q(\\sqrt[n]{a + b\\sqrt{d}}), n >= 3, computed with PARI (via cypari2).

These generalise the imprimitive quartics Q(\\sqrt{A + B\\sqrt{D}}) (the case n = 2, see
imprimitive_quartic.py). For integers a, b != 0 and squarefree d, alpha = \\sqrt[n]{a + b\\sqrt{d}}
is a root of (x^n - a)^2 - b^2 d, and when that is irreducible Q(alpha) has degree 2n and
contains Q(\\sqrt{d}) = Q(alpha^n). We run over |a|, |b|, |d| <= bound and identify each field
by its polredabs polynomial. Of the (a, b, d) giving the same field we keep the simplest:

  * b > 0 (b and -b give conjugate, so isomorphic, fields);
  * no prime p with p^n dividing both a and b (dividing by p^n gives the same field);
  * smallest max(|a|, |b|), then |b|, then |a|, then |d|, then a > 0, then d > 0.

Pure fields are left out, as generate_pure.py names them: those with a = 0, where
Q(\\sqrt[n]{b\\sqrt{d}}) = Q(\\sqrt[2n]{b^2 d}), and those like Q(\\sqrt[3]{5 + 5i}) = Q(\\sqrt[6]{-400})
that only have a = 0 beyond the bound (found with recognise.pure_radicand).
As for the quartics, d = -1 is written i, e.g. \\Q(\\sqrt[3]{1 + 2i}). Where this family overlaps
another one (e.g. Q(\\sqrt{2}, \\sqrt[3]{3}) when a + b\\sqrt{d} is 3 times a cube in
Q(\\sqrt{d})), the name with the smaller name_key wins (see nf_names.py).
"""

from nf_names import is_squarefree
from pari_fields import make_field, pari, polredabs_coeffs


def nested_name(n, a, b, d):
    sqrt_d = "i" if d == -1 else r"\sqrt{%d}" % d
    return r"\Q(\sqrt[%d]{%d + %s%s})" % (n, a, "" if b == 1 else b, sqrt_d)


def simplicity(abd):
    a, b, d = abd
    return (max(abs(a), abs(b)), abs(b), abs(a), abs(d), a < 0, d < 0)


def nested_poly(n, a, b, d):
    return pari(f"(x^{n} - ({a}))^2 - ({b})^2 * ({d})")


def has_nth_power_factor(n, a, b):
    """Whether some prime power p^n divides both a and b."""
    g = abs(pari.gcd(a, b))
    return any(e >= n for e in pari.factor(g)[1]) if g > 1 else False


def nested_radical_fields(n, bound, progress=None):
    """Return the distinct fields Q(\\sqrt[n]{a + b\\sqrt{d}}) of degree 2n with |a|, |b|, |d| <= bound
    and a != 0 (see above), each named by its simplest (a, b, d), as Fields.

    If given, progress(done, total, extra) is called after each d.
    """
    ds = [d for d in range(-bound, bound + 1) if d not in (0, 1) and is_squarefree(d)]
    reps = {}  # polredabs coeffs -> list of (a, b, d)
    for i, d in enumerate(ds, start=1):
        for a in range(-bound, bound + 1):
            for b in range(1, bound + 1):
                if has_nth_power_factor(n, a, b):
                    continue
                f = nested_poly(n, a, b, d)
                if pari.polisirreducible(f):
                    reps.setdefault(polredabs_coeffs(f), []).append((a, b, d))
        if progress:
            progress(i, len(ds), f", {len(reps)} fields so far")
    from recognise import from_coeffs, pure_radicand
    return [make_field(nested_name(n, *min(abds, key=simplicity)), coeffs)
            for coeffs, abds in reps.items()
            if all(a != 0 for a, _, _ in abds) and pure_radicand(from_coeffs(coeffs)) is None]
