import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import FAMILIES, family_rank, name_key

try:
    from cyclotomic import cyclotomic_fields, cyclotomic_name, real_cyclotomic_fields, real_cyclotomic_name
    from imprimitive_quartic import imprimitive_quartic_fields, quartic_name
    from multiquadratic import greedy_bases, multiquadratic_coeffs, multiquadratic_name
    from nested_radical import nested_name, nested_radical_fields
    from named_composita import named_composita_fields
    from pari_fields import make_field, pari, polredabs_coeffs
    from pure import pure_fields, pure_name
    from radical import composita_fields, kummer_fields, kummer_name
    from recognise import (all_names, best_name, composita_names, cyclotomic_conductor, from_coeffs,
                           from_label, kummer_params, multiquadratic_basis, nested_params,
                           pure_radicand, quartic_params)
except SystemExit:  # cypari2 not installed
    best_name = None


def family(name):
    rank = family_rank(name)
    return FAMILIES[rank][0] if rank < len(FAMILIES) else None


def field(poly):
    return from_coeffs(polredabs_coeffs(pari(poly)))


@unittest.skipIf(best_name is None, "needs cypari2")
class TestRecognise(unittest.TestCase):
    """The backward tests give the names the forward generators give (whenever the forward
    bounds are big enough to reach the best name)."""

    def check(self, fields, recognise):
        for F in fields:
            self.assertEqual(recognise(from_coeffs(F.coeffs)), F.name)

    def test_pure(self):
        for n, bound in [(3, 30), (4, 30), (5, 10), (6, 20), (8, 5)]:
            self.check(pure_fields(n, bound), lambda K: pure_name(K.degree, pure_radicand(K)))

    def test_cyclotomic(self):
        self.check(cyclotomic_fields(16), lambda K: cyclotomic_name(cyclotomic_conductor(K)))
        self.check(real_cyclotomic_fields(16),
                   lambda K: real_cyclotomic_name(cyclotomic_conductor(K, real=True)))

    def test_multiquadratic(self):
        fields = [make_field(multiquadratic_name(b), multiquadratic_coeffs(b))
                  for k, bound in [(2, 10), (3, 5)] for b in greedy_bases(k, bound)]
        self.check(fields, lambda K: multiquadratic_name(multiquadratic_basis(K)))

    def test_imprimitive_quartic(self):
        fields = [F for F in imprimitive_quartic_fields(6) if r"\sqrt[4]" not in F.name]
        self.check(fields, lambda K: quartic_name(*quartic_params(K)))

    def test_nested_radical(self):
        for n, bound in [(3, 4), (4, 3), (6, 2)]:
            self.check(nested_radical_fields(n, bound), lambda K: nested_name(*nested_params(K)))

    def test_kummer(self):
        self.check(kummer_fields([3, 4], 10), lambda K: kummer_name(*kummer_params(K)))

    def test_radical_composita(self):
        def best(K):
            return min((x for x in composita_names(K) if family(x) == "radical compositum"),
                       key=name_key)
        self.check(composita_fields([(2, 3), (3, 4)], 5), best)

    def test_named_composita(self):
        # With small forward bounds, the backward name can be better (e.g. \sqrt{-5} before
        # \sqrt{6}) but never worse.
        for F in named_composita_fields(8, 3, 3):
            names = all_names(from_coeffs(F.coeffs))
            self.assertIn(F.name, names)
            self.assertLessEqual(name_key(min(names, key=name_key)), name_key(F.name))

    def test_best_names(self):
        examples = {
            "x^6 + 3": r"\Q(\sqrt[6]{-3})",                                 # not \Q(\zeta_{3}, \sqrt[3]{3})
            "polcompositum(x^2 + 3, x^3 - 2)[1]": r"\Q(\zeta_{3}, \sqrt[3]{2})",  # = \Q(\sqrt[6]{-108})
            "polcompositum(x^2 - 2, x^3 - 3)[1]": r"\Q(\sqrt{2}, \sqrt[3]{3})",   # = \Q(\sqrt[6]{72})
            "x^4 - 10*x^2 + 1": r"\Q(\zeta_{24})^+",                         # = \Q(\sqrt{2}, \sqrt{3})
            "x^4 + 1": r"\Q(\zeta_{8})",
            "x^4 - 2*x^2 - 1": r"\Q(\sqrt{1 + \sqrt{2}})",
            "polcompositum(polcyclo(5), x^2 - 2)[1]": r"\Q(\zeta_{5}, \sqrt{2})",
            "polcompositum(x^2 + 1, x^4 - 2*x^2 - 1)[1]": r"\Q(i, \sqrt{1 + \sqrt{2}})",
            "(x^3 - 1)^2 + 24": r"\Q(\sqrt[3]{1 + 2\sqrt{-6}})",
            "x^4 - 12": r"\Q(\sqrt[4]{12})",                              # not \Q(\sqrt{3 + 2\sqrt{3}})
            "polcompositum(polcyclo(8), x^3 - 3)[1]": r"\Q(\zeta_{8}, \sqrt[3]{3})",  # not \Q(i, \sqrt{2}, ...)
            "polcompositum(polcompositum(x^2 - 2, x^2 - 3)[1], x^3 - 5)[1]":
                r"\Q(\sqrt{2}, \sqrt{3}, \sqrt[3]{5})",
            "polcompositum(polcompositum(x^2 + 3, x^3 - 2)[1], x^3 - 3)[1]":
                r"\Q(\zeta_{3}, \sqrt[3]{2}, \sqrt[3]{3})",
            "polcompositum(polcompositum(polcyclo(5), x^2 - 2)[1], x^2 - 3)[1]":
                r"\Q(\zeta_{5}, \sqrt{2}, \sqrt{3})",
            "polcompositum(x^3 - 2, x^5 - 3)[1]": r"\Q(\sqrt[3]{2}, \sqrt[5]{3})",
        }
        for poly, name in examples.items():
            self.assertEqual(best_name(field(poly)), name, poly)
        self.assertIsNone(best_name(field("x^4 - x - 1")))  # S4 quartic

    def test_from_label(self):
        K = from_label("3.1.108.1", [-2, 0, 0, 1])
        self.assertEqual((K.degree, K.r1, K.disc_abs, K.ramps), (3, 1, 108, (2, 3)))
        self.assertEqual(best_name(K), r"\Q(\sqrt[3]{2})")


if __name__ == "__main__":
    unittest.main()
