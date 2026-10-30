import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import FAMILIES, family_rank, parse_nested_pure, pure_basis

try:
    from nested_pure import nested_pure_fields, nested_pure_name, nested_pure_poly
    from pari_fields import is_isomorphic, polredabs_coeffs
    from pure import pure_fields
except SystemExit:  # cypari2 not installed
    nested_pure_fields = None


class TestNames(unittest.TestCase):
    def test_basis(self):
        self.assertEqual(pure_basis(3, 2), [(1, 1), (1, 2), (1, 4)])
        self.assertEqual(pure_basis(3, 12), [(1, 1), (1, 12), (2, 18)])  # \sqrt[3]{144} = 2\sqrt[3]{18}

    def test_parse(self):
        self.assertEqual(parse_nested_pure(r"\Q(\sqrt{1 + \sqrt[3]{2} + \sqrt[3]{4}})"), (3, 2, (1, 1, 1)))
        self.assertEqual(parse_nested_pure(r"\Q(\sqrt{-1 + 2\sqrt[3]{18}})"), (3, 12, (-1, 0, 2)))
        self.assertEqual(parse_nested_pure(r"\Q(\sqrt{2 + \sqrt[5]{4}})"), (5, 2, (2, 0, 1, 0, 0)))
        for name in [r"\Q(\sqrt{\sqrt[3]{4} + \sqrt[3]{2}})",   # not in the order of the basis
                     r"\Q(\sqrt{1 + \sqrt[3]{3} + \sqrt[3]{4}})",  # not one pure field
                     r"\Q(\sqrt{1 + \sqrt[4]{2}})"]:               # n not an odd prime
            self.assertIsNone(parse_nested_pure(name), name)

    def test_family(self):
        name = r"\Q(\sqrt{-1 + \sqrt[3]{2}})"
        self.assertEqual(FAMILIES[family_rank(name)][0], "nested pure")
        # every other name wins (the nested biquadratic and quartic fields have degree 8, never 2n)
        self.assertEqual([f for f, _ in FAMILIES][-4:-1], ["nested biquadratic", "nested pure", "nested quartic"])


@unittest.skipIf(nested_pure_fields is None, "needs cypari2")
class TestNestedPure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fields = nested_pure_fields(3, 1, 12)
        cls.by_name = {K.name: K for K in cls.fields}

    def test_names(self):
        self.assertEqual(nested_pure_name(3, 2, (1, 1, 1)), r"\Q(\sqrt{1 + \sqrt[3]{2} + \sqrt[3]{4}})")
        self.assertEqual(nested_pure_name(3, 12, (0, 1, -3)), r"\Q(\sqrt{\sqrt[3]{12} - 3\sqrt[3]{18}})")

    def test_simplest_name_kept(self):
        # 1 + \sqrt[3]{2} + \sqrt[3]{4} = 1 / (\sqrt[3]{2} - 1), the same up to squares
        self.assertIn(r"\Q(\sqrt{-1 + \sqrt[3]{2}})", self.by_name)
        self.assertNotIn(r"\Q(\sqrt{1 + \sqrt[3]{2} + \sqrt[3]{4}})", self.by_name)

    def test_names_are_certified(self):
        for K in self.fields:
            p = parse_nested_pure(K.name)
            self.assertIsNotNone(p, K.name)
            self.assertEqual(nested_pure_name(*p), K.name)
            self.assertEqual(K.degree, 6)
            self.assertTrue(is_isomorphic(polredabs_coeffs(nested_pure_poly(*p)), K.coeffs), K.name)

    def test_single_terms_are_pure(self):
        # \sqrt{\sqrt[3]{2}} = \sqrt[6]{2}: left to generate_pure.py
        self.assertNotIn(polredabs_coeffs(nested_pure_poly(3, 2, (0, 1, 0))), {K.coeffs for K in self.fields})
        self.assertIn(polredabs_coeffs(nested_pure_poly(3, 2, (0, 1, 0))), {K.coeffs for K in pure_fields(6, 2)})


if __name__ == "__main__":
    unittest.main()
