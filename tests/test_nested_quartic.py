import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import FAMILIES, family_rank, parse_nested_quartic

try:
    from nested_quartic import nested_quartic_fields, nested_quartic_name, nested_quartic_poly
    from pari_fields import is_isomorphic, polredabs_coeffs
except SystemExit:  # cypari2 not installed
    nested_quartic_fields = None


class TestNames(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_nested_quartic(r"\Q(\sqrt{1 + \sqrt{1 + \sqrt{2}}})"), (1, 0, 1, 0, 1, 1, 2))
        self.assertEqual(parse_nested_quartic(
            r"\Q(\sqrt{-1 + 2\sqrt{2} - \sqrt{1 + \sqrt{2}} + 3\sqrt{2}\sqrt{1 + \sqrt{2}}})"),
            (-1, 2, -1, 3, 1, 1, 2))
        self.assertEqual(parse_nested_quartic(r"\Q(\sqrt{1 - 2i\sqrt{1 + 2i}})"), (1, 0, 0, -2, 1, 2, -1))
        for name in [r"\Q(\sqrt{1 + \sqrt{2}})",                    # imprimitive quartic
                     r"\Q(\sqrt{\sqrt{1 + \sqrt{2}}})",              # a single term
                     r"\Q(\sqrt{\sqrt{1 + \sqrt{2}} + \sqrt{2}})",   # not in the order of the basis
                     r"\Q(\sqrt{1 + \sqrt{3}\sqrt{1 + \sqrt{2}}})"]:  # \sqrt{3} is not in the field
            self.assertIsNone(parse_nested_quartic(name), name)

    def test_family(self):
        self.assertEqual(FAMILIES[family_rank(r"\Q(\sqrt{1 + \sqrt{1 + \sqrt{2}}})")][0], "nested quartic")
        self.assertEqual(FAMILIES[-2][0], "nested quartic")  # every other name wins (but Cardano: cubic)


@unittest.skipIf(nested_quartic_fields is None, "needs cypari2")
class TestNestedQuartic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fields = nested_quartic_fields(1, 2)
        cls.by_name = {K.name: K for K in cls.fields}

    def test_names(self):
        self.assertEqual(nested_quartic_name(1, 0, 1, 0, 1, 1, 2), r"\Q(\sqrt{1 + \sqrt{1 + \sqrt{2}}})")
        self.assertEqual(nested_quartic_name(0, 1, 0, -2, 1, 2, -1), r"\Q(\sqrt{i - 2i\sqrt{1 + 2i}})")

    def test_names_are_certified(self):
        for K in self.fields:
            p = parse_nested_quartic(K.name)
            self.assertIsNotNone(p, K.name)
            self.assertEqual(nested_quartic_name(*p), K.name)
            self.assertEqual(K.degree, 8)
            self.assertTrue(is_isomorphic(polredabs_coeffs(nested_quartic_poly(*p)), K.coeffs), K.name)

    def test_simplest_name_kept(self):
        # gamma -> -gamma flips the signs of c and d, so the first of them is positive
        self.assertIn(r"\Q(\sqrt{1 + \sqrt{1 + \sqrt{2}}})", self.by_name)
        self.assertNotIn(r"\Q(\sqrt{1 - \sqrt{1 + \sqrt{2}}})", self.by_name)
        for K in self.fields:
            a, b, c, d = parse_nested_quartic(K.name)[:4]
            self.assertTrue(c > 0 or (c == 0 and d > 0), K.name)


if __name__ == "__main__":
    unittest.main()
