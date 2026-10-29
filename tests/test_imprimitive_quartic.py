import os
import re
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
try:
    from imprimitive_quartic import imprimitive_quartic_fields, quartic_name
    from pari_fields import is_isomorphic, pari, polredabs_coeffs
    from pure import pure_fields
except SystemExit:  # cypari2 not installed
    imprimitive_quartic_fields = None

NAME_RE = re.compile(r"^\\Q\(\\sqrt\{(-?\d+) \+ (\d*)(i|\\sqrt\{(-?\d+)\})\}\)$")


def coeffs_from_name(name):
    """The polynomial (constant term first) that a quartic name describes."""
    m = NAME_RE.match(name)
    if m:
        A, B = int(m[1]), int(m[2] or 1)
        D = -1 if m[3] == "i" else int(m[4])
        return [A * A - B * B * D, 0, -2 * A, 0, 1]
    m = re.match(r"^\\Q\(\\sqrt\[4\]\{(-?\d+)\}\)$", name)
    return [-int(m[1]), 0, 0, 0, 1]


@unittest.skipIf(imprimitive_quartic_fields is None, "needs cypari2")
class TestImprimitiveQuartic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fields = imprimitive_quartic_fields(8)
        cls.by_name = {K.name: K for K in cls.fields}

    def test_names(self):
        self.assertEqual(quartic_name(2, 1, 5), r"\Q(\sqrt{2 + \sqrt{5}})")
        self.assertEqual(quartic_name(-1, 2, -3), r"\Q(\sqrt{-1 + 2\sqrt{-3}})")
        self.assertEqual(quartic_name(1, 2, -1), r"\Q(\sqrt{1 + 2i})")
        self.assertEqual(quartic_name(0, 2, 3), r"\Q(\sqrt[4]{12})")

    def test_examples(self):
        K = self.by_name[r"\Q(\sqrt{2 + \sqrt{5}})"]
        self.assertEqual((K.disc, K.r1), (-400, 2))
        self.assertEqual(self.by_name[r"\Q(\sqrt{-5 + 2\sqrt{5}})"].disc, 125)    # Q(zeta_5)

    def test_names_are_certified(self):
        # Every name describes the field it is attached to.
        for K in self.fields:
            self.assertTrue(is_isomorphic(coeffs_from_name(K.name), K.coeffs), K.name)

    def test_no_biquadratic(self):
        # Q(sqrt(5 + 2 sqrt 6)) = Q(sqrt 2, sqrt 3) is left for the biquadratic names.
        v4 = polredabs_coeffs(pari("x^4 - 10*x^2 + 1"))
        self.assertNotIn(v4, {K.coeffs for K in self.fields})
        for K in self.fields:
            self.assertEqual(len(pari.nfsubfields(pari.Polrev(list(K.coeffs)), 2)), 1, K.name)

    def test_no_pure_fields(self):
        # Pure quartics are left to generate_pure.py, also Q(\sqrt[4]{12}) = Q(\sqrt{3 + 2\sqrt{3}}).
        pure = {K.coeffs for K in pure_fields(4, 8 ** 3)}
        self.assertIn(polredabs_coeffs(pari("x^4 - 12")), pure)
        self.assertFalse(pure & {K.coeffs for K in self.fields})


if __name__ == "__main__":
    unittest.main()
