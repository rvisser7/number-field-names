import os
import re
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import family_rank, FAMILIES

try:
    from named_composita import named_composita_fields
    from pari_fields import pari, polredabs_coeffs
except SystemExit:  # cypari2 not installed
    named_composita_fields = None


def generator_poly(s):
    """A polynomial for one generator, read off its LaTeX."""
    if s == "i":
        return pari("x^2 + 1")
    m = re.match(r"^\\zeta_\{(\d+)\}$", s)
    if m:
        return pari.polcyclo(int(m[1]))
    m = re.match(r"^\\sqrt\[(\d+)\]\{(-?\d+)\}$", s)
    if m:
        return pari(f"x^{m[1]} - ({m[2]})")
    m = re.match(r"^\\sqrt\{(-?\d+) \+ (\d*)(i|\\sqrt\{(-?\d+)\})\}$", s)
    if m:
        A, B, D = int(m[1]), int(m[2] or 1), (-1 if m[3] == "i" else int(m[4]))
        return pari(f"(x^2 - ({A}))^2 - ({B})^2 * ({D})")
    return pari(f"x^2 - ({re.match(r'^.sqrt.(-?[0-9]+).$', s)[1]})")


def coeffs_from_name(name):
    g1, g2 = re.match(r"^\\Q\((\\sqrt\{[^,]*\}|[^,]*), (.*)\)$", name).groups()
    f = pari.polcompositum(generator_poly(g1), generator_poly(g2))[0]
    return polredabs_coeffs(pari.polredbest(f))


@unittest.skipIf(named_composita_fields is None, "needs cypari2")
class TestNamedComposita(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fields = named_composita_fields(12, 6, 2)
        cls.names = {K.name for K in cls.fields}

    def test_examples(self):
        for name in [r"\Q(\zeta_{5}, \sqrt{2})", r"\Q(i, \sqrt{1 + \sqrt{2}})",
                     r"\Q(\sqrt{2}, \sqrt[4]{3})", r"\Q(i, \sqrt[4]{3})"]:
            self.assertIn(name, self.names)
        # Not linearly disjoint (sqrt 5 is in Q(zeta_5)), or covered by other families.
        for name in [r"\Q(\zeta_{5}, \sqrt{5})", r"\Q(i, \sqrt{2})", r"\Q(\sqrt{2}, \sqrt[3]{3})"]:
            self.assertNotIn(name, self.names)

    def test_names_are_certified(self):
        for K in self.fields:
            self.assertEqual(coeffs_from_name(K.name), K.coeffs, K.name)
            self.assertEqual(K.degree, pari.poldegree(pari.Polrev(list(K.coeffs))))

    def test_ranked_below_other_families(self):
        names = [name for name, _ in FAMILIES]
        self.assertEqual(family_rank(r"\Q(\zeta_{5}, \sqrt{2})"), names.index("compositum"))
        self.assertGreater(family_rank(r"\Q(\zeta_{5}, \sqrt{2})"), names.index("radical compositum"))
        self.assertLess(family_rank(r"\Q(\zeta_{15})"), family_rank(r"\Q(\zeta_{5}, \sqrt{-3})"))


if __name__ == "__main__":
    unittest.main()
