import os
import re
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import family_rank

try:
    from pari_fields import pari, polredabs_coeffs
    from radical import composita_fields, compositum_name, kummer_fields, kummer_name
except SystemExit:  # cypari2 not installed
    kummer_fields = None


def coeffs_from_name(name):
    """polredabs polynomial of the field a name describes, built without radical.py."""
    def poly(s):
        if s == "i":
            return pari("x^2 + 1")
        m = re.match(r"^\\zeta_\{(\d+)\}$", s)
        if m:
            return pari.polcyclo(int(m[1]))
        m = re.match(r"^\\sqrt\[(\d+)\]\{(-?\d+)\}$", s)
        if m:
            return pari(f"x^{m[1]} - ({m[2]})")
        return pari(f"x^2 - ({re.match(r'^.sqrt.(-?[0-9]+).$', s)[1]})")
    parts = re.match(r"^\\Q\((.*)\)$", name)[1].split(", ")
    f = poly(parts[0])
    for s in parts[1:]:
        f = pari.polredbest(pari.polcompositum(f, poly(s))[0])
    return polredabs_coeffs(f)


@unittest.skipIf(kummer_fields is None, "needs cypari2")
class TestRadical(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.kummer = kummer_fields([3, 4, 6], 8)
        cls.composita = composita_fields([(2, 3), (2, 5)], 6)

    def test_names(self):
        self.assertEqual(kummer_name(3, 2), r"\Q(\zeta_{3}, \sqrt[3]{2})")
        self.assertEqual(kummer_name(4, 2), r"\Q(i, \sqrt[4]{2})")
        self.assertEqual(kummer_name(6, 2), r"\Q(\zeta_{3}, \sqrt[6]{2})")
        self.assertEqual(compositum_name(2, -1, 3, 2), r"\Q(i, \sqrt[3]{2})")
        self.assertEqual(compositum_name(2, 5, 3, 2), r"\Q(\sqrt{5}, \sqrt[3]{2})")

    def test_names_are_certified(self):
        for K in self.kummer + self.composita:
            self.assertEqual(coeffs_from_name(K.name), K.coeffs, K.name)

    def test_kummer(self):
        by_name = {K.name: K for K in self.kummer}
        K = by_name[r"\Q(\zeta_{3}, \sqrt[3]{2})"]
        self.assertEqual((K.degree, K.r1, K.disc), (6, 0, -34992))
        self.assertEqual(by_name[r"\Q(i, \sqrt[4]{2})"].degree, 8)
        # Q(\sqrt[3]{4}) = Q(\sqrt[3]{2}), and Q(\sqrt[6]{-3}) already contains zeta_3.
        self.assertNotIn(r"\Q(\zeta_{3}, \sqrt[3]{4})", by_name)
        self.assertNotIn(r"\Q(\zeta_{3}, \sqrt[6]{-3})", by_name)

    def test_priorities(self):
        # Q(sqrt(-3), cbrt(2)) is the Kummer field Q(zeta_3, cbrt(2)): the Kummer name wins.
        kummer = {K.coeffs: K.name for K in self.kummer}
        overlap = [(kummer[K.coeffs], K.name) for K in self.composita if K.coeffs in kummer]
        self.assertIn((r"\Q(\zeta_{3}, \sqrt[3]{2})", r"\Q(\sqrt{-3}, \sqrt[3]{2})"), overlap)
        for k, c in overlap:
            self.assertLess(family_rank(k), family_rank(c))
        # Q(sqrt 2, cbrt 2) = Q(\sqrt[6]{2}): the pure name wins.
        c = {K.name for K in self.composita}
        self.assertIn(r"\Q(\sqrt{2}, \sqrt[3]{2})", c)
        self.assertLess(family_rank(r"\Q(\sqrt[6]{2})"), family_rank(r"\Q(\sqrt{2}, \sqrt[3]{2})"))
        self.assertEqual({family_rank(K.name) for K in self.composita}, {8})


if __name__ == "__main__":
    unittest.main()
