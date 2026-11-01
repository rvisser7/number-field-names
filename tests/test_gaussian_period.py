import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import FAMILIES, family_rank, name_key, parse_period

try:
    from gaussian_period import (conductors, is_conductor, period_fields, period_name, period_poly,
                                 periods, small_subgroups)
    from pari_fields import is_isomorphic, pari, polredabs_coeffs
except SystemExit:  # cypari2 not installed
    periods = None

NAME_13 = r"\Q(\zeta_{13} + \zeta_{13}^{5} + \zeta_{13}^{8} + \zeta_{13}^{12})"


class TestNames(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_period(NAME_13), (13, (1, 5, 8, 12)))
        self.assertEqual(parse_period(r"\Q(\zeta_{21} + \zeta_{21}^{13})"), (21, (1, 13)))
        for name in [r"\Q(\zeta_{13} + \zeta_{13}^{5})",              # not a subgroup
                     r"\Q(\zeta_{13} + \zeta_{12}^{5} + \zeta_{13}^{8} + \zeta_{13}^{12})",
                     r"\Q(\zeta_{13} + \zeta_{13}^{8} + \zeta_{13}^{5} + \zeta_{13}^{12})",  # order
                     r"\Q(\zeta_{13})"]:
            self.assertIsNone(parse_period(name), name)

    def test_family(self):
        self.assertEqual(FAMILIES[family_rank(NAME_13)][0], "gaussian period")
        # radical names win, the nested families and Cardano cubics lose
        self.assertLess(name_key(r"\Q(\sqrt[3]{28 + 84\sqrt{-3}})"), name_key(r"\Q(\zeta_{21} + \zeta_{21}^{13})"))
        self.assertLess(name_key(r"\Q(\zeta_{5}, \sqrt{-7})"), name_key(r"\Q(\zeta_{35} + \zeta_{35}^{11} + \zeta_{35}^{16})"))
        self.assertLess(name_key(NAME_13), name_key(r"\Q(\sqrt[3]{2 + \sqrt{3}} + \sqrt[3]{2 - \sqrt{3}})"))
        self.assertLess(name_key(NAME_13), name_key(r"\Q(\sqrt{1 + \sqrt{1 + \sqrt{2}}})"))


@unittest.skipIf(periods is None, "needs cypari2")
class TestPeriods(unittest.TestCase):
    def test_subgroups(self):
        self.assertEqual(small_subgroups(13, 4), [(1, 5, 8, 12)])
        self.assertEqual(small_subgroups(15, 4), [(1, 2, 4, 8), (1, 4, 7, 13), (1, 4, 11, 14)])
        self.assertTrue(is_conductor(13, (1, 5, 8, 12)))
        self.assertFalse(is_conductor(21, (1, 8)))  # = Q(\zeta_7)
        self.assertEqual(conductors(3), [7, 9, 13, 21, 28, 36])

    def test_poly(self):
        self.assertEqual(period_poly(13, (1, 5, 8, 12)), pari("x^3 + x^2 - 4*x + 1"))
        self.assertIsNone(period_poly(9, (1, 4, 7)))  # zeta_9 (1 + zeta_3 + zeta_3^2) = 0

    def test_cubics(self):
        # Only 3.3.169.1: the other cyclic cubics need more terms or are real cyclotomic.
        self.assertEqual([(P.f, P.H, P.disc_abs, P.r1) for P in periods(3)], [(13, (1, 5, 8, 12), 169, 3)])

    def test_fields_are_right(self):
        # Each name describes its field: compare with PARI's galoissubcyclo, which also gives
        # the conductor.
        for n in range(3, 9):
            for P, K in zip(periods(n), period_fields([n])):
                self.assertEqual(K.name, period_name(P.f, P.H))
                g, f = pari.galoissubcyclo(P.f, list(P.H), 2)
                self.assertEqual(int(f), P.f)
                self.assertTrue(is_isomorphic(polredabs_coeffs(g), K.coeffs), K.name)
                self.assertEqual((abs(K.disc), K.r1), (P.disc_abs, P.r1))
            # and no two names give the same field
            self.assertEqual(len({K.coeffs for K in period_fields([n])}), len(periods(n)))


if __name__ == "__main__":
    unittest.main()
