import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import FAMILIES, family_rank, parse_cardano

try:
    from cardano import cardano_fields, cardano_k, cardano_name, cardano_params, cardano_poly
    from pari_fields import is_isomorphic, pari, polredabs_coeffs
except SystemExit:  # cypari2 not installed
    cardano_fields = None


class TestNames(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_cardano(r"\Q(\sqrt[3]{2 + \sqrt{3}} + \sqrt[3]{2 - \sqrt{3}})"), (2, 1, 3))
        self.assertEqual(parse_cardano(r"\Q(\sqrt[3]{68 + 17i} + \sqrt[3]{68 - 17i})"), (68, 17, -1))
        for name in [r"\Q(\sqrt[3]{2 + \sqrt{3}} + \sqrt[3]{2 - \sqrt{5}})",
                     r"\Q(\sqrt[3]{2 + 1\sqrt{3}} + \sqrt[3]{2 - 1\sqrt{3}})"]:
            self.assertIsNone(parse_cardano(name), name)

    def test_family(self):
        name = r"\Q(\sqrt[3]{2 + \sqrt{3}} + \sqrt[3]{2 - \sqrt{3}})"
        self.assertEqual(FAMILIES[family_rank(name)][0], "cardano cubic")
        self.assertEqual(FAMILIES[-1][0], "cardano cubic")  # every other name wins


@unittest.skipIf(cardano_fields is None, "needs cypari2")
class TestCardano(unittest.TestCase):
    def test_params(self):
        self.assertEqual(cardano_params(2, 1), (2, 1, 3))  # x^3 - 3x - 4
        self.assertEqual(cardano_params(-2, 1), (2, 1, 3))  # alpha -> -alpha
        self.assertEqual(cardano_params(16, 4), (2, 1, 3))  # 8 (2 + \sqrt{3}): alpha -> 2 alpha
        self.assertIsNone(cardano_params(3, 2))             # 9 - 8 = 1: pure
        self.assertEqual(cardano_k(2, 1, 3), 1)

    def test_names_are_certified(self):
        for K in cardano_fields(15):
            u, e, d = parse_cardano(K.name)
            self.assertEqual(cardano_name(u, e, d), K.name)
            self.assertEqual(int(pari.core(-3 * K.disc)), d)  # the quadratic resolvent
            self.assertTrue(is_isomorphic(polredabs_coeffs(cardano_poly(u, cardano_k(u, e, d))), K.coeffs))


if __name__ == "__main__":
    unittest.main()
