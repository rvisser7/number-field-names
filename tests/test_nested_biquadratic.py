import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import FAMILIES, family_rank, name_key, parse_nested_biquadratic

try:
    from nested_biquadratic import (nested_biquadratic_fields, nested_biquadratic_name,
                                    nested_biquadratic_poly)
    from pari_fields import is_isomorphic, pari, polredabs_coeffs
except SystemExit:  # cypari2 not installed
    nested_biquadratic_fields = None

Q8 = r"\Q(\sqrt{6 + 2\sqrt{3} + \sqrt{6}})"  # = Q(\sqrt{(2 + \sqrt{2})(3 + \sqrt{3})})


class TestNames(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_nested_biquadratic(r"\Q(\sqrt{6 + 3\sqrt{2} + 2\sqrt{3} + \sqrt{6}})"),
                         (6, 3, 2, 1, 2, 3))
        self.assertEqual(parse_nested_biquadratic(r"\Q(\sqrt{i - \sqrt{2}})"), (0, 1, -1, 0, -1, 2))
        self.assertEqual(parse_nested_biquadratic(Q8), (6, 0, 2, 1, 2, 3))
        for name in [r"\Q(\sqrt{1 + 2\sqrt{3}})",            # imprimitive quartic
                     r"\Q(\sqrt{\sqrt{3} + \sqrt{2}})",       # not in the order d_1, d_2, d_3
                     r"\Q(\sqrt{1 + \sqrt{2} + \sqrt{3} + \sqrt{5}})"]:  # not biquadratic
            self.assertIsNone(parse_nested_biquadratic(name), name)

    def test_family(self):
        self.assertEqual(FAMILIES[family_rank(Q8)][0], "nested biquadratic")
        # every other name wins, except nested pure (degree 2n, n an odd prime, so never 8) and
        # nested quartic, which this one wins over
        self.assertEqual([f for f, _ in FAMILIES][-3:], ["nested biquadratic", "nested pure", "nested quartic"])
        self.assertLess(name_key(r"\Q(i, \sqrt{1 + \sqrt{2}})"), name_key(r"\Q(\sqrt{i + \sqrt{2}})"))


@unittest.skipIf(nested_biquadratic_fields is None, "needs cypari2")
class TestNestedBiquadratic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fields = nested_biquadratic_fields(1, 3)
        cls.by_name = {K.name: K for K in cls.fields}

    def test_names(self):
        self.assertEqual(nested_biquadratic_name(6, 3, 2, 1, 2, 3),
                         r"\Q(\sqrt{6 + 3\sqrt{2} + 2\sqrt{3} + \sqrt{6}})")
        self.assertEqual(nested_biquadratic_name(-1, 0, 2, -3, -1, 2), r"\Q(\sqrt{-1 + 2\sqrt{2} - 3\sqrt{-2}})")
        self.assertEqual(nested_biquadratic_name(0, 1, 1, 0, -1, 2), r"\Q(\sqrt{i + \sqrt{2}})")

    def test_quaternion_field(self):
        f = nested_biquadratic_poly(6, 3, 2, 1, 2, 3)
        self.assertEqual(abs(int(pari.nfdisc(f))), 2 ** 24 * 3 ** 6)  # 8.0.12230590464.1
        # its simplest name (test_recognise checks it is the simplest)
        self.assertTrue(is_isomorphic(polredabs_coeffs(f), polredabs_coeffs(
            nested_biquadratic_poly(*parse_nested_biquadratic(Q8)))))

    def test_names_are_certified(self):
        for K in self.fields:
            p = parse_nested_biquadratic(K.name)
            self.assertIsNotNone(p, K.name)
            self.assertEqual(nested_biquadratic_name(*p), K.name)
            self.assertTrue(is_isomorphic(polredabs_coeffs(nested_biquadratic_poly(*p)), K.coeffs), K.name)

    def test_simplest_name_kept(self):
        # conjugates flip the signs of two of B, C, D; C >= 0 when B = 0, and so on
        for K in self.fields:
            A, B, C, D, _, _ = parse_nested_biquadratic(K.name)
            self.assertTrue(B >= 0 and C >= 0, K.name)


if __name__ == "__main__":
    unittest.main()
