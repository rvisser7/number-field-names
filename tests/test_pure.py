import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
try:
    from pure import is_isomorphic, pure_fields
    from lmfdb_lookup import find_label
except ImportError:  # cypari2 not installed
    pure_fields = None


@unittest.skipIf(pure_fields is None, "needs cypari2")
class TestPure(unittest.TestCase):
    def test_cubic(self):
        # 4 ~ 2, 9 ~ 3 and -m ~ m; 1 and 8 give reducible polynomials.
        fields = {K.m: K for K in pure_fields(3, 10)}
        self.assertEqual(sorted(fields), [2, 3, 5, 6, 7, 10])
        K = fields[2]
        self.assertEqual(K.name, r"\Q(\sqrt[3]{2})")
        self.assertEqual((K.coeffs, K.disc, K.r1), ((-2, 0, 0, 1), -108, 1))

    def test_quartic(self):
        fields = {K.m: K for K in pure_fields(4, 4)}
        # x^4 - 1, x^4 - 4 and x^4 + 4 are reducible.
        self.assertEqual(sorted(fields), [-3, -2, -1, 2, 3])
        self.assertEqual((fields[-1].disc, fields[-1].r1), (256, 0))    # Q(zeta_8)
        self.assertEqual((fields[2].disc, fields[2].r1), (-2048, 2))

    def test_find_label(self):
        candidates = [("3.1.108.1", [-2, 0, 0, 1]), ("3.3.108.1", [1, -3, 0, 1])]
        self.assertEqual(find_label((-2, 0, 0, 1), 1, candidates, is_isomorphic), "3.1.108.1")
        # Not the same polynomial, but the same field: found by the isomorphism test.
        self.assertEqual(find_label((-4, 0, 0, 1), 1, candidates, is_isomorphic), "3.1.108.1")
        self.assertIsNone(find_label((-3, 0, 0, 1), 1, candidates, is_isomorphic))


if __name__ == "__main__":
    unittest.main()
