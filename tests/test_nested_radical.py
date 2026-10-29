import os
import re
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import FAMILIES, family_rank

try:
    from nested_radical import nested_name, nested_radical_fields
    from pari_fields import is_isomorphic, pari, polredabs_coeffs
    from pure import pure_fields
except SystemExit:  # cypari2 not installed
    nested_radical_fields = None

NAME_RE = re.compile(r"^\\Q\(\\sqrt\[(\d+)\]\{(-?\d+) \+ (\d*)(i|\\sqrt\{(-?\d+)\})\}\)$")


def poly_from_name(name):
    """The polynomial (x^n - a)^2 - b^2 d that a name describes."""
    m = NAME_RE.match(name)
    n, a, b = int(m[1]), int(m[2]), int(m[3] or 1)
    d = -1 if m[4] == "i" else int(m[5])
    return pari(f"(x^{n} - ({a}))^2 - ({b})^2 * ({d})")


@unittest.skipIf(nested_radical_fields is None, "needs cypari2")
class TestNestedRadical(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fields = {n: nested_radical_fields(n, 8) for n in (3, 4)}
        cls.by_name = {K.name: K for fs in cls.fields.values() for K in fs}

    def test_names(self):
        self.assertEqual(nested_name(3, 1, 1, 2), r"\Q(\sqrt[3]{1 + \sqrt{2}})")
        self.assertEqual(nested_name(3, -1, 2, -1), r"\Q(\sqrt[3]{-1 + 2i})")
        self.assertEqual(nested_name(5, 2, 3, -7), r"\Q(\sqrt[5]{2 + 3\sqrt{-7}})")

    def test_examples(self):
        K = self.by_name[r"\Q(\sqrt[3]{1 + \sqrt{2}})"]
        self.assertEqual((K.degree, K.r1), (6, 2))
        self.assertEqual(self.by_name[r"\Q(\sqrt[4]{1 + 4i})"].degree, 8)

    def test_simplest_name_kept(self):
        # -1 + sqrt 2 = 1 / (1 + sqrt 2), so its cube root gives the same field;
        # 8 + 8 sqrt 2 is 2^3 (1 + sqrt 2).
        for name in [r"\Q(\sqrt[3]{-1 + \sqrt{2}})", r"\Q(\sqrt[3]{8 + 8\sqrt{2}})"]:
            self.assertNotIn(name, self.by_name)

    def test_names_are_certified(self):
        for n, fs in self.fields.items():
            for K in fs:
                self.assertEqual(K.degree, 2 * n, K.name)
                self.assertTrue(is_isomorphic(polredabs_coeffs(poly_from_name(K.name)), K.coeffs),
                                K.name)

    def test_pure_fields_left_out(self):
        # Q(\sqrt[6]{m}) with |m| <= 8 is reached with a = 0, so it is left to generate_pure.py.
        pure = {K.coeffs for K in pure_fields(6, 8)}
        self.assertFalse(pure & {K.coeffs for K in self.fields[3]})

    def test_family(self):
        rank = [name for name, _ in FAMILIES].index("nested radical")
        for name in self.by_name:
            self.assertEqual(family_rank(name), rank, name)


if __name__ == "__main__":
    unittest.main()
