import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import is_fundamental_discriminant, label_key, quadratic_name


class TestQuadratic(unittest.TestCase):
    def test_names(self):
        # These match field_pretty on the LMFDB.
        self.assertEqual(quadratic_name("2.0.3.1"), r"\Q(\sqrt{-3})")
        self.assertEqual(quadratic_name("2.0.4.1"), r"\Q(\sqrt{-1})")
        self.assertEqual(quadratic_name("2.0.8.1"), r"\Q(\sqrt{-2})")
        self.assertEqual(quadratic_name("2.2.5.1"), r"\Q(\sqrt{5})")
        self.assertEqual(quadratic_name("2.2.8.1"), r"\Q(\sqrt{2})")
        self.assertEqual(quadratic_name("2.2.12.1"), r"\Q(\sqrt{3})")

    def test_invalid(self):
        for label in ["2.2.4.1", "2.1.5.1", "3.1.23.1", "2.2.9.1"]:
            with self.assertRaises(ValueError):
                quadratic_name(label)

    def test_fundamental(self):
        fund = [D for D in range(-30, 31) if is_fundamental_discriminant(D)]
        self.assertEqual(fund, [-24, -23, -20, -19, -15, -11, -8, -7, -4, -3,
                                5, 8, 12, 13, 17, 21, 24, 28, 29])


class TestLabelOrder(unittest.TestCase):
    def test_numeric(self):
        self.assertLess(label_key("2.0.8.1"), label_key("2.0.11.1"))
        self.assertLess(label_key("2.0.1000.1"), label_key("2.2.5.1"))


if __name__ == "__main__":
    unittest.main()
