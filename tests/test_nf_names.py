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
    def test_lmfdb_order(self):
        # The order the LMFDB lists quadratic fields in.
        labels = ["2.0.3.1", "2.0.4.1", "2.2.5.1", "2.0.7.1", "2.0.8.1", "2.2.8.1",
                  "2.0.11.1", "2.2.12.1", "2.2.13.1", "2.0.15.1", "2.2.17.1", "2.0.19.1",
                  "2.0.20.1", "2.2.21.1", "2.0.23.1", "2.0.24.1", "2.2.24.1", "2.2.28.1",
                  "2.2.29.1", "2.0.31.1", "2.2.33.1", "2.0.35.1", "2.2.37.1", "2.0.39.1",
                  "2.0.40.1", "2.2.40.1", "2.2.41.1", "2.0.43.1", "2.2.44.1", "2.0.47.1",
                  "2.0.51.1"]
        self.assertEqual(sorted(reversed(labels), key=label_key), labels)

    def test_sign(self):
        # disc sign is (-1)^s: 3.1 is negative, 3.3 positive; 4.2 negative, 4.0 and 4.4 positive.
        self.assertLess(label_key("3.1.49.1"), label_key("3.3.49.1"))
        self.assertLess(label_key("4.2.400.1"), label_key("4.0.400.1"))
        self.assertLess(label_key("4.0.400.1"), label_key("4.4.400.1"))
        self.assertLess(label_key("2.2.5.1"), label_key("2.0.7.1"))


if __name__ == "__main__":
    unittest.main()
