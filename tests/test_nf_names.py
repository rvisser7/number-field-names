import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
import tempfile

from nf_names import (family_rank, is_fundamental_discriminant, label_key, merge_data,
                      quadratic_name, read_data, write_data)


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


class TestMerge(unittest.TestCase):
    def test_family_rank(self):
        names = [r"\Q", r"\Q(\sqrt{-3})", r"\Q(\zeta_{8})", r"\Q(\zeta_{7})^+", r"\Q(\sqrt[4]{-1})",
                 r"\Q(\sqrt{2 + \sqrt{5}})", r"\Q(\sqrt{1 + 2i})", r"\Q(\sqrt{2}, \sqrt{3})"]
        self.assertEqual([family_rank(s) for s in names], [0, 1, 2, 3, 4, 5, 5, 6])

    def test_merge(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "degree4.txt")
            write_data(path, [("4.0.256.1", r"\Q(\sqrt[4]{-1})"), ("4.4.1125.1", "custom")])
            # Cyclotomic names win over pure names and over names of no known family...
            self.assertEqual(merge_data(path, {"4.0.256.1": r"\Q(\zeta_{8})",
                                               "4.4.1125.1": r"\Q(\zeta_{15})^+",
                                               "4.0.125.1": r"\Q(\zeta_{5})"}), (1, 2, 0))
            self.assertEqual(dict(read_data(path)), {"4.0.125.1": r"\Q(\zeta_{5})",
                                                     "4.0.256.1": r"\Q(\zeta_{8})",
                                                     "4.4.1125.1": r"\Q(\zeta_{15})^+"})
            # ...but a pure name does not replace a cyclotomic one,
            self.assertEqual(merge_data(path, {"4.0.256.1": r"\Q(\sqrt[4]{-1})"}), (0, 0, 1))
            # unless asked to.
            self.assertEqual(merge_data(path, {"4.0.256.1": r"\Q(\sqrt[4]{-1})"}, replace=True),
                             (0, 1, 0))


if __name__ == "__main__":
    unittest.main()
