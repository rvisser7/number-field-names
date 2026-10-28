import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
try:
    from cyclotomic import conductors, cyclotomic_fields, real_cyclotomic_fields
except SystemExit:  # cypari2 not installed
    cyclotomic_fields = None

# A few entries of cycloinfo and rcycloinfo in the LMFDB's web_number_field.py.
LMFDB_CYCLO = {5: "4.0.125.1", 7: "6.0.16807.1", 8: "4.0.256.1", 12: "4.0.144.1",
               15: "8.0.1265625.1", 16: "8.0.16777216.1", 24: "8.0.5308416.1",
               132: "40.0.118511797886229481159007653491590053243629014721874976833536.1"}
LMFDB_RCYCLO = {7: "3.3.49.1", 9: "3.3.81.1", 15: "4.4.1125.1", 16: "4.4.2048.1",
                20: "4.4.2000.1", 24: "4.4.2304.1", 300: "40.40.3247321025468409061431884765625"
                "0000000000000000000000000000000000000000.1"}


def label_prefix(K):
    return f"{K.degree}.{K.r1}.{abs(K.disc)}."


@unittest.skipIf(cyclotomic_fields is None, "needs cypari2")
class TestCyclotomic(unittest.TestCase):
    def test_conductors(self):
        # Degree 4: Q(zeta_n) for n = 5, 8, 12 (10 = 2*5 gives the same field as 5).
        self.assertEqual(conductors(4, real=False), [5, 8, 12])
        self.assertEqual(conductors(6, real=False), [5, 7, 8, 9, 12])
        self.assertEqual(conductors(3, real=True), [7, 9])

    def test_against_lmfdb(self):
        cyc = {K.name: K for K in cyclotomic_fields(40)}
        for n, label in LMFDB_CYCLO.items():
            self.assertTrue(label.startswith(label_prefix(cyc[r"\Q(\zeta_{%d})" % n])), n)
        rcyc = {K.name: K for K in real_cyclotomic_fields(40)}
        for n, label in LMFDB_RCYCLO.items():
            self.assertTrue(label.startswith(label_prefix(rcyc[r"\Q(\zeta_{%d})^+" % n])), n)

    def test_no_small_degree(self):
        names = {K.name for K in cyclotomic_fields(6) + real_cyclotomic_fields(6)}
        for n in [3, 4]:
            self.assertNotIn(r"\Q(\zeta_{%d})" % n, names)      # quadratic
        for n in [5, 8, 12]:
            self.assertNotIn(r"\Q(\zeta_{%d})^+" % n, names)    # quadratic


if __name__ == "__main__":
    unittest.main()
