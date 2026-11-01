import os
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)
from nf_names import data_path, read_data, write_data

try:
    from cardano import cardano_fields
    from cyclotomic import cyclotomic_fields, real_cyclotomic_fields
    from gaussian_period import period_fields
    from imprimitive_quartic import imprimitive_quartic_fields
    from multiquadratic import greedy_bases, multiquadratic_coeffs, multiquadratic_name
    from name_field import NameError_, check_name, name_poly, parse_name
    from named_composita import named_composita_fields
    from nested_biquadratic import nested_biquadratic_fields
    from nested_pure import nested_pure_fields
    from nested_quartic import nested_quartic_fields
    from nested_radical import nested_radical_fields
    from pari_fields import is_isomorphic, make_field, pari, polredabs_coeffs
    from pure import pure_fields
    from radical import composita_fields, kummer_fields
    from recognise import best_name, from_coeffs
except SystemExit:  # cypari2 not installed
    name_poly = None


def coeffs(P):
    return tuple(int(c) for c in pari.Vecrev(P))


def label(F, i=1):
    return f"{F.degree}.{F.r1}.{abs(F.disc)}.{i}"


@unittest.skipIf(name_poly is None, "needs cypari2")
class TestNamePoly(unittest.TestCase):
    def test_examples(self):
        for name, poly in [(r"\Q", "x"),
                           (r"\Q(\sqrt{-3})", "x^2 - x + 1"),
                           (r"\Q(i, \sqrt{2})", "x^4 + 1"),
                           (r"\Q(\zeta_{7})^+", "x^3 - x^2 - 2*x + 1"),
                           (r"\Q(\sqrt[3]{-2})", "x^3 - 2"),
                           (r"\Q(\sqrt{1 + 2i})", "x^4 - 2*x^3 + 2"),
                           (r"\Q(\sqrt[3]{2 + \sqrt{3}} + \sqrt[3]{2 - \sqrt{3}})", "x^3 - 3*x - 4"),
                           (r"\Q(\zeta_{13} + \zeta_{13}^{5} + \zeta_{13}^{8} + \zeta_{13}^{12})", "x^3 - x^2 - 4*x - 1"),
                           (r"\Q(\sqrt{2}\sqrt{1 + \sqrt{2}})", "x^4 - 2*x^2 - 1")]:
            self.assertEqual(pari.polredabs(name_poly(name)), pari(poly), name)

    def test_principal_roots(self):
        # \sqrt{-2}\sqrt{-3} = -\sqrt{6}: these two names are different fields
        a = name_poly(r"\Q(\sqrt{1 + \sqrt{-2} + \sqrt{-3} + \sqrt{6}})")
        b = name_poly(r"\Q(\sqrt{1 + \sqrt{-2} + \sqrt{-3} - \sqrt{6}})")
        self.assertFalse(is_isomorphic(coeffs(a), coeffs(b)))
        self.assertTrue(is_isomorphic(coeffs(a), coeffs(name_poly(r"\Q(\sqrt{1 + \sqrt{-2} + \sqrt{-3} - \sqrt{-2}\sqrt{-3}})"))))

    def test_parse(self):
        self.assertEqual(parse_name(r"\Q(i, \sqrt[3]{2})"),
                         [("add", [(1, ("mul", [("i",)]))]),
                          ("add", [(1, ("mul", [("root", 3, ("add", [(1, ("mul", [("int", 2)]))]))]))])])
        for name in [r"\Q(\alpha)", r"\Q(\sqrt{2}", "Q(S{2})", r"\Q(\sqrt{})"]:
            with self.assertRaises(NameError_, msg=name):
                name_poly(name)

    def test_every_family(self):
        # Each generator's names, read as LaTeX, give the generator's fields.
        families = [cyclotomic_fields(12), real_cyclotomic_fields(12), pure_fields(3, 12), pure_fields(4, 6),
                    imprimitive_quartic_fields(3),
                    [make_field(multiquadratic_name(b), multiquadratic_coeffs(b)) for b in greedy_bases(2, 8)],
                    kummer_fields([3, 4], 4), composita_fields([(2, 3)], 3), named_composita_fields(8, 2, 2),
                    nested_radical_fields(3, 3), nested_biquadratic_fields(1, 6)[::3], nested_pure_fields(3, 1, 3),
                    nested_quartic_fields(1, 2)[::5], cardano_fields(6), period_fields(range(3, 7))]
        for fields in families:
            for F in fields:
                self.assertTrue(is_isomorphic(coeffs(name_poly(F.name)), F.coeffs), F.name)

    def test_backward_names(self):
        for poly in ["polcompositum(polcompositum(x^2 - 2, x^2 - 3)[1], x^3 - 5)[1]",
                     "x^8 - 24*x^6 + 144*x^4 - 288*x^2 + 144", "((x^2 - 1)^2 - 1)^2 - 2"]:
            c = polredabs_coeffs(pari(poly))
            self.assertTrue(is_isomorphic(coeffs(name_poly(best_name(from_coeffs(c)))), c), poly)


@unittest.skipIf(name_poly is None, "needs cypari2")
class TestCheckName(unittest.TestCase):
    def test_check(self):
        good = r"\Q(\sqrt[3]{14 + 7\sqrt{-3}} + \sqrt[3]{14 - 7\sqrt{-3}})"  # x^3 - 21x + 28
        self.assertIsNone(check_name("3.3.3969.1", good))
        self.assertIsNone(check_name("3.3.3969.1", good, [-28, -21, 0, 1]))
        # the other cyclic cubic field of discriminant 3969
        self.assertIn("not isomorphic", check_name("3.3.3969.2", good, [-35, -21, 0, 1]))
        self.assertIn("|disc| 3969, not 49", check_name("3.3.49.1", good))
        self.assertIn("real places", check_name("3.1.3969.1", good))
        self.assertIn("degree 3", check_name("4.4.3969.1", good))
        self.assertIn("cannot read", check_name("3.3.3969.1", "custom"))
        self.assertIsNone(check_name("1.1.1.1", r"\Q"))
        self.assertIsNone(check_name("2.0.4.1", r"\Q(\sqrt{-1})"))

    def test_check_data(self):
        with tempfile.TemporaryDirectory() as d:
            write_data(data_path(3, d), [("3.1.108.1", r"\Q(\sqrt[3]{2})"), ("3.3.49.1", r"\Q(\zeta_{7})^+"),
                                         ("3.3.81.1", r"\Q(\zeta_{7})^+")])  # wrong
            cmd = [sys.executable, os.path.join(SCRIPTS, "check_data.py"), "-d", d]
            self.assertEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
            r = subprocess.run(cmd + ["--fields"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 1)
            self.assertIn("3.3.81.1|\\Q(\\zeta_{7})^+: the name gives a field of |disc| 49, not 81", r.stdout)
            self.assertIn("checked 3 names, 1 wrong", r.stdout)
            r = subprocess.run(cmd + ["--fields", "--drop-wrong"], capture_output=True, text=True)
            self.assertEqual([lab for lab, _ in read_data(data_path(3, d))], ["3.3.49.1", "3.1.108.1"])
            self.assertEqual(subprocess.run(cmd + ["--fields"], capture_output=True).returncode, 0)


if __name__ == "__main__":
    unittest.main()
