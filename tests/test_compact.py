import glob
import itertools
import os
import re
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)

from cardano import cardano_name
from cyclotomic import cyclotomic_name, real_cyclotomic_name
from imprimitive_quartic import quartic_name
from multiquadratic import multiquadratic_name
from nested_biquadratic import nested_biquadratic_name
from nested_pure import nested_pure_name
from nested_quartic import nested_quartic_name
from nested_radical import nested_name
from nf_names import (compress_name, data_path, expand_name, is_compact, merge_data, name_key,
                      read_data, write_data)
from pure import pure_name
from radical import compositum_name, kummer_name

R = range(-3, 4)


def sample_names():
    """Names of every family, over small parameters (including negative and zero ones)."""
    yield r"\Q"
    for n in [3, 4, 5, 7, 8, 12, 15]:
        yield cyclotomic_name(n)
        yield real_cyclotomic_name(n)
    for n, m in itertools.product([2, 3, 4, 6], [-12, -2, -1, 2, 3, 10]):
        yield pure_name(n, m)
        yield kummer_name(n, m)
    yield compositum_name(2, -3, 3, 2)
    yield compositum_name(3, 2, 4, -5)
    for basis in [(-1, 2), (2, 3), (-1, 2, -3), (-2, 5, 7)]:
        yield multiquadratic_name(basis)
    for A, B, D in itertools.product(R, [1, 2], [-1, 2, -5]):
        yield quartic_name(A, B, D)
        yield nested_name(3, A, B, D)
        yield nested_name(5, A, -B, D)
    for A, B, C, D in itertools.product(R, R, R, R):
        if sum(1 for c in (B, C, D) if c) >= 2:
            yield nested_biquadratic_name(A, B, C, D, -1, 2)
            yield nested_biquadratic_name(A, B, C, D, 2, -3)
    for coeffs in itertools.product(R, R, R):
        if any(coeffs[1:]):
            yield nested_pure_name(3, 2, coeffs)
            yield nested_pure_name(3, 12, coeffs)
    for a, b, c, d in itertools.product(R, R, R, R):
        if c or d:
            yield nested_quartic_name(a, b, c, d, 1, 1, 2)
            yield nested_quartic_name(a, b, c, d, 1, 2, -1)
    for u, e, d in [(1, 1, -1), (2, 1, 3), (5, 2, -3), (-2, 3, 5)]:
        yield cardano_name(u, e, d)
    # Named composita: two generators of other families.
    yield r"\Q(\zeta_{5}, \sqrt{2})"
    yield r"\Q(i, \sqrt[4]{3})"
    yield r"\Q(\sqrt{-3}, \sqrt{1 + 2i})"


def test_literals():
    """Every name written out in the other tests and the scripts."""
    root = os.path.dirname(SCRIPTS)
    for path in glob.glob(os.path.join(root, "tests", "*.py")) + glob.glob(os.path.join(SCRIPTS, "*.py")):
        if os.path.basename(path) == "test_compact.py":
            continue
        with open(path) as f:
            for name in re.findall(r'r"(\\Q[^"%]*)"', f.read()):
                yield name


class TestCompact(unittest.TestCase):
    def test_examples(self):
        for full, compact in [(r"\Q", "Q"),
                              (r"\Q(\sqrt{-3})", "Q(S{-3})"),
                              (r"\Q(\sqrt{1 + 2i})", "Q(S{1+2i})"),
                              (r"\Q(\sqrt[3]{2})", "Q(S[3]{2})"),
                              (r"\Q(\zeta_{7})", "Q(Z_{7})"),
                              (r"\Q(\zeta_{7})^+", "Q(Z_{7})^+"),
                              (r"\Q(i, \sqrt{2}, \sqrt{-3})", "Q(i,S{2},S{-3})"),
                              (r"\Q(\sqrt{-3 - 2\sqrt{2}})", "Q(S{-3-2S{2}})"),
                              (r"\Q(\sqrt{1 - i})", "Q(S{1-i})"),
                              (r"\Q(\sqrt[3]{2 + \sqrt{3}} + \sqrt[3]{2 - \sqrt{3}})",
                               "Q(S[3]{2+S{3}}+S[3]{2-S{3}})")]:
            self.assertEqual(compress_name(full), compact)
            self.assertEqual(expand_name(compact), full)

    def test_round_trip(self):
        names = set(sample_names()) | set(test_literals())
        self.assertGreater(len(names), 1000)
        for name in names:
            c = compress_name(name)
            self.assertTrue(is_compact(c), name)
            self.assertEqual(expand_name(c), name)

    def test_full_names_unchanged(self):
        # Names in data files from before the compact form read back as they are.
        for name in [r"\Q(\sqrt{-3})", r"\Q(\zeta_{7})^+", r"\Q(i, \sqrt{2})"]:
            self.assertEqual(expand_name(name), name)

    def test_other_names(self):
        # Names that the compact form cannot hold are stored as they are...
        for name in ["custom", r"\Q(\alpha)", r"\mathbb{Q}(\sqrt{2})"]:
            self.assertEqual(expand_name(compress_name(name)), name)
        # ...unless that would not read back the same.
        with self.assertRaises(ValueError):
            compress_name("Some field")

    def test_files(self):
        with tempfile.TemporaryDirectory() as d:
            path = data_path(4, d)
            # A data file in the old full form, as written before the compact form.
            with open(path, "w") as f:
                f.write("label|pretty_name\ntext|text\n\n4.0.125.1|\\Q(\\zeta_{5})\n"
                        "4.2.1088.1|\\Q(\\sqrt{1 + 2i})\n")
            self.assertEqual(read_data(path), [("4.0.125.1", r"\Q(\zeta_{5})"),
                                               ("4.2.1088.1", r"\Q(\sqrt{1 + 2i})")])
            # A merge stores every name compact, the old ones too.
            self.assertEqual(merge_data(path, {"4.4.1125.1": r"\Q(\zeta_{15})^+",
                                               "4.2.1088.1": r"\Q(\sqrt{1 + 2i})"}), (1, 0, 0))
            with open(path) as f:
                self.assertEqual(f.read(), "label|pretty_name\ntext|text\n\n4.0.125.1|Q(Z_{5})\n"
                                           "4.2.1088.1|Q(S{1+2i})\n4.4.1125.1|Q(Z_{15})^+\n")
            # Comparing a new name with a stored one uses the full names.
            self.assertLess(name_key(r"\Q(\zeta_{5})"), name_key(r"\Q(\sqrt[4]{-5})"))
            self.assertEqual(merge_data(path, {"4.0.125.1": r"\Q(\sqrt[4]{-5})"}), (0, 0, 1))
            write_data(path, [("4.0.125.1", r"\Q(\zeta_{5})")], compact=False)
            self.assertEqual(read_data(path, expand=False), [("4.0.125.1", r"\Q(\zeta_{5})")])

    def test_scripts(self):
        with tempfile.TemporaryDirectory() as d:
            data, out = os.path.join(d, "data"), os.path.join(d, "export")
            os.mkdir(data)
            with open(data_path(3, data), "w") as f:
                f.write("label|pretty_name\ntext|text\n\n3.3.49.1|\\Q(\\zeta_{7})^+\n"
                        "3.1.108.1|\\Q(\\sqrt[3]{2})\n")
            write_data(data_path(2, data), [("2.0.3.1", r"\Q(\sqrt{-3})")])

            def run(*args):
                return subprocess.run([sys.executable, os.path.join(SCRIPTS, args[0])] + list(args[1:]),
                                      capture_output=True, text=True)

            # check_data.py flags names not in the compact form,
            self.assertEqual(run("check_data.py", "-d", data).returncode, 1)
            # compact_data.py converts them,
            r = run("compact_data.py", "-d", data)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("2 of 2 names compressed", r.stdout)
            self.assertEqual(read_data(data_path(3, data), expand=False),
                             [("3.3.49.1", "Q(Z_{7})^+"), ("3.1.108.1", "Q(S[3]{2})")])
            self.assertEqual(run("check_data.py", "-d", data).returncode, 0)
            # and export_lmfdb.py writes the full names, per degree and in one file.
            r = run("export_lmfdb.py", "-d", data, "-o", out, "--combined", os.path.join(out, "all.txt"))
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(read_data(data_path(3, out), expand=False),
                             [("3.3.49.1", r"\Q(\zeta_{7})^+"), ("3.1.108.1", r"\Q(\sqrt[3]{2})")])
            self.assertEqual(read_data(os.path.join(out, "all.txt"), expand=False),
                             [("2.0.3.1", r"\Q(\sqrt{-3})"), ("3.3.49.1", r"\Q(\zeta_{7})^+"),
                              ("3.1.108.1", r"\Q(\sqrt[3]{2})")])


if __name__ == "__main__":
    unittest.main()
