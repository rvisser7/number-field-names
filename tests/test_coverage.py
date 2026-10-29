import contextlib
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from coverage import galois_counts, parse_int, print_galois_report
from nf_names import data_path, write_data

FIELDS = [("3.1.23.1", "3T2"), ("3.3.49.1", "3T1"), ("3.1.31.1", "3T2"), ("3.1.44.1", "3T2"),
          ("3.3.81.1", "3T1"), ("3.1.108.1", "3T2"), ("3.3.148.1", "3T2"), ("3.3.169.1", "3T1"),
          ("3.3.169.2", "3T1")]


class FakeTable:
    def __init__(self, rows):
        self.rows = rows

    def search(self, query, cols):
        out = []
        for r in self.rows:
            disc = query.get("disc_abs", {}).get("$lte")
            if r["degree"] == query["degree"] and (disc is None or r["disc_abs"] <= disc):
                out.append({c: r[c] for c in cols})
        return out


class FakeDB:
    nf_fields = FakeTable([{"label": lab, "galois_label": g, "degree": 3,
                            "disc_abs": int(lab.split(".")[2])} for lab, g in FIELDS])
    gps_transitive = None  # a search on it fails, so the report leaves out group names


class TestCoverage(unittest.TestCase):
    def test_galois_counts(self):
        named = {"3.3.49.1", "3.3.81.1", "3.1.108.1"}
        counts = galois_counts(FIELDS, named, examples=2)
        self.assertEqual(counts["3T1"], (4, 2, ["3.3.169.1", "3.3.169.2"]))
        # smallest unnamed in LMFDB order, whatever order the rows come in
        self.assertEqual(galois_counts(FIELDS[::-1], named, examples=2)["3T2"],
                         (5, 1, ["3.1.23.1", "3.1.31.1"]))

    def test_parse_int(self):
        self.assertEqual([parse_int(s) for s in ["10^8", "1e6", "42"]], [10 ** 8, 10 ** 6, 42])

    def test_report(self):
        with tempfile.TemporaryDirectory() as d:
            write_data(data_path(3, d), [("3.3.49.1", r"\Q(\zeta_{7})^+"), ("3.1.108.1", r"\Q(\sqrt[3]{2})")])
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                print_galois_report(FakeDB, 3, d, max_disc=100, top=1)
        lines = out.getvalue().splitlines()
        self.assertEqual(lines[0], "degree 3, |disc| <= 100: 5 fields, 1 named (20.00%)")
        self.assertTrue(lines[2].split()[:3] == ["3T2", "3", "0"], lines[2])
        self.assertIn("3.1.23.1, 3.1.31.1, 3.1.44.1", lines[2])
        self.assertEqual(lines[3], "  (1 more groups, 1 unnamed fields)")


if __name__ == "__main__":
    unittest.main()
