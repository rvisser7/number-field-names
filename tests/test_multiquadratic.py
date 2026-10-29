import itertools
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from multiquadratic import (disc_and_r1, greedy_bases, multiquadratic_name, order_key, span,
                            squarefree_part)
from nf_names import is_squarefree

try:
    from multiquadratic import multiquadratic_coeffs
    from pari_fields import make_field
except SystemExit:  # cypari2 not installed
    make_field = None


def lmfdb_greedy(ds):
    """The basis field_pretty picks from the list of quadratic subfields ds."""
    out, spanned = [], set()
    for d in sorted(ds, key=order_key):
        if d not in spanned:
            out.append(d)
            spanned |= {d} | {squarefree_part(d * a) for a in spanned}
    return tuple(out)


class TestMultiquadratic(unittest.TestCase):
    def test_names(self):
        self.assertEqual(multiquadratic_name((-1, 2)), r"\Q(i, \sqrt{2})")
        self.assertEqual(multiquadratic_name((2, 3, 5)), r"\Q(\sqrt{2}, \sqrt{3}, \sqrt{5})")

    def test_greedy_bases_match_lmfdb(self):
        # Brute force: every independent set of small d, reduced to the basis field_pretty uses.
        for rank, N in [(2, 20), (3, 10), (4, 7)]:
            cands = [d for d in range(-N, N + 1) if d not in (0, 1) and is_squarefree(d)]
            want = set()
            for sub in itertools.combinations(cands, rank):
                ds = span(sub)
                if 1 not in ds and len(set(ds)) == 2 ** rank - 1:
                    g = lmfdb_greedy(ds)
                    if max(abs(d) for d in g) <= N:
                        want.add(g)
            got = list(greedy_bases(rank, N))
            self.assertEqual(len(got), len(set(got)))
            self.assertEqual(set(got), want)

    def test_examples(self):
        bases = set(greedy_bases(2, 6))
        self.assertIn((2, 3), bases)
        self.assertNotIn((3, 6), bases)    # same field as (2, 3)
        # Labels from the LMFDB.
        for basis, label in [((-1, 2), "4.0.256.1"), ((-1, 3), "4.0.144.1"),
                             ((2, 3), "4.4.2304.1"), ((-3, 5), "4.0.225.1")]:
            disc, r1 = disc_and_r1(basis)
            self.assertTrue(label.startswith(f"4.{r1}.{abs(disc)}."), basis)

    @unittest.skipIf(make_field is None, "needs cypari2")
    def test_disc_formula(self):
        for basis in [(-1, 2), (2, 3), (-3, 5), (2, 3, 5), (-1, 2, 3), (-3, 5, -7), (-1, 2, 3, 5)]:
            K = make_field("", multiquadratic_coeffs(basis))
            self.assertEqual((K.disc, K.r1), disc_and_r1(basis), basis)


if __name__ == "__main__":
    unittest.main()
