import os
import sys
import tempfile
import unittest
from multiprocessing import Process

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from nf_names import merge_data, read_data


def worker(path, start):
    # Many small merges, each adding new labels, from several processes at once.
    for i in range(start, start + 200):
        merge_data(path, {f"6.0.{i}.1": r"\Q(\sqrt[6]{%d})" % i})


class TestParallelMerge(unittest.TestCase):
    def test_no_rows_lost(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "degree6.txt")
            procs = [Process(target=worker, args=(path, 1000 * k)) for k in range(1, 5)]
            for p in procs:
                p.start()
            for p in procs:
                p.join()
            self.assertEqual(len(read_data(path)), 800)
            self.assertEqual(sorted(os.listdir(d)), [".merge.lock", "degree6.txt"])


if __name__ == "__main__":
    unittest.main()
