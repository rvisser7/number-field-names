import io
import os
import sys
import tempfile
import time
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from run_jobs import read_jobs, run_jobs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestRunJobs(unittest.TestCase):
    def test_read_jobs(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "jobs.txt")
            with open(path, "w") as f:
                f.write("# comment\n\ngenerate_pure.py --degrees 5 --bound 10  # pure\n"
                        "generate_pure.py --degrees 6 --bound 10 --source api\n")
            self.assertEqual(read_jobs(path, "db"), [
                ["generate_pure.py", "--degrees", "5", "--bound", "10", "--source", "db"],
                ["generate_pure.py", "--degrees", "6", "--bound", "10", "--source", "api"]])

    def test_the_example_jobs_file(self):
        for job in read_jobs(os.path.join(ROOT, "jobs.txt")):
            self.assertTrue(os.path.isfile(os.path.join(ROOT, "scripts", job[0])), job[0])

    def test_run_in_parallel(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "job.py"), "w") as f:
                f.write("import sys, time\ntime.sleep(1)\nprint('hello', sys.argv[1])\n"
                        "sys.exit(int(sys.argv[1]))\n")
            jobs = [["job.py", "0"], ["job.py", "3"], ["job.py", "0"]]
            out = io.StringIO()
            t = time.time()
            codes = run_jobs(jobs, 3, os.path.join(d, "logs"), d, out)
            self.assertLess(time.time() - t, 2.5)  # all three at once, not one after another
            self.assertEqual(codes, [0, 3, 0])
            self.assertEqual(sorted(os.listdir(os.path.join(d, "logs"))),
                             ["01-job.log", "02-job.log", "03-job.log"])
            with open(os.path.join(d, "logs", "02-job.log")) as f:
                self.assertEqual(f.read(), "hello 3\n")
            self.assertIn("3 jobs finished, 1 failed: 2", out.getvalue())


if __name__ == "__main__":
    unittest.main()
