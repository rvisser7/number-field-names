#!/usr/bin/env python3
"""Run a list of generate_*.py jobs in parallel, each in its own process.

The jobs file has one job per line: a script in scripts/ and its arguments, e.g.

    generate_pure.py --degrees 5 --bound 100
    generate_multiquadratic.py --degrees 8 --bound 40

Blank lines and lines starting with # are skipped. At most --parallel jobs run at once, and the
others wait for a free slot. Each job writes its output (progress, fields not in the LMFDB,
errors) to logs/<n>-<script>.log, and a line is printed here when it starts and finishes.
Merges into the data files are locked (see merge_data in nf_names.py), so jobs writing to the
same degree are safe.

Jobs run with the same Python as this script, so run it with sage -python if that is how you
get cypari2. --source is added to every job that does not give its own.

Usage:
    sage -python scripts/run_jobs.py jobs.txt --parallel 4 --source db
    sage -python scripts/run_jobs.py jobs.txt --dry-run
"""

import argparse
import os
import shlex
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from progress import fmt_time

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(os.path.dirname(SCRIPTS_DIR), "logs")


def read_jobs(path, source=None):
    """The jobs in a jobs file, as argument lists [script, arg, ...]."""
    jobs = []
    with open(path) as f:
        for line in f:
            args = shlex.split(line, comments=True)
            if not args:
                continue
            if source and not any(a == "--source" or a.startswith("--source=") for a in args):
                args += ["--source", source]
            jobs.append(args)
    return jobs


def log_name(n, args):
    return "%02d-%s.log" % (n, os.path.splitext(os.path.basename(args[0]))[0])


def run_jobs(jobs, parallel, log_dir=LOG_DIR, scripts_dir=SCRIPTS_DIR, out=sys.stdout):
    """Run the jobs, at most `parallel` at a time. Returns the exit code of each job."""
    for args in jobs:
        if not os.path.isfile(os.path.join(scripts_dir, args[0])):
            raise SystemExit(f"no such script: {os.path.join(scripts_dir, args[0])}")
    os.makedirs(log_dir, exist_ok=True)
    waiting = list(enumerate(jobs, start=1))
    running = {}  # n -> (process, log file, start time)
    codes = {}
    t0 = time.time()

    def say(msg):
        print(f"[{fmt_time(time.time() - t0):>9}] {msg}", file=out, flush=True)

    try:
        while waiting or running:
            while waiting and len(running) < parallel:
                n, args = waiting.pop(0)
                log = open(os.path.join(log_dir, log_name(n, args)), "w")
                cmd = [sys.executable, "-u", os.path.join(scripts_dir, args[0])] + args[1:]
                running[n] = (subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                               stdin=subprocess.DEVNULL), log, time.time())
                say(f"started  {n}/{len(jobs)}: {shlex.join(args)}  (log: {log.name})")
            time.sleep(0.2)
            for n, (proc, log, start) in list(running.items()):
                if proc.poll() is None:
                    continue
                log.close()
                del running[n]
                codes[n] = proc.returncode
                status = "done" if proc.returncode == 0 else f"FAILED (exit {proc.returncode})"
                say(f"{status:8} {n}/{len(jobs)}: {shlex.join(jobs[n - 1])}  "
                    f"in {fmt_time(time.time() - start)}")
    except KeyboardInterrupt:
        for proc, log, _ in running.values():
            proc.terminate()
        for proc, log, _ in running.values():
            proc.wait()
            log.close()
        say(f"interrupted: stopped {len(running)} running job(s), {len(waiting)} not started")
        raise SystemExit(130)
    failed = [n for n in sorted(codes) if codes[n]]
    say(f"all {len(jobs)} jobs finished, {len(failed)} failed"
        + (": " + ", ".join(map(str, failed)) + f" (see the logs in {log_dir})" if failed else ""))
    return [codes[n] for n in sorted(codes)]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("jobs", help="the jobs file")
    parser.add_argument("-p", "--parallel", type=int, default=4,
                        help="how many jobs to run at once (default 4)")
    parser.add_argument("--source", choices=["api", "db"],
                        help="add --source to every job that does not give one")
    parser.add_argument("--log-dir", default=LOG_DIR, help="where the logs go (default logs/)")
    parser.add_argument("--dry-run", action="store_true", help="only list the jobs")
    args = parser.parse_args()

    jobs = read_jobs(args.jobs, args.source)
    if args.dry_run:
        for n, job in enumerate(jobs, start=1):
            print(f"{n}: {shlex.join(job)}  (log: {log_name(n, job)})")
        return
    codes = run_jobs(jobs, args.parallel, args.log_dir)
    sys.exit(1 if any(codes) else 0)


if __name__ == "__main__":
    main()
