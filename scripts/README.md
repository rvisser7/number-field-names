# Scripts

Small Python utilities for working with the data in `../data/`. They need only the standard library,
except for the scripts that compute fields, which need PARI through cypari2:
either `pip install cypari2`, or run them with `sage -python` since SageMath ships cypari2.

- `nf_names.py`: shared helpers (reading, writing and merging data files, LMFDB sort order, quadratic names).
  `FAMILIES` sets which name wins when a field is in more than one family.
- `pari_fields.py`: PARI setup, and the `Field` record (name, polredabs polynomial, discriminant, signature).
- `progress.py`: the progress line the long-running scripts print on stderr.
- `pipeline.py`: the shared end of the newer `generate_*.py` scripts (LMFDB lookup, merge, summary).
- `check_data.py`: checks every `data/degree<n>.txt` is well formed: right header, one row per label,
  labels of the right degree, in LMFDB order, no duplicates, and balanced brackets in each name.
- `generate_degree2.py`: writes `data/degree2.txt`, either for all fundamental discriminants up to
  `--bound`, or for the labels listed in a file given with `--labels` (e.g. exported from the LMFDB).
- `pure.py`: finds the distinct pure fields $\mathbb{Q}(\sqrt[n]{m})$ with $|m| \leq B$. Each field is identified by its
  polredabs polynomial and named by the $m$ of smallest $|m|$ that generates it, preferring $m > 0$
  (so $\mathbb{Q}(\sqrt[3]{4})$ and $\mathbb{Q}(\sqrt[3]{-2})$ are both named `\Q(\sqrt[3]{2})`).
- `lmfdb_lookup.py`: finds the LMFDB label of a field, from the fields with the same degree and $|\mathrm{disc}|$,
  by matching the polredabs polynomial, falling back to an isomorphism test (`nfisisom`). Uses either the public
  LMFDB API or a database connection (`from lmfdb import db`, or lmfdb-lite). The database is faster and
  more reliable; to get lmfdb-lite inside Sage:
  `sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git"`.
- `generate_pure.py`: puts the two together and merges the names into `data/degree<n>.txt`, keeping any name
  already there (unless pure names win over it, or `--replace` is given), and lists pure fields that are not in the LMFDB. For example
  `python3 scripts/generate_pure.py --degrees 3-10 --bound 100 --source db`.
- `cyclotomic.py`: the fields $\mathbb{Q}(\zeta_n)$ and $\mathbb{Q}(\zeta_n)^+$ of degree 3 to $D$, for $n \not\equiv 2 \bmod 4$
  (degree 1 and 2 fields keep their rational or quadratic names, as on the LMFDB).
- `generate_cyclotomic.py`: looks these up in the LMFDB and merges them into the data files, e.g.
  `python3 scripts/generate_cyclotomic.py --max-degree 47 --source db`.
- `imprimitive_quartic.py`: the quartic fields $\mathbb{Q}(\sqrt{A + B\sqrt{D}})$ with $|A|, |B|, |D| \leq N$ and $D$ squarefree,
  leaving out the biquadratic (V4) ones. Their quadratic subfield $\mathbb{Q}(\sqrt{D})$ is unique, and $(A, B)$ is chosen
  with $B > 0$ and smallest $\max(|A|, |B|)$, then $|B|$, then $|A|$, then $A > 0$. If that gives $A = 0$ the field is
  named as a pure quartic $\mathbb{Q}(\sqrt[4]{B^2 D})$, the same name `generate_pure.py` gives. $D = -1$ is written $i$,
  as in the LMFDB, e.g. `\Q(\sqrt{1 + 2i})`.
- `generate_imprimitive_quartic.py`: looks these up in the LMFDB and merges them into `data/degree4.txt`, e.g.
  `python3 scripts/generate_imprimitive_quartic.py --bound 50 --source db`.
- `multiquadratic.py`: the fields $\mathbb{Q}(\sqrt{d_1}, \dots, \sqrt{d_k})$, $k \geq 2$, named as in the LMFDB's `field_pretty`:
  the $d_i$ are the greedy basis of the quadratic subfields ordered by $|d|$ (positive first), and $\sqrt{-1}$ is
  written $i$, e.g. `\Q(i, \sqrt{2})`. Runs over the fields whose basis has every $|d_i| \leq N$; the discriminant
  comes from the conductor-discriminant formula, so no polynomial is needed to decide whether to look a field up.
- `generate_multiquadratic.py`: looks these up in the LMFDB and merges them into the data files, e.g.
  `python3 scripts/generate_multiquadratic.py --degrees 4 --bound 200 --source db`. The number of fields grows fast with
  the degree (bound 200, 40, 20, 15 give about 27000, 12000, 3200, 500 fields in degree 4, 8, 16, 32).
- `radical.py`: Kummer fields $\mathbb{Q}(\zeta_n, \sqrt[n]{a})$ and radical composita $\mathbb{Q}(\sqrt[m]{a}, \sqrt[n]{b})$ with
  $\gcd(m, n) = 1$. Both are found by running over small $a$ (and $b$), building the compositum with PARI and grouping
  by polredabs polynomial; a field reached several ways keeps the first name, i.e. the smallest $n$ (or pair), then the
  smallest $|a|$, positive first. $\zeta_n$ is written as for cyclotomic fields ($\zeta_{3}$ for $n = 6$, $i$ for $n = 4$).
  Kummer fields with $\zeta_n$ already in $\mathbb{Q}(\sqrt[n]{a})$ are left out, and where the families overlap the
  simpler name wins: `\Q(\sqrt[6]{-3})` over `\Q(\zeta_{3}, \sqrt[3]{3})` over `\Q(\sqrt{-3}, \sqrt[3]{3})`.
- `generate_kummer.py`, `generate_composita.py`: look these up in the LMFDB and merge them into the data files, e.g.
  `python3 scripts/generate_kummer.py --n 3,4,6 --bound 100 --source db` and
  `python3 scripts/generate_composita.py --pairs 2x3 --bound 50 --source db`.
- `named_composita.py`: composita $K_1 K_2$ of two fields we already name, where $K_1$ is cyclotomic or quadratic
  (so Galois, and the compositum is well defined) and $K_2$ is cyclotomic, quadratic, pure or an imprimitive quartic
  $\mathbb{Q}(\sqrt{A + B\sqrt{D}})$. Only linearly disjoint pairs are kept, named `\Q(g_1, g_2)` from the two generators, e.g.
  `\Q(\zeta_{5}, \sqrt{2})` or `\Q(i, \sqrt[4]{3})`. Pairs another family covers are skipped (two quadratics, two cyclotomics,
  quadratic with odd degree pure), and a field reached several ways keeps the first name (cyclotomic blocks first, smallest
  degree first, then quadratics by $|d|$, then pure, then quartic). Other families' names win over these.
- `generate_named_composita.py`: looks these up in the LMFDB and merges them into the data files, e.g.
  `python3 scripts/generate_named_composita.py --max-degree 16 --bound 20 --source db`.
- `coverage.py`: counts the named fields per degree and per family; with `--db` also the LMFDB fields of each degree,
  to show what fraction has a name.
- `run_jobs.py`: runs the jobs in a jobs file (one `generate_*.py` script and its arguments per line, like `../jobs.txt`)
  as separate processes, at most `--parallel` at a time, each logging to `logs/<n>-<script>.log`, e.g.
  `sage -python scripts/run_jobs.py jobs.txt --parallel 4 --source db`. `--source` is added to each job that does not
  give one, and `--dry-run` only lists the jobs.
