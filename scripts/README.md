# Scripts

Small Python utilities for working with the data in `../data/`. They need only the standard library,
except for the pure field scripts, which need PARI through cypari2 (`pip install cypari2`).

- `nf_names.py`: shared helpers (reading and writing data files, LMFDB sort order, quadratic names).
- `check_data.py`: checks every `data/degree<n>.txt` is well formed: right header, one row per label,
  labels of the right degree, in LMFDB order, no duplicates, and balanced brackets in each name.
- `generate_degree2.py`: writes `data/degree2.txt`, either for all fundamental discriminants up to
  `--bound`, or for the labels listed in a file given with `--labels` (e.g. exported from the LMFDB).
- `pure.py`: finds the distinct pure fields $\mathbb{Q}(\sqrt[n]{m})$ with $|m| \leq B$. Each field is identified by its
  polredabs polynomial and named by the $m$ of smallest $|m|$ that generates it, preferring $m > 0$
  (so $\mathbb{Q}(\sqrt[3]{4})$ and $\mathbb{Q}(\sqrt[3]{-2})$ are both named `\Q(\sqrt[3]{2})`).
- `lmfdb_lookup.py`: finds the LMFDB label of a field, from the fields with the same degree and $|\mathrm{disc}|$,
  by matching the polredabs polynomial, falling back to an isomorphism test (`nfisisom`). Uses either the public
  LMFDB API or a database connection (`from lmfdb import db`, or lmfdb-lite).
- `generate_pure.py`: puts the two together and merges the names into `data/degree<n>.txt`, keeping any name
  already there unless `--replace` is given, and lists pure fields that are not in the LMFDB. For example
  `python3 scripts/generate_pure.py --degrees 3-10 --bound 100 --source db`.
