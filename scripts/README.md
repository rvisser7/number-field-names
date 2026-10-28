# Scripts

Small Python utilities for working with the data in `../data/`. They need only the standard library.

- `nf_names.py`: shared helpers (reading and writing data files, label sort order, quadratic names).
- `check_data.py`: checks every `data/degree<n>.txt` is well formed: right header, one row per label,
  labels of the right degree, sorted, no duplicates, and balanced brackets in each name.
- `generate_degree2.py`: writes `data/degree2.txt`, either for all fundamental discriminants up to
  `--bound`, or for the labels listed in a file given with `--labels` (e.g. exported from the LMFDB).
