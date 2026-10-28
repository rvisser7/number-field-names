# Pretty LaTeX names for number fields

This repository computes "pretty" LaTeX names for number fields in the [LMFDB](https://www.lmfdb.org/NumberField/), such as $\mathbb{Q}(\sqrt{-3})$, $\mathbb{Q}(\zeta_{7})$ or $\mathbb{Q}(\sqrt[5]{2})$.

At the moment the LMFDB computes these names on the fly in `field_pretty` (`lmfdb/number_fields/web_number_field.py`), which only covers a few families.
The aim here is to precompute names offline, where expensive checks are affordable, so that they can be stored in a `pretty_name` column of `nf_fields`.
See [LMFDB issue #6991](https://github.com/LMFDB/lmfdb/issues/6991).

## Data

The data lives in [`data/`](data/), with one file per degree: `degree2.txt`, `degree3.txt`, ...
Each file is in psycodict upload format:

```
label|pretty_name
text|text

2.0.3.1|\Q(\sqrt{-3})
2.0.4.1|\Q(\sqrt{-1})
```

Names are stored without the surrounding `\(` `\)` and use the LMFDB's macros (such as `\Q`).
A field only appears in a file if we have a name for it.

## Families

| Family | Name | Status |
|---|---|---|
| Rationals | `\Q` | done |
| Quadratic | `\Q(\sqrt{D})` | read off the label; discriminants up to 1000 in absolute value so far |
| Cyclotomic | `\Q(\zeta_n)` | todo |
| Real cyclotomic | `\Q(\zeta_n)^+` | todo |
| Multiquadratic | `\Q(\sqrt{a}, \sqrt{b}, \dots)` | todo |
| Pure cubic | `\Q(\sqrt[3]{n})` | todo |
| Imprimitive quartic | `\Q(\sqrt{A + B\sqrt{D}})` | todo |
| Pure, degree $m \geq 4$ | `\Q(\sqrt[m]{n})` | todo |
| Composita | e.g. `\Q(\sqrt{n}, \sqrt[3]{m})` | todo |
| Kummer | `\Q(\zeta_n, \sqrt[n]{a_1}, \dots)` | todo |

Every name other than those read off the label should be certified by checking that the
polynomial it describes defines the same field as the LMFDB polynomial (e.g. with PARI's `nfisisom`).

## Scripts

See [`scripts/README.md`](scripts/README.md).

```
python3 scripts/check_data.py          # validate every data/degree<n>.txt
python3 -m unittest discover tests     # run the tests
```
