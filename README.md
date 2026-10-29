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
When a field has names from more than one family, the family higher up in the table below wins
(as in the LMFDB's `field_pretty`), e.g. `4.0.256.1` is `\Q(\zeta_{8})` rather than `\Q(\sqrt[4]{-1})`.

## Families

| Family | Name | Status |
|---|---|---|
| Rationals | `\Q` | done |
| Quadratic | `\Q(\sqrt{D})` | read off the label; discriminants up to 1000 in absolute value so far |
| Cyclotomic | `\Q(\zeta_{n})` | script ready (`generate_cyclotomic.py`), degree 3 to 47 |
| Real cyclotomic | `\Q(\zeta_{n})^+` | script ready (`generate_cyclotomic.py`), degree 3 to 47 |
| Multiquadratic | `\Q(\sqrt{a}, \sqrt{b}, \dots)` | script ready (`generate_multiquadratic.py`), degree 4 to 32 |
| Pure, degree $n \geq 3$ | `\Q(\sqrt[n]{m})` | script ready (`generate_pure.py`) |
| Imprimitive quartic (C4, D4) | `\Q(\sqrt{A + B\sqrt{D}})` | script ready (`generate_imprimitive_quartic.py`) |
| Kummer | `\Q(\zeta_{n}, \sqrt[n]{a})`, e.g. `\Q(\zeta_{3}, \sqrt[3]{2})`, `\Q(i, \sqrt[4]{2})` | script ready (`generate_kummer.py`) |
| Radical composita | `\Q(\sqrt[m]{a}, \sqrt[n]{b})`, $\gcd(m, n) = 1$, e.g. `\Q(\sqrt{2}, \sqrt[3]{3})` | script ready (`generate_composita.py`) |

Every name other than those read off the label should be certified by checking that the
polynomial it describes defines the same field as the LMFDB polynomial (e.g. with PARI's `nfisisom`).

## Scripts

See [`scripts/README.md`](scripts/README.md).

```
python3 scripts/check_data.py          # validate every data/degree<n>.txt
python3 -m unittest discover tests     # run the tests (most tests need cypari2)
```
