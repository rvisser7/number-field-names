# Scripts

Small Python utilities for working with the data in `../data/`. They need only the standard library,
except for the scripts that compute fields, which need PARI through cypari2:
either `pip install cypari2`, or run them with `sage -python` since SageMath ships cypari2.

- `nf_names.py`: shared helpers (reading, writing and merging data files, LMFDB sort order, quadratic names).
  `FAMILIES` lists the families, and `name_key` orders the names of one field (best first): by family, except
  that among the radical families the smallest radicands come first, then by the simplest parameters.
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
  left out as pure, as are pure fields that only have $A = 0$ beyond the bound (`generate_pure.py` names them). $D = -1$ is written $i$,
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
- `nested_radical.py`: the fields $\mathbb{Q}(\sqrt[n]{a + b\sqrt{d}})$ of degree $2n$, $n \geq 3$, with $|a|, |b|, |d| \leq N$ and $d$
  squarefree, the analogue of the imprimitive quartics for higher $n$. Each field gets its simplest $(a, b, d)$: $b > 0$, no
  $p^n$ dividing both $a$ and $b$, smallest $\max(|a|, |b|)$, then $|b|$, $|a|$, $|d|$, then $a > 0$, $d > 0$. Pure fields (those
  that $a = 0$ also gives, and the likes of $\mathbb{Q}(\sqrt[3]{5 + 5i}) = \mathbb{Q}(\sqrt[6]{-400})$) are left to `generate_pure.py`;
  $d = -1$ is written $i$.
- `generate_nested_radical.py`: looks these up in the LMFDB and merges them into the data files, e.g.
  `python3 scripts/generate_nested_radical.py --n 3 --bound 20 --source db` (about 13 s of PARI for $n = 3$, bound 20).
- `nested_biquadratic.py`: the octic fields $\mathbb{Q}(\sqrt{A + B\sqrt{d_1} + C\sqrt{d_2} + D\sqrt{d_3}})$, quadratic
  extensions of a biquadratic field $L = \mathbb{Q}(\sqrt{d_1}, \sqrt{d_2})$ (every octic with a biquadratic subfield is
  one), with $|A|, |B|, |C|, |D| \leq N$ and $|d_1|, |d_2| \leq N'$. $d_1, d_2$ is the greedy basis of $L$ and $d_3$ the
  squarefree part of $d_1 d_2$; at least two of $B, C, D$ are nonzero. Each field gets its simplest parameters: smallest
  $\max(|A|, |B|, |C|, |D|)$, then $|D|$, $|C|$, $|B|$, $|A|$, then the smallest basis, then $B, C \geq 0$ (conjugating flips
  the signs of two of $B, C, D$), $A \geq 0$, $D \geq 0$. Every other family's name wins over these.
- `generate_nested_biquadratic.py`: looks these up in the LMFDB and merges them into `data/degree8.txt`, e.g.
  `python3 scripts/generate_nested_biquadratic.py --bound 3 --basis-bound 10 --source db` (about 2 s per
  biquadratic field for bound 3).
- `nested_pure.py`: the fields $\mathbb{Q}(\sqrt{A_0 + A_1\sqrt[n]{m_1} + \dots + A_{n-1}\sqrt[n]{m_{n-1}}})$ of degree $2n$, $n$ an
  odd prime: quadratic extensions of a pure field $F = \mathbb{Q}(\sqrt[n]{m})$ (every field of degree $2n$ with a pure
  subfield of degree $n$ is one), with $|A_k| \leq N$ and $m \leq N'$. $m$ is $F$'s pure radicand, and $m_k$ is the $n$-th power
  free part of $m^k$, so for $m = 12$ the basis is $1, \sqrt[3]{12}, \sqrt[3]{18}$. $\beta$ is not rational and has at least
  two terms (a single term gives a pure field, e.g. $\sqrt{\sqrt[3]{2}} = \sqrt[6]{2}$). Each field gets its simplest
  coefficients: smallest $\max |A_k|$, then $|A_{n-1}|, \dots, |A_0|$, then positive signs. Every other family's name wins.
- `generate_nested_pure.py`: looks these up in the LMFDB and merges them into the data files, e.g.
  `python3 scripts/generate_nested_pure.py --n 3 --bound 5 --m-bound 100 --source db` (about 1 ms per field tried,
  $(2N + 1)^n$ for each $m$).
- `nested_quartic.py`: the octic fields $\mathbb{Q}(\sqrt{a + b\sqrt{D} + c\gamma + d\sqrt{D}\gamma})$, $\gamma = \sqrt{A + B\sqrt{D}}$:
  quadratic extensions of an imprimitive quartic $F = \mathbb{Q}(\gamma)$ (every octic with a C4 or D4 quartic subfield is
  one), with $|a|, |b|, |c|, |d| \leq N$ and $F$ from `imprimitive_quartic_fields(N')`. $(A, B, D)$ is $F$'s own name, $c$ or $d$
  is nonzero, and there are at least two terms. Each field gets its simplest coefficients: smallest
  $\max(|a|, |b|, |c|, |d|)$, then $|d|, |c|, |b|, |a|$, then the simplest quartic, then $c \geq 0$ ($\gamma \to -\gamma$ flips $c$
  and $d$), $d, b, a \geq 0$. Every other family's name wins, also the octics over a biquadratic field.
- `generate_nested_quartic.py`: looks these up in the LMFDB and merges them into `data/degree8.txt`, e.g.
  `python3 scripts/generate_nested_quartic.py --bound 2 --quartic-bound 5 --source db` (about 0.5 s per quartic
  field for bound 2).
- `coverage.py`: counts the named fields per degree and per family; with `--db` also the LMFDB fields of each degree,
  to show what fraction has a name. With `--by-galois --degrees 3-8` it goes through the LMFDB fields instead and
  lists, per Galois group, how many fields there are and how many are named, most unnamed first, with the smallest
  unnamed fields as examples (`--max-disc 10^8` for a quicker look, `--top` and `--examples` for more rows). This
  shows which family to add next.
- `run_jobs.py`: runs the jobs in a jobs file (one `generate_*.py` script and its arguments per line, like `../jobs.txt`)
  as separate processes, at most `--parallel` at a time, each logging to `logs/<n>-<script>.log`, e.g.
  `sage -python scripts/run_jobs.py jobs.txt --parallel 4 --source db`. `--source` is added to each job that does not
  give one, and `--dry-run` only lists the jobs.
- `recognise.py`: the backward direction: given a field (polynomial and ramified primes), test each family and
  return the certified names, and the best one. Pure and Kummer fields: the radicand is a product of ramified
  primes with exponents below $n$ (fixed by the discriminant for tame primes), sieved with small primes and
  checked with `nfroots`. Imprimitive quartics: from $K = L(\sqrt{\beta})$, the simplest $A + B\sqrt{D}$ in
  $\beta L^{*2}$ is found by running over the ideals of small norm in the right class and the units of $L$.
  Nested radicals: for each quadratic subfield $L$, the $\beta$ with $K = L(\sqrt[n]{\beta})$ are $S$-units times $n$-th
  powers ($S$ the primes above the ramified primes and $n$), cut down by $n$-th power conditions at small primes and
  checked with `nffactor` and `nfroots`, then simplified as for the quartics. Fields with more than 20000 candidates
  (e.g. $\mathbb{Q}(\zeta_5, \sqrt{2}, \sqrt{3})$, which has a compositum name) are skipped.
  Octics over a biquadratic field $L$: $K = L(\gamma)$ with $\gamma^2 = \beta \in L$ (from the automorphism of $K/L$), and the
  simplest $A + B\sqrt{d_1} + C\sqrt{d_2} + D\sqrt{d_3}$ in $\beta L^{*2}$ is found as for the quartics, over the ideals of $L$
  of small norm (at most 5000) and its units; only run when $K$ has no other name. Quadratic extensions of a pure
  field $F$, and of an imprimitive quartic $F$: the same, over the ideals and units of $F$ (and its automorphisms).
  Multiquadratic: the quadratic subfields. Composita: pairs of named subfields, and a Galois named subfield
  (quadratic, cyclotomic, multiquadratic or Kummer) with another named subfield, giving names like
  $\mathbb{Q}(i, \sqrt[3]{3}, \sqrt[4]{2})$; fewer generators win.
- `name_lmfdb_fields.py`: runs `recognise.py` over the LMFDB fields of the given degrees and merges the best names,
  every `--batch` fields. Fields come from the database (`--max-disc`, `--min-disc`, `--galois` narrow the query)
  or from a `label|coeffs` file (`--input`); `--shard i/k` splits a degree between jobs, `--families` picks the tests,
  `--only-unnamed` skips fields that already have a name, and `--print` prints the names instead of merging.
  It takes about 1 ms per field in degree 3 to 5, 2 to 10 ms in degree 6 to 10 and 20 to 40 ms in degree 12,
  mostly computing subfields for composita and nested radicals (leave those out with `--families` for a faster pass). E.g.
  `sage -python scripts/name_lmfdb_fields.py --degrees 6 --max-disc 10^8 --shard 0/4`.
