# Wannier90 v3.1.0 — Test Suite & Validation Targets

Reference tree: `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90`
All paths below are relative to that root unless absolute.

This note maps the ctest/testcode validation harness, extracts the exact
numerical tolerances and comparison rule, inventories the `testw90_*` cases and
the shipped tutorials, and picks the M0/M1/M2 target for the Julia
reimplementation (isolated bands, no disentanglement, all inputs + reference
numbers shipped).

---

## 1. How the harness works

The `wannier90.x` tests use **testcode2** (vendored under
`test-suite/testcode/`), driven by two config files:

- `test-suite/tests/jobconfig` — one INI `[section]` per test case:
  which `program` (extraction+tolerance profile) to run, the input file(s) and
  CLI args (`inputs_args`), and the primary `output` file to parse.
- `test-suite/tests/userconfig` — defines each `program`: the executable
  (`exe`), the Python `extract_fn` parser, and the **tolerances**.

Runner: `test-suite/run_tests` (or `ctest` via `test-suite/CMakeLists.txt`).
`clean_tests` wipes generated files.

### 1.1 Where the reference numbers live

There is **no separate `benchmark/` directory**. Each test folder contains a
file literally named:

```
benchmark.out.default.inp=<winfile>
```

e.g. `test-suite/tests/testw90_example01/benchmark.out.default.inp=gaas.win`.
This file **is a full reference `.wout`** (the golden run output). The parser is
run on *both* this benchmark file and the freshly-produced `output`, and the
extracted dictionaries are compared field-by-field. The `default` token is the
`benchmark` name set in `userconfig` `[user] benchmark = default`.

### 1.2 The extraction parsers

`extract_fn = tools parsers.parse_wout.parse` → `test-suite/tools/parsers/parse_wout.py`.
For the main W90 minimizer tests (`WANNIER90_WOUT_OK`) it regex-scrapes the
`.wout` and returns lists keyed by label. Key regexes
(`parse_wout.py`):

- WF centre/spread line (used for BOTH `Initial State` and `Final State`
  blocks; only the block after `"Final State"` is stored as
  `final_centres_{x,y,z}` / `final_spreads`):
  `^\s*WF centre and spread\s+(\d+)\s+\(\s*([0-9\.-]+)\s*,\s*([0-9\.-]+)\s*,\s*([0-9\.-]+)\s*\)\s*([0-9\.-]+)\s*$`
- `Omega I  = ...`, `Omega D = ...`, `Omega OD = ...`, `Omega Total = ...`
  → keys `omegaI`, `omegaD`, `omegaOD`, `omegaTotal`.
- Nearest-neighbour shell table → `near_neigh_dist`, `near_neigh_mult`.
- b_k vectors/weights and directions, completeness relation → `completeness_*`.

**Gotcha:** the centre/spread regex only accepts characters `[0-9.-]`, i.e. it
does **not** match exponential notation (`E`/`D`) or `NaN`. Reference outputs
are fixed-decimal (`f`-format), so our reimplementation's `.wout` must print
centres/spreads in the same `F` format for the harness to parse them, but for
*our own* numeric comparison we simply compare the floating values.

### 1.3 The comparison rule (exact)

`test-suite/testcode/lib/testcode2/validation.py`:

- `validate_absolute`: `err = abs(test - benchmark); passed = err < absolute`
  (strict `<`, not `<=`).
- `validate_relative`: `err = abs((test-benchmark)/benchmark); passed = err < relative`.
  If `benchmark == 0`, relative check is skipped (treated as passing).
- Default `strict = True` ⇒ a value passes only if it is **within BOTH the
  absolute AND the relative tolerance** (`status_relative + status_absolute`,
  logical AND). When a tolerance tuple sets one of the two to `None`, that check
  is skipped.

Tolerance tuple format in `userconfig` is **`(absolute, relative, 'label')`**
(absolute first).

### 1.4 The tolerances that matter for the minimizer (`WANNIER90_WOUT_OK`)

From `test-suite/tests/userconfig`, section `[WANNIER90_WOUT_OK]`
(`abs`, `rel`):

| label                | absolute | relative |
|----------------------|----------|----------|
| `final_centres_x/y/z`| 1.0e-5   | 1.0e-5   |
| `final_spreads`      | 3.0e-6   | 3.0e-6   |
| `omegaI`             | 1.0e-6   | 1.0e-6   |
| `omegaD`             | 1.0e-6   | 5.0e-6   |
| `omegaOD`            | 1.0e-6   | 1.0e-6   |
| `omegaTotal`         | 1.0e-6   | 1.0e-6   |
| `near_neigh_dist`    | 1.0e-6   | 1.0e-6   |
| `near_neigh_mult`    | 1.0e-6   | 1.0e-6   |
| `completeness_{x,y,z}`| 1.0e-6  | 1.0e-6   |
| `completeness_weight`| 1.0e-6   | 1.0e-6   |

So for M0/M1/M2 we must reproduce **Omega I/D/OD/Total to ~1e-6 absolute** and
**Wannier centres to ~1e-5** — a genuinely tight numerical match, not a loose
smoke test. (Selectively-localized `SLWFC` and `SAWFS` variants use slightly
different Omega labels: `omegaIOD_C`, `omegaRest`, `penaltyfunc`,
`omegaTotal_C`; irrelevant for plain isolated-band targets.)

Other program profiles (not minimizer targets): `WANNIER90_NNKP_OK`
(b-vectors only, 1e-6), `WANNIER90_BVEC` (`.bvec` file, 1e-10),
`WANNIER90_RMN_OK` (r-matrix, abs 2e-6 / rel 1.0), `WANNIER90_BANDS_PLOT`,
`WANNIER90_LABELINFO`, `WANNIER90_WERR_FAIL` (expected-crash), and the
`POSTW90_*` post-processing profiles.

---

## 2. Inventory of W90 minimizer cases + tutorials

"Isolated" = `num_bands == num_wann` and no disentanglement. When `num_bands`
is absent from the `.win`, it defaults to `num_wann` ⇒ isolated. "Benchmark?"
= whether a `benchmark.out.default.inp=*` golden `.wout` is shipped in the test
folder. Tutorials ship **inputs only, no benchmark output**.

### 2.1 `testw90_*` main-code cases (subset that runs the minimizer / relevant)

| Case (dir under `test-suite/tests/`) | system | num_wann | num_bands | disentangle? | k-mesh | inputs shipped | benchmark? |
|---|---|---|---|---|---|---|---|
| `testw90_example01` | GaAs | 4 | 4 | no | 2×2×2 full | `gaas.amn`,`gaas.mmn`,`gaas.win`,`UNK*` (no `.eig`) | yes |
| `testw90_example02` | Pb (lead) | 4 | 4 | no | 4×4×4 full | `lead.amn/.mmn/.eig/.win` | yes |
| `testw90_example02_restart` | Pb | 4 | 4 | no | 4×4×4 | same as ex02 | yes |
| `testw90_example05` | diamond C | 4 | 4 | no | 4×4×4 full | `diamond.amn/.mmn/.eig/.win` | yes |
| `testw90_example07` | silane SiH4 | 4 | 4 | no | 1×1×1 **Γ-only** | `silane.amn/.mmn/.eig/.win` | yes |
| `testw90_na_chain_gamma` | Na chain | 10 | 30 | **yes** | 1×1×1 Γ | `Na_chain.*` | yes |
| `testw90_benzene_gamma_val` | benzene (valence) | 15 | 15 | no | 1×1×1 **Γ-only** | `benzene.amn/.mmn/.eig/.win` | yes |
| `testw90_benzene_gamma_valcond` | benzene val+cond | — | — | yes (dis) | Γ | shipped | yes |
| `testw90_basic1` | GaAs-like | 4 | 4 (dflt) | no* | 2×2×2 | `wannier.amn/.mmn/.win` (no `.eig`) | yes |
| `testw90_basic2` | GaAs-like | 4 | 8 | **yes** | 4×4×4 | `wannier.amn/.mmn/.eig/.win` | yes |
| `testw90_bvec` | Pb | 4 | 4 | n/a | 4×4×4 | `lead.amn/.mmn/.eig/.win` | yes (**`.bvec` only**, not spread) |
| `testw90_example03*` | Si val+cond | 8 | 12 | **yes** | 4×4×4 | shipped | band/kpt/labelinfo checks |
| `testw90_example04` | Cu | 7 | 12 | **yes** | 4×4×4 | shipped | yes |
| `testw90_example11_1/2`, `_36` | Si | 4/8 | varies | some | full | shipped | yes |
| `testw90_example21_As_sp` | GaAs SAWF | — | — | sym-adapted | full | shipped | yes |
| `testw90_example26` | GaAs SLWF | — | — | selective-loc | full | shipped | yes (SLWFC profile) |
| `testw90_disentanglement_sawfs` | H3S | — | — | dis+SAWF | full | shipped | yes (only `omegaI`) |
| `testw90_knbo3_higher*` | KNbO3 | — | — | higher-order FD | full | shipped | yes |

\* `basic1` lists `dis_*` keywords but ships **no `.eig`** and `num_bands`
defaults to `num_wann=4`; the benchmark shows no `Disentanglement` block, so it
runs as an isolated-band case (it is deliberately a "messy input formatting"
parser test: comments after values, `!` inline comments, `Bohr` units,
`num_print_cycles=13`). It is NOT a clean physics target.

### 2.2 Tutorials (ship inputs; NO benchmark output — not standalone targets)

| dir | system | num_wann | num_bands | disentangle? | k-mesh | inputs |
|---|---|---|---|---|---|---|
| `tutorials/tutorial03` | silicon | 8 | 12 | **yes** (`dis_win_max=17`, `dis_froz_max=6.4`) | 4×4×4 | `silicon.amn/.mmn/.eig/.win` |
| `tutorials/tutorial04` | copper | 7 | 12 | **yes** (`dis_win_max=38`, `dis_froz_max=13`) | 4×4×4 | `copper.amn/.mmn/.eig/.win` |

Because tutorials ship no golden `.wout`, they are useful as *additional inputs*
but the reference numbers to match live in the `test-suite` mirror cases
(e.g. `example03/example04` mirror the silicon/copper tutorials with a benchmark).

### 2.3 Library-mode reference

`test-suite/library-mode-test/` ships `gaas.amn/.mmn/.eig` and a golden
`ref/gaas_ref.wout` — but it is **8 WF / 12 bands with disentanglement**
(`Number of Wannier Functions: 8`, `input Bloch states: 12`), and the `.win`
is *generated at runtime by `demo.f90`* (no static `gaas.win`). Final numbers:
Omega I = 13.567803346, Omega D = 0.314383305, Omega OD = 5.186173138,
Omega Total = 19.068359789. Not an isolated target; keep for a later
disentanglement milestone.

---

## 3. Ranked shortlist — ISOLATED-BANDS targets (all inputs + ref numbers)

Ranking criteria: (a) all inputs + benchmark shipped, (b) isolated
(`num_bands==num_wann`, no dis), (c) simplest first — fewest WFs, then Γ-only
before full mesh, then a full-k case that exercises the k→b machinery.

1. **`testw90_example07` (silane, Γ-only, 4 WF)** — SIMPLEST. 1×1×1 mesh,
   `gamma_only=true`, isolated 4 valence WFs. Ships `.amn/.mmn/.eig/.win` + a
   benchmark. Good **M0** smoke target: at Γ the b-vectors/overlaps are minimal
   and the Γ-only real-gauge code path is exercised. *Caveat:* the Γ-only branch
   in W90 uses a special real-valued algorithm (`wann_main_gamma`), a different
   code path from the general complex minimizer. If you implement the general
   complex path first, Γ-only is NOT the same code and may mismatch — so it is a
   good *first-light* target but not a validator of the general path.

2. **`testw90_example01` (GaAs, 2×2×2 full mesh, 4 WF)** — BEST M0/M1/M2
   general-path target. Small full k-mesh (8 k-points), isolated 4 valence WFs
   (`As:sp3` projections), converges essentially in ONE iteration, and all four
   Omega components are non-trivial (nonzero `Omega D`), exercising the full
   k→b finite-difference + gauge-gradient machinery. Ships
   `gaas.amn`,`gaas.mmn`,`gaas.win`. **Gotcha:** ships NO `gaas.eig` (the
   `UNK*.1` wavefunctions are shipped instead for plotting; the eig-dependent
   features are off). If your reimplementation requires eigenvalues, note this
   case runs without them (`Omega D` still computed from mmn/amn gauge).

3. **`testw90_example05` (diamond C, 4×4×4 full mesh, 4 WF)** — clean
   validator. Isolated 4 WFs, ships `.amn/.mmn/.eig/.win` + benchmark. Larger
   mesh (64 k-points), `Omega D == 0` exactly by symmetry (nice invariant to
   check), converges fully (20 iters, Delta ~1e-15). Good **M2** convergence
   test with `.eig` present.

4. **`testw90_example02` (Pb, 4×4×4, 4 WF)** — isolated 4/4, full mesh, ships
   all inputs + benchmark. Metallic ⇒ slightly harder convergence; secondary.

5. **`testw90_benzene_gamma_val` (benzene, Γ-only, 15 WF)** — isolated but 15
   WFs and Γ-only path; more moving parts. Use later.

**Decision:** M0 = example07 (Γ first-light) *and/or* example01 (general path,
1 iter). M1/M2 = example01 → example05 (full convergence with `.eig`).
Primary numeric target = **`testw90_example01` (GaAs)**.

---

## 4. Reference numbers to match — top targets

### 4.1 TARGET #1 — `testw90_example01` (GaAs, 4 WF, 2×2×2)

Files:
- Inputs: `test-suite/tests/testw90_example01/{gaas.win,gaas.amn,gaas.mmn}`
  (+ `UNK0000{1..8}.1`; **no `gaas.eig`**).
- Golden: `test-suite/tests/testw90_example01/benchmark.out.default.inp=gaas.win`

Key `.win` settings: `num_wann=4`, `num_iter=20`, `use_ws_distance=.false.`,
`search_shells=12`, projections `As:sp3`, unit cell in **bohr** (FCC, a=5.367),
`mp_grid: 2 2 2`, `wvfn_formatted=.true.`.

**Initial State** (benchmark lines 221–226):
```
WF 1 ( -0.866632,  1.973462,  1.973462 )  1.11720303
WF 2 ( -0.866632,  0.866632,  0.866632 )  1.11720303
WF 3 ( -1.973462,  1.973462,  0.866632 )  1.11720303
WF 4 ( -1.973462,  0.866632,  1.973462 )  1.11720303
Sum   ( -5.680188,  5.680188,  5.680188 )  4.46881212   (initial Omega_Total)
Iter 0:  O_D = 0.0083198  O_OD = 0.5036294  O_TOT = 4.4688121
```

**Final State** (benchmark lines 451–461), the numbers our M-milestones must
reproduce to the §1.4 tolerances:
```
WF 1 ( -0.866253,  1.973841,  1.973841 )  1.11672024
WF 2 ( -0.866253,  0.866253,  0.866253 )  1.11672024
WF 3 ( -1.973841,  1.973841,  0.866253 )  1.11672024
WF 4 ( -1.973841,  0.866253,  1.973841 )  1.11672024
Sum   ( -5.680188,  5.680188,  5.680188 )  4.46688098

Omega I     = 3.956862958
Omega D     = 0.008030049
Omega OD    = 0.501987969
Omega Total = 4.466880976
```
Centres/spreads are in **Ångström / Ångström²** (`.wout` convention), even
though the `.win` cell is in bohr. Converges after 1 cycle (num_iter=20 but
already at minimum).

### 4.2 TARGET #2 — `testw90_example05` (diamond, 4 WF, 4×4×4)

Files: `test-suite/tests/testw90_example05/{diamond.win,diamond.amn,diamond.mmn,diamond.eig}`
+ `benchmark.out.default.inp=diamond.win`.
`num_wann=4`, `num_iter=20`, projections = 4 explicit `s` centres, cell in Å
(a-vectors ±1.613990), `mp_grid: 4 4 4`.

**Initial State** (lines 229–234):
```
WF 1 ( -0.000000, -0.000000,  0.000000 )  0.58137959
WF 2 ( -0.806995,  0.806995, -0.000000 )  0.58137959
WF 3 (  0.000000,  0.806995,  0.806995 )  0.58137959
WF 4 ( -0.806995,  0.000000,  0.806995 )  0.58137959
Iter 0:  O_TOT = 2.3255183666
```
**Final State** (lines 459–469):
```
WF 1 ( -0.000000, -0.000000,  0.000000 )  0.58022623
WF 2 ( -0.806995,  0.806995, -0.000000 )  0.58022623
WF 3 (  0.000000,  0.806995,  0.806995 )  0.58022623
WF 4 ( -0.806995,  0.000000,  0.806995 )  0.58022623
Sum   ( -1.613990,  1.613990,  1.613990 )  2.32090491

Omega I     = 1.954619860
Omega D     = 0.000000000   (exactly zero by symmetry — good invariant)
Omega OD    = 0.366285055
Omega Total = 2.320904915
```
Fully converged (iter 20, Delta O_TOT ≈ -1.3e-15).

### 4.3 (Reference-only) `testw90_example07` silane Γ

4 WF, Γ-only; ships `.amn/.mmn/.eig/.win` + benchmark
(`benchmark.out.default.inp=silane.win`). Number of Wannier functions = 4,
input Bloch states = 4. Uses the **Γ-only real-gauge** minimizer path (distinct
from the general complex path) — first-light only.

---

## 5. Gotchas that would cause a silent numerical mismatch

- **Units in `.wout`:** WF centres and spreads are always printed in
  **Ångström / Ångström²**, regardless of whether the `.win` cell/atoms are in
  `bohr` (example01) or `ang`/`Å` (example05). The internal minimizer works in
  Å⁻¹ b-vectors (see the `.wout` header "b_k Vectors (Ang^-1) and Weights
  (Ang^2)"). Get the Bohr↔Å conversion (`bohr = 0.52917721...`) and its exact
  value matching W90's `constants` module, or centres drift past the 1e-5 tol.
- **Comparison is AND of abs & rel (strict `<`).** Omega components must match
  to ~1e-6 absolute; this is a *precision* test, not a smoke test.
- **`Omega D` sensitivity:** its relative tol is looser (5e-6) precisely because
  it is the small, gauge-dependent, numerically delicate piece
  (`Omega D = Σ N |Im ln M + b·r̄|²` type term). Branch cut of `Im ln` of the
  diagonal overlaps and the choice `r̄ = -(1/N) Σ w_b b Im ln M_kb` must match
  W90 exactly; a wrong branch flips O_D silently while O_I/O_Total look fine.
- **b-vector shell selection / weights** must reproduce W90's finite-difference
  shells (`search_shells=12` in these inputs) and the `w_b` least-squares
  weights, else O_I and centres shift. `use_ws_distance=.false.` in the target
  inputs — do NOT apply Wigner-Seitz distance remapping for these.
- **example01 ships no `.eig`** — do not assume eigenvalues are always present.
- **Γ-only cases (07, benzene) use a different algorithm** than the general
  complex minimizer; don't validate the general path against them.
- **Parser format constraint:** `parse_wout.py` centre/spread regex rejects
  exponential notation; emit `F`-format centres/spreads in the `.wout` for the
  ctest harness (our own Julia numeric comparison can bypass this).
- **Tolerance tuple order is `(absolute, relative, label)`** — absolute first.
