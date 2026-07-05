# Guiding centres + select_projections (branch-cut control in wannierise)

Reference: Wannier90 v3.1.0. Line numbers are for `src/wannierise.F90` unless another
file is named. This note is a companion to `localization.md` (the MV gauge optimizer) and
`disentanglement.md`. It covers **only** the standard non-gamma path (`wann_main`, line 61).
The gamma-only twin (`wann_main_gamma`, line 2836+) is out of scope; its one extra rule is
noted in §6 (trap G).

**One-line summary.** Guiding centres do **not** change the spread/gradient formulas. They
change **only the branch of the Im-ln** used to define the diagonal-overlap phase
`q_n^{k,b}`, by setting a per-(n,nn,k) real offset `sheet(n,nn,k) = b·rguide_n` (and its phase
`csheet = exp(i·sheet)`), where `rguide_n` is a per-WF guiding centre initialised to the
projection centre and periodically re-estimated by a least-squares fit to the average
diagonal-overlap phases. `select_projections` is an orthogonal, purely-I/O feature: it picks
a subset of `.amn` columns to become the WFs.

---

## 0. Where guiding centres plug into the MV machinery

Everything in `localization.md` §1–§9 is unchanged. The *only* coupling is through the
pair `(sheet, csheet)` that appears inside the Im-ln, in exactly two places:

- **`wann_omega`** (spread), line 2116-2117:
  `ln_tmp(n,nn,k) = aimag(log(csheet(n,nn,k)·M_nn)) − sheet(n,nn,k)` = `q_n^{k,b}` (no `w_b`).
- **`wann_domega`** (gradient), line 2580-2581:
  `ln_tmp(n,nn,k) = w_b·(aimag(log(csheet·M_nn)) − sheet)` = `w_b·q_n^{k,b}`.

With guiding **off** (default): `csheet ≡ 1`, `sheet ≡ 0` (init line 310-311), so
`q_n^{k,b} = Im ln M_nn = atan2(Im,Re) ∈ (−π,π]` (principal value). This is the
`localization.md` baseline.

With guiding **on**: `csheet`, `sheet` are overwritten by `wann_phases` (§2) so that
```
q_n^{k,b} = Im ln( e^{i·sheet}·M_nn ) − sheet
          = wrap( φ_nn + b·rguide_n ) − b·rguide_n ,   φ_nn := Im ln M_nn
```
where `wrap(x) = atan2(sin x, cos x) ∈ (−π,π]` is the principal Im-ln applied to the
guide-rotated overlap. So the guide **shifts which 2π branch** of `φ_nn` is taken: instead of
the principal value of `φ_nn`, you take the value nearest to `−b·rguide_n`. When
`|φ_nn + b·rguide_n| < π` for all shells (a well-localized, well-guessed WF) the guided and
principal branches **coincide** and guiding is a no-op at the fixed point (this is exactly the
situation in the reference test — see §5, trap A).

**Modified centre / spread / gradient (structurally identical, only `q_n` re-branched).**
Substitute the guided `q_n^{k,b}` above into `localization.md` §2b/§2d/§2e/§3:
```
r_n        = −(1/N_k) Σ_{k,b} w_b · b · q_n^{k,b}                       (rave, minus sign, /N_k)
<r²>_n     =  (1/N_k) Σ_{k,b} w_b [ 1 − |M_nn|² + (q_n^{k,b})² ]        (r2ave)
Ω_D        =  (1/N_k) Σ_{k,b} w_b ( q_n^{k,b} + b·r_n )²                (om_d)
Ω_I, Ω_OD  =  unchanged (independent of sheet — pure |M|² sums)
G = dΩ/dW  =  (4/N_k) Σ_{k,b} w_b ( A[R] − S[T] ),  T_mn = R̃_mn·(q_n + b·r_n)
```
No new terms; `Ω_I` and `Ω_OD` do not depend on `sheet` at all (they are `|M|²` sums). At a
**fixed** M, changing the branch changes only `q_n`, hence only `Ω_D` and the centre
positions `r_n` — never `Ω_I` or `Ω_OD`. (Across a full guided-vs-unguided *run* the converged
M differs, so the endpoint `Ω_OD` can legitimately differ; the branch itself never touches it.)

---

## 1. Data structures (kmesh half-shell), and the schedule

### 1a. Half-shell quantities used by the guide (`kmesh.F90`)
- `nnh = nntot/2` (line 678): number of **b-directions** ignoring inversion (one per ±b pair).
- `bka(1:3, 1:nnh)` (line 693-695): the `nnh` distinct b-directions, taken from the neighbours
  of **k-point 1**, deleting any vector that is the negative of one already collected
  (`utility_compar` `ifneg`, lines 682-697). Cartesian, Å⁻¹ internally.
- `neigh(nkp, na)` (lines 735-743): for k-point `nkp` and half-direction `na∈[1,nnh]`, the
  full-neighbour index `nn∈[1,nntot]` whose b-vector equals `+bka(:,na)`. Used to sum the
  diagonal overlap along one consistent direction across all k.
- For the reference test: shell 1, `nntot = 8`, so `nnh = 4`, `bka` = the 4 distinct
  ⟨111⟩-type directions.

### 1b. Activation schedule (driver)
Keyword defaults (`wannier90_types.F90:178-180`):
`guiding_centres%enable = .false.`, `num_guide_cycles = 1`, `num_no_guide_iter = 0`.

- **Init `rguide` = projection centres** (lines 426-431), only if enabled:
  ```
  do n = 1, num_wann
      rguide(:,n) = utility_frac_to_cart( guiding_centres%centres(:,n), real_lattice )
  ! i.e. rguide(i,n) = Σ_a real_lattice(a,i)·centres(a,n)   (frac·real_lattice, Å)
  ```
- **First `wann_phases` call, before iteration** (lines 446-454), only if
  `num_no_guide_iter ≤ 0`. Passed `irguide = 0` (see §2, this suppresses the least-squares
  re-estimate, so this call just sets `sheet = b·(projection centres)`). Then `irguide := 1`.
- **In-loop call** (lines 557-566), every iteration for which
  `iter > num_no_guide_iter  .and.  mod(iter, num_guide_cycles) == 0`. Passed the current
  `irguide = 1`, so it **does** re-estimate `rguide`.
- **Final call, after the loop** (lines 895-900), to fix the branch for the final
  centre/spread report and the written `.chk`.

With the defaults (`num_no_guide_iter = 0`, `num_guide_cycles = 1`): the init call fires
(`0 ≤ 0`), and `mod(iter,1) == 0` is always true with `iter > 0`, so **guiding runs on every
single iteration** — maximally active. This is the reference-test configuration.

Schedule semantics to reproduce exactly:
- `num_no_guide_iter = K` ⇒ guide is silent for iters `1..K`, active from `K+1`. The
  pre-loop init call is skipped when `K > 0` (`K ≤ 0` guard fails), so early iters use the
  principal branch (`sheet ≡ 0`).
- `num_guide_cycles = C` ⇒ guide re-estimates only when `mod(iter,C) == 0`. Between guide
  calls `sheet`/`csheet` are frozen (they are `intent(out)` of `wann_phases`, untouched
  elsewhere).

---

## 2. `wann_phases` — the guide estimator (line 1776)

Computes `rguide` (per-WF Cartesian centre, Å) and from it `sheet`/`csheet`. Signature args
of interest: `csheet(out)`, `sheet(out)`, `rguide(inout)`, `irguide(in)` (0 on first call).

### 2a. Average diagonal overlap per b-direction (lines 1830-1854)
For each WF `n` (`loop_wann`) and each half-direction `na∈[1,nnh]`:
```
csum(na) = Σ_k  M^{k, neigh(k,na)}_{nn}          (sum of diagonal overlap over all k, along +bka(na))
```
(MPI: local partial sum then `comms_allreduce SUM`; serial: full k-sum. `m_w` optional
argument is the Gamma real-storage variant — ignore for a non-gamma port.)

### 2b. Least-squares fit for `rguide` (lines 1864-1929)
Model: `phase of csum(na) ≈ phase of exp(−i·bka(na)·rguide_n)`, i.e. want
`bka(na)·rguide_n ≈ xx(na)` with `xx(na) = −Im ln csum(na)`. Solve the normal equations
```
Σ_i smat(j,i)·rguide(i,n) = svec(j),   smat(j,i)=Σ_na bka(j,na)bka(i,na),  svec(j)=Σ_na bka(j,na)xx(na)
rguide(:,n) = (sinv/det)·svec        via utility_inv3(smat,sinv,det), only if |det|>eps6 (1e-6)
```
**Branch handling of `xx`, incremental over `na=1..nnh` (lines 1888-1911):**
- For `na ≤ 3`: `xx(na) = −aimag(log(csum(na)))` — arbitrary (principal) branch.
- For `na > 3`: predict `xx0 = Σ_j bka(j,na)·rguide(j,n)` from the current `rguide`, then
  ```
  xx(na) = xx0 − aimag( log( csum(na)·exp(i·xx0) ) )
  ```
  i.e. choose the 2π branch of `xx(na)` nearest the current guide prediction (consistency).
- The `rguide` update at line 1920-1926 runs **only if `irguide ≠ 0`**. On the pre-loop call
  (`irguide = 0`) the fit is computed but the assignment is skipped, so `rguide` stays at the
  projection centres. From the first in-loop call on (`irguide = 1`) it is re-estimated.
- `smat`/`svec` accumulate as `na` grows; `rguide` is re-solved once `na ≥ 3` and updated each
  further `na` (nested inside the `na` loop, lines 1913-1929).

### 2c. Build `sheet`/`csheet` (lines 1934-1957)
```
sheet(n,nn,nkp) = Σ_j bk(j,nn,nkp)·rguide(j,n)      (full nntot neighbours, all k)   [non-gamma]
csheet          = exp(i·sheet)
```
(`use_ss_functional` variant collapses to `nkp=1` only — out of scope.) `bk` is the full
per-(nn,k) Cartesian b-vector (Å⁻¹); note this uses `bk`, not the half-shell `bka`.

**Net effect:** `q_n^{k,b} = Im ln(csheet·M_nn) − sheet = wrap(φ_nn + b·rguide_n) − b·rguide_n`.
The dead commented block at 1959-1987 documents the intent: pick the sheet so that
`q_n^{k,b} ≈ 0` for a good solution and `≈ 2π·integer` for a bad one, avoiding a wrong wrap in
poorly-localized starting WFs.

---

## 3. `select_projections` — subset of `.amn` columns

Purpose: when the `projections` block (and hence the `.amn`) has **more** projections
`num_proj` than target Wannier functions `num_wann`, pick which `num_wann` columns of the
`.amn` become the WFs. Type: `select_projection_type` (`wannier90_types.F90:262-272`),
`lselproj=.false.` default, `proj2wann_map(:)` the column→WF map.

### 3a. Parsing (`wannier90_readwrite.F90:1622-1696`)
- Read as a range vector `select_projections` of length `num_select_projections`
  (`w90_readwrite_get_range_vector`, supports `1 2 3 4` and `5-12` syntax).
- Validation (lines 1649-1668): all entries `≥1`; `num_select_projections == num_wann`
  exactly (else "too few"/"too many projections selected"); requires `lhasproj` (a
  `projections` block must exist); `max(select_projections) ≤ num_proj`.
- Build the map (lines 1686-1696):
  ```
  if lselproj:  proj2wann_map(i) = j   where select_projections(j) == i   (else −1)   for i=1..num_proj
  else:         proj2wann_map(i) = i                                       for i=1..num_wann
  ```
  So `proj2wann_map(input_col) = target_WF_index`, `−1` for unselected columns.

### 3b. Effect on `.amn` read (`overlap.F90:361-377`)
```
if num_proj > num_wann and not lselproj:  ERROR (too many projections to use without selecting)
read each amn line (m, n, nkp, a_real, a_imag):
    if proj2wann_map(n) < 0: cycle                              ! drop unselected column
    au_matrix(m, proj2wann_map(n), nkp) = cmplx(a_real,a_imag)  ! place into WF slot
```
So the `.amn` header must still declare `num_proj` (all columns), but only the selected
`num_wann` columns land in `au_matrix` (the `A^{k}` used to build the starting gauge / disentangle).

### 3c. Interaction with `num_wann` / `num_bands`
- `num_wann` = number of WFs = number of selected columns (enforced equal).
- `num_bands` (=12 here) is independent; it is the number of Bloch states / `.amn` rows `m`.
  `num_bands > num_wann` triggers disentanglement (§4), orthogonal to select_projections.
- `num_proj` (=12 here, the row count of the `projections` block) is the `.amn` column count;
  `select_projections` reduces the *used* columns from `num_proj` to `num_wann`.

### 3d. Interaction with guiding-centre `centres` array (load-bearing under non-identity select)
`guiding_centres%centres` is filled from the **reordered** `proj` array
(`wannier90_readwrite.F90:169-178`): `centres(:,ip) = proj(ip)%site(:)`, and `proj` was built
by `proj(proj2wann_map(loop)) = proj_input(loop)` (lines 1708-1716). Therefore
`centres(:,n)` for `n=1..num_wann` holds the **site of the selected projection that maps to WF
`n`**, i.e. it is already in **WF order over the selected subset** — not raw input order. The
driver reads only `centres(:,1..num_wann)` (line 427), so entries `num_wann+1..num_proj`
(allocated but junk) are never used. In the reference test `select_projections = 1 2 3 4` is
the identity on the first 4, so `centres(:,1..4)` = the 4 bond-centred s-orbital sites
(fractional `(-0.125,-0.125,0.375)` etc.) — exactly the guide seeds.

---

## 4. Interplay with disentanglement (this test sets both)

The test runs disentanglement (`num_bands=12 > num_wann=4`, `dis_num_iter=100`) **and**
guiding centres. They are sequential stages, coupled only through the `M`/`U` that
disentanglement hands to `wann_main`:

1. **Disentanglement first** (`dis_main`), independent of guiding centres. Here it is
   **trivial/fully frozen**: `dis_froz_max = dis_win_max = 6.5 eV` captures exactly the 4
   valence bands, so every k-point has `Ndimwin = Ndimfroz = 4` (benchmark lines 369-436).
   The DIS iteration table (benchmark 447-454) shows `Omega_I = 20.71710600` **constant from
   iteration 1** — nothing to optimize, the optimal subspace is the frozen one.
   `Final Omega_I = 20.71710600 Bohr²`. Guiding centres play **no role** in disentanglement
   (`wann_phases` is called only inside `wann_main`).
2. **Wannierise second** (`wann_main`), with guiding active every iteration (§1b). It takes the
   `num_wann × num_wann` gauge and relaxes it. Because `Ω_I` is fixed by disentanglement and
   independent of `sheet`, the branch choice directly affects only `Ω_D` and the centres at
   fixed M; `Ω_I` and `Ω_OD` are `|M|²` sums independent of `sheet` (the converged `Ω_OD` moves
   only through the gauge relaxation, not the branch).

**Why the two features are safe to port independently:** `Ω_I` is a `|M|²` sum, invariant to
both the gauge rotation and the branch choice; disentanglement fixes it; guiding centres and
the MV rotation then act purely on `Ω_D + Ω_OD` at fixed `Ω_I`. A gauge-invariant Julia
validation (Ω decomposition + centres + interpolated bands) will match as long as the port
(a) reproduces the fully-frozen disentanglement (Ω_I) and (b) reaches the same symmetric MV
minimum — the guided vs principal branch coincide there (§5).

---

## 5. Reference test — `testw90_guidingcentre_selectproj`

### 5a. `silicon.win` keywords (all relevant lines quoted)
```
num_bands        =   12
num_iter         =  100
dis_num_iter     =  100
length_unit      =  bohr            ! ⇒ all wout numbers below are in Bohr / Bohr²
guiding_centres  = .true.
search_shells    =  12
iprint           =    2             ! verbose enough to print SELECTED PROJECTIONS block
begin projections                   ! 12 projections (num_proj = 12)
  f=-0.125,-0.125, 0.375:s          ! 4 bond-centred s-orbitals ...
  f= 0.375,-0.125,-0.125:s
  f=-0.125, 0.375,-0.125:s
  f=-0.125,-0.125,-0.125:s
  Si:sp3                            ! ... + 8 Si sp3  (2 atoms × 4)
end projections
num_wann          =   4
select_projections = 1 2 3 4        ! pick the 4 bond-centred s-orbitals as the WFs
dis_froz_max     =   6.5
dis_win_max      =   6.5            ! froz==win ⇒ fully frozen, trivial disentanglement
mp_grid          =  4 4 4           ! 64 k-points
```
`num_no_guide_iter`, `num_guide_cycles` are **not** set ⇒ defaults 0 and 1 ⇒ guide runs
every iteration from the pre-loop init on.

### 5b. What is checked and tolerances (jobconfig / userconfig)
- `jobconfig` (line 581-584): program `WANNIER90_WOUT_OK`, input `silicon.win`, output
  `silicon.wout`, parsed by `tools/parsers/parse_wout.parse`.
- `userconfig` `[WANNIER90_WOUT_OK]` tolerances `(abs, rel, name)` — the guiding/select test is
  validated on the standard `.wout` quantities:
  ```
  final_centres_x/y/z : 1.0e-5 / 1.0e-5
  final_spreads       : 3.0e-6 / 3.0e-6
  omegaI              : 1.0e-6 / 1.0e-6
  omegaD              : 1.0e-6 / 5.0e-6
  omegaOD             : 1.0e-6 / 1.0e-6
  omegaTotal          : 1.0e-6 / 1.0e-6
  near_neigh_dist/mult, completeness_x/y/z/weight : 1.0e-6
  ```
  These are exactly the gauge-invariant quantities (centres, per-WF spreads, Ω
  decomposition, shell geometry) — **not** `U(k)` or `H(R)`. Matches the WannierFunctions.jl
  validation policy.

### 5c. Anchor numbers from `benchmark.out.default.inp=silicon.win` (units = Bohr / Bohr²)
Final State (benchmark lines 630-640), all 4 WFs symmetric with identical spread:
```
Final State
  WF centre and spread  1  ( -1.275000,  1.275000, -1.275000 )     5.68462413
  WF centre and spread  2  ( -1.275000, -1.275000,  1.275000 )     5.68462413
  WF centre and spread  3  (  1.275000,  1.275000,  1.275000 )     5.68462413
  WF centre and spread  4  (  1.275000, -1.275000, -1.275000 )     5.68462413
  Sum of centres and spreads ( -0.000000, -0.000000,  0.000000 )   22.73849650

  Omega I      =    20.717105998
  Omega D      =     0.000000000
  Omega OD     =     2.021390507
  Omega Total  =    22.738496505
```
Cross-checks in the benchmark:
- Init (iter 0): `Ω_OD = 2.0262387`, `Ω_TOT = 22.7433447`, spreads `5.68583618` — guiding is
  active from the start yet the relaxation only moves `22.7433 → 22.7385` (tiny), because the
  projection guess is already the symmetric bond-centred solution.
- `Ω_D = 0.000000000` at **every** printed iteration (0,1,10,20,…,100). The symmetric
  solution has each WF centre exactly at a bond centre, `q_n + b·r_n = 0` with no wrap.
- Converges by iter ~10 (`Delta ~1e-13`); runs full 100 iters (`conv_window = -1`, off).

**Unit note (CODATA2006).** These are **Bohr** because `length_unit = bohr`. Internal math is
Å²; `lenconfac = 1/bohr` is a display multiplier (centres ×lenconfac, spreads ×lenconfac²).
With CODATA2006 `bohr = 0.52917720859 Å` (so `bohr² = 0.28002834 Å²`): centres scale ×`bohr`
(`1.275 Bohr = 0.674701 Å`) and spreads/Ω scale ×`bohr²` (`Ω_TOT 22.738496505 Bohr² = 6.3674 Å²`,
`Ω_I 20.717106 Bohr² = 5.8014 Å²`, per-WF spread `5.68462413 Bohr² = 1.59185 Å²`,
`4 × 1.59185 = 6.3674 = Ω_TOT ✓`). A Julia port computing in Å² must apply `bohr²` before
comparing to these Bohr² numbers (or convert these to Å²). See `localization.md` §0.

---

## 6. Traps for the Julia port

**A. The benchmark header prints `Use guiding centre to control phases : F`** (benchmark line
259) even though `guiding_centres = .true.` and the run is non-gamma with `nntot = 8`. In the
**current** reference source the only place `guiding_centres%enable` is forced off is
`wannierise.F90:3060` — inside `wann_main_gamma`, gated on `nntot == 3` — which does **not**
apply here. So current source runs guiding **actively** and would echo `T`; the benchmark's
`F` is a version artifact of the header echo and must **not** be trusted as behavior. Do **not**
key the port off that echo. The Ω/centre anchors are unaffected because the guided and
principal branches coincide at this symmetric fixed point (`Ω_D = 0`, integer-wrap-free), so
on-vs-off give the same endpoint — which is also why the `WOUT_OK` extractor (Ω/centres/
spreads, not the header) passes regardless.

**B. Two `ln_tmp` conventions still apply** (inherited from `localization.md` §2a/§3a): `q_n`
(no `w_b`) in `wann_omega`, `w_b·q_n` in `wann_domega`. Guiding only changes what `q_n` is
(re-branched), not this asymmetry.

**C. `sheet` uses full `bk`, the guide fit uses half-shell `bka`.** `csum`/`xx`/`rguide` are
built from `bka(:,1:nnh)` and `neigh` (§2a-b); the final `sheet(n,nn,k) = Σ_j bk(j,nn,k)·
rguide(j,n)` uses the full `bk(:,nntot,nkpts)` (§2c). Don't conflate the two b-arrays.

**D. `irguide` gate.** First (pre-loop) call passes `irguide = 0` ⇒ `rguide` stays at the
projection centres (fit computed, assignment skipped). All later calls pass `irguide = 1` ⇒
`rguide` is re-estimated from `csum`. Reproduce this or the very first branch choice differs.

**E. `rguide` init is `frac · real_lattice`** (row convention, `utility_frac_to_cart`,
`utility.F90:459-461`): `rguide(i,n) = Σ_a real_lattice(a,i)·centres(a,n)`, in Å. `centres` is
the selected projection sites in WF order (§3d).

**F. `select_projections` reorders `proj` before `centres` is filled.** `centres(:,n)` is the
site of the projection mapped to WF `n` (via `proj2wann_map`), not input row `n`. Identity in
this test, but a non-identity `select_projections` (e.g. the commented `5-12` alternative in
the win file) would seed the guide from the *selected* subset.

**G. Gamma-only extra rule (only if you ever port `wann_main_gamma`):** `wannierise.F90:3060`
sets `guiding_centres%enable = .false.` when `nntot == 3` (orthorhombic 3-shell case). Not hit
by the reference test.

**H. `Ω_I` is independent of `sheet`.** Guiding centres and disentanglement both leave `Ω_I`
untouched by the branch choice; only `Ω_D` (and centre positions) respond. A regression in the
port that shifts `Ω_I` when toggling guiding is therefore definitely a bug.

**I. select_projections requires an explicit `projections` block** (`lhasproj`), the `.amn`
header must declare the full `num_proj` (12) column count, and `num_select == num_wann`
exactly. `num_proj > num_wann` without `lselproj` is a fatal error (`overlap.F90:361`).
