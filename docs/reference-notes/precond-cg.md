# Preconditioned CG for maximal localization (`precond = .true.`)

Reference: Wannier90 v3.1.0, `src/wannierise.F90`. All line numbers are for that file
unless another file is named. This note is a **delta on top of
`docs/reference-notes/localization.md`** — read that first. It covers ONLY what changes
when `precond = .true.`; the spread `wann_omega`, gradient `wann_domega`, line search,
and unitary update are byte-for-byte identical to the plain path and are not repeated here.

The one-sentence answer: **the preconditioner is a real-space Lorentzian low-pass filter
applied to the gradient** (Fourier k→R, multiply each R-component by
`1/(1 + |R_cart|²/α)`, transform back). It is **not** a diagonal scaling by b-vector
weights and **not** a k-independent scalar — it is diagonal in the Wigner–Seitz R-index,
hence k-*dependent*. In preconditioned-CG language it is the operator `M⁻¹` applied to the
gradient before the Fletcher–Reeves direction update. Everything else (FR ratio, `>3`
reset, `num_cg_steps` cap, parabolic line search) is structurally unchanged; only the
inner product and the steepest-descent component of the direction see the filtered
gradient.

---

## 0. Keyword and defaults

- `precond` (logical) default **`.false.`** — `wannier90_types.F90:202`
  (`logical :: precond = .false.`), read in
  `wannier90_readwrite.F90:741-742` (`w90_readwrite_get_keyword(settings,'precond', … l_value=wann_control%precond)`).
- `optimisation` (integer) default **3** — `library_interface.F90:129`
  (`integer :: optimisation = 3`), read in `readwrite.F90:149-150`. Selects the transform
  implementation only (see §2); `optimisation >= 3` = BLAS `zgemm` path,
  `optimisation <= 2` = explicit double DFT loop. **Both paths are mathematically identical**
  — that is exactly why `testw90_precond_2` sets `optimisation = 2`.
- Unchanged CG/step defaults inherited from the plain path:
  `num_cg_steps = 5` (`wannier90_types.F90:194`), `trial_step = 2.0`
  (`wannier90_types.F90:201`), `num_iter = 100` default (overridden to 40 in the tests).
- Extra memory allocated only when `precond`: `cdodq_r(num_wann,num_wann,nrpts)`,
  `cdodq_precond(num_wann,num_wann,num_kpts)`, `cdodq_precond_loc(…,nkrank)`, and (for
  `optimisation>=3`) `k_to_r(num_kpts,nrpts)` (lines 282-307, 348-354). A
  `hamiltonian_setup` call at line 276-280 exists ONLY to build the Wigner–Seitz R-grid
  (`irvec`, `ndegen`, `nrpts`); the port needs the R-grid, **not** the `ham_k`/`ham_r`
  interpolation machinery it also allocates.

---

## 1. The preconditioner operator — `precond_search_direction` (line 1264)

Called at line 591-594 in the driver, immediately **before** `internal_search_direction`,
only if `precond`. Inputs the true k-space gradient `cdodq(m,n,k)` (from `wann_domega`);
outputs `cdodq_precond(m,n,k) = M⁻¹·cdodq` and its per-rank copy `cdodq_precond_loc`.

Let `g^k = cdodq(:,:,k)` be the anti-Hermitian gradient block at k-point k (§3 of
localization.md). The operator is three steps:

### 1a. Forward Fourier transform k → R (lines 1322-1334)

```
g_R(:,:,R) = (1/N_k) Σ_k exp(-i·2π k·R) · g^k          (R over WS supercell vectors)
```
- `optimisation >= 3` (lines 1323-1325): `zgemm('N','N', …)` contracting the pre-tabulated
  `k_to_r(k,R) = exp(-i·2π k·R)` (built once at driver lines 301-306), then
  `cdodq_r = cdodq_r/num_kpts`.
- `optimisation <= 2` (lines 1327-1333): explicit double loop,
  `rdotk = 2π·kpt_latt(:,k)·irvec(:,R)`, `fac = exp(-i·rdotk)/num_kpts`, accumulate.
- `kpt_latt` is the k-point in **fractional** (crystal) coordinates; `irvec` are integer
  lattice-vector triples. Their dot product times `2π` is the Bloch phase — the standard
  W90 k↔R convention (same as H(R) interpolation).

### 1b. Real-space Lorentzian filter (lines 1336-1351) — THE PRECONDITIONER

```
alpha_precond = 10 · om_tot / num_wann                         (line 1346)
for each R:
    R_cart = real_lattice · irvec(:,R)                         (Cartesian Å, line 1348)
    g_R(:,:,R) *= 1 / ( 1 + (R_cart·R_cart)/alpha_precond )    (line 1349-1350)
```
- **`alpha_precond` is state-dependent and recomputed every iteration** from the *current*
  total spread `wann_spread%om_tot`. It is a scalar (same for all R, all m,n), NOT tied to
  b-vectors. Units: `om_tot` is Å², `num_wann` dimensionless ⇒ `α` in Å²; `R_cart` in Å ⇒
  `|R_cart|²/α` dimensionless. Consistent with the Å/Å²-internal convention
  (localization.md §0). The code comment (lines 1338-1345) states the value is "more or
  less arbitrary", chosen only to have the right units and neither over- nor under-filter.
- `real_lattice` is the lattice in Å (columns are lattice vectors); `matmul(real_lattice,
  irvec)` = Cartesian R.
- The filter is a **low-pass**: the `R = 0` component is unattenuated (`1/(1+0)=1`); large-|R|
  Fourier components of the gradient (short-wavelength / high-"frequency" in k) are damped.
  This suppresses the far-neighbour couplings that make plain SD/CG take small steps.

### 1c. Backward Fourier transform R → k (lines 1353-1369)

```
g̃^k = cdodq_precond(:,:,k) = Σ_R (1/ndegen(R)) exp(+i·2π k·R) · g_R(:,:,R)
```
- `optimisation >= 3` (lines 1354-1359): first divide each R-slice by `ndegen(R)` (line 1356),
  then `zgemm('N','C', …)` against `conjg(k_to_r)` = `exp(+i·2π k·R)`.
- `optimisation <= 2` (lines 1360-1368): explicit loop, `fac = exp(+i·rdotk)/ndegen(R)`.
- `ndegen(R)` is the Wigner–Seitz degeneracy weight of R-vector R (same array as H(R)
  interpolation). Dividing by it makes the forward+backward pair the **identity** when the
  filter is 1 — so the preconditioner is *purely* the `1/(1+R²/α)` weighting and nothing
  else. Confirm this in a port by setting the filter to 1 and checking `g̃^k == g^k`.
- Then copy to per-rank storage (lines 1370-1373): `cdodq_precond_loc(:,:,nkp_loc) =
  cdodq_precond(:,:,global_k(nkp_loc))`. For a serial port `nkp_loc == nkp`.

### 1d. Provenance of `irvec`, `ndegen`, `nrpts`

From `hamiltonian_wigner_seitz` (via `hamiltonian_setup`, lines 276-280;
`hamiltonian.F90:113-165`): the Wigner–Seitz supercell of the `mp_grid` real-space lattice,
with degeneracy weights `ndegen`. **The R-set is symmetric under R ↔ −R** and the filter
`1/(1+|R|²/α)` is even in R. This is why applying the filter preserves the
anti-Hermiticity of `g^k` in the (m,n) WF indices — see traps.

---

## 2. Modified formulas vs plain Fletcher–Reeves CG

Notation: `g = cdodq` (true gradient), `g̃ = cdodq_precond = M⁻¹g` (filtered gradient),
`d = cdq` (search direction), `d_prev = cdqkeep` (previous direction). All four changes are
localized in `internal_search_direction` (line 1378); the branch is `if
(wann_control%precond)`.

| Quantity | Plain FR-CG | Preconditioned (`precond`) | Code |
|---|---|---|---|
| **FR norm `gcnorm1`** | `Re⟨g, g⟩ = Σ_k Σ_{mn}\|g_{mn}^k\|²` | `Re⟨g̃, g⟩ = Re Σ_k Σ_{mn} conj(g̃_{mn}^k)·g_{mn}^k` (mixed) | 1435 vs 1439 |
| **FR coefficient `gcfac = β`** | `gcnorm1/gcnorm0` | same ratio, but numerator/denominator are the mixed norm | 1451 |
| **Direction `d`** | `d = g + β·d_prev` | `d = g̃ + β·d_prev` | 1472 vs 1474 |
| **Slope `doda0`** | `-Re⟨g, d⟩/(4·wbtot)` | **identical** — uses the **TRUE** gradient `g`, not `g̃` | 1488-1494 |
| **Uphill reset direction** | `d = g` | `d = g` (**TRUE** gradient, not `g̃`) | 1502 |

Precise statements:

**(1) Steepest-descent component of the direction** (line 1472):
```
cdq_loc = cdodq_precond_loc + cdqkeep_loc · gcfac          ! precond
cdq_loc = cdodq_loc        + cdqkeep_loc · gcfac           ! plain (line 1474)
```
Only the leading (steepest-descent) term is replaced by the filtered gradient; the CG
memory term `β·d_prev` is untouched.

**(2) Fletcher–Reeves numerator uses the mixed inner product** (lines 1433-1441):
```
precond:  gcnorm1 = Re( zdotc(cdodq_precond_loc, cdodq_loc) ) = Re⟨M⁻¹g, g⟩
plain:    gcnorm1 = Re( zdotc(cdodq_loc,        cdodq_loc) ) = Re⟨g, g⟩
```
Then `gcfac = gcnorm1/gcnorm0` (line 1451) with `gcnorm0` = previous iter's `gcnorm1`
(line 1468). The `>3.0` reset-to-SD guard (lines 1453-1457), the `iter==1 ||
ncg>=num_cg_steps` SD reset (line 1446), and the `ncg` counter are **unchanged**. This is
the textbook preconditioned-CG modification: β becomes `⟨g_i, M⁻¹g_i⟩ / ⟨g_{i-1},
M⁻¹g_{i-1}⟩` (here as `Re⟨M⁻¹g,g⟩`, the anti-Hermitian gradient making the real part the
relevant scalar). It is still the **Fletcher–Reeves** form (a ratio), NOT
Polak–Ribière — preconditioning changes the inner product, not FR→PR.

**(3) The line-search slope `doda0` still uses the TRUE gradient** (lines 1485-1494):
```
doda0 = -Re( zdotc(cdodq_loc, cdq_loc) ) / (4·wbtot)      ! cdodq_loc = TRUE gradient g
```
This is **correct and deliberate**, not an oversight: `doda0` is the directional derivative
of the objective Ω along the search direction `d`, which is `⟨∇Ω, d⟩` and must use the true
gradient `g = ∇Ω`, regardless of how `d` was preconditioned. A port must NOT "fix" this to
use `g̃`. Formula and `4·wbtot` normalization are identical to the plain path
(localization.md §4).

**(4) The uphill-reset fallback resets to the TRUE gradient** (lines 1497-1531): if `doda0
> 0`, a CG step resets `cdq_loc = cdodq_loc` (line 1502) — the *unfiltered* gradient — then
recomputes `doda0`, and if still uphill reverses the direction. An SD step just reverses.
Same logic as plain; the only subtlety is it reverts to `g`, not `g̃`.

**Unchanged downstream:** the search direction `d = cdq_loc` then flows into the **identical**
trial-step / parabolic line search / `internal_new_u_and_m` machinery (driver lines
607-702). The generator exponentiated is still `(α_step/(4·wbtot))·d` with `α_step =
trial_step` then `alphamin` (localization.md §5-6). The `4·wbtot` step normalization and the
`×4/N_k` gradient prefactor are untouched.

---

## 3. Driver wiring (`wann_main`, differences only)

```
if precond: hamiltonian_setup → irvec, ndegen, nrpts        (276)    [WS R-grid only]
            allocate cdodq_r, cdodq_precond, cdodq_precond_loc, k_to_r  (282-307,348-354)
for iter = 1 .. num_iter:
    [guiding centres if enabled]
    wann_domega → cdodq_loc  (and full-BZ cdodq when precond|sitesym)  (569-585)
    if precond: precond_search_direction(cdodq → cdodq_precond[_loc])  (591-594)
    internal_search_direction(cdodq_precond_loc, …, cdodq_loc, …)      (596)
        → cdq_loc = g̃ + β·d_prev ; doda0 = -⟨g,d⟩/(4 wbtot)
    … identical trial step, line search, U/M update …                 (607-702)
```
Note at line 569: when `precond` (or `lsitesymmetry`) the `wann_domega` call requests the
full-BZ gradient array `cdodq` (extra actual arg) in addition to `cdodq_loc`, because the
k→R transform needs all k-points on every rank. For a serial port `cdodq == cdodq_loc` over
all k.

---

## 4. Reference tests

Both tests are the GaAs Example-1 system (`num_wann=4`, `sp3` on As, 2×2×2 k-mesh,
`num_iter=40`, `search_shells=12`, `use_ws_distance=.false.`, `wvfn_formatted=.true.`).
They exist to prove preconditioning reaches the **same minimum** as the plain path, just via
a different (faster) trajectory — so the anchor numbers below are the *converged* spreads,
identical between the two tests.

### 4a. `testw90_precond_1/gaas1.win`
Keywords that differ from a plain MV run:
```
precond = true
```
(plus `num_wann=4`, `num_iter=40`, `use_ws_distance=.false.`, `search_shells=12`;
`optimisation` defaults to 3 → GEMM transform path.) Unit cell in `bohr`.

### 4b. `testw90_precond_2/gaas2.win`
Identical to test 1 plus:
```
precond = true
optimisation = 2
```
`optimisation = 2` forces the explicit double-DFT-loop transform (§1a/1c, `else` branches)
instead of `zgemm`. Same math ⇒ same result.

### 4c. Benchmark converged state (IDENTICAL for both tests)
From `benchmark.out.default.inp=gaas1.win:687-697` (and `…gaas2.win:687-697`, bit-identical):
```
 Final State
  WF centre and spread    1  (  0.866253, -0.866253,  0.866254 )     1.11671997
  WF centre and spread    2  (  0.866254,  0.866254, -0.866253 )     1.11672036
  WF centre and spread    3  ( -0.866254,  0.866253,  0.866253 )     1.11672024
  WF centre and spread    4  ( -0.866253, -0.866254, -0.866254 )     1.11672041
  Sum of centres and spreads (  0.000000, -0.000000,  0.000000 )     4.46688098

         Spreads (Ang^2)       Omega I      =     3.956862958
        ================       Omega D      =     0.008030049
                               Omega OD     =     0.501987969
    Final Spread (Ang^2)       Omega Total  =     4.466880976
```
Anchor numbers (Å², Å):
- Ω_I  = **3.956862958**   (gauge-invariant; the disentanglement floor)
- Ω_D  = **0.008030049**
- Ω_OD = **0.501987969**
- Ω_Total = **4.466880976**  (= Ω_I + Ω_D + Ω_OD; check: 3.956862958 + 0.008030049 +
  0.501987969 = 4.466880976 ✓)
- Four WF centres at the four As–Ga bond midpoints, each a permutation of
  (±0.866253, ±0.866253, ±0.866254) Å; per-WF spread ≈ **1.11672** Å²; sum of spreads =
  **4.46688098** Å²; sum of centres ≈ 0.

`|centre| = 0.866253 Å = (√3/4)·a_cubic·... ` — the sp³ WFs sit on the Ga–As bonds, the
expected physical result. All four spreads equal to ~5 sig figs (symmetry-equivalent bonds).

### 4d. Tolerances (`WANNIER90_WOUT_OK`, `tests/userconfig:6-20`)
Both tests use `program = WANNIER90_WOUT_OK` (`tests/jobconfig:184-194`). Compared
quantities and (abs, rel) tolerances:
```
final_centres_x/y/z : 1.0e-5 / 1.0e-5
final_spreads       : 3.0e-6 / 3.0e-6
omegaI              : 1.0e-6 / 1.0e-6
omegaD              : 1.0e-6 / 5.0e-6
omegaOD             : 1.0e-6 / 1.0e-6
omegaTotal          : 1.0e-6 / 1.0e-6
near_neigh_dist/mult, completeness_x/y/z/weight : 1.0e-6 (b-vector/mesh metadata)
```
These are all **gauge-invariant** quantities (spreads, centres, Ω decomposition, b-vector
mesh) — consistent with the project rule: validate the port by Ω/centres, never raw
U(k)/H(R). A Julia port passes if it reaches this same minimum within tolerance; it need NOT
reproduce the preconditioned iteration trajectory.

---

## 5. Traps for the port

1. **The preconditioner is a real-space filter, not a b-weight/diagonal scaling.** Do not
   look for a `w_b`-based or per-WF scaling — it is `1/(1+|R_cart|²/α)` diagonal in the WS
   R-index, applied via k→R→k Fourier round-trip (§1). This is `M⁻¹g` in preconditioned-CG
   terms.
2. **`alpha_precond = 10·om_tot/num_wann` is recomputed every iteration** from the current
   total spread. Not a fixed constant. Getting the factor 10 or the `/num_wann` wrong shifts
   the filter strength (but converges to the same minimum, just slower/faster).
3. **Three inner-product sites, three different gradients** (the single most dangerous
   confusion):
   - FR norm `gcnorm1` = `Re⟨g̃, g⟩` — the **mixed** product (line 1435).
   - CG direction SD term = `g̃` (filtered, line 1472).
   - Slope `doda0` and uphill-reset = **`g`** (true, unfiltered, lines 1488, 1502).
   Using `g̃` for `doda0` or `g` for the direction both break faithfulness (though may still
   converge — the line search self-corrects). Match all three.
4. **Fourier round-trip must use `1/N_k` forward and `1/ndegen(R)` backward.** Omitting
   `ndegen` breaks the identity property (filter=1 must give `g̃=g`) and mis-weights
   degenerate WS R-vectors.
5. **R-grid must be the symmetric Wigner–Seitz supercell** (R ↔ −R) from
   `hamiltonian_wigner_seitz`. Because the filter is even in R, this symmetry is what keeps
   the filtered generator `g̃^k` anti-Hermitian. An asymmetric or wrong R-set silently
   destroys anti-Hermiticity and the exp-map step stops being unitary.
6. **`trial_step = 2.0` is NOT "consistent" under preconditioning** — the code comment
   (lines 1343-1345) warns the preconditioned direction has a different magnitude than the
   plain one, so the trial step probes a different physical distance. The parabolic line
   search (`doda0` + trial spread → `alphamin`) self-corrects, so convergence is fine; do
   not hand-tune `trial_step` for the preconditioned path.
7. **Still Fletcher–Reeves, not Polak–Ribière.** The prompt's "PR+" label is wrong for this
   code: β is the FR ratio of mixed norms; preconditioning changes the inner product, not
   the β formula family. Keep terminology consistent with localization.md §4.
8. **`optimisation` selects only the transform implementation** (GEMM vs loop) — identical
   results. Do not treat `optimisation=2` (test 2) as a different physics path; it is a
   coverage test of the non-GEMM branch.
9. **Validate on Ω/centres, not trajectory.** The path-level subtleties (items 3, 6) are for
   *faithfulness*; the tests only check the gauge-invariant converged minimum (§4d). A port
   that reaches Ω_Total=4.466880976 with the centres above is correct even if its
   iteration-by-iteration spreads differ from the benchmark's intermediate values.
