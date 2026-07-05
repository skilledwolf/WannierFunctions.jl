# Projectability-Disentangled Wannier Functions (PDWF)

Implementation-grade notes for the projectability-based frozen-window variant of
disentanglement (`dis_froz_proj`), extracted from the Wannier90 v3.1.0 reference
source. Reference paper: **Qiao, Pizzi, Marzari, *npj Comput. Mater.* 9, 208
(2023)** ("Projectability disentanglement for accurate and automated
electronic-structure Hamiltonians"). File:line citations are to
`reference/wannier90/src/…` as named. Read this **alongside**
`docs/reference-notes/disentanglement.md` — PDWF changes ONLY the *window/frozen
selection* step; everything downstream (`dis_project`, `dis_proj_froz`,
`dis_extract` Z-matrix iteration, `internal_find_u`, handoff to MV) is the
**identical** energy-window disentanglement code and is documented there. This
note covers the delta.

Constants: reference build defaults to **CODATA2006** (`constants.F90:92–96`,
`#define CODATA2006` when no other year is set). No energies/lengths in the PDWF
selection depend on the constant set (projectabilities are pure numbers,
eigenvalues come verbatim from `.eig` in eV); CODATA only matters for the
Bohr↔Å conversions in kmesh/lattice, unchanged from the energy path.

--------------------------------------------------------------------------------
## 0. What PDWF is, in one paragraph

Standard SMV disentanglement freezes states inside an **energy** window
(`dis_froz_min/max`). PDWF instead freezes states by their **projectability**
`p_{nk}` onto the target atomic-orbital subspace defined by the `.amn`
projections: high-projectability states (`p ≥ dis_proj_max`) are frozen (locked
into the optimal subspace unchanged), mid-projectability states
(`dis_proj_min ≤ p < dis_proj_max`) are the disentanglement pool that the
Z-matrix iteration optimizes over, and low-projectability states
(`p < dis_proj_min`) are **discarded entirely** — removed from the window even if
they lie in the middle of the outer energy window. The `.amn` is expected to come
from pseudo-atomic-orbital (or SCDM) projections generated *externally*
(`auto_projections`); w90 itself does **not** compute SCDM.

--------------------------------------------------------------------------------
## 1. The projectability `p_{nk}` and the frozen/disentangle/discard criterion

### 1a. Definition of `p_{nk}` — EXACT

Computed in `dis_windows_proj` (`disentangle.F90:1335–1349`) on the **raw
`a_matrix`** (the `.amn` projections `A_{n,w}(k)`), **before** any SVD
orthonormalization (`dis_project` runs *after*, at `disentangle.F90:166`):

```
projs(i) = Σ_{j=1..num_wann}  [ Re A(i,j,k)^2 + Im A(i,j,k)^2 ]     ! line 1338
         = Σ_{w=1..num_wann} |A_{i w}(k)|^2
```

with `i = 1..num_bands` the band index. So:

**`p_{nk} = Σ_{w=1..num_wann} |A_{nw}(k)|²`** — the **row-norm² of the `.amn`
matrix**, i.e. the diagonal of `A A^†`:  `p_{nk} = (A A^†)_{nn}`.

- It is the **sum over Wannier/projection columns of |A|²**, NOT SVD singular
  values, NOT a diagonal of `A^† A`. (Column count here = `num_wann`; see §2 on
  why `num_proj == num_wann` for this to be a valid probability.)
- `a_matrix` at this point is the un-orthonormalized projection read from disk;
  `a_matrix = u_matrix_opt` on entry (`disentangle.F90:120`), and
  `u_matrix_opt` was loaded with the raw `.amn` A (`overlap.F90:376`).
- **Hard bound check:** `projs(i)` must lie in `[0, 1]`, else fatal error
  (`disentangle.F90:1340–1348`, `'projectability < 0.0 or > 1.0'`). This holds
  only if the `.amn` columns are (quasi-)orthonormal, i.e. `A^† A ≈ I` — true for
  SCDM/atomic-orbital projections with `num_proj = num_wann`. A port must
  replicate the check; if it trips, the `.amn` is not a valid projectability
  source.

### 1b. Selection criterion — EXACT (`disentangle.F90:1366–1390`)

Per k, loop `i = 1..num_bands` over ALL bands (eigenvalues assumed ascending):

```
! 1. DISCARD if outside outer ENERGY window (inclusive test via <, >):
if (eig(i) < win_min  .or.  eig(i) > win_max)  cycle          ! line 1368-1369

! 2. FROZEN  (union of high-projectability AND energy-inner-window):
if ( projs(i) >= proj_max                                     ! line 1371
     .or.  ( frozen_states .and.
             froz_min <= eig(i) .and. eig(i) <= froz_max ) )  ! line 1372-1373
  → mark band i FROZEN

! 3. else DISENTANGLE:
else if ( projs(i) >= proj_min  .and.  projs(i) < proj_max )  ! line 1382
  → mark band i as non-frozen (disentanglement pool)

! 4. else (implicit): projs(i) < proj_min  → DISCARDED (dropped from window)
```

Three-way partition of the bands inside the outer energy window, by
projectability `p ≡ projs(i)`:

| region        | condition                                | fate                                   |
|---------------|------------------------------------------|----------------------------------------|
| **Frozen**    | `p ≥ dis_proj_max`                        | locked into subspace, always reproduced |
| **Disentangle** | `dis_proj_min ≤ p < dis_proj_max`       | Z-matrix optimization pool             |
| **Discarded** | `p < dis_proj_min`                        | removed from window entirely           |

Boundary conventions (get these exactly right — off-by-epsilon changes the count):
- Frozen lower bound **inclusive**: `p >= proj_max` (`>=`, line 1371).
- Disentangle band **`[proj_min, proj_max)`**: lower inclusive `>=`, upper
  strict `<` (line 1382). A state exactly at `proj_max` is FROZEN, not
  disentangled.
- Discard is the strict complement `p < proj_min`.
- Outer-window energy test is **strict-outside** (`<`/`>`, line 1368) → a state
  with `eig == win_min` or `eig == win_max` is KEPT. (Contrast the energy-path
  `dis_windows`, which uses `.ge./.le.` — same inclusive effect.)

**The frozen set is a UNION** (`.or.`, line 1371): a state is frozen if it has
high projectability **OR** falls in the energy inner window
(`dis_froz_min/max`, only when `dis_froz_max` was supplied so
`frozen_states == .true.`). Consequence a port must reproduce: a state with
`p < proj_min` (would otherwise be discarded) that lies in the energy frozen
window is **frozen, not discarded** — the first `if` catches it before the
`else if`. In the graphene test `frozen_states == .false.` (the
`dis_froz_max = 0.5` line is commented out), so the union degenerates to pure
projectability freezing.

### 1c. Counts and bookkeeping (`disentangle.F90:1363–1478`)

Counters accumulated in the same loop:
- `j` = kept states (frozen + disentangle) → `ndimwin(nkp) = j` (line 1391).
- `k` = frozen states → `ndimfroz(nkp) = k` (line 1392).
- `l` = disentangle states; require `j == k + l` else fatal (line 1394).
- `linner = .true.` if any k has `ndimfroz > 0` (line 1435) — this is what turns
  on `dis_proj_froz` downstream (identical trigger to the energy path).
- `indxkeep(1..j, nkp) = i` — original band index of the m-th kept state
  (line 1376, 1384). Used to slim arrays and to compact out discarded
  mid-window states.
- `indxfroz`, `indxnfroz` — filled with **original** band indices first
  (lines 1378, 1386), then **remapped to window-local (slimmed) indices** via
  `invindxkeep` (lines 1461–1478). After remap they index the compacted
  1..ndimwin arrays, exactly as the energy path expects.
- `dis_manifold%lwindow(i, nkp) = .true.` for kept bands, indexed by **original**
  band index (lines 1381, 1388) — used by checkpoint/symmetry.

Error guards (all fatal): `ndimwin == 0` (1402), `ndimwin < num_wann` (1412),
`ndimfroz > num_wann` (1424) — same as energy path but with projectability-aware
error messages suggesting to lower `dis_proj_min` / raise `dis_win_max`.

### 1d. `nfirstwin` and the array slimming — DIFFERENT from energy path

Because discarded states can sit **in the middle** of the outer window, PDWF sets
`nfirstwin(nkp) = 1` for all k (line 1357) and does the slimming with the general
`indxkeep` gather map **inside `dis_windows_proj` itself**, rather than the
`nfirstwin`-offset compaction the energy path uses. Specifically:

- **eigval slim** (lines 1442–1448): `eigval_opt(i,k) ← eigval_opt(indxkeep(i,k),k)`
  for `i=1..ndimwin`, then zero above `ndimwin`.
- **a_matrix slim** (lines 1449–1457):
  `a_matrix(i,j,k) ← a_matrix(indxkeep(i,k),j,k)`, zero rows above `ndimwin`.
- **M-matrix slim** (lines 1490–1504, in a **separate k-loop** so neighbor
  k2's `indxkeep` is available):
  `M(i,j,nn,k) ← M(indxkeep(i,k_global), indxkeep(j,k2), nn, k)`.

Because of this, `dis_main` **skips** the usual `internal_slim_m` and the
`lwindow` fill for the PDWF branch (`disentangle.F90:193–208`,
`if (.not. frozen_proj)`): the slimming and `lwindow` are already done inside
`dis_windows_proj`. A port must NOT double-slim.

**Trap for the Julia port:** the energy path's window-local index is
`global_band = nfirstwin + i − 1` (a contiguous offset). The PDWF path's
window-local index is `global_band = indxkeep(i)` (a **gather list**, possibly
non-contiguous). Use the gather list; do not assume contiguity.

--------------------------------------------------------------------------------
## 2. `auto_projections` — what it does and does NOT do

`auto_projections = .true.` (`select_proj%auto_projections`,
`wannier90_types.F90:270`, default `.false.`) is parsed in
`wannier90_readwrite.F90:1574–1585`. In a **wannier90.x post-processing run**
(reading `.amn/.mmn/.eig`) it has exactly TWO effects:

1. **Mutual exclusion with a `projections` block** (`:1581–1585`): giving both is
   a fatal input error. With `auto_projections`, there is **no** projection
   specification inside w90; `lhasproj` stays `.false.`, `num_proj = num_wann`
   (`:1596–1598` for the CLI path).

2. In the **`-pp` preprocessing** step it writes a `begin auto_projections … end
   auto_projections` block into the `.nnkp` file (`kmesh.F90:1092–1097`,
   emitting `num_proj` then `0`). This tells the *interface code*
   (pw2wannier90 / QE) to generate the `.amn` automatically — via SCDM or
   atomic-orbital projection.

**w90 does NOT compute SCDM, atomic projections, or any initial gauge.** There is
no SCDM machinery anywhere in `src/` (grep confirms: `scdm` appears in zero
files; `overlap.F90:335` simply `Read A_matrix from file wannier.amn`). The
initial gauge is **whatever `A` is in the `.amn` on disk**, orthonormalized by
the standard SVD/polar step in `dis_project` (`A → Z V^†`, see
disentanglement.md §3a). So for the graphene test:

- The test **ships `graphene.amn`** (header `20 9 8`: `num_bands=20`,
  `num_kpts=9`, `num_proj=8`; and `num_wann=8`, so **`num_proj == num_wann`** —
  required for `p_{nk} = Σ_w |A_{nw}|²` to be a valid ∈[0,1] probability and for
  `overlap.F90:361` not to reject it).
- `auto_projections = .true.` here is **just a flag asserting the `.amn` came
  from an automatic (pseudo-atomic-orbital) method**; w90 reads it verbatim.
  Nothing is computed by w90 from `auto_projections` in the wout run — the
  reproduction only needs the `.amn` file. **A Julia port replicating this test
  reads the shipped `.amn` and does NOT need to implement SCDM or
  `auto_projections` at all.** (If you want to *generate* PDWF `.amn` yourself,
  that is the SCDM/atomic-projection step done by the DFT interface, out of scope
  for w90 and for this note.)

--------------------------------------------------------------------------------
## 3. Interaction with the SMV Z-matrix disentanglement

**None beyond the window selection.** After `dis_windows_proj` returns, control
flow in `dis_main` (`disentangle.F90:151–231`) is *identical* to the energy path:

1. `dis_project` (line 166) — SVD/polar orthonormalize `a_matrix → u_matrix_opt`
   (SMV §III.D). Uses the **slimmed** `a_matrix` (window-local rows). Identical
   code; see disentanglement.md §3a. Output header confirms:
   `A_mn = <psi_m|g_n> → S = A.A^+ → U = S^-1/2.A`.
2. `dis_proj_froz` (line 180) — run iff `linner` (i.e. any frozen states). Locks
   the `ndimfroz` frozen bands and builds the complementary `num_wann − ndimfroz`
   subspace from the `Q P_s Q` diagonalization with the ortho-fix. **Identical
   code**; see disentanglement.md §3b. In the graphene run this IS triggered
   (`Using an inner window (linner = T)` in the wout) because high-projectability
   states are frozen.
3. `dis_extract` (line 226) — the SMV Z-matrix eigenvalue iteration over the
   **disentangle pool** (`indxnfroz` rows), holding the frozen `indxfroz` rows
   fixed, minimizing `Ω_I`. **Byte-for-byte identical** to the energy path;
   `ndimfroz`/`indxnfroz` are the only inputs that changed, and they now come
   from projectability rather than energy. See disentanglement.md §4 for the
   Z-matrix formula, mixing, convergence.
4. Handoff (`internal_find_u`, M-matrix recomputation, zeroing) — identical;
   disentanglement.md §6.

So **the entire Z-matrix machinery is reused unchanged**: PDWF is purely a
different rule for populating `ndimfroz`, `indxfroz`, `indxnfroz`, `ndimwin`,
and the slimmed `eigval_opt`/`a_matrix`/M. The gauge-invariant quantities
(`Ω_I`, `Ω_D`, `Ω_OD`, centres, interpolated bands) are then produced by exactly
the same code as any disentanglement + MV run.

--------------------------------------------------------------------------------
## 4. Keywords: exact names, defaults, file:line

Parsed in `w90_readwrite_read_dis_manifold` (`readwrite.F90:848–876`):

| keyword          | type    | default            | file:line        | meaning |
|------------------|---------|--------------------|------------------|---------|
| `dis_froz_proj`  | logical | `.false.`          | `readwrite.F90:848–850` → `dis_manifold%frozen_proj` | master switch: use projectability freezing (`dis_windows_proj` instead of `dis_windows`) |
| `dis_proj_min`   | real    | **`0.01`**         | `readwrite.F90:853–863` → `dis_manifold%proj_min` | below this `p` → discard; must be ∈[0,1] else fatal |
| `dis_proj_max`   | real    | **`0.95`**         | `readwrite.F90:853,864–872` → `dis_manifold%proj_max` | at/above this `p` → freeze; must be ∈[0,1] else fatal |
| `auto_projections` | logical | `.false.`        | `wannier90_readwrite.F90:1574` → `select_proj%auto_projections` | assert `.amn` is auto-generated; no projections block allowed |

Validation:
- `dis_proj_max < dis_proj_min` → fatal (`readwrite.F90:873–875`).
- `dis_proj_min`/`dis_proj_max` outside `[0,1]` → fatal (859–862, 867–871).
- `auto_projections` + `projections` block → fatal
  (`wannier90_readwrite.F90:1582–1585`).

**Note the DEFAULT `dis_proj_max = 0.95` differs from the graphene test's
`0.85`.** The disentanglement.md note (§2) cites the defaults 0.01/0.95; those
are correct — the test *overrides* `dis_proj_max = 0.85`.

Dispatch: `dis_main` calls `dis_windows_proj` iff
`dis_manifold%frozen_proj == .true.` (`disentangle.F90:151`), else the energy
`dis_windows`. `dis_froz_max`/`dis_froz_min` (energy inner window) are still read
and still contribute to the frozen UNION (§1b) even under `dis_froz_proj`.

### 4a. `.win` keywords in the graphene_pdwf test (all quoted verbatim)

`testw90_graphene_pdwf/graphene.win`:
```
num_wann = 8
num_bands = 20
auto_projections = .true.        ! pseudo-atomic-orbital projection (amn shipped)
dis_froz_proj = .true.           ! enable projectability disentanglement
dis_proj_max =   0.85            ! freeze  p >= 0.85
dis_proj_min =   0.01            ! discard p <  0.01; disentangle 0.01 <= p < 0.85
fermi_energy =  -2.3043
! dis_froz_max =  0.5            ! (COMMENTED OUT → frozen_states = .false.)
num_iter = 10                    ! MV wannierise iterations
dis_num_iter = 10                ! disentanglement Z-matrix iterations
mp_grid = 3 3 1                  ! → num_kpts = 9
bands_plot = .true.
```
No `dis_win_min/max` given → they default to ∓huge, then are set to
`minval/maxval(eigval)` (`library_interface.F90:576–577`), giving the
Outer window **[-21.91167, 15.67273] eV** seen in the wout. The `fermi_energy`
line does not affect disentanglement (only used for plotting/DOS).

--------------------------------------------------------------------------------
## 5. Exact output format (`dis_windows_proj` header) — for wout matching

The projectability header block, `disentangle.F90:1293–1331` (`f10.5` fields):
```
 +----------------------------------------------------------------------------+
 |                              Energy  Windows                               |
 |                              ---------------                               |
 |                   Outer:  -21.91167  to   15.67273  (eV)                   |
 |                   No frozen states were specified                          |
 |----------------------------------------------------------------------------|
 |                          Projectability  Windows                           |
 |                          -----------------------                           |
 |               Discarded:    0.00000  to    0.01000                         |
 |            Disentangled:    0.01000  to    1.00000                         |
 |                  Frozen:    0.85000  to    1.00000                         |
 +----------------------------------------------------------------------------+
   Number of target bands to extract:    8
```
- "Discarded" row prints `0.0` → `proj_min`; "Disentangled" prints
  `proj_min` → `1.0` (NOT `proj_max` — the printed range is cosmetic and
  overlaps the frozen range); "Frozen" prints `proj_max` → `1.0`.
- The "No frozen states were specified" line for the Energy block appears because
  `frozen_states == .false.` (`:1306–1308`); with `dis_froz_max` set it would
  print the Inner energy window instead.

Then the standard disentanglement iteration table and final `Ω_I` follow,
formatted identically to the energy path.

--------------------------------------------------------------------------------
## 6. Benchmark quantities, tolerances, and anchor values

Test: `testw90_graphene_pdwf`, driver `WANNIER90_WOUT_OK`
(`test-suite/tests/jobconfig:256–259`, output `graphene.wout`). Tolerances
(`test-suite/tests/userconfig:6–22`), format `(abs, rel, quantity)`:

| quantity          | abs tol | rel tol |
|-------------------|---------|---------|
| `final_centres_x` | 1.0e-5  | 1.0e-5  |
| `final_centres_y` | 1.0e-5  | 1.0e-5  |
| `final_centres_z` | 1.0e-5  | 1.0e-5  |
| `final_spreads`   | 3.0e-6  | 3.0e-6  |
| `omegaI`          | 1.0e-6  | 1.0e-6  |
| `omegaD`          | 1.0e-6  | 5.0e-6  |
| `omegaOD`         | 1.0e-6  | 1.0e-6  |
| `omegaTotal`      | 1.0e-6  | 1.0e-6  |
| `near_neigh_dist/mult`, `completeness_*` | 1.0e-6 | 1e-6/5e-6 |

**Anchor spread decomposition** (final state, `graphene.wout:487–490`,
`benchmark.out.default.inp=graphene.win`):
```
Omega I      =     7.962090079   (Ang^2)   ← gauge-invariant, fixed after disentangle
Omega D      =     0.048018588
Omega OD     =     7.793241243
Omega Total  =    15.803349910
```
`Final Omega_I 7.96209008` also appears at the end of the disentanglement
(`:297`). Note: `dis_num_iter = 10` is **not converged** — the wout carries
`Maximum number of disentanglement iterations reached` (`:294–295`), so `Ω_I` is
the value after exactly 10 Z-matrix iterations, not the fixed point. A port MUST
run exactly `dis_num_iter` iterations (no early exit) to hit `7.962090079`.

**Anchor final WF centres/spreads** (Å; `graphene.wout:477–485`):
```
 WF 1  (  0.000117,  1.421688, -0.000000 )   1.56254284
 WF 2  (  0.000000,  1.420282,  0.000000 )   0.95002823
 WF 3  ( -0.001555,  1.505457,  0.000000 )   3.47068055
 WF 4  ( -0.000143,  1.342765, -0.000000 )   1.99009911
 WF 5  (  1.230211,  0.708826, -0.000000 )   1.56042110
 WF 6  (  1.230000,  0.710141,  0.000000 )   0.95002823
 WF 7  (  1.233495,  0.641676, -0.000000 )   2.38026068
 WF 8  (  1.231300,  0.829338, -0.000000 )   2.93928886
 Sum   (  4.923426,  8.580174, -0.000000 )  15.80334958
```
**Initial-state** spreads/centres (`:313–321`, after disentangle+project, before
MV) for cross-checking the disentanglement handoff:
`Ω_total(initial) = 16.47594638`, `O_D = 0.0616041`, `O_OD = 8.4522525`
(`:324`); sum of centres `(4.922461, 8.537652, 0.000000)`.

**Disentanglement iteration table anchors** (`graphene.wout:283–292`), columns
`Iter  Omega_I(i-1)  Omega_I(i)  Delta(frac.)`:
```
 1   9.08771181   8.61701617   5.462E-02
 2   8.72831044   8.03144412   8.677E-02
 …
10   7.96729451   7.96209008   6.537E-04
```
These validate the Z-matrix iteration gauge-invariantly (`Ω_I(i)` per iter).

--------------------------------------------------------------------------------
## 7. Traps (silent-mismatch axes specific to PDWF)

1. **`p_{nk}` is the row-norm² of the raw `.amn` A**, `Σ_w |A_{nw}|²` =
   `(A A^†)_{nn}`, computed **before** SVD orthonormalization. Do NOT use singular
   values, `A^† A`, or the post-`dis_project` orthonormal U. The sum runs
   `w = 1..num_wann` (= `num_proj` here).
2. **Freeze boundary `p ≥ proj_max` inclusive; disentangle `[proj_min, proj_max)`
   upper-strict; discard `p < proj_min`.** A state exactly at `proj_max` is
   frozen. Off-by-epsilon here changes `ndimfroz` and the whole subspace.
3. **Frozen set is a UNION** of high-projectability and the energy inner window.
   Even a `p < proj_min` state is frozen (not discarded) if it lies in
   `[froz_min, froz_max]` and `dis_froz_max` was supplied. In the graphene test
   `frozen_states == .false.`, so only projectability freezes — but a general
   port must implement the union.
4. **Discarded states can be mid-window** → slimming uses a gather list
   `indxkeep` (possibly non-contiguous), NOT an `nfirstwin` offset.
   `nfirstwin ≡ 1` for all k in the PDWF path.
5. **PDWF slims `eigval_opt`, `a_matrix`, and M *inside* `dis_windows_proj`; the
   M-slim is a separate k-loop** (needs neighbor k2's `indxkeep`). `dis_main`
   then SKIPS `internal_slim_m` and the `lwindow` fill (`if .not. frozen_proj`).
   Do not double-slim.
6. **`indxfroz`/`indxnfroz` are remapped from original band indices to slimmed
   window-local indices** via `invindxkeep` (lines 1461–1478) at the END of the
   per-k loop. Downstream code (Z-matrix) expects window-local indices.
7. **`p_{nk} ∈ [0,1]` is enforced** (fatal otherwise). Requires
   `num_proj == num_wann` and quasi-orthonormal `.amn` columns. A `.amn` with
   `num_proj > num_wann` cannot be used for projectability freezing without
   `select_projections` (which reduces it to `num_wann` columns).
8. **w90 computes no SCDM / no initial gauge.** `auto_projections` is a flag +
   an `.nnkp` directive to the DFT interface. The reproduction reads the shipped
   `.amn` verbatim; the initial gauge is `SVD-orthonormalize(A_from_amn)`.
9. **`dis_num_iter = 10` is not converged** — the reference stops at the iter cap
   with a warning; run exactly 10 Z-matrix iterations (no `dis_conv_tol` early
   exit) to reproduce `Ω_I = 7.962090079`.
10. **Outer energy window defaults to `[min,max](eigval)`** when
    `dis_win_min/max` are absent (`library_interface.F90:576–577`), NOT ±∞ at the
    comparison. Here that is `[-21.91167, 15.67273]` eV; with PDWF and no energy
    window every band passes the energy gate and projectability alone selects.
