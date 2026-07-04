# Disentanglement (Souza–Marzari–Vanderbilt subspace selection)

Implementation-grade notes extracted from the Wannier90 v3.1.0 reference source.
All file:line citations are to
`reference/wannier90/src/disentangle.F90` unless another file is named.
"SMV" = Souza, Marzari, Vanderbilt, PRB 65, 035109 (2001). Equation numbers below
are the SMV equation numbers quoted in the reference comments.

The task is: given `num_bands > num_wann` Bloch states inside an outer energy
window at each k, pick a `num_wann`-dimensional subspace at every k that minimizes
the gauge-invariant spread `Omega_I` (equivalently, maximizes k-to-neighbor
subspace overlap). Output is a per-k rectangular matrix `u_matrix_opt`
(`num_bands × num_wann`) that embeds the optimal subspace into the window states,
plus a square `u_matrix` initial guess for the subsequent MV localization.

--------------------------------------------------------------------------------
## 0. Index / array conventions (critical for a Julia port)

- Fortran is 1-based and column-major. All quoted formulas use 1-based indices.
  Array shapes (as declared in `dis_main`, lines 71–78):
  - `eigval(num_bands, num_kpts)` — DFT eigenvalues, **eV**, read verbatim from
    the `.eig` file (`readwrite.F90:767`), no unit conversion.
  - `u_matrix(num_wann, num_wann, num_kpts)` — square rotation (output).
  - `u_matrix_opt(num_bands, num_wann, num_kpts)` — rectangular embedding (I/O).
    On *input* it carries the initial projections A (`a_matrix = u_matrix_opt`,
    line 120).
  - `m_matrix_orig_local(num_bands, num_bands, nntot, nkrank_local)` — the Mmn
    overlaps `<u_{m,k}|u_{n,k+b}>`. Note the **b-index (nn) comes before the
    local-k index**. This is the ONLY M-matrix used here; it is overwritten in
    place (slimmed, then finally replaced by the `num_wann × num_wann` overlaps).
- Eigenvalues are assumed **sorted ascending in band index** (`dis_windows`
  comment, line 998). A port must sort per-k before windowing.
- MPI: `nkrank`/`ranknk` = number of k on this rank; `global_k(nkp_loc)` maps
  local→global k index. For a serial port `ranknk = num_kpts`,
  `global_k(i) = i`. All `comms_allreduce(..., 'SUM')` reduce partial sums over
  k; serially they are no-ops. `womegai` etc. are summed over local k then
  allreduced, then divided by `num_kpts`.

--------------------------------------------------------------------------------
## 1. Quantity minimized: Omega_I

`Omega_I` (called `womegai` in code, public name `omega_invariant`) is computed
from subspace-to-neighbor overlaps. Per-k contribution (lines 2885–2906):

```
wkomegai = sum_nn  wb(nn) * sum_{m,n=1..num_wann} |cww(m,n)|^2
wkomegai = num_wann * wbtot  -  wkomegai            ! line 2904
womegai  = sum_k wkomegai                            ! line 2905
```
then (line 2910)
```
womegai = womegai / num_kpts
```
Here (lines 2891–2895):
```
cwb = u_matrix_opt(:,:,k)^H  ·  M(:,:,nn,k)          ! zgemm 'C','N'
cww = cwb  ·  u_matrix_opt(:,:,k2)                    ! zgemm 'N','N'
```
so `cww(m,n) = sum_{a,b} conj(U_opt(a,m,k)) M(a,b,nn,k) U_opt(b,n,k2)`,
= `<w_m,k | w_n,k+b>` between the two optimal subspaces. `k2 = nnlist(k,nn)`.

`rsum = sum_{m,n} |cww|^2` uses `real(cww)^2 + aimag(cww)^2` (line 2899) — i.e.
`|.|^2`, all `num_wann×num_wann` elements.

**This is SMV Eq. (12):**
`Omega_I = (1/N_k) Σ_k Σ_b w_b [ N_wann − Σ_{m,n} |<u_mk|u_n,k+b>|^2 ]`.

Key numerical conventions:
- **Normalization is 1/N_k (`num_kpts`), NOT 1/(N_k N_b).** The b-sum is weighted
  by `wb(nn)`; `wbtot = Σ_nn wb(nn)` (kmesh.F90:668–674). The term
  `num_wann * wbtot` is the "full" value; overlap reduces it.
- `wb` weights have units **1/length²** (Å⁻² by default). `Omega_I` therefore has
  units length² (Å²). Output multiplies by `lenconfac**2` (line 2916, 3036),
  `lenconfac = 1` when length_unit = Ang.
- b-vectors: the shell weights `wb(nn)` and neighbor list `nnlist`, `nntot` come
  from `kmesh` and MUST match the b-vector ordering used to build M. Any
  reordering of shells silently changes `Omega_I`.

`womegai1` (`womegai_(i-1)` in the printed table) is the SAME quantity but with
the current-iteration subspace at k paired against the **previous** iteration's
subspace at neighbors; at self-consistency `womegai1 == womegai`. It is computed
incrementally from Z-matrix eigenvalues (Section 3) rather than by brute force.

--------------------------------------------------------------------------------
## 2. Energy windows and the frozen (inner) window — `dis_windows` (886–1207)

Input keys (all **eV**), parsed in `readwrite.F90:809–840`:
- `dis_win_min` / `dis_win_max` → `dis_manifold%win_min/max`
  (default win_min=−huge, win_max=+huge; `types.F90:213–215`). Outer window.
- `dis_froz_min` / `dis_froz_max` → `dis_manifold%froz_min/max`
  (default ∓huge). Supplying `dis_froz_max` sets `frozen_states=.true.`
  (line 827). It is an **error** to give `dis_froz_min` without `dis_froz_max`
  (line 837). `froz_max < froz_min` is an error (833).

Per k-point (loop `do nkp`, line 980):

Outer window selection (lines 999–1011), inclusive comparisons with `.ge./.le.`:
- `imin` = first band with `win_min ≤ eig ≤ win_max`.
- `imax` = last band with `eig ≤ win_max`.
- `ndimwin(nkp) = imax − imin + 1`, `nfirstwin(nkp) = imin`.
- Error if window empty (982) or `ndimwin < num_wann` (1037).

Frozen window (lines 1048–1071), only if `frozen_states`:
- Scan `i = imin..imax`; `kifroz_min/max` are indices **relative to bottom of
  outer window** (`i − imin + 1`), for bands with `froz_min ≤ eig ≤ froz_max`
  (inclusive).
- `ndimfroz(nkp) = kifroz_max − kifroz_min + 1` (0 if none: init
  `kifroz_min=0, kifroz_max=−1`, line 1049–1050).
- Error if `ndimfroz > num_wann` (line 1073).
- `linner = .true.` if any k has `ndimfroz > 0` (line 1084).

Index bookkeeping:
- `lfrozen(i,nkp)` true for frozen bands (line 1094), indexed 1..ndimwin.
- `indxfroz(i,nkp) = kifroz_min + i − 1` — window-relative index of i-th frozen
  band (line 1093).
- `indxnfroz(i,nkp)` — window-relative index of i-th NON-frozen band (1109–1116).
  Must satisfy `count = ndimwin − ndimfroz` (checked line 1118).

**eigval slimming (lines 1126–1134):** eigenvalues are compacted so that band
`i=1..ndimwin` holds `eigval(nfirstwin+i−1)`; entries above ndimwin set to 0.
After this, `eigval_opt` index 1..ndimwin is *window-local*. A port must keep the
same convention because `u_matrix_opt` rows are window-local from here on.
`internal_slim_m` (488–564) does the analogous row/column compaction of M so that
its indices are window-local (`M(i,j) ← M(nfirstwin_k+i−1, nfirstwin_k2+j−1)`).

`dis_manifold%lwindow(band, k)` is set true for the ndimwin bands (dis_main
lines 202–207) — used by symmetry code and checkpoint.

**Frozen-by-projectability variant** `dis_windows_proj` (1209–1598, key
`dis_froz_proj`, defaults proj_min=0.01, proj_max=0.95, `readwrite.F90:848–854`):
selects the frozen set by projectability `Σ_n |A(b,n)|^2` thresholds instead of
energy. Same downstream bookkeeping. Off by default — document but a first port
can implement only the energy path.

--------------------------------------------------------------------------------
## 3. Initial subspace: projection + frozen locking

### 3a. `dis_project` (1600–1856) — unitarized projection (SMV Sec. III.D)
For each k, `a_matrix` (the A_mn projections, window-slimmed at lines 1730–1740
to `a_matrix(i,j) = A(nfirstwin+i−1, j)`) is SVD'd:
```
call zgesvd('A','A', ndimwin, num_wann, A, ..., svals, cz, ..., cvdag, ...)   ! 1744
```
LAPACK returns `A = cz · Σ · cvdag` with `cvdag = V^H`. Then (lines 1776–1780):
```
U_opt(i,j) = Σ_{l=1..num_wann} cz(i,l) · cvdag(l,j)
```
i.e. **`U_opt = Z · V^H` (drop the singular values)** = the SVD/polar
orthonormalization `A (A^H A)^{-1/2}`. This is **SVD-based orthonormalization, NOT
explicit Löwdin `S^{-1/2}A`** (though mathematically equal). The inner sum runs
only `l=1..num_wann` (rectangular). Columns of U_opt are orthonormal:
`U_opt^H U_opt = I_{num_wann}` (checked to `eps5`, lines 1793–1820). Note:
`U_opt U_opt^H ≠ I`.

### 3b. `dis_proj_froz` (1858–2318) — lock frozen states (SMV Eq. 27, Sec. III.G)
Only run if `linner` (dis_main line 173). For each k with `num_wann > ndimfroz`:
- Build `P_s = Σ_{l=1..num_wann} U_opt(:,l) U_opt(:,l)^H` (projector onto
  projected-gaussian subspace, lines 2021–2028).
- Build `Q_froz = 1 − P_froz`: diagonal, `Q_froz(n,n)=1` iff band n NOT frozen
  (line 2027).
- Form `CQPQ = Q_froz · P_s · Q_froz` (lines 2030–2046), Hermitian.
- Diagonalize with `zhpevx('V','A','U', ndimwin, ...)` (line 2070). Take the
  `num_wann − ndimfroz` **leading (largest-eigenvalue) eigenvectors**
  (`il = ndimwin − (num_wann−ndimfroz) + 1 .. iu = ndimwin`, lines 2068–2069),
  placed in columns `ndimfroz+1 .. num_wann` of U_opt (lines 2217–2220).
- Eigenvalues must lie in [−eps8, 1+eps8] (line 2116).
- **Ortho-fix (default on, lines 2138–2227):** if a required eigenvalue is
  `< eps8` (degenerate with frozen states, sign-ambiguous), instead of blindly
  taking `il..iu`, it re-selects eigenvectors by checking orthogonality
  (`|<frozen|v>| ≤ eps8`) to the frozen states. **Must be reproduced** — omitting
  it can silently pull a frozen state into U_opt and fail later orthonormality.
- Finally the frozen columns 1..ndimfroz of U_opt are set to unit vectors on the
  frozen bands: `U_opt(:,l)=0; U_opt(indxfroz(l),l)=1` (lines 2247–2253).

--------------------------------------------------------------------------------
## 4. Z-matrix / eigenvalue subspace iteration — `dis_extract` (2320–3295)

### 4a. Z-matrix `internal_zmatrix` (3353–3419) — SMV Eq. (21)
For a k with non-frozen states (`num_wann > ndimfroz`), over neighbors nn:
```
cbw = M(:,:,nn,k) · u_matrix_opt(:,:,k2)             ! zgemm, line 3398; k2=nnlist(k,nn)
for n=1..ndimk (ndimk = ndimwin−ndimfroz):  q = indxnfroz(n,k)
  for m=1..n:                                p = indxnfroz(m,k)
     csum = Σ_{l=1..num_wann} cbw(p,l) · conj(cbw(q,l))     ! line 3407
     Z(m,n) += wb(nn) * csum                                ! line 3409
     Z(n,m)  = conj(Z(m,n))                                 ! line 3410 (Hermitian)
```
So `Z_{mn} = Σ_b w_b Σ_l [M U_opt(k+b)]_{p,l} conj([M U_opt(k+b)]_{q,l})`,
restricted to the **non-frozen** window rows p,q. This is
`Σ_b w_b (P_{k+b})` projected into the non-frozen subspace, where
`P_{k+b} = Σ_l |M U_opt(k+b)>_l <...|`. **Complex-conjugate is on the second
factor `cbw(q,l)`; the frozen bands are excluded from Z rows/cols but the
`l`-sum over neighbor subspace runs full `1..num_wann`.**

### 4b. Iteration loop (2605–2961)
```
for iter = 1 .. dis_control%num_iter:        ! dis_num_iter, default 200
```
- **iter 1:** build `czmat_in` from initial U_opt via `internal_zmatrix`
  (2611–2620).
- **iter > 1:** mixing (2646–2657):
  ```
  Z_in(j,i) = mix_ratio*Z_out(j,i) + (1−mix_ratio)*Z_in(j,i)   ! upper triangle j≤i
  Z_in(i,j) = conj(Z_in(j,i))                                   ! keep Hermitian
  ```
  `mix_ratio = dis_control%mix_ratio` (dis_mix_ratio, default 0.5). Mixing is on
  the **Z matrix**, not on U.

- **wkomegai1 (SMV Eq. 18)** initialized to `num_wann * wbtot` (line 2669).
  Frozen contribution subtracted first, for all k, BEFORE any k is updated
  (comment 2666–2668), so neighbor overlaps use previous-iteration non-frozen
  states:
  ```
  cww = U_opt(k)^H M(nn,k) U_opt(k2)   restricted to ndimfroz(k) rows   ! 2691–2695
  wkomegai1(k) -= wb(nn) * Σ_{m=1..ndimfroz,n=1..num_wann} |cww(m,n)|^2  ! 2696–2702
  ```
- **Diagonalize Z** (non-frozen k), lines 2732–2739:
  ```
  pack Z_in upper triangle → cap;  ZHPEVX('V','A','U', ndiff, ...)   ! ndiff=ndimwin−ndimfroz
  ```
  ZHPEVX returns eigenvalues `w` **ascending**. Take the `num_wann − ndimfroz`
  **largest** eigenvectors: loop `j = ndimwin−num_wann+1 .. ndimwin−ndimfroz`
  (line 2759). For each, place eigenvector into U_opt column `m` (m runs
  `ndimfroz+1..num_wann`), scattering back to full window rows via
  `U_opt(indxnfroz(i,k), m) = cz(i,j)` (lines 2762–2767) and
  ```
  wkomegai1(k) -= w(j)      ! largest eigenvalues, line 2761
  ```
- `womegai1 += wkomegai1(k)` (line 2775), allreduce, `/num_kpts` (2822,2845).
- Then recompute `womegai` (Section 1) from the UPDATED subspaces (2885–2910).
- `delta_womegai = womegai1/womegai − 1` (line 2911) — **fractional**, signed.
- Build `czmat_out` from updated U_opt (2922–2931) → feeds next-iter mixing.

### 4c. Convergence `internal_test_convergence` (3297–3351)
- `dis_conv_window` (default 3) rolling history of `delta_womegai`
  (`eoshift`, line 3334).
- Converged when `iter ≥ conv_window` AND `all(|history| < dis_conv_tol)`
  (default `dis_conv_tol = 1e-10`), line 3339–3340. Note **|delta|** — the
  fractional change over ALL of the last `conv_window` iterations must be below
  tol. If never met, loop runs full `num_iter` with a warning (3000–3006).

--------------------------------------------------------------------------------
## 5. Post-iteration: diagonalize H in the optimal subspace (SMV Sec. III.E)

After convergence (root only, 3043–3117; then broadcast):
```
cham(i,j,k) = Σ_{l=1..ndimwin} conj(U_opt(l,i,k)) U_opt(l,j,k) eigval_opt(l,k)  ! 3050–3053
```
i.e. `H_sub = U_opt^H diag(eig) U_opt` (num_wann × num_wann, eV). Diagonalize
with ZHPEVX (3063). Store eigenvalues into `eigval_opt(1:num_wann,k)` (line 3080).
Rotate U_opt by the eigenvectors:
```
ceamp(i,j,k) = Σ_{l=1..num_wann} cz(l,j) U_opt(i,l,k)          ! 3084–3091
U_opt(:,j,k) = ceamp(:,j,k)                                     ! 3106–3111
```
So on output, U_opt columns are the eigenstates of the subspace Hamiltonian
(convenient for interpolation), **unless** `lsitesymmetry` (then this rotation is
skipped, line 3106). `omega_invariant = womegai` (line 3040) is stored and must
stay fixed through the later MV localization (sanity-checked in checkpoint).

--------------------------------------------------------------------------------
## 6. Handoff to MV localization

Back in `dis_main` (247–280):
1. Recompute the `num_wann × num_wann` neighbor overlaps from the final U_opt and
   overwrite the top-left block of `m_matrix_orig_local`:
   ```
   cwb = U_opt(k)^H · M(nn,k)                          ! zgemm 'C','N', 251
   cww = cwb · U_opt(k2)                                ! zgemm 'N','N', 254
   M(1:num_wann,1:num_wann,nn,k) = cww                  ! 256
   ```
   These are the M-matrices the wannierise step consumes.
2. `internal_find_u` (566–733) builds the initial square `u_matrix` (SMV Sec.
   III.D, square case):
   ```
   caa = U_opt(k)^H · A(k)                              ! num_wann×num_wann, 665
   zgesvd(caa) → cz, cv (=V^H)                          ! 668
   u_matrix(:,:,k) = cz · cv                            ! = Z V^H, polar factor, 681
   ```
   i.e. the SVD orthonormalization of `<psi_tilde | g>`.
3. Unused rows of U_opt (rows > ndimwin) are zeroed (275–279).

The wannierise (MV) step then minimizes `Omega_tilde = Omega_D + Omega_OD` over
the square `u_matrix`, holding U_opt (and hence `Omega_I`) fixed.

--------------------------------------------------------------------------------
## 7. Gotchas that cause silent numerical mismatch

1. **Normalization 1/N_k only** (not 1/(N_k N_b)); b-weighting via `wb(nn)`,
   `wbtot`. Units Å² by default.
2. **b-vector / shell ordering** must match the M-matrix build and the `wb`,
   `nnlist`, `nntot` from kmesh. Reordering shells changes results silently.
3. **Window-local indexing** after `dis_windows`/`slim_m`: `eigval_opt`,
   `u_matrix_opt` rows, and M rows/cols are compacted to 1..ndimwin. `indxfroz`,
   `indxnfroz` are window-local indices; `nfirstwin` maps back to global bands.
4. **Complex-conjugate placement:** in Z-matrix `csum += cbw(p,l)·conj(cbw(q,l))`
   (row p conjugated on the SECOND factor); in overlaps `cww` the bra U_opt(k) is
   conjugated (zgemm 'C'). `|cww|^2 = real^2 + aimag^2`.
5. **Eigenvector selection = LARGEST eigenvalues.** ZHPEVX returns ascending; the
   leading `num_wann−ndimfroz` are the TOP indices. Off-by-one here inverts the
   algorithm.
6. **Orthonormalization is SVD/polar (`Z·V^H`), not literal `S^{-1/2}·A`.**
   Numerically equivalent but a naive `S^{-1/2}` can differ in degenerate/near-
   singular cases; use SVD to match the reference bit-for-bit-ish.
7. **Ortho-fix in `dis_proj_froz`** (default on) for near-zero QPQ eigenvalues —
   must be reproduced when frozen window is used.
8. **Mixing is on Z (`mix_ratio·Z_out + (1−mix_ratio)·Z_in`), Hermitized after
   mixing.** iter 1 has no mixing. Default mix 0.5.
9. **Convergence tests |delta| over a full window of `conv_window` iters**, delta
   is the signed fractional `womegai1/womegai − 1`. Default tol 1e-10, window 3,
   num_iter 200.
10. **Units:** eigenvalues eV throughout (from `.eig`, no conversion). `wb` in
    Å⁻²; `Omega_I` in Å² (multiply by `lenconfac**2` only for printing).
    Lattice handling (Bohr↔Å) is confined to kmesh; not in this module.
11. **Frozen contribution to wkomegai1 computed before any k is updated**, so
    neighbor states are from the previous iteration — ordering matters for the
    printed `Omega_I(i-1)`.
12. **num_wann == ndimfroz(k) at some k:** Z-matrix/diagonalization skipped there
    (`if num_wann > ndimfroz` guards), the frozen states alone fill the subspace.
