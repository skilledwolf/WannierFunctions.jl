# Maximal Localization (Marzari–Vanderbilt gauge optimization)

Reference: Wannier90 v3.1.0, `src/wannierise.F90`. All line numbers below are for
that file unless another file is named. This note covers the
**`num_bands == num_wann`, `use_ss_functional = .false.`, `selective_loc = .false.`**
branch — i.e. the standard MV localization. The Stengel–Spalding (`use_ss_functional`,
`nnord`/`nnrev` reordering) and SLWF+C selective-localization / constrained-centre
(`selective_loc`, `constrain`, `lambda`) code paths are read below only to isolate them
as **out of scope**; do NOT port those `cdodq` / `om_*` blocks.

Driver subroutine: `wann_main` (line 61). Key inner subroutines:
`wann_omega` (1994, spread), `wann_domega` (2400, gradient),
`internal_search_direction` (1378, CG), `internal_optimal_step` (1543, line search),
`internal_new_u_and_m` (1614, unitary update), `internal_test_convergence` (1084),
`wann_phases` (1776, guiding centres).

---

## 0. Units, conventions, index base (CRITICAL — silent-mismatch axes)

- **Length unit is Ångström internally.** `real_lattice` is stored in Å.
  `recip_lattice = 2π·inv(real_lattice)` (`utility.F90:346-352`,
  `recip_lat = twopi*recip_lat/volume`), so **`bk` (b-vectors) are in Å⁻¹** and
  **`wb` (shell weights) are in Å²**. Consequently the spread `Omega` and all
  `om_*` are in **Å²** internally.
- `print_output%lenconfac` (`types.F90:60`) defaults to **1.0** (`length_unit='ang'`,
  `types.F90:58`). It is set to `1.0_dp/bohr ≈ 1.889726` **only** when the user requests
  `length_unit = bohr` (`readwrite.F90:181`). `lenconfac` is a *display-only* multiplier
  applied to centres (×lenconfac) and spreads (×lenconfac²) in write statements — it never
  enters the math. The "Distance (Bohr⁻¹)" header at `kmesh.F90:209` is the `else` branch
  of an `if length_unit=='Ang'` conditional (`kmesh.F90:205`); the Ang branch (:206) prints
  "Ang⁻¹". So do NOT conclude the internal unit is Bohr.
- **Energies are irrelevant here** — localization is purely geometric (M-matrix overlaps).
  No Ha/eV appears.
- **Fortran arrays are 1-based.** All loops `do n = 1, num_wann`, `do nn = 1, nntot`,
  `do nkp = 1, num_kpts`, `do ind = 1, 3`. A Julia port keeps 1-based indexing naturally.
- `M^{k,b}_{mn}` is stored as `m_matrix(m, n, nn, nkp)` (band m = row, band n = column;
  neighbour `nn`, k-point `nkp`). In the MPI code it is `m_matrix_loc(:,:,nn,nkp_loc)`;
  for a serial port treat `nkp_loc == nkp`, `nkrank == num_kpts`, `global_k(nkp)=nkp`,
  and drop all `comms_allreduce` calls (they only sum partial k-sums across ranks).
- `N_k = num_kpts`, `N_b = nntot` (total neighbours per k), `wbtot = Σ_nn wb(nn)`
  (`kmesh.F90:668-676`).
- The shell weights satisfy the **B1 condition** `Σ_b w_b b_α b_β = δ_αβ`
  (`kmesh.F90:256`, "Get the shell weights to satisfy the B1 condition"). This is assumed
  throughout the spread/gradient algebra.

---

## 1. Branch cut of Im ln — `sheet` / `csheet`

The imaginary part of a complex log is used for the centres and the spread. To keep a
consistent branch, W90 carries per-(n,nn,k) real `sheet(n,nn,nkp)` and its phase
`csheet = exp(i·sheet)`.

- **Default initialization (no guiding centres):** `csheet = cmplx_1`, `sheet = 0.0`
  (line 310). Then `aimag(log(csheet·M_nn)) = aimag(log(M_nn))` = principal-value
  `atan2(Im, Re) ∈ (-π, π]`. **A from-scratch port that starts from unitary U close to a
  reasonable projection will match this as long as it uses the principal branch and
  `csheet ≡ 1`, `sheet ≡ 0`.**
- `q_n^{k,b}` (see §3) is defined as `Im ln(csheet·M_nn) - sheet`. With defaults this is
  just `Im ln M_nn` (principal value).
- Guiding centres (§7) modify `sheet`/`csheet` to shift the branch consistently; only
  active if `guiding_centres%enable` (default `.false.`, `wannier90_types.F90:178`).

---

## 2. Wannier centres and spread — `wann_omega` (line 1994)

### 2a. `ln_tmp` in wann_omega (NO wb factor)

Line 2116-2117 (the code explicitly warns this differs from wann_domega):
```
ln_tmp_loc(n,nn,nkp) = aimag(log(csheet(n,nn,nkp)*m_matrix_loc(n,n,nn,nkp))) - sheet(n,nn,nkp)
```
So in **wann_omega**, `ln_tmp(n,nn,k) = q_n^{k,b} = Im ln(csheet·M_nn) − sheet`
(principal-value `Im ln M_nn` with defaults). **No `w_b`.**

### 2b. Wannier centres `rave = r_n` (lines 2122-2138)

```
rave(ind,iw) = Σ_k Σ_nn  wb(nn) * bk(ind,nn,nkp) * ln_tmp(iw,nn,nkp)     ! accumulate
rave = -rave / num_kpts                                                   ! sign + 1/N_k
```
i.e.
```
r_n = -(1/N_k) Σ_{k,b} w_b · b · Im ln M^{k,b}_{nn}
```
**Sign convention: leading minus. Normalization: 1/N_k. `w_b` present.**
`bk(ind,nn,nkp)` is the Cartesian b-vector (Å⁻¹) for neighbour `nn` at k-point `nkp`.

### 2c. `rave2(iw) = |r_n|²` (lines 2140-2143)
```
rave2(iw) = Σ_ind rave(ind,iw)²  =  |r_n|²
```

### 2d. `r2ave(iw) = <r²>_n` (lines 2145-2160)
```
mnn2 = Re( M_nn · conj(M_nn) ) = |M^{k,b}_nn|²
r2ave(iw) = Σ_k Σ_nn wb(nn) * ( 1 - |M_nn|² + ln_tmp(iw,nn,nkp)² )
r2ave = r2ave / num_kpts
```
i.e.
```
<r²>_n = (1/N_k) Σ_{k,b} w_b [ 1 − |M^{k,b}_nn|² + (Im ln M^{k,b}_nn)² ]
```
**The per-WF spread reported** is `spread_n = r2ave(n) − rave2(n) = <r²>_n − |r_n|²`
(`wannier_data%spreads = r2ave - rave2`, line 482/773).

### 2e. Spread decomposition `Omega = Omega_I + Omega_D + Omega_OD`

**Omega_I (`om_i`, invariant)** — computed only on the FIRST pass (`first_pass`), then cached
in `omega%invariant` and reused every subsequent call (lines 2293-2318). It is gauge-invariant
so recomputing is wasteful:
```
om_i = (1/N_k) Σ_k Σ_nn wb(nn) * ( num_wann − Σ_m Σ_n |M^{k,b}_{nm}|² )
```
(lines 2299-2314; `summ = Σ_{m,n} Re(M_{nm}·conj(M_{nm}))`, then
`om_i += wb·(num_wann − summ)`, `om_i /= num_kpts`).
NOTE index order in the source: `m_matrix_loc(n,m,...)` summed over both m,n, so it is
`Σ_{mn}|M_{nm}|²` = full Frobenius norm² of the M block (order irrelevant for the sum).

**Omega_OD (`om_od`, off-diagonal)** (lines 2320-2336):
```
om_od = (1/N_k) Σ_k Σ_nn wb(nn) * Σ_{m≠n} |M^{k,b}_{nm}|²
```
(loop skips `m==n`; `m_matrix_loc(n,m,...)` with the m≠n guard).

**Omega_D (`om_d`, diagonal)** (lines 2372-2388):
```
brn = Σ_ind bk(ind,nn,nkp) * rave(ind,n)          ! = b · r_n
om_d = (1/N_k) Σ_k Σ_nn wb(nn) * ( ln_tmp(n,nn,nkp) + brn )²
     = (1/N_k) Σ_{k,b} w_b ( Im ln M^{k,b}_nn + b·r_n )²
```

**Total** (line 2391):
```
om_tot = om_i + om_d + om_od
```
Public exports (lines 471-473 / 777-778):
`omega%total = om_tot`, `omega%invariant = om_i`, `omega%tilde = om_d + om_od`.

---

## 3. Gradient G = dΩ/dW — `wann_domega` (line 2400)

Result is `cdodq(m,n,nkp)` = the gradient matrix (anti-Hermitian) at each k.

### 3a. `ln_tmp` in wann_domega (WITH wb factor — DIFFERENT from wann_omega!)

Lines 2580-2581 (code warns: "this ln_tmp is defined differently wrt the one in wann_omega"):
```
ln_tmp_loc(n,nn,nkp) = wb(nn) * ( aimag(log(csheet·M_nn)) − sheet )
                     = w_b · q_n^{k,b}
```
**So here `ln_tmp` carries an extra `w_b`.** This is the single most dangerous porting
trap: the same-named quantity is `q_n` in wann_omega (§2a) and `w_b·q_n` in wann_domega.

### 3b. Recompute rave (lines 2587-2599)
`rave` is recomputed here from the wb-weighted `ln_tmp`, so the wb cancels the missing
wb consistently:
```
rave(ind,iw) = Σ_k Σ_nn bk(ind,nn,nkp) * ln_tmp(iw,nn,nkp)   ! ln_tmp already has wb
rave = -rave / num_kpts
```
Same `r_n` as §2b (leading minus, 1/N_k).

### 3c. `rnkb = b · r_n` (lines 2618-2626)
```
rnkb_loc(n,nn,nkp) = Σ_ind bk(ind,nn,nkp) * rave(ind,n)   ! = b · r_n
```

### 3d. R and T matrices, then G (the `selective_loc=.false.` block, lines 2632-2739)

For each k, nn, build per-column-n auxiliary matrices (lines 2635-2639):
```
mnn      = M^{k,b}_{nn}                       (diagonal element, column n)
crt(:,n) = M^{k,b}_{:,n} / mnn                = R̃^{k,b}  (Rt, "tilde")
cr(:,n)  = M^{k,b}_{:,n} * conj(mnn)          = R^{k,b}
```
So `R_{mn} = M_{mn} · conj(M_nn)` and `R̃_{mn} = M_{mn} / M_nn` (MV notation).

Gradient accumulation (lines 2720-2736), the standard MV formula
`G = 4 Σ_b w_b ( A[R] − S[T] )` with `T_{mn} = R̃_{mn} · q_n`,
`q_n = Im ln M_nn + b·r_n`:
```
! A[R] = (R − R†)/2
cdodq(m,n) += wb(nn) * 0.5 * ( cr(m,n) − conj(cr(n,m)) )

! −S[T],  S[T]=(T+T†)/2i , here written with cmplx(0,-0.5)=-1/(2i)
! T split into two pieces: q_n = ln_tmp/wb (already ×wb) + rnkb
cdodq(m,n) −= ( crt(m,n)*ln_tmp(n,nn) + conj(crt(n,m)*ln_tmp(m,nn)) ) * cmplx(0,-0.5)
cdodq(m,n) −= wb(nn) * ( crt(m,n)*rnkb(n,nn) + conj(crt(n,m)*rnkb(m,nn)) ) * cmplx(0,-0.5)
```
Notes on factors that must be reproduced exactly:
- The **first −S[T] term uses `ln_tmp` which already contains `w_b`** (§3a), so it is NOT
  multiplied by `wb(nn)` again.
- The **second −S[T] term (the `rnkb = b·r_n` part) IS multiplied by `wb(nn)`** explicitly,
  because `rnkb` has no wb.
- `cmplx(0.0, -0.5)` = `-1/(2i)`. So `−S[T] = −(T+T†)/(2i)` is encoded as
  `−(T+T†)·cmplx(0,-0.5) = +(T+T†)·(1/(2i))` — check sign carefully:
  the term written is `− (Rt·q_n + conj(Rt·q_m)) · cmplx(0,-0.5)`. With
  `A[R]=(R−R†)/2` and `S[T]=(T+T†)/2i`, the total spatial gradient sign convention is
  `G_mn = 4 Σ_b w_b ( A[R]_mn − S[T]_mn )`.

### 3e. Global prefactor (line 2740)
```
cdodq_loc = cdodq_loc / num_kpts * 4.0
```
**So the full gradient is `G = (4/N_k) Σ_{k,b} w_b ( A[R] − S[T] )`.** The `×4` and `/N_k`
are applied once, globally, at the very end (NOT inside the nn-loop). `w_b` is applied
per-term inside (partly via the wb-carrying `ln_tmp`, partly explicit).

`cdodq` is **anti-Hermitian** in the (m,n) WF indices at each k (used as the generator of
the unitary rotation). If `lsitesymmetry`, it is symmetrized (`sitesym_symmetrize_gradient`,
lines 2752-2760) — out of scope for a first port.

---

## 4. Search direction (CG) — `internal_search_direction` (line 1378)

Fletcher–Reeves conjugate gradient over the anti-Hermitian gradient `cdodq`.

- `gcnorm1 = Re Tr[ G† · G ] = Σ_k Σ_{mn} |cdodq(m,n,k)|²` (via `zgemv 'c'`, line 1439;
  `zdotc` of the flattened array). NB gradient is anti-Hermitian so this is `Tr[G·G]` up to
  sign but the code takes the real part of the conjugate dot product.
- CG coefficient (lines 1446-1465):
  ```
  if (iter == 1 .or. ncg >= num_cg_steps):  gcfac = 0   (steepest descent), ncg = 0
  else if gcnorm0 > eps:                    gcfac = gcnorm1 / gcnorm0   (Fletcher-Reeves)
       if gcfac > 3.0:                       gcfac = 0, ncg = 0  (reset — too large)
       else                                  ncg = ncg + 1
  ```
  `gcnorm0` is `gcnorm1` from the previous iter (saved line 1468).
- Search direction (line 1474): `cdq = cdodq + gcfac · cdqkeep` where `cdqkeep` is the
  previous search direction (saved at line 607).
- Slope along search direction (lines 1488-1494):
  ```
  doda0 = − Re Tr[ G† · cdq ] / (4 · wbtot)
  ```
  (`zdotc(cdodq, cdq)`, negated, allreduced, then `/(4·wbtot)`).
- Uphill guard (lines 1497-1532): if `doda0 > 0` reset CG to steepest descent (or reverse
  the direction and flip `doda0`).
- **Random noise** (`internal_random_noise`, line 1175): only if `lrandom` (driven by
  `conv_noise_amp > 0`, default off). Adds an anti-Hermitian random matrix scaled by
  `conv_noise_amp` to `cdq`. Uses `random_seed()`/`random_number()` — **not reproducible**;
  only enabled when the user sets `conv_noise_amp`.

`num_cg_steps` default = **5** (`wannier90_types.F90:194`).

---

## 5. Line search / step length — `internal_optimal_step` (line 1543)

Parabolic (quadratic) line search using spread at the current point and at a trial step.

**Trial step taken first** (line 618, driver):
```
cdq = cdqkeep * ( trial_step / (4·wbtot) )
```
Then `internal_new_u_and_m` rotates U and M by this trial `cdq`, and `wann_omega` gives
`trial_spread`. Parabola fit (lines 1574-1595):
```
fac   = 1 / (trial_spread%om_tot − wann_spread%om_tot)      (or 1e6 fallback)
shift = 1
eqb   = fac · doda0
eqa   = shift − eqb · trial_step
alphamin  = −0.5 · eqb/eqa · trial_step²
falphamin = wann_spread%om_tot − 0.25 · eqb²/(fac·eqa) · trial_step²
```
Stability guards fall back to `alphamin = trial_step`, `falphamin = trial_spread%om_tot`
(`lquad=.false.`) if the parabola is unstable or predicts an uphill step (lines 1589-1603).

**Optimal step taken** (line 678, driver), after restoring the original U0/M0 (lines 681-688):
```
cdq = cdqkeep * ( alphamin / (4·wbtot) )
```
then `internal_new_u_and_m` again.

- `trial_step` (α trial) default = **2.0** (`wannier90_types.F90:201`).
- Fixed-step mode: if `lfixstep` (i.e. user set `fixed_step`, default `-999.0`), skip the
  line search and use `alphamin = fixed_step` directly (lines 610-612). `lquad` forced
  `.false.` (line 484).

**Net effective step:** the generator actually exponentiated is
`(alphamin / (4·wbtot)) · cdq`, and `cdq` = search direction built from the gradient
`G = (4/N_k)Σ w_b(...)`. The `4·wbtot` in the denominator and the `4` in the gradient are
the MV `Δ W = (α / 4Σw_b) · G` normalization — reproduce both.

---

## 6. Unitary update — `internal_new_u_and_m` (line 1614)

Given anti-Hermitian generator `cdq(:,:,k)` (already scaled by `α/(4·wbtot)`):

1. Form Hermitian `tmp_cdq = i · cdq` (line 1677).
2. Diagonalize with `zheev('V','U', ...)` → eigenvalues `evals`, eigenvectors in `tmp_cdq`
   (line 1679). (Schur fallback via `zgees` if `zheev` fails, lines 1684-1698.)
3. Build `exp(cdq)` (lines 1700-1704):
   ```
   cmtmp(:,i) = tmp_cdq(:,i) * exp(-i·evals(i))
   cdq        = cmtmp · tmp_cdq†                 ! = exp(cdq), unitary
   ```
   i.e. `exp(cdq) = V · diag(exp(-i·evals)) · V†` where `i·cdq = V·diag(evals)·V†`.
   Since `evals` are eigenvalues of `i·cdq`, `exp(-i·evals)` are eigenvalues of `exp(cdq)`.
4. Rotate the gauge (lines 1744-1748):
   ```
   U(k) ← U(k) · exp(cdq(k))
   ```
5. Update M (lines 1752-1768) using the rotation at k and at each neighbour k2:
   ```
   M^{k,b}_new = exp(cdq(k))† · M^{k,b}_old · exp(cdq(k2))
   ```
   (`nkp2 = nnlist(nkp,nn)`; `tmp_cdq = cdq(k)† · M`, then `cmtmp = tmp_cdq · cdq(k2)`.)
   Note it reuses the array name `cdq` for the exponentiated rotation after step 3.

**The rotation is `U ← U · exp(ΔW)` with `ΔW = (α/(4·wbtot)) · (search direction built from
G)`.** Orthonormality is exact because the update is a matrix exponential of an
anti-Hermitian matrix (NOT Löwdin/SVD re-orthonormalization — there is no explicit
re-orthonormalization step in the loop).

---

## 7. Guiding centres — `wann_phases` (line 1776)

Optional (`guiding_centres%enable`, default `.false.`). Chooses `sheet`/`csheet` to fix a
consistent branch cut, especially at the start (poorly-localized WFs).

- `rguide` initialized to projection centres in **Cartesian Å** via
  `utility_frac_to_cart` (driver lines 426-431).
- Called at start if `num_no_guide_iter <= 0` (lines 447-454), and during the loop every
  `num_guide_cycles` iterations once `iter > num_no_guide_iter` (lines 557-566).
- Algorithm (lines 1830-1932): for each WF, average phase over each unique b-direction
  (`csum(na) = Σ_k M_nn` over neighbour pairs via `neigh`/`nnh`), set
  `xx(nn) = −Im ln csum` (arbitrary branch for first 3), then least-squares solve
  `Σ_i smat(j,i)·rguide(i) = svec(j)` with `smat(j,i)=Σ bka_j bka_i`,
  `svec(j)=Σ bka_j·xx(nn)` (lines 1881-1926, via `utility_inv3`). For nn>3 the branch is
  chosen consistent with the current `rguide` (lines 1892-1902).
- Then `sheet(n,nn,k) = Σ_j bk(j,nn,k)·rguide(j,n)` and `csheet = exp(i·sheet)`
  (lines 1946-1957). This effectively picks `Im ln` branch near `−b·rguide`.
- Defaults: `num_guide_cycles = 1`, `num_no_guide_iter = 0` (`wannier90_types.F90:179-180`).

**Restarts:** there is no in-loop restart of U beyond the CG-reset and random-noise escape
(§4). The driver reads the starting `u_matrix` (from projection/`.chk`) into `u_matrix_loc`
(lines 358-361) and writes it back at the end (lines 833-838). `optimisation`/`page_unit`
scratch file (lines 541, 624, 684) just caches M0 to restore after a failed trial step —
not a physics restart.

---

## 8. Convergence — `internal_test_convergence` (line 1084)

Only active if `conv_window > 1` (driver line 812). **Default `conv_window = -1`
(`wannier90_readwrite.F90:681`), i.e. convergence checking is OFF by default and the loop
runs the full `num_iter` iterations.** If `conv_noise_amp > 0`, `conv_window` defaults to 5
(`wannier90_readwrite.F90:682`). Named-type default 3 (`wannier90_types.F90:149`) is
overridden by the readwrite logic.

- Convergence measure: `delta_omega = wann_spread%om_tot − old_spread%om_tot` (line 1122),
  i.e. change in **total** spread between consecutive iterations.
- A ring buffer `history(1:conv_window)` of the last `conv_window` deltas is kept
  (lines 1124-1129, via `eoshift`).
- Converged when **every** entry in the window satisfies `|delta| ≤ conv_tol`
  (lines 1133-1139): loop `j=1..conv_window`, `if abs(history(j)) > conv_tol return`
  (not yet converged); if all pass → converged.
- If `conv_noise_amp > 0`, on reaching the window it adds random noise and requires the
  spread to be stable across `conv_noise_num` noise injections before declaring convergence
  (lines 1141-1160).
- Defaults: `conv_tol = 1e-10` (Å², `wannier90_types.F90:196`), `num_iter = 100`
  (`wannier90_types.F90:192`), `conv_noise_num = 3`, `conv_noise_amp = -1.0` (off).

---

## 9. Driver iteration outline — `wann_main` (line 545 loop)

```
initialize: csheet=1, sheet=0 (310); U_loc from U (358); lambda=0
wann_omega  → initial rave, r2ave, om_* ; cache omega%invariant (463)
for iter = 1 .. num_iter:
    [guiding: wann_phases if enabled]                                    (557)
    wann_domega → cdodq (gradient)                                       (570/578)
    internal_search_direction → cdq (CG dir), doda0, gcnorm1            (596)
    cdqkeep = cdq                                                        (607)
    if lfixstep: alphamin = fixed_step                                  (612)
    else:
        cdq = cdqkeep * trial_step/(4 wbtot)                            (618)
        save U0, M0                                                     (621)
        internal_new_u_and_m (trial rotate)                            (631)
        wann_omega → trial_spread                                       (639)
        internal_optimal_step → alphamin, falphamin, lquad             (646)
    if lfixstep or lquad:
        cdq = cdqkeep * alphamin/(4 wbtot)                             (678)
        restore U0, M0                                                  (682)
        internal_new_u_and_m (real rotate)                             (692)
        old_spread = wann_spread ; wann_omega → wann_spread            (699/702)
    else:  (parabola failed → keep trial step)
        old_spread = wann_spread ; wann_spread = trial_spread          (711/712)
    export omega%total, omega%tilde ; centres, spreads                 (772/777)
    if conv_window>1: internal_test_convergence → lconverged           (813)
    if lconverged: exit                                                (820)
copy U_loc → U, allreduce                                             (833)
```

---

## 10. Silent-mismatch checklist for the Julia port

1. **Two `ln_tmp` conventions**: `q_n` (no wb) in `wann_omega`; `w_b·q_n` in `wann_domega`
   (§2a vs §3a). Do not conflate.
2. **Centre sign**: `r_n = −(1/N_k) Σ w_b b Im ln M_nn`. Leading minus, divide by N_k, wb inside.
3. **Global gradient prefactor** `4/N_k` applied ONCE at the end (§3e); wb applied per-term
   inside, half of it hidden in the wb-carrying `ln_tmp`.
4. **Step normalization**: exponentiate `(α/(4·wbtot))·cdq`; the `4` and `wbtot` here
   pair with the `4` in the gradient (§5, §6).
5. **Branch cut**: principal-value `atan2` for `Im ln`, `csheet=1`/`sheet=0` unless guiding
   centres are on (§1, §7).
6. **Omega_I computed once** and cached (`first_pass`, `omega%invariant`) — recomputing each
   step is harmless numerically but note the code reuses the cached value (§2e).
7. **Units are Å / Å⁻¹ / Å²** internally; `lenconfac` is display-only, default 1 (§0).
8. **B1 condition** `Σ_b w_b b_α b_β = δ_αβ` is assumed by the algebra (§0).
9. **Unitary update is exp of anti-Hermitian** (matrix exponential), no Löwdin/SVD
   re-orthonormalization (§6).
10. **Convergence off by default** (`conv_window=-1`); loop runs full `num_iter=100` with
    `trial_step=2.0`, `num_cg_steps=5`, `conv_tol=1e-10` (§5, §8).
11. **M update at trial and optimal step**: `M_new = exp(cdq_k)† M_old exp(cdq_{k2})`
    with `k2 = nnlist(k,nn)` — the neighbour rotation matters (§6).
12. Out of scope (do not port): `use_ss_functional`, `selective_loc`, `constrain`/`lambda`,
    `lsitesymmetry`, `precond`, random-noise escape.
