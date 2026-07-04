# Wannier Interpolation — Reference Notes (Wannier90 v3.1.0)

Scope: construction of the real-space Hamiltonian `H(R)` from `H(k)`, the
Wigner-Seitz `R`-vector selection with degeneracy weights, the
`use_ws_distance` minimal-image improvement, band interpolation `H(k')` and
diagonalization, and the exact `seedname_hr.dat` file layout.

All file:line citations refer to the reference tree at
`reference/wannier90/src/`. Fortran arrays are **1-based** and **column-major**;
matrices are declared `(row, col, ...)` so `A(:, :, irpt)` is a full
`num_wann × num_wann` slice.

---

## 0. Key constants and defaults

- `twopi` and `cmplx_i = (0,1)` from `w90_constants` (imported e.g.
  `hamiltonian.F90:249`). `rdotk` is always built as `twopi * dot_product(k_frac, R_int)`.
- Reciprocal lattice **includes the 2π factor**:
  `utility_recip_lattice_base` sets `recip_lat = twopi*recip_lat/volume`
  (`utility.F90:349`). Consequently k-points are stored in **fractional
  (crystal) coordinates** `kpt_latt(:,ik)` and R-vectors `irvec(:,irpt)` are
  **integer** lattice-vector components; the phase is `exp(±i·2π·k_frac·R_int)`.
  No explicit Cartesian dot product enters the Fourier phase.
- Defaults (`types.F90:88-91`):
  - `use_ws_distance = .true.`  ← **ON by default** (major gotcha; see §3)
  - `ws_distance_tol = 1.0e-5`
  - `ws_search_size = 2` (a 3-vector, all components 2 by default)
- Energies: `eigval` are the DFT eigenvalues in **eV** as read from
  `seedname.eig`; `H(R)` and interpolated bands are therefore in **eV**. No
  Ha↔eV conversion happens in the interpolation path.

---

## 1. Building H(k) from eigenvalues and U matrices

Routine: `hamiltonian_get_hr` (`hamiltonian.F90:237-628`).

### 1a. Assemble the Wannier-gauge H(k)

`H(k) = U†(k) · diag(eps_k) · U(k)`, computed element-wise. The general
(no site-symmetry) branch, `hamiltonian.F90:393-405`:

```fortran
do loop_kpt = 1, num_kpts
  do j = 1, num_wann
    do i = 1, j                       ! upper triangle only
      do m = 1, num_wann
        ham_k(i, j, loop_kpt) = ham_k(i, j, loop_kpt) + eigval2(m, loop_kpt)* &
                                conjg(u_matrix(m, i, loop_kpt))*u_matrix(m, j, loop_kpt)
      end do
      if (i .lt. j) ham_k(j, i, loop_kpt) = conjg(ham_k(i, j, loop_kpt))  ! Hermitize
    end do
  end do
end do
```

Convention to reproduce exactly:
- Sum index `m` runs over bands/Wannier states; the **complex conjugate sits on
  the left factor** `conjg(u_matrix(m,i))`, i.e. `H_ij = Σ_m U*_{mi} eps_m U_{mj}`.
  This is `U†` with `(U†)_{im} = conj(U_{mi})`.
- Only the **upper triangle `i ≤ j`** is computed; the lower is set by explicit
  Hermitian conjugation. Reproduce this to get bit-identical Hermiticity (avoids
  tiny asymmetries from finite-precision sums).

### 1b. Disentanglement case (`have_disentangled`)

`eigval2(m,k)` is not simply `eigval`; the outer-window states are first slimmed
(`hamiltonian.F90:339-347`, using `dis_manifold%lwindow`), then rotated into the
optimal subspace. Non-sitesym path (`hamiltonian.F90:354-363`):

```fortran
eigval2(j,k) = Σ_m eigval_opt(m,k) * real( conjg(u_matrix_opt(m,j,k)) * u_matrix_opt(m,j,k) )
```

i.e. the diagonal of `U_opt† eps U_opt` (valid because `U_opt` is chosen to
diagonalize the disentangled Hamiltonian). With `lsitesymmetry` the full
`ham_k` is instead built from `utmp = U_opt · U` (`hamiltonian.F90:367-380`).
For a first Julia reimplementation targeting the standard (non-sitesym) path,
the two-step reduces to: `eigval2 = diag(U_opt† eps U_opt)` then
`H(k) = U† diag(eigval2) U`.

`postw90/get_oper.F90` (`get_HH_R`) confirms the same convention via
`v_matrix = u_matrix_opt·u_matrix` and `HH_q(m,n) = Σ conjg(v(i,m)) eps(i) v(i,n)`
(`get_oper.F90:242-245`), Hermitized `HH_q(m,n)=conjg(HH_q(n,m))`.

---

## 2. Fourier transform H(k) → H(R)

Routine body: `hamiltonian.F90:409-461`. The comment states the formula
(`hamiltonian.F90:412`):

`H_ij(R) = (1/N_kpts) Σ_k e^{-ikR} H_ij(k)`

Non-translation (default) branch (`hamiltonian.F90:420-428`):

```fortran
ham_r = cmplx_0
do irpt = 1, nrpts
  do loop_kpt = 1, num_kpts
    rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec(:, irpt), dp))
    fac   = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
    ham_r(:, :, irpt) = ham_r(:, :, irpt) + fac*ham_k(:, :, loop_kpt)
  end do
end do
```

Conventions to reproduce:
- **Phase sign: `exp(-i·2π·k·R)`** for the forward (k→R) transform.
- **Normalization: divide by `N_kpts` = `num_kpts`** (the number of k-points in
  the full MP grid), applied per-term inside `fac`.
- `kpt_latt` fractional, `irvec` integer, cast to real for the dot product.
- The sum runs over **all** `num_kpts`; each `irpt` gets an independent full
  matrix slice.

`postw90` `fourier_q_to_R` is byte-identical in convention
(`get_oper.F90:3130`, `:3147-3155`): `O(R)=(1/N) Σ_q e^{-iqR} O(q)`, phase
`exp(-cmplx_i*rdotq)`, then divide by `num_kpts` at the end.

### Optional WF-centre translation branch (`use_translation`)

`hamiltonian.F90:444-457`: when `translate_home_cell`/`use_translation` is set,
each WF is shifted into the home cell and the phase uses
`irvec_tmp = irvec + shift_vec(:,i) - shift_vec(:,j)`
(`hamiltonian.F90:450-451`). `shift_vec` are integer lattice shifts from
`internal_translate_centres`. Note the index order `ham_r(j,i,irpt)` here.
This is **not** the default and can be ignored unless `translate_home_cell=.true.`.

---

## 3. Wigner-Seitz R-vector set and degeneracies (`ndegen`)

Routine: `hamiltonian_wigner_seitz` (`hamiltonian.F90:695-859`). Called **twice**
from `hamiltonian_setup` (`hamiltonian.F90:113` with `count_pts=.true.` to count,
then `:147` with `count_pts=.false.` to fill) so arrays can be sized first.

Algorithm:
- Metric `real_metric = A·Aᵀ` from `utility_metric` (`utility.F90:416-442`,
  `metric(i,j)=Σ_l lattice(i,l)·lattice(j,l)`).
- Outer loop over candidate lattice points `(n1,n2,n3)` in a supercell
  `±ws_search_size(i)*mp_grid(i)` (`hamiltonian.F90:772-774`).
- Inner loop over Born–von-Karman supercell translations
  `(i1,i2,i3) ∈ [-(ws_search_size+1), ws_search_size+1]`
  (`hamiltonian.F90:780-782`). `ndiff = (n1,n2,n3) - (i1,i2,i3)*mp_grid`
  (`:785-787`), and `dist(icnt) = ndiffᵀ · real_metric · ndiff` (`:788-794`).
- The candidate point is **in the WS cell** iff its distance to R=0 equals the
  minimum over all supercell images, within tolerance:
  `abs(dist(center) - dist_min) < ws_distance_tol²`  (`:800`), where
  `center = (dist_dim+1)/2` is the `R=0` supercell image (`dist_dim` is the total
  count of `(i1,i2,i3)` images, `:745-748`).
- **Degeneracy** `ndegen(nrpts)` = number of supercell images tied for the
  minimum distance within `ws_distance_tol²` (`:803-807`). Weight of point is
  `1/ndegen`.
- `irvec(:,nrpts) = (n1,n2,n3)` (`:808-810`); `rpt_origin` records the index of
  `R=0` (`:813`).
- Tolerance is squared because distances are squared (`ws_distance_tol**2`).

**Sum rule** (`hamiltonian.F90:836-853`): `Σ_irpt 1/ndegen(irpt)` must equal
`mp_grid(1)*mp_grid(2)*mp_grid(3)` within `eps8`, else fatal error. Use this as
a self-check in the reimplementation.

Gotcha: the `(i1,i2,i3)` search range is `ws_search_size+1` (one larger than the
`(n1,n2,n3)` range's `ws_search_size`), on purpose, so a candidate can be
compared against images slightly outside its own search box.

---

## 4. `use_ws_distance` minimal-image improvement

Module: `w90_ws_distance` (`ws_distance.F90`). Parameter `ndegenx = 8`
(`ws_distance.F90:62`) — max degeneracy (cube vertex).

Purpose: for each pair `(i,j)` of Wannier functions and each `R` (`irvec(:,ir)`),
translate WF `j` (at `R + center_j`) by integer multiples of the supercell so it
lands in the Wigner-Seitz cell around WF `i` (at `center_i`). This yields, per
`(i,j,R)`, a set of `ndeg(i,j,R)` equivalent minimal-image R-vectors
`irdist(:,ideg,i,j,R)` with a **pair-specific degeneracy** `ndeg`.

`ws_translate_dist` (`ws_distance.F90:70-176`):
- For each `(ir, jw, iw)`: `irvec_cart = frac_to_cart(irvec(:,ir))`
  (`:150`), then calls `R_wz_sc` with input vector
  `-centre(iw) + (irvec_cart + centre(jw))`, target `R0 = 0`
  (`:160-164`).
- Stores `irdist(:,ideg,iw,jw,ir) = irvec(:,ir) + shifts(:,ideg)` (`:168`) and
  its Cartesian form `crdist` (`:169-171`).

`R_wz_sc` (`ws_distance.F90:179-302`):
- First loop finds the single shortest image over
  `(i,j,k) ∈ [-(ws_search_size+1), ws_search_size+1]`, shift being
  `(i*mp_grid(1), j*mp_grid(2), k*mp_grid(3))` (`:227-252`). Note: shifts are
  **integer multiples of `mp_grid`** (supercell displacements).
- If the minimal vector is essentially zero (`< ws_distance_tol²`), degeneracy
  is 1 and returns (`:264-269`).
- Second loop collects **all** images whose Cartesian distance equals the
  minimum within `ws_distance_tol` (linear tol on `sqrt`, not squared here —
  `:282`), incrementing `ndeg` and **summing** the extra shift onto the
  first-loop shift (`:283-296`). Fatal if `ndeg > ndegenx=8` (`:284-286`).

### How the WS-distance shift enters interpolation

In `plot_interpolate_bands` s-k branch, `use_ws_distance=.true.`
(`plot.F90:786-797`):

```fortran
do j = 1, num_wann
do i = 1, num_wann
  do ideg = 1, ws_distance%ndeg(i, j, irpt)
    rdotk = twopi*dot_product(plot_kpoint(:, loop_kpt), &
                              real(ws_distance%irdist(:, ideg, i, j, irpt), dp))
    fac = cmplx(cos(rdotk), sin(rdotk), dp) &
          /real(ndegen(irpt)*ws_distance%ndeg(i, j, irpt), dp)
    ham_kprm(i, j) = ham_kprm(i, j) + fac*ham_r(i, j, irpt)
  end do
end do
end do
```

Critical points:
- The phase uses the **per-pair shifted vector** `irdist(:,ideg,i,j,irpt)`
  instead of `irvec(:,irpt)`.
- The weight is `1 / (ndegen(irpt) · ndeg(i,j,irpt))` — **both** the WS
  degeneracy from §3 and the pair minimal-image degeneracy divide.
- `ham_r(i,j,irpt)` (the same stored H(R)) is reused; only the phase/weight per
  matrix element changes. This is element-wise, so each `(i,j)` gets its own
  phase — you cannot do a single matrix-level `Σ_R e^{ikR} H(R)` when
  `use_ws_distance` is on.

`ws_distance` is (re)computed via `ws_translate_dist(..., force_recompute=.true.)`
right before the k-loop (`plot.F90:754-759`). The module caches with `done`
flag; `force_recompute` clears it (`ws_distance.F90:111-118`).

`_wsvec.dat` is written by `ws_write_vec` (`ws_distance.F90:306-376`): per
`(irpt, iw, jw)` a header line `irvec(:), iw, jw`, then `ndeg`, then `ndeg`
lines of `irdist - irvec` (the shift). When `use_ws_distance=.false.` it writes
`ndeg=1` and shift `0 0 0`. Note the **iw,jw loop order** in the header vs the
`iw` (outer) / `jw` (inner) order used elsewhere — see §7 index caveat.

---

## 5. Band interpolation H(k') and diagonalization

Routine: `plot_interpolate_bands` (`plot.F90:337`+). Core k-loop
`plot.F90:779-860`.

### 5a. Forward Fourier (R→k'), no ws_distance (`plot.F90:799-802`)

```fortran
rdotk = twopi*dot_product(plot_kpoint(:, loop_kpt), irvec(:, irpt))
fac   = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(irpt), dp)
ham_kprm = ham_kprm + fac*ham_r(:, :, irpt)
```

- **Phase sign: `+i·2π·k'·R`** (opposite sign to the k→R transform in §2) —
  built via `cmplx(cos(rdotk), sin(rdotk)) = exp(+i·rdotk)`.
- **Weight `1/ndegen(irpt)`** — the WS degeneracy divides here, at
  interpolation, **not** when H(R) is built. So `H(R)` stored in `_hr.dat`
  is the *undivided* `(1/N_k)Σ_k e^{-ikR}H(k)`; the reader must divide by
  `ndegen` when reconstructing H(k'). This matches the comment formula
  `H(k')=Σ_R e^{ik'R} H(R)/ndegen(R)`.
- Full matrix slice `ham_r(:,:,irpt)` used at once (matrix-level) in this branch.

### 5b. Diagonalization (`plot.F90:830-847`)

Upper-triangle packed storage, then LAPACK `ZHPEVX`:

```fortran
do j = 1, num_wann
  do i = 1, j
    ham_pack(i + ((j - 1)*j)/2) = ham_kprm(i, j)   ! column-major upper-triangle pack
  end do
end do
call ZHPEVX('V','A','U', num_wann, ham_pack, 0,0,0,0, -1.0_dp, &
            nfound, eig_int(1,loop_kpt), U_int, num_wann, ...)
```

- Only the **upper triangle `i ≤ j`** of `ham_kprm` is packed (index
  `i + j(j-1)/2`). Any non-Hermiticity in the lower triangle is silently
  discarded → Hermitize consistently or rely on the packing.
- `'V'` eigenvectors, `'A'` all, `'U'` upper. Eigenvalues `eig_int(:,loop_kpt)`
  ascending (LAPACK convention) → interpolated band energies in eV.
- Optional projection onto selected WFs: `|U_int(p,w)|²` summed
  (`plot.F90:849-858`).

### 5c. k-path sampling (`plot.F90:527-707`)

Two modes:
- **Explicit** (`bands_kpt_explicit`): `plot_kpoint = bands_kpt_frac`
  (`:640`); `xval` is cumulative path length using `recip_metric`
  (`sqrt(vecᵀ·recip_metric·vec)`, `:646`), where
  `recip_metric = utility_metric(recip_lattice)` (`plot.F90:457`) and
  `recip_lattice` includes 2π.
- **Segment mode** (default, `bands_num_spec_points` special points → `num_paths
  = bands_num_spec_points/2` segments, `:540`). Number of points per segment is
  proportional to segment length: `kpath_pts(seg) = nint(num_points_first_segment
  * kpath_len(seg)/kpath_len(1))` (`plot.F90:571-572`); the first segment gets
  exactly `num_points_first_segment` (`:569`). Points are linearly interpolated
  in fractional coords between the two special points
  (`plot.F90:697-700`):
  `k = P(2s-1) + (P(2s)-P(2s-1)) * (loop_i/kpath_pts(s))`, `loop_i=1..kpath_pts`.
  First point of a segment is printed only when there is a discontinuity
  (`kpath_print_first_point`, `:542-559`). `total_pts` includes these extra
  first points (`:577-580`). The last point is set explicitly to the final
  special point (`:706`).
- `_band.kpt` lists `total_pts` then `k_frac  1.0` per line (`:712-716`);
  `_band.labelinfo.dat` gives label, point index, `xval`, and coords for each
  special point (`:721-741`).

---

## 6. `seedname_hr.dat` file layout

Writer: `hamiltonian_write_hr` (`hamiltonian.F90:631-692`). Exact sequence
(`hamiltonian.F90:677-688`):

1. `write(unit,*) header` — free-format string
   `'written on '//cdate//' at '//ctime` (`:674-677`).
2. `write(unit,*) num_wann` — free-format integer.
3. `write(unit,*) nrpts` — free-format integer.
4. `ndegen` list: `write(unit,'(15I5)') (ndegen(i), i=1,nrpts)` — **15 integers
   per line, width 5** (`:680`). Continues wrapping until all `nrpts` values
   written.
5. Data block, triple loop **`irpt` (outer) → `i` → `j` (inner)**
   (`hamiltonian.F90:681-688`):

```fortran
do irpt = 1, nrpts
  do i = 1, num_wann
    do j = 1, num_wann
      write(unit,'(5I5,2F12.6)') irvec(:, irpt), j, i, ham_r(j, i, irpt)
    end do
  end do
end do
```

Per data line, format `'(5I5,2F12.6)'`:
`Rx Ry Rz  j  i   Re(H)  Im(H)`, where each integer is width-5 and each real is
`F12.6`.

**Critical index conventions:**
- The **fastest-varying printed index is `j`** (the inner loop). Columns 4,5 are
  `j` then `i` (in that print order).
- The written matrix element is `ham_r(j, i, irpt)` — i.e. **row `j`, column
  `i`**. So the value on a line labeled `... j i ...` is `H_{j,i}(R)`. Reading
  back: element at printed columns `(m,n)` = `H_{m,n}(R)`, since the loop prints
  `j=m` (row), `i=n` (col) and value `ham_r(m,n)`. In practice standard readers
  index `H(row=col4, col=col5, R)` and it round-trips.
- `Re/Im` are the real and imaginary parts of the **undivided** H(R)
  (before `/ndegen`). Consumers must apply `1/ndegen(irpt)` at interpolation.

Related outputs (same directory of code):
- `_tb.dat` (`hamiltonian_write_tb`, `:862`+): header, 3 lattice-vector lines
  `real_lattice(1,:)`,`(2,:)`,`(3,:)` (rows are a1,a2,a3), `num_wann`, `nrpts`,
  `ndegen` `(15I5)`, then per-R block `irvec` on its own line then
  `j i Re Im` with `'(2I5,3x,2(E15.8,1x))'` and value `ham_r(j,i,irpt)`
  (`:941-948`), then the `<0n|r|Rm>` position-operator blocks.
- `_wsvec.dat` — see §4.

---

## 7. Silent-mismatch gotchas checklist

1. **Fourier sign asymmetry**: k→R uses `e^{-i2πkR}` (§2); R→k' uses
   `e^{+i2πkR}` (§5a). Getting either sign wrong flips band structure phases but
   may still look plausible.
2. **`ndegen` division timing**: NOT divided when building/storing H(R); divided
   `1/ndegen(irpt)` at interpolation. `_hr.dat` holds undivided values.
3. **`use_ws_distance=.true.` is the DEFAULT.** Reproducing reference `.dat`
   bands requires the per-element minimal-image shifts and the **double** weight
   `1/(ndegen(irpt)·ndeg(i,j,irpt))` (§4). A plain matrix-level
   `Σ_R e^{ikR}H(R)/ndegen` will NOT match reference bands unless
   `use_ws_distance=.false.`.
4. **WS tolerance is squared** in `hamiltonian_wigner_seitz` (`ws_distance_tol²`,
   compared against squared distances) but **linear** in the second loop of
   `R_wz_sc` (compared against `sqrt`). Mixing these changes degeneracy counts at
   cell edges.
5. **`ws_search_size+1`**: both WS routines inflate the image search range by one
   over `ws_search_size` (`hamiltonian.F90:780`, `ws_distance.F90:227`). Off-by-one
   here changes which edge points are found.
6. **Metric = A·Aᵀ** with lattice rows as vectors (`utility.F90:416`); reciprocal
   lattice carries the 2π (`utility.F90:349`), so k dot R uses only integer R and
   fractional k times 2π — never Cartesian.
7. **Hermitization / upper-triangle only**: H(k) built for `i≤j` then
   conjugate-filled; diagonalization packs only `i≤j`. Fill the lower triangle
   consistently or you diagonalize a different matrix.
8. **`_hr.dat` index order**: data prints `irvec, j, i, ham_r(j,i,irpt)` with
   `j` innermost; the `ndegen` line is `(15I5)` wrapped. Both easy to transpose
   or mis-wrap.
9. **Energy units eV** throughout; no Ha conversion. Lattice in the `_tb.dat`
   is in Å (whatever `real_lattice` holds), rows = a1,a2,a3.
10. **N_k = num_kpts** = full MP grid size (`mp_grid(1)*mp_grid(2)*mp_grid(3)`),
    used both as Fourier normalization and in the WS sum-rule check.
