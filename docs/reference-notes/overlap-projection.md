# Overlap & Projection — Implementation Notes (Wannier90 v3.1.0)

Scope: `src/overlap.F90` plus the disentanglement counterparts in `src/disentangle.F90`
and the driver glue in `src/library_interface.F90` / `src/library_extra.F90`.
All file:line citations are into `reference/wannier90/src/`.

This stage does three things:
1. Read the `.mmn` overlaps → `M^{k,b}_{mn}` (`overlap_read`).
2. Read the `.amn` projections → `A^k_{mn}` (`overlap_read`).
3. Build the **initial Wannier gauge** `U^k` from `A` by SVD-based Löwdin
   orthonormalization, and rotate `M` into that gauge (`overlap_project`,
   `overlap_project_gamma`; disentanglement analogues `dis_project`,
   `internal_find_u`, `internal_find_u_gamma`).

All arrays are **1-based, column-major Fortran**. `dp = kind(1.0d0)` (constants.F90:52).
`cmplx_0=(0,0)`, `cmplx_1=(1,0)`, `cmplx_i=(0,1)` (constants.F90:58-62).
`eps5=1e-5` (unitarity check tol), `eps8=1e-8` (constants.F90:70,76).

---

## 0. What is NOT in this stage (avoid hunting for missing factors)

Confirmed absent from `overlap.F90` and the projection routines:
- **No `Im ln` / branch cut.** Centres/spreads (which contain the `Im ln` with its
  branch subtlety) are computed later in `wannierise.F90`. Not here.
- **No b-vector weights `w_b`** and **no `N_k` / `N_b` normalization factors.**
  This stage only reads/rotates matrices; weights live in `kmesh.F90` and are applied
  in `wannierise.F90`/`disentangle.F90` spread functionals. There is *no* `1/N_k` or
  `1/N_b` anywhere in the read or the gauge construction.
- **No unit conversions** (Bohr/Ang, Ha/eV) in this stage. `M` and `A` are
  dimensionless overlaps; conversions happen at I/O of positions/energies elsewhere.

---

## 1. The M-matrix (Mmn overlaps)

### 1.1 Physical definition and index convention (LOAD-BEARING)

```
M^{k,b}_{mn} = < u_{m,k} | u_{n,k+b} >
```
- `m` = **bra** state, at k, and is the state that gets **complex-conjugated**.
- `n` = **ket** state, at the neighbour k+b.

This is fixed by the read loop `overlap.F90:282-285`:
```fortran
do n = 1, num_bands          ! slow / column
  do m = 1, num_bands        ! fast / row
    read (mmn_in, *, ...) m_real, m_imag
    mmn_tmp(m, n) = cmplx(m_real, m_imag, kind=dp)
```
`m` is the fast (row) index, `n` the column. Storage: `M(m, n, nn, k)`.
Getting the bra/ket backwards transposes every subsequent `U^dag M U`, so this must
match the reimplementation exactly.

### 1.2 `.mmn` file format (overlap_read, overlap.F90:246-320)

- Line 1: comment (echoed).
- Line 2: `num_bands  num_kpts  nntot` — checked against expected (overlap.F90:256-271).
- Then `num_kpts*nntot` blocks. Each block header (overlap.F90:281):
  `nkp  nkp2  nnl nnm nnn` = (this k index, neighbour k index, 3-int G-vector shift
  `nncell` bringing k+b back into the BZ).
- Each block body: `num_bands*num_bands` complex numbers, in order `do n; do m`
  (m fastest), one `(real imag)` pair per line (overlap.F90:282-287).

### 1.3 b-vector ordering is CANONICALIZED, not file order (GOTCHA)

The stored neighbour slot `nn` is **not** the order rows appear in `.mmn`. For each
file block, W90 searches `kmesh_info%nnlist`/`nncell` for the unique match
(overlap.F90:288-312):
```fortran
do inn = 1, kmesh_info%nntot
  if ((nkp2 == nnlist(nkp,inn)) .and. (nnl == nncell(1,nkp,inn)) .and.
      (nnm == nncell(2,nkp,inn)) .and. (nnn == nncell(3,nkp,inn))) then
      nn = inn                          ! must be unique, else fatal error
m_matrix_local(:, :, nn, map_kpts(nkp)) = mmn_tmp(:, :)
```
Consequence: to reproduce W90's `nn` axis you must reproduce **kmesh's** neighbour
construction (`kmesh.F90`) exactly. If your neighbour list is ordered differently,
every b-indexed array is silently permuted — no error, wrong numbers. `map_kpts`
(overlap.F90:233-239) is just the global→local k compaction for MPI; single-rank it
is the identity.

### 1.4 Storage arrays and the disentanglement split (overlap_allocate, overlap.F90:56-163)

`disentanglement = (num_bands > num_wann)` (overlap.F90:95).

- **Disentanglement path:** `.mmn` (which is `num_bands x num_bands`) is read into
  `m_matrix_orig` / `m_matrix_orig_local`, both dimensioned
  `(num_bands, num_bands, nntot, num_kpts|nkl)` (overlap.F90:106-118).
  The `num_wann`-sized `m_matrix`/`m_matrix_local` is allocated but filled later by
  `setup_m_loc` after the optimal subspace is found.
- **No-disentanglement path:** `m_matrix`/`m_matrix_local` dimensioned
  `(num_wann, num_wann, nntot, ...)`, and since `num_bands == num_wann` the read fills
  them directly (overlap.F90:127-141).

NB: `overlap_read` itself is **size-agnostic** — its `m_matrix_local(:,:,:,:)` dummy
receives whatever the caller passes, and it always writes a `num_bands x num_bands`
block (`mmn_tmp`). So the *caller* is responsible for passing the num_bands-sized array
when disentangling. In the newer library API the array is a caller-owned pointer
(`lib_common_type%m_matrix_local`, library_interface.F90:100,850-852); the classic
sizing above is what `overlap_allocate` provides.

---

## 2. The A-matrix (Amn projections)

### 2.1 Definition
```
A^k_{mn} = < psi_{m,k} | g_n >
```
Bloch state `m` (bra, conjugated) projected onto trial orbital `g_n`. The `.amn` file
is produced **externally** (e.g. pw2wannier90); W90 only reads it. There is no analytic
projection *generator* inside W90 core — see §6.

### 2.2 `.amn` read (overlap_read, overlap.F90:335-378)

- Header: comment, then `num_bands  num_kpts  num_proj` (checked, overlap.F90:345-359).
- Body: `num_bands*num_proj*num_kpts` lines, each:
  `m  n  nkp  a_real  a_imag`   (overlap.F90:374)
  ```fortran
  read (amn_in, *, ...) m, n, nkp, a_real, a_imag
  if (select_projection%proj2wann_map(n) < 0) cycle
  au_matrix(m, proj2wann_map(n), nkp) = cmplx(a_real, a_imag, kind=dp)
  ```
  - `m` = band (row, 1..num_bands), `n` = projection (1..num_proj), `nkp` = k.
  - `proj2wann_map` handles `select_projections`: `num_proj >= num_wann`; projections
    mapped to `<0` are dropped, others remap the column into `1..num_wann`. Without
    selection this is the identity and `num_proj == num_wann`.

### 2.3 CRITICAL aliasing: A is read into `u_matrix_opt` (`au_matrix`)

The dummy `au_matrix` in `overlap_read` is the **`u_matrix_opt`** array
(library_extra.F90:243-245: "projections are stored in u_opt"). The projection matrix
`A` therefore *lives in `u_matrix_opt`* until the gauge step overwrites it. The name
`a_matrix` inside `overlap_project` / `dis_project` is likewise the same storage passed
as `u_matrix`/`u_matrix_opt`; ZGESVD overwrites its input, so `A` is destroyed by the
SVD. A reimplementation should keep a separate copy of `A` if it needs it afterward.

### 2.4 `use_bloch_phases` (overlap.F90:384-392)

If `use_bloch_phases=.true.`, `.amn` is **not** read; instead `A` is set to the
identity: `au_matrix(m,m,n)=cmplx_1` for `m=1..num_wann`, all k. The Bloch states are
used directly as the initial gauge (trivial projection).

---

## 3. Initial-gauge construction: SVD Löwdin (NO disentanglement)

Routine `overlap_project` (overlap.F90:838-1032). Here `num_wann == num_bands`
(disentanglement is FALSE). Precondition from the driver
(library_interface.F90:658-665): `u_matrix` is set to `u_matrix_opt` (i.e. holds `A`),
and `u_matrix_opt` is reset to the identity.

### 3.1 The math actually computed (GOTCHA — SVD form, not S^{-1/2})

Comment (overlap.F90:919-920): `CU = CS^{-1/2}.CA, CS = CA.CA^dagger`. But the code
**never forms `S^{-1/2}`.** It takes the SVD of `A` and multiplies the unitary factors,
discarding the singular values:

```fortran
call zgesvd('A','A', num_bands, num_bands, u_matrix(1,1,nkp), num_bands,
            svals, cz, num_bands, cvdag, num_bands, cwork, 4*num_bands, rwork, info)   ! :927
call utility_zgemm(u_matrix(:,:,nkp), cz, 'N', cvdag, 'N', num_wann)                   ! U = Z . Vdag  :940
```
So with `A = Z S V^dag` (LAPACK `zgesvd` returns `V^dag` directly, not `V`):
```
U^k = Z V^dag
```
This equals the Löwdin `U = A (A^dag A)^{-1/2}` (and, for square `A`, also
`(A A^dag)^{-1/2} A`) because `A(A^dag A)^{-1/2} = Z S V^dag (V S^2 V^dag)^{-1/2}
= Z S V^dag V S^{-1} V^dag = Z V^dag`. **Reproduce numerically by doing the same
Z·V^dag from an SVD.** An eigendecomposition-based `S^{-1/2}` is mathematically equal
but a *different* code path (different rounding); prefer the SVD path to match bit-for-bit.

`utility_zgemm(c,a,transa,b,transb,n)` computes `C = op(A) op(B)` via `zgemm`
(utility.F90:78-108): `'N'`→as-is, `'T'`→transpose, `'C'`→conjugate-transpose.

### 3.2 Unitarity check
After building `U`, checks `sum_m U(m,j) conj(U(m,i)) == delta_ij` to `eps5`
(overlap.F90:945-972) — i.e. columns of `U` orthonormal. Fatal error otherwise.
MPI: non-owned k are zeroed then `comms_allreduce`'d (overlap.F90:974-979); single-rank
is a no-op.

### 3.3 Rotate M into the initial gauge (EXACT ORDER)

overlap.F90:989-1003:
```fortran
do nn = 1, nntot
  nkp2 = nnlist(nkp, nn)                                    ! neighbour k+b
  call utility_zgemm(cvdag, u_matrix(:,:,nkp),  'C', m_matrix_local(:,:,nn,nkp_loc), 'N', num_wann) ! cvdag = U(k)^dag . M
  call utility_zgemm(cz,    cvdag,              'N', u_matrix(:,:,nkp2),             'N', num_wann) ! cz = cvdag . U(k+b)
  m_matrix_local(:,:,nn,nkp_loc) = cz(:,:)
```
Result:
```
M^{k,b}_new = U(k)^dag . M^{k,b} . U(k+b)
```
- Left factor: **conjugate-transpose of `U` at k** (the `'C'`).
- Right factor: **`U` at the neighbour k+b** (plain `'N'`), where `k+b = nnlist(k,nn)`.
- Matrix product order is exactly `Udag · M · U(neighbour)`. Using `U(k)` on the right
  instead of `U(k+b)`, or dropping the conjugate, silently gives wrong overlaps.

---

## 4. Gamma-only variant (`overlap_project_gamma`, overlap.F90:1036-1267)

Single k-point, real arithmetic. Take the real part of `A` (`u_matrix_r =
real(u_matrix(:,:,1))`, :1130), then real SVD:
```fortran
call dgesvd('A','A', num_wann, num_wann, u_matrix_r, ..., svals, rz, ..., rv, ...)   ! :1171
call dgemm('N','N', ..., rz, ..., rv, ..., u_matrix_r, ...)                          ! U_r = Z . V^T  :1182
u_matrix(:,:,1) = cmplx(u_matrix_r, 0.0_dp)                                          ! :1214
```
Same `Z·V^dag` structure, real. Unitarity checked with `u_r(m,j)*u_r(m,i)`
(no conjugation, real) to `eps5` (:1187-1212). M rotation identical form,
single k (:1219-1225): `M = U^dag M U`.

---

## 5. Disentanglement path (num_bands > num_wann)

Called instead of `overlap_project` when disentangling. Three relevant routines in
`disentangle.F90`:

### 5.1 `dis_project` — rectangular Löwdin of A within the window (:1600-1856)

Analogous "unitarised projection" `A_mn=<psi_m|g_n> -> S=A A^+ -> U = S^{-1/2} A`
(printed :1686), but rectangular. Only `A(A^dag A)^{-1/2}` is well-defined here
(`A A^dag` is rank-deficient / `ndimwin x ndimwin` but rank `num_wann`).

- First "slims" `A` to the outer-window rows `1..ndimwin(nkp)` (:1730-1740), rows
  offset by `nfirstwin(nkp)`.
- SVD of the `ndimwin(nkp) x num_wann` block:
  `zgesvd('A','A', ndimwin(nkp), num_wann, a_matrix, num_bands, svals, cz, ..., cvdag, ...)` (:1744)
- Build `U_opt` as `Z·V^dag` with the inner sum **truncated to `num_wann`** (rectangular
  truncation, :1776-1783):
  ```fortran
  do j=1,num_wann; do i=1,ndimwin(nkp); do l=1,num_wann
     u_matrix_opt(i,j,nkp) += cz(i,l)*cvdag(l,j)
     a_matrix(i,j,nkp)     += cz(i,l)*svals(l)*cvdag(l,j)   ! reconstruct A only because zgesvd destroyed it
  ```
  The `l=1..num_wann` cap is the `Z S S^{-1} V^dag = Z_{:,1:nw} V^dag` reduction
  (comment :1759-1769). `U_opt` is `ndimwin(k) x num_wann` with orthonormal columns
  (checked to `eps5`, :1793-1820). This is the initial `u_matrix_opt` (the subspace).

### 5.2 `internal_find_u` — square U from optimized subspace (:566-685)

After the subspace optimization, build the square `num_wann x num_wann` gauge:
```fortran
call zgemm('C','N', num_wann, num_wann, ndimwin(nkp), cmplx_1, u_matrix_opt(:,:,nkp), num_bands,
           a_matrix(:,:,nkp), num_bands, cmplx_0, caa(:,:,nkp), num_wann)      ! caa = U_opt^dag . A   :665
call zgesvd('A','A', num_wann, num_wann, caa(:,:,nkp), ..., svals, cz, ..., cv, ...)  ! :668
call zgemm('N','N', ..., cz, ..., cv, ..., u_matrix(:,:,nkp), ...)             ! U = Z . Vdag   :681
```
`caa = U_opt^dag A` (the projection of trial orbitals onto the optimal subspace), then
`U = Z·V^dag` (same Löwdin). Serial-on-root then `comms_bcast` (:685).

### 5.3 `setup_m_loc` — rotate M into the subspace (:310-405)

Maps the num_bands-sized `m_matrix_orig_local` into the num_wann-sized
`m_matrix_local` using the subspace `u_matrix` (which here holds `u_matrix_opt`
contracted appropriately):
```fortran
do nn: nkp2 = nnlist(nkp_global, nn)
  zgemm('C','N', ..., u_matrix(:,:,nkp_global), ..., m_matrix_orig_local(:,:,nn,nkp), num_bands, ..., cwb) ! cwb = U(k)^dag M
  zgemm('N','N', ..., cwb, ..., u_matrix(:,:,nkp2), ..., cww)                                             ! cww = cwb U(k+b)
  m_matrix_local(1:num_wann,1:num_wann,nn,nkp) = cww
```
Same convention as §3.3: `M_new = U(k)^dag M U(k+b)`, neighbour `k+b = nnlist(k,nn)`.
Note the **leading dimension `num_bands`** on the M argument (it is stored num_bands-wide
but only the top-left num_wann block is used after slimming).

### 5.4 Gamma disentanglement: `internal_find_u_gamma` (:735-884)

Real analogue of `internal_find_u`: `raa = U_opt_r^T A_r` via `dgemm('T','N',...)`
(:820), `dgesvd` (:823), `U = Z·V^T` via `dgemm('N','N',...)` (:836),
`u_matrix = cmplx(raa,0)` (:839).

---

## 6. Analytic projection functions (high level only)

Wannier90 v3.1.0 core does **not** generate `A_mn` from analytic orbitals — that is the
job of the interface code (pw2wannier90 etc.), which evaluates
`<psi_{m,k}|g_n>` for trial orbitals `g_n` (Gaussian×real spherical harmonics, with
`l,mr,r` quantum numbers, position, axes/z-axis, radial index, zona) parsed from the
`projections` block. In this reference tree the only in-core "projection" logic is the
**orthonormalization** of an already-read `A` (§3, §5) — the "unitarised projection"
(`A -> S^{-1/2} A`). For the Julia reimplementation of *this* stage, `A` is an input
read from `.amn`; the analytic generator is out of scope and would need to match the
interface code's orbital conventions if reimplemented.

---

## 7. Reimplementation checklist (silent-mismatch traps)

1. `M(m,n,nn,k) = <u_{m,k}|u_{n,k+b}>`; `m`=row=bra=conjugated. `m` fastest in file.
2. Assign the b-slot `nn` by **matching `(nkp2, G)` to the kmesh neighbour list**, not by
   file order. Reproduce `kmesh.F90` neighbour ordering or everything permutes.
3. `A` is read into `u_matrix_opt`; ZGESVD **destroys** it — copy first if needed.
4. Initial gauge = `U = Z·V^dag` from `zgesvd(A)` (LAPACK returns `V^dag`). Do **not**
   form `S^{-1/2}` explicitly — match the SVD path.
5. Disentanglement rectangular case: SVD gives `ndimwin x num_wann`; truncate inner
   sum to `num_wann` (`U_opt = Z_{:,1:nw} V^dag`).
6. M rotation is exactly `U(k)^dag · M^{k,b} · U(k+b)`, `k+b = nnlist(k,nn)`. Conjugate
   the k factor (`'C'`), plain neighbour factor (`'N'`), in that product order.
7. Gamma-only paths use the real part and real SVD/GEMM; unitarity check without conjugation.
8. No `w_b`, no `1/N_k`, no `Im ln`, no unit conversions in this stage.
9. Unitarity tolerance `eps5=1e-5`; column-orthonormality is what is checked.
