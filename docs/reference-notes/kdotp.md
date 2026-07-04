# postw90 reference notes: k·p expansion coefficients (berry_task = kdotp)

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Key reference (`src/postw90/berry.F90:45`):
- IAdJS19 = Ibañez-Azpiroz, de Juan, Souza, arXiv:1910.06172 (2019) — quasi-degenerate
  (Löwdin) k·p from Wannier functions. User-guide derivation:
  `docs/docs/user_guide/postw90/berry.md:481-634` (cites Winkler App. B for the
  quasi-degenerate PT formulas and WYSV06 Sec. III.B for the Wannier-gauge notation).

No physical constants enter this task (no CODATA dependence): everything is in eV (from
`seedname.eig`) and Å (from `real_lattice`); the only numeric constant used is `twopi`.

---

## A. Keywords and parsing

All parsed in `w90_wannier90_readwrite_read_berry` (`src/postw90/postw90_readwrite.F90`):

| keyword | parse site | type / default | notes |
|---|---|---|---|
| `berry` | (calc flags) | logical, F | must be `true` to enter `berry_main` |
| `berry_task` | postw90_readwrite.F90:865-879 | string, required when `berry=T` | kdotp active iff `index(task,'kdotp')>0` (substring match; combinable, e.g. `ahc-kdotp`) |
| `kdotp_kpoint` | postw90_readwrite.F90:970-972 | real(3), default `0.0 0.0 0.0` (`src/postw90/postw90_types.F90:178`) | **fractional (reduced) coordinates** of the reciprocal lattice — used directly as `kpt` in `twopi*dot_product(kpt, irvec)` phases. (The docs CSV `docs/docs/parameters/postw90-berry-parameters.csv:29` says "2π/a units", which is misleading; the code treats it as fractional.) |
| `kdotp_num_bands` | postw90_readwrite.F90:998-1006 | integer, no default (0 sentinel) | must be ≥ 1 if present; if absent, `kdotp_bands` is never allocated (see Traps) |
| `kdotp_bands` | postw90_readwrite.F90:1012-1018 | integer list, no default | read via `w90_readwrite_get_range_vector` (accepts `4,5` or `4-5`), must have exactly `kdotp_num_bands` entries, all ≥ 1. Indices refer to **Wannier-interpolated band indices at kdotp_kpoint** (ascending-eigenvalue order of the num_wann interpolated bands), NOT ab-initio band indices |
| `fermi_energy` | (fermi block) | required | not used by kdotp itself, but `berry_main` errors with `'Must specify one or more Fermi levels when berry=true'` if `fermi_n==0` (`src/postw90/berry.F90:254-257`) |
| `berry_kmesh` | (kmesh block) | — | the berry interpolation k-loop still runs (does nothing for kdotp-only); set `1 1 1` to make it trivial (as the test does) |

`kdotp_bands` is broadcast to all ranks in `postw90_param_dist`
(`src/postw90/postw90_common.F90:711-724`), but the actual kdotp computation runs **on root
only** (`berry.F90:579-590`, inside `if (print_output%iprint > 0)`; see comment at
berry.F90:580-582 about the historical segfault).

Parameter-summary stdout block (`src/postw90/postw90_readwrite.F90:2341-2345, 2371-2378`):

```
|  Compute k.p expansion coefficients        :        T                     |     ! '(1x,a46,10x,a8,13x,a1)'
|  Chosen k-point kdotp_kpoint                 :   0.000    0.500    0.000  |     ! '(1x,a46,10x,f8.3,1x,f8.3,1x,f8.3,1x,13x,a1)'
|  kdotp_num_bands                             :     2                      |     ! '(1x,a46,10x,i4,13x,a1)'
|  kdotp_bands                                 :     4   5                  |     ! '(1x,a46,10x,*(i4))'  (no trailing '|')
```

---

## B. Pipeline (call graph)

`berry_main` (`src/postw90/berry.F90`):

1. `eval_kdotp = index(pw90_berry%task,'kdotp')>0` (berry.F90:280, 287).
2. All ranks: `get_HH_R` (berry.F90:502-507) — builds the real-space Hamiltonian
   (see berry-ahc.md §A.1; formulas repeated in §C below). **No AA_R / .mmn needed.**
3. Root allocates `kdotp(kdotp_nbands, kdotp_nbands, 3, 3, 3) = cmplx_0` (berry.F90:508-510);
   `kdotp_nbands = size(pw90_berry%kdotp_bands)`.
4. Root: `berry_get_kdotp` (berry.F90:583-588 → subroutine at berry.F90:3340-3498).
5. Root: writes the three output files (berry.F90:1731-1770).

Inside `berry_get_kdotp` (berry.F90:3421-3455):

1. `wham_get_eig_UU_HH_AA_sc` (`src/postw90/wan_ham.F90:770-839`): one Fourier call
   `pw90common_fourier_R_to_k_new_second_d` giving `HH`, `HH_da(:,:,3)`, `HH_dadb(:,:,3,3)`
   at `kpt = pw90_berry%kdotp_kpoint`, then `utility_diagonalize(HH) → eig, UU`.
2. `wham_get_eig_deleig` (wan_ham.F90:442-543): **recomputes and overwrites** `HH`
   (alpha=0 Fourier), re-diagonalizes → same `eig, UU` (same matrix, same ZHPEVX call),
   recomputes `HH_da` via three alpha=1,2,3 Fourier calls (identical formula), and computes
   `eig_da` (band velocities) — `eig_da` is **never used** afterwards.
3. `wham_get_D_h_P_value` (wan_ham.F90:145-193): computes
   `D_h(n,m,i) = H̄_i(n,m) · ΔE/(ΔE² + sc_eta²)`, `ΔE = eig(m)−eig(n)` — `D_h` is **never
   used** afterwards (but the call reads `pw90_berry%sc_eta`, default 0.04,
   postw90_types.F90:171). Dead code with no output effect.
4. Rotation to the Hamiltonian gauge (berry.F90:3446-3455) with `utility_rotate`
   (`src/utility.F90:699-716`, `rot† · mat · rot`):

```
HH_bar        = U† HH U            ! computed but NEVER used (eig used instead)
HH_da_bar(a)  = U† HH_da(a) U      ! a = 1..3 Cartesian
HH_dadb_bar(a,b) = U† HH_dadb(a,b) U
```

`utility_diagonalize` (`src/utility.F90:652-696`) packs the upper triangle and calls LAPACK
`ZHPEVX('V','A','U', ...)` → eigenvalues ascending in `eig`, eigenvectors as columns of `UU`.

---

## C. Matrix elements needed and their conventions

### C.1 HH_R (only real-space operator required)

`get_HH_R` (`src/postw90/get_oper.F90:64-291`), identical to the AHC path:

```
HH_q(n,m,ik) = Σ_{i=1..num_states(ik)} conjg(v_matrix(i,n,ik)) · eigval(winmin_q+i-1,ik) · v_matrix(i,m,ik)   ! get_oper.F90:237-247, hermitized
HH_R_temp(:,:,ir) = (1/N_kpts) Σ_q e^{−i 2π q·R} HH_q(:,:,q)        ! fourier_q_to_R, get_oper.F90:3124-3157
HH_R = operator_wigner_setup(HH_R_temp)                              ! /ndegen(ir) (+ ws_distance remap), get_oper.F90:3275-3327
```

Optional scissors shift applies here too (get_oper.F90:255-278). After
`operator_wigner_setup`, interpolation uses `nrpts_pw90/irvec_pw90/crvec_pw90` with degeneracy
weights already baked in (irvec_pw90 = irvec when `use_ws_distance=F`).

### C.2 Interpolation to kdotp_kpoint

`pw90common_fourier_R_to_k_new_second_d` (`src/postw90/postw90_common.F90:1161-1232`), with
`rdotk = twopi·dot_product(kpt, irvec_pw90(:,ir))`, `phase = cos(rdotk) + i·sin(rdotk)`
(= e^{+i 2π k·R}):

```
HH(k)        = Σ_R e^{+ik·R} H(R)                                   ! line 1214
HH_da(k)_a   = Σ_R  i·Rc_a        · e^{+ik·R} H(R)                  ! lines 1215-1220, factor +i
HH_dadb(k)_ab = − Σ_R Rc_a · Rc_b · e^{+ik·R} H(R)                  ! lines 1221-1229, factor −1 (from i²)
```

`Rc = crvec_pw90(:,ir) = matmul(transpose(real_lattice), irvec(:,ir))` — **Cartesian, in Å**
(real_lattice rows = lattice vectors in Å). Hence derivatives are w.r.t. Cartesian k in Å⁻¹:
HH in eV, HH_da in eV·Å, HH_dadb in eV·Å². The first-derivative that actually survives into
the output is recomputed by `pw90common_fourier_R_to_k` (postw90_common.F90:1032-1093,
alpha=1,2,3) — numerically the identical sum. Neither R→k routine applies any per-(i,j)
ws_distance shift at interpolation time in this refactored source; `use_ws_distance` acts
solely through the HH_R remap of §C.1.

---

## D. Löwdin / quasi-degenerate expansion as implemented (berry.F90:3457-3496)

Set A = `kdotp_bands` (size N = kdotp_num_bands); set B = all other interpolated bands
r ∈ {1..num_wann} \ A. All formulas use the H-gauge (barred) matrices of §B step 4 and the
ascending eigenvalues `eig`. Below `bn = kdotp_bands(n)`, `bm = kdotp_bands(m)` are the
global band indices; n,m = 1..N index the output matrix.

**Order 0** (berry.F90:3463) — diagonal only, off-diagonals stay 0 from the cmplx_0 init:

```
kdotp(n,m,1,1,1) = eig(bn)          if n == m, else 0
```

**Order 1** (berry.F90:3465-3467) — for a = 1..3 (Cartesian x,y,z):

```
kdotp(n,m,2,a,1) = HH_da_bar(bn, bm, a)                 ! = [U† ∂_a H^W U]_{bn,bm},  eV·Å
```

Row index bn first, column bm second — no conjugation beyond the U†…U rotation.

**Order 2** (berry.F90:3469-3493) — for a,b = 1..3:

```
kdotp(n,m,3,a,b) = 0.5 · HH_dadb_bar(bn, bm, a, b)
                 + 0.5 · Σ_{r ∉ A}  HH_da_bar(bn, r, a) · HH_da_bar(r, bm, b)
                        · ( 1/(eig(bn) − eig(r)) + 1/(eig(bm) − eig(r)) )        ! eV·Å²
```

i.e. the stored order-2 tensor is `½ [ (H̄_ab)_{nm} + (T_ab)_{nm} ]` with the
virtual-transition matrix `T_ab` of IAdJS19 / user guide Eq. (Tab). The r-loop runs over all
`r = 1..num_wann` and `cycle`s when r equals any entry of `kdotp_bands` (berry.F90:3476-3483).
There is **no regularization** of the energy denominators — if a B-state is degenerate with an
A-state at kdotp_kpoint the coefficient diverges.

**Reconstruction of the effective Hamiltonian** (contraction convention): with
κ = Cartesian (k − k0) in Å⁻¹,

```
H_kp(κ)_{nm} = kdotp(n,m,1,1,1) + Σ_a kdotp(n,m,2,a,1)·κ_a + Σ_{a,b} kdotp(n,m,3,a,b)·κ_a·κ_b
```

— no extra ½ on the second-order contraction (the ½ is already stored). The full double sum
over (a,b) must be kept: for fixed (a,b) the T-block satisfies `T_ab(n,m)* = T_ba(m,n)`, so a
single (a,b) block is NOT Hermitian; Hermiticity holds only after the symmetric κ_a κ_b sum.

---

## E. Output files (root only, berry.F90:1731-1770)

stdout banner (berry.F90:1736-1741, then `'(/,3x,a)'` per file):

```
 ----------------------------------------------------------
 Output data files related to k.p:                         
 ----------------------------------------------------------

   * gaas-kdotp_0.dat
```

All three files: `open(FILE=..., STATUS='UNKNOWN', FORM='FORMATTED')`, **no header lines, no
comment lines** — pure data. Every record is one complex number written with Fortran format
`'(2E18.8E3)'` (format reversion: the 2-field format consumes exactly one complex = re, im per
line). Each field is width 18, 8 mantissa digits, **3-digit exponent**, Fortran leading-zero
form, e.g. `   0.64827435E+001` (positive: 3 leading blanks + 15 chars; negative uses one blank
less: `  -0.64827435E+001`).

Array sections are written in Fortran column-major element order, so **n (the row/bra index)
varies fastest**: (n,m) = (1,1), (2,1), …, (N,1), (1,2), …, (N,N).

1. `seedname-kdotp_0.dat` (berry.F90:1743-1748): single write of `kdotp(:,:,1,1,1)`
   → N² lines. Only diagonal entries are nonzero (the eigenvalues, eV).
2. `seedname-kdotp_1.dat` (berry.F90:1751-1757): `do i = 1,3: write kdotp(:,:,2,i,1)`
   → 3·N² lines; three consecutive N²-line blocks for a = x, y, z (eV·Å).
3. `seedname-kdotp_2.dat` (berry.F90:1760-1768): `do i = 1,3: do j = 1,3: write kdotp(:,:,3,i,j)`
   → 9·N² lines; nine N²-line blocks ordered (a,b) = (x,x),(x,y),(x,z),(y,x),(y,y),(y,z),(z,x),(z,y),(z,z)
   — **b varies fastest** (eV·Å²).

---

## F. Test: test-suite/tests/testpostw90_gaas_kdotp

Files in the test dir: `gaas.win`, `gaas.eig`, `gaas.mmn`, `gaas.amn`, `gaas.chk.fmt.bz2`
(+ `Makefile` that bunzips and runs `w90chk2chk.x -f2u` to produce `gaas.chk`), and the
benchmark `benchmark.out.default.inp=gaas.win`. Only `.chk` + `.eig` + `.win` are actually
consumed by the kdotp path (.mmn/.amn unused).

### F.1 .win keywords (`testpostw90_gaas_kdotp/gaas.win`)

```
berry = true
berry_task=kdotp
fermi_energy=7.7414
use_ws_distance = .false.
search_shells=12
kdotp_kpoint  =  0.0000 0.5000 0.0000       ! the L point in this cell's reduced coords
kdotp_num_bands = 2
kdotp_bands =  4,5
berry_kmesh = 1 1 1
sc_eta=0.040
num_bands = 12,  num_wann = 8,  exclude_bands : 1-5
dis_win_max = 24.0d0, dis_froz_max = 14.0d0, dis_num_iter = 1200, dis_mix_ratio = 1.d0
unit_cell_cart (bohr): (-5.34 0 5.34 / 0 5.34 5.34 / -5.34 5.34 0)
projections: As s,p at (¼,¼,¼); Ga p,s at (0,0,0)
mp_grid : 4 4 4  (+ explicit 64-point kpoints block)
```

(plus `num_iter=1000`, bands_plot keywords irrelevant to kdotp.)

### F.2 Test harness config

`test-suite/tests/jobconfig:527-530`:

```
[testpostw90_gaas_kdotp/]
program = POSTW90_KDOTP_OK
inputs_args = ('gaas.win', '')
output = gaas-kdotp_0.dat
```

→ **only the order-0 file is regression-checked**; kdotp_1/kdotp_2 are produced but not
compared.

`test-suite/tests/userconfig:177-181`:

```
[POSTW90_KDOTP_OK]
exe = ../../postw90.x
extract_fn = tools parsers.parse_kdotp_dat.parse
tolerance = (  (1.0e-4, 1.0e+2, 'real_part'),
               (1.0e-4, 1.0e+2, 'imag_part'))
```

testcode tolerance tuple order is `(absolute, relative, key)`
(`test-suite/testcode/lib/testcode2/config.py:40-53`, `validation.py:104-121`; strict=True).
With relative = 1e+2 the relative check is vacuous → effective criterion is
**|Δ| ≤ 1e-4 absolute on every real and imaginary part**. The parser
(`test-suite/tools/parsers/parse_kdotp_dat.py`) reads every non-blank, non-`#` line as two
floats into flat lists `real_part` / `imag_part` (order-preserving).

### F.3 Benchmark (`benchmark.out.default.inp=gaas.win`, complete file)

```
   0.64827435E+001   0.00000000E+000
   0.00000000E+000   0.00000000E+000
   0.00000000E+000   0.00000000E+000
   0.86209080E+001   0.00000000E+000
```

Interpretation (N=2, column-major): line 1 = (n=1,m=1) = E_band4(L) = 6.4827435 eV;
lines 2,3 = off-diagonal zeros; line 4 = (2,2) = E_band5(L) = 8.6209080 eV.

---

## G. Traps

1. **kdotp_bands indices are interpolated-band indices** (ascending eigenvalue order at
   kdotp_kpoint, 1..num_wann), not ab-initio band indices; `exclude_bands`/disentanglement act
   upstream.
2. **`kdotp_num_bands` missing → crash, not error**: `kdotp_bands` stays unallocated and
   `size(pw90_berry%kdotp_bands)` at berry.F90:508 is invoked on an unallocated array. There is
   also no check that `kdotp_bands ≤ num_wann`.
3. **Column-major file layout**: within each N² block the FIRST (bra) index varies fastest.
   For Hermitian order-1 blocks this is easy to get transposed silently — only detectable in
   the imaginary parts.
4. **Order-2 block ordering**: (a,b) blocks are written with b fastest
   (xx,xy,xz,yx,…); and a single (a,b) block is not Hermitian — `T_ab(n,m)* = T_ba(m,n)`.
   Symmetrize only via the κ_a κ_b contraction, never per-block, when validating.
5. **The ½ is stored**: kdotp_2 already contains ½(H̄_ab + T_ab); contract with κ_a κ_b
   without another ½.
6. **Gauge dependence**: orders 1 and 2 depend on the eigenvector phases/degenerate-subspace
   rotation chosen by ZHPEVX at kdotp_kpoint. Only order 0 (and |elements| in non-degenerate
   cases) is robustly comparable — consistent with the project's gauge-invariance validation
   policy; the reference test only compares order 0 with abs tol 1e-4.
7. **Unregularized denominators**: B-states degenerate with A-states at kdotp_kpoint make
   order 2 blow up; the code does not guard against it (`sc_eta` is NOT used here — it only
   feeds the dead `wham_get_D_h_P_value` call).
8. **Dead computations that must not be mistaken for inputs**: `HH_bar`, `D_h`, `eig_da` are
   computed inside `berry_get_kdotp` but never used for the output. A reimplementation can skip
   them (they only cost time; no file/stdout effect).
9. **fermi_energy is required** by `berry_main` even though kdotp never uses it
   (berry.F90:254-257).
10. **kdotp_kpoint is fractional**, despite the docs CSV claiming 2π/a units.
11. **Number format**: `E18.8E3` uses Fortran leading-zero mantissa `0.XXXXXXXXE±eee` with a
    3-digit exponent — Julia's default `@sprintf` scientific form (`6.48e0`) does not match;
    a custom formatter is needed for byte-identical files.
12. **Substring task matching**: `index(task,'kdotp')` — any berry_task string containing
    `kdotp` enables it, and it composes with other tasks in one run.
13. `use_ws_distance` defaults to TRUE in postw90; the test disables it. Its entire effect on
    kdotp comes through the HH_R remap in `operator_wigner_setup` (expanded R list
    `irvec_pw90`, weights 1/(ndegen·ndeg) folded into HH_R); the R→k routines themselves apply
    no per-orbital shifts.
