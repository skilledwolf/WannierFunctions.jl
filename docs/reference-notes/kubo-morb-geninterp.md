# postw90 reference notes: adaptive-kmesh/Fermi scan, Kubo optical conductivity, orbital magnetisation, geninterp

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Shared machinery (get_HH_R, get_AA_R, JJ± matrices, occupation matrices f/g, the imf
J0/J1/J2 trace formula, and the AHC assembly/units) is **already documented in
`docs/reference-notes/berry-ahc.md`** (§A get_HH_R/get_AA_R, §B wan_ham k-space machinery,
§C berry_get_imf(gh)_klist, §D berry_main AHC) and is only referenced here, not repeated.

Key papers cited by the code (`src/postw90/berry.F90:38-46`): WYSV06 = PRB 74, 195118 (2006);
YWVS07 = PRB 75, 195121 (2007) (Kubo); LVTS12 = PRB 85, 014435 (2012) (M_orb); CTVR06 = PRB 74,
024408 (2006).

Harness tolerance semantics (used throughout below): testcode tuples are
`(abs_tol, rel_tol, name)` — `test-suite/testcode/lib/testcode2/config.py:35-53`
(`'''Parse (abs_tol,rel_tol,name,strict).'''`). With both set and `strict=True` (default),
**both** absolute and relative must pass: `validation.py:114-…` Tolerance and
`Status.__add__` returns `max(status)` (`validation.py`, "Return the maximum level (ie most
"failed") status"). `validate_absolute`: `passed = |test-bench| < absolute`;
`validate_relative`: `passed = |diff/benchmark| < relative` (benchmark==0 & diff/=0 → Inf →
fail; benchmark==0 & diff==0 → 0). NOTE: berry-ahc.md §E states the tuple order as
"(relative, absolute)" — that is inverted; the correct order is (absolute, relative).

---

## 1. Adaptive refinement + Fermi-energy scan (berry_task=ahc, berry_main)

### 1.1 Input parameters and defaults

`src/postw90/postw90_types.F90:160-189` (`pw90_berry_mod_type`):

```fortran
integer :: curv_adpt_kmesh = 1
real(kind=dp) :: curv_adpt_kmesh_thresh = 100.0_dp
character(len=20) :: curv_unit = 'ang2' ! postw90/kpath, kslice as well
```

Parsed in `src/postw90/postw90_readwrite.F90:881-899`: keyword `berry_curv_adpt_kmesh`
(single integer, must be >= 1, error `'Error:  berry_curv_adpt_kmesh must be a positive
integer'`), `berry_curv_adpt_kmesh_thresh` (real), `berry_curv_unit` (must be `'ang2'` or
`'bohr2'`). So defaults: **adpt_kmesh = 1, thresh = 100.0, unit = ang2** (thresh is in the
units selected by `berry_curv_unit`, i.e. Å² by default — the trigger value is Ω in Å² and is
divided by bohr² only when `curv_unit == 'bohr2'`, see below).

### 1.2 Fermi-energy list (fermi_energy_min/max/step → fermi_n)

`src/readwrite.F90:609-692` `w90_readwrite_read_fermi_energy`:
- `fermi_energy` alone → n=1 list `[fermi_energy]` (lines 639-646). Default value 0.0 if
  nothing given (lines 634-637 → n=1, list=[0.0]).
- `fermi_energy_min` present (mutually exclusive with `fermi_energy`, lines 651-656):
  - `fermi_energy_max` default `fermi_energy_min + 1.0_dp` (line 662); must be > min.
  - `fermi_energy_step` default `0.01_dp` (line 673); must be positive.
  - `n = nint(abs((fermi_energy_max - fermi_energy_min)/fermi_energy_step)) + 1` (line 679)
  - `fermi_energy_step = (fermi_energy_max - fermi_energy_min)/real(n - 1, dp)` (line 680)
  - `fermi_energy_list(i) = fermi_energy_min + (i - 1)*fermi_energy_step` (lines 689-691).

`fermi_n = size(fermi_energy_list)`; berry_main errors if 0
(`berry.F90:254-257`: `'Must specify one or more Fermi levels when berry=true'`).
`berry_task=kubo` (and shc freq_scan) demand fermi_n == 1 (`berry.F90:340-347`,
`not_scannable = eval_kubo .or. (eval_shc .and. pw90_spin_hall%freq_scan)`); **ahc and morb
are scannable**.

### 1.3 Refinement sub-mesh placement

`berry.F90:593-626`. `adkpt(3, curv_adpt_kmesh**3)` offsets in reduced coordinates, where
`db1 = 1/mesh(1)` etc. (`berry.F90:271-273`):

```fortran
do i = 0, pw90_berry%curv_adpt_kmesh - 1
  do j = 0, pw90_berry%curv_adpt_kmesh - 1
    do k = 0, pw90_berry%curv_adpt_kmesh - 1
      ikpt = ikpt + 1
      adkpt(1, ikpt) = db1*((i + 0.5_dp)/pw90_berry%curv_adpt_kmesh - 0.5_dp)
      adkpt(2, ikpt) = db2*((j + 0.5_dp)/pw90_berry%curv_adpt_kmesh - 0.5_dp)
      adkpt(3, ikpt) = db3*((k + 0.5_dp)/pw90_berry%curv_adpt_kmesh - 0.5_dp)
```

i.e. an n³ *cell-centred* sub-grid spanning `[-db/2, +db/2)` **around** the triggering k
(works for even and odd n; the coarse point itself is a sub-point only for odd n). Sub-point
k-vectors are `kpt(:) + adkpt(:, loop_adpt)` (line 888).

Weights (regular-grid branch, `berry.F90:843-844`):

```fortran
kweight = db1*db2*db3
kweight_adpt = kweight/pw90_berry%curv_adpt_kmesh**3
```

(kpoint-file branch analogously `kweight_adpt = kweight/curv_adpt_kmesh**3`, line 649.)

### 1.4 Trigger norm and per-Fermi-energy logic (ahc)

Regular-grid loop `berry.F90:846-903` (identical logic in the wanint_kpoint_file branch,
lines 653-698). After `berry_get_imf_klist` for the coarse point:

```fortran
ladpt = .false.
do if = 1, fermi_n
  vdum(1) = sum(imf_k_list(:, 1, if))
  vdum(2) = sum(imf_k_list(:, 2, if))
  vdum(3) = sum(imf_k_list(:, 3, if))
  if (pw90_berry%curv_unit == 'bohr2') vdum = vdum/physics%bohr**2
  rdum = sqrt(dot_product(vdum, vdum))
  if (rdum > pw90_berry%curv_adpt_kmesh_thresh) then
    adpt_counter_list(if) = adpt_counter_list(if) + 1
    ladpt(if) = .true.
  else
    imf_list(:, :, if) = imf_list(:, :, if) + imf_k_list(:, :, if)*kweight
  end if
end do
```

(`berry.F90:867-880`). So the norm is: **sum J0+J1+J2** (first index 1:3 of
`imf_k_list(J, α, ife)`) per Cartesian pseudovector component α, then the **Euclidean norm of
the resulting 3-vector** (this is |Ω(k)| in Å², or bohr² after division by
`physics%bohr**2`, `bohr = 0.52917720859` for the default CODATA2006 build,
`src/constants.F90:182,217-224`).

- The **coarse value is discarded** for every Fermi index that triggers (only non-triggering
  `if` accumulate the coarse `imf_k_list*kweight`); it is **replaced** by the sub-mesh sum.
- Refinement runs once if ANY Fermi energy triggered (`if (any(ladpt)) then`, line 881), but
  the sub-mesh contributions are accumulated **only into the triggered Fermi indices**:

```fortran
do loop_adpt = 1, pw90_berry%curv_adpt_kmesh**3
  call berry_get_imf_klist(..., kpt(:) + adkpt(:, loop_adpt), ..., imf_k_list_dummy, ...,
                           ladpt=ladpt)
  do if = 1, fermi_n
    if (ladpt(if)) then
      imf_list(:, :, if) = imf_list(:, :, if) + imf_k_list_dummy(:, :, if)*kweight_adpt
    end if
  end do
end do
```

(`berry.F90:881-902`; a dummy array is used because "Using imf_k_list here would corrupt
values for other frequencies", lines 883-884). The `ladpt` mask propagates into
`berry_get_imfgh_klist` as `todo` so the J0/J1/J2 traces are only evaluated for triggered
Fermi indices (`berry.F90:1974-1978` `todo = ladpt`, and `if (todo(ife))` at 2029).

Note: for ahc there is **no** `curv_adpt_kmesh > 1` guard (unlike shc, lines 791-803, which
skips re-evaluation when adpt_kmesh==1 and triggers on `abs(shc_k_fermi(if))` with an `exit`
after the first triggering Fermi energy, updating all Fermi energies at once, lines 989-1031).

After the k-loop: `comms_reduce(imf_list…,'SUM')` and `comms_reduce(adpt_counter_list…,'SUM')`
(`berry.F90:1190-1195`).

### 1.5 Output: .wpout block and Fe-ahc-fermiscan.dat

Conversion (see berry-ahc.md §D for derivation, `berry.F90:1329-1361`):
`fac = -1.0e8_dp*physics%elem_charge_SI**2/(physics%hbar_SI*cell_volume)` (line 1360);
`ahc_list = imf_list*fac` (line 1361). CODATA2006 defaults:
`elem_charge_SI = 1.602176487e-19`, `hbar_SI = 1.054571628e-34` (`src/constants.F90:164-168`).

Stdout header for adaptive runs (`berry.F90:1244-1271`), printed when
`eval_ahc .and. curv_adpt_kmesh /= 1`:

```
Regular interpolation grid: <'(1x,a28,3(i0,1x))'>
Adaptive refinement grid:   n n n
Refinement threshold:  'Berry curvature >' F6.2 ' Ang^2'   (or ' bohr^2')
```

and if fermi_n == 1 also `' Points triggering refinement: ', I5, '(', F5.2, '%)'`
(count/product(mesh)*100).

Fermi-scan file (`berry.F90:1362-1416`), only when `fermi_n > 1`:

```fortran
file_name = trim(seedname)//'-ahc-fermiscan.dat'
...
do if = 1, fermi_n
  if (fermi_n > 1) write (file_unit, '(4(F12.6,1x))') &
    fermi_energy_list(if), sum(ahc_list(:, 1, if)), &
    sum(ahc_list(:, 2, if)), sum(ahc_list(:, 3, if))
```

i.e. one line per Fermi energy: `E_F  σ_x  σ_y  σ_z` in `4(F12.6,1x)` (σ = J0+J1+J2 sum in
S/cm). Simultaneously, per Fermi energy, the .wpout gets (lines 1377-1414):
`'(/,1x,a18,F10.4)') 'Fermi energy (ev):'`, then (fermi_n>1) the per-E_F
`' Points triggering refinement: '` line (format `'(1x,a30,i5,a,f5.2,a)'`), then
`'AHC (S/cm)       x          y          z'` and, for iprint<=1, the one-liner
`write (stdout, '(1x,a10,1x,3(f10.4,1x),/)') '==========', sum(ahc_list(:,1,if)), ...`
(iprint>1 additionally prints `'J0 term :'/'J1 term :'/'J2 term :'` rows, format
`'(1x,a9,2x,3(f10.4,1x))'`).

### 1.6 Test testpostw90_fe_ahc_adaptandfermi

`test-suite/tests/jobconfig:286-289`:

```
[testpostw90_fe_ahc_adaptandfermi/]
program = POSTW90_FERMISCAN_OK
inputs_args = ('Fe.win', '')
output = Fe-ahc-fermiscan.dat
```

Parser `tools/parsers/parse_fermiscan_dat.py`: skips `#` and blank lines, requires exactly 4
whitespace-separated fields per line → keys `fermienergy, ahc_x, ahc_y, ahc_z`.

Tolerances `tests/userconfig` `[POSTW90_FERMISCAN_OK]`:

```
tolerance = (  (1.0e-6, 5.0e-6, 'fermienergy'),
               (1.0e-4, 2.0e-4, 'ahc_x'),
               (1.0e-4, 2.0e-4, 'ahc_y'),
               (1.0e-4, 2.0e-4, 'ahc_z'))
```

(abs, rel; strict — both must hold.)

`Fe.win` deltas vs the fe_ahc test (same bcc-Fe 2.71175-bohr cell, num_bands=28, num_wann=18,
spinors, `use_ws_distance=.false.`, `search_shells=12`, mp_grid 2 2 2, `berry=true`,
`berry_task = ahc`, `berry_kmesh = 10 10 10`):

```
fermi_energy_min = 11.6279
fermi_energy_max = 13.6279
fermi_energy_step = 0.2
berry_curv_adpt_kmesh = 5 5 5        ! keyword read as single integer -> 5
berry_curv_adpt_kmesh_thresh = 10
uHu_formatted = .true.               ! irrelevant for ahc (no uHu read)
```

fermi_n = nint(2.0/0.2)+1 = 11. Full benchmark
(`benchmark.out.default.inp=Fe.win`, 11 lines):

```
   11.627900    50.845203   -50.619458   157.962756
   11.827900    20.541495   -20.578261   568.011017
   12.027900    73.684613   -73.701303   533.909220
   12.227900    32.332682   -32.350845   210.635328
   12.427900     1.291984    -1.290424    11.323568
   12.627900    41.926963   -41.913963   722.267035
   12.827900    35.739944   -35.798745  2804.956845
   13.027900    -5.553881     5.561349   124.550233
   13.227900     4.113229    -4.105007   142.046112
   13.427900    -1.363268     1.404604   -10.104289
   13.627900     0.001088     0.006733    90.262635
```

---

## 2. Kubo optical conductivity (berry_task=kubo)

### 2.1 Parameters

- Frequency grid (`src/postw90/postw90_readwrite.F90:1684-1724`): `kubo_freq_min` default
  `0.0_dp`, `kubo_freq_step` default `0.01_dp` (`postw90_readwrite.F90:56-58`);
  `kubo_freq_max` default: `froz_max − fermi_energy_list(1) + 0.6667` if frozen states, else
  `maxval(eigval) − minval(eigval) + 0.6667`, else `win_max − win_min + 0.6667` (1688-1694).
  Then

  ```fortran
  pw90_berry%kubo_nfreq = nint((kubo_freq_max - kubo_freq_min)/kubo_freq_step) + 1
  if (pw90_berry%kubo_nfreq <= 1) pw90_berry%kubo_nfreq = 2
  kubo_freq_step = (kubo_freq_max - kubo_freq_min)/(kubo_nfreq - 1)
  kubo_freq_list(i) = kubo_freq_min + (i-1)*(kubo_freq_max - kubo_freq_min)/(kubo_nfreq-1)
  ```

  (1708-1724; the list is complex but built real here).
- `kubo_eigval_max` default (`postw90_readwrite.F90:1758-1769`): `froz_max + 0.6667` if
  frozen states, else `maxval(eigval) + 0.6667`, else `win_max + 0.6667`.
- Smearing: `pw90_smearing_type` defaults (`postw90_types.F90:129-141`):
  `use_adaptive = .true.`, `adaptive_prefactor = sqrt(2.0_dp)`,
  `type_index = 0 ! Gaussian default`, `fixed_width = 0.0_dp !none`,
  `adaptive_max_width = 1.0_dp ! 1 eV`. Kubo-specific overrides
  (`postw90_readwrite.F90:908-960`): `kubo_adpt_smr`, `kubo_adpt_smr_fac`,
  `kubo_adpt_smr_max`, `kubo_smr_fixed_en_width`, `kubo_smr_type` — each defaulting to the
  global (`smr_type`/`adpt_smr`/`adpt_smr_fac`/`adpt_smr_max`/`smr_fixed_en_width`,
  readwrite 580-636). Smearing-type indices (`src/readwrite.F90:1864-1914`,
  `w90_readwrite_get_smearing_index`): `'m-v'`/`'cold'` → −1 (Marzari-Vanderbilt),
  `'m-p'`/`'m-pN'` → N (Methfessel-Paxton), `'f-d'` → −99 (Fermi-Dirac), `'gauss'` → 0.
  **Default kubo_smr_index = 0 = Gaussian.**
- **The Fe kubo tests use adaptive smearing** (their Fe.win sets none of the smearing
  keywords → adaptive, prefactor √2, max 1 eV, Gaussian type; see §2.5).

### 2.2 berry_main plumbing

`berry.F90:349-390`: kubo needs get_HH_R + get_AA_R only (plus get_SS_R when
`spin_decomp`); allocates `kubo_H_k/kubo_H/kubo_AH_k/kubo_AH(3,3,kubo_nfreq)` and
`jdos_k/jdos(kubo_nfreq)`, zeroed. Accumulation per k (regular grid branch):
`kubo_H = kubo_H + kubo_H_k*kweight` etc. (`berry.F90:943-945`), MPI-reduced at 1206-1212.
`fermi_n` must be 1 (lines 340-347).

### 2.3 berry_get_kubo_k (berry.F90:2145-2363)

Setup:
- Adaptive path (`2241-2250`): `wham_get_eig_deleig` (H, ∂H, eig, del_eig; see berry-ahc.md
  §B) and `Delta_k = pw90common_kmesh_spacing(pw90_berry%kmesh%mesh, recip_lattice)`
  where `kmesh_spacing_mesh` = `maxval(|b_i|/mesh(i))` (`postw90_common.F90:1013-1030`,
  b_i are the rows of recip_lattice, cf. 1002-1010).
- Fixed path (`2252-2259`): plain Fourier of HH_R (+derivatives) and diagonalisation.
- `call pw90common_get_occ(fermi_energy_list(1), eig, occ, num_wann)` (2261) — **T=0 step
  function**: `occ=1 if eig<ef else 0` (`postw90_common.F90:942-985`).
- `wham_get_D_h(delHH, D_h, UU, eig, num_wann)` (2263):
  `D_h(n,m,i) = (U†·∂H_i·U)(n,m)/(eig(m)−eig(n))`, zero for n==m or |ΔE|<1e-7
  (`wan_ham.F90:102-142`, Eq.(24) WYSV06).
- Band-gauge position matrix (2265-2273):

  ```fortran
  call pw90common_fourier_R_to_k_vec(..., AA_R, ..., OO_true=AA)
  do i = 1, 3
    AA(:, :, i) = utility_rotate(AA(:, :, i), UU, num_wann)
  end do
  AA = AA + cmplx_i*D_h ! Eq.(25) WYSV06
  ```

- Fixed-smearing trick (2275-2278): if not adaptive and fixed_width /= 0,
  `kubo_freq_list = real(kubo_freq_list) + cmplx_i*kubo_smearing%fixed_width` — the
  broadening enters the AH part as Im(ω).
- `spin_decomp` (2283-2304): `spin_get_nk` (src/postw90/spin.F90) gives band spin projection
  `spn_nk`; transitions classified ispn = 1 (up→up: both >= 0), 2 (dn→dn: both < 0),
  3 (spin-flip) — used only to accumulate the same integrand into `*_spn(:,:,ispn,:)`.

Main double band loop (2292-2361), skipping `n == m` and any pair with
`eig(m) > kubo_eigval_max .or. eig(n) > kubo_eigval_max` (2295):

```fortran
if (pw90_berry%kubo_smearing%use_adaptive) then
  ! Eq.(35) YWVS07
  vdum(:) = del_eig(m, :) - del_eig(n, :)
  joint_level_spacing = sqrt(dot_product(vdum(:), vdum(:)))*Delta_k
  eta_smr = min(joint_level_spacing*pw90_berry%kubo_smearing%adaptive_prefactor, &
                pw90_berry%kubo_smearing%adaptive_max_width)
else
  eta_smr = pw90_berry%kubo_smearing%fixed_width
end if
rfac1 = (occ(m) - occ(n))*(eig(m) - eig(n))
occ_prod = occ(n)*(1.0_dp - occ(m))
do ifreq = 1, pw90_berry%kubo_nfreq
  if (pw90_berry%kubo_smearing%use_adaptive) then
    omega = real(pw90_berry%kubo_freq_list(ifreq), dp) + cmplx_i*eta_smr
  else
    omega = pw90_berry%kubo_freq_list(ifreq)
  end if
  arg = (eig(m) - eig(n) - real(omega, dp))/eta_smr
  delta = utility_w0gauss(arg, pw90_berry%kubo_smearing%type_index, error, comm)/eta_smr
  jdos_k(ifreq) = jdos_k(ifreq) + occ_prod*delta
  cfac = cmplx_i*rfac1/(eig(m) - eig(n) - omega)
  rfac2 = -pi*rfac1*delta
  do j = 1, 3
    do i = 1, 3
      kubo_H_k(i, j, ifreq)  = kubo_H_k(i, j, ifreq)  + rfac2*AA(n, m, i)*AA(m, n, j)
      kubo_AH_k(i, j, ifreq) = kubo_AH_k(i, j, ifreq) + cfac*AA(n, m, i)*AA(m, n, j)
```

(`berry.F90:2305-2359`). So per (n→m) transition:
- σ^H integrand: `−π (f_m−f_n)(ε_m−ε_n) δ_η(ε_m−ε_n−ħω) A_nm,α A_mn,β` (real prefactor
  times complex dyad — kubo_H_k is complex).
- σ^AH integrand: `i (f_m−f_n)(ε_m−ε_n)/(ε_m−ε_n−ħω_c) A_nm,α A_mn,β`, with
  ħω_c = ω + iη (η = eta_smr adaptive, or Im added to the freq list for fixed smearing).
- JDOS integrand (jdos_k): `f_n (1−f_m) δ_η(ε_m−ε_n−ħω)` — **occupation product**, not
  (f_n−f_m); units eV⁻¹; no unit conversion is ever applied to jdos.

δ-representation: `utility_w0gauss(x, n)/eta_smr` (`src/utility.F90:1008-1091`):
n=0 Gaussian `exp(−min(200,x²))/√π`; n=−1 M-V cold
`1/√π · exp(−(x−1/√2)²)(2−√2 x)`; n=−99 Fermi-Dirac `1/(2+e^x+e^−x)` (|x|<=36);
n=1..10 Methfessel-Paxton recursion; n>10 or other negative → error.

### 2.4 Units, symmetrisation, and output files

`berry.F90:1501-1510`: `fac = 1.0e8_dp*physics%elem_charge_SI**2/(physics%hbar_SI*cell_volume)`
(**positive**, unlike AHC's −fac), `kubo_H = kubo_H*fac`, `kubo_AH = kubo_AH*fac` → S/cm
(same reasoning chain as AHC: Å⁻¹ × 10⁸ → cm⁻¹, times e²/ħ in Siemens).

Output (on root only). Symmetric components, `berry.F90:1519-1553` — "real (imaginary) part
is Hermitean (anti-Hermitean)":

```fortran
do n = 1, 6
  i = alpha_S(n); j = beta_S(n)                    ! (xx,yy,zz,xy,xz,yz), berry.F90:84-85
  file_name = trim(seedname)//'-kubo_S_'//achar(119 + i)//achar(119 + j)//'.dat'
  ...
  write (file_unit, '(3E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
    real(0.5_dp*(kubo_H(i, j, ifreq) + kubo_H(j, i, ifreq)), dp), &
    aimag(0.5_dp*(kubo_AH(i, j, ifreq) + kubo_AH(j, i, ifreq)))
```

Antisymmetric components, `berry.F90:1555-1589` — "real (imaginary) part is anti-Hermitean
(Hermitean)":

```fortran
do n = 1, 3
  i = alpha_A(n); j = beta_A(n)                    ! (yz,zx,xy), berry.F90:72-73
  file_name = trim(seedname)//'-kubo_A_'//achar(119 + i)//achar(119 + j)//'.dat'
  ...
  write (file_unit, '(3E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
    real(0.5_dp*(kubo_AH(i, j, ifreq) - kubo_AH(j, i, ifreq)), dp), &
    aimag(0.5_dp*(kubo_H(i, j, ifreq) - kubo_H(j, i, ifreq)))
```

(with spin_decomp, 9 columns `'(9E16.8)'`: the total pair followed by up-up, down-down,
spin-flip pairs). JDOS, `berry.F90:1591-1605`:

```fortran
file_name = trim(seedname)//'-jdos.dat'
...
write (file_unit, '(2E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), jdos(ifreq)
```

(spin_decomp: `'(5E16.8)'` with jdos_spn(1:3)). None of these files has a header line. The
.wpout merely lists the produced file names under
`'Output data files related to complex optical conductivity:'` (1512-1517).

### 2.5 Fe kubo tests

The three tests share a byte-identical `Fe.win` (verified by diff): bcc Fe, num_bands=28,
num_wann=18, spinors=T, `use_ws_distance = .false.`, `search_shells=12`,
`fermi_energy = 12.6279`, `uHu_formatted = .true.` (unused for kubo), `berry = true`,
`berry_task = kubo`, `kubo_freq_max = 7.0`, `kubo_freq_step = 0.5`,
`berry_kmesh = 10 10 10`, mp_grid 2 2 2. **No kubo_smr/adpt keywords → adaptive smearing
(default), prefactor √2, max width 1 eV, Gaussian.** nfreq = nint(7.0/0.5)+1 = 15
(frequencies 0.0, 0.5, …, 7.0); benchmarks have 15 rows.

`tests/jobconfig:346-361`:

```
[testpostw90_fe_kubo_jdos/]   program = POSTW90_JDOS_OK   output = Fe-jdos.dat
[testpostw90_fe_kubo_Axy/]    program = POSTW90_KUBO_OK   output = Fe-kubo_A_xy.dat
[testpostw90_fe_kubo_Szz/]    program = POSTW90_KUBO_OK   output = Fe-kubo_S_zz.dat
```

(all `inputs_args = ('Fe.win', '')`).

Benchmark first data lines:
- `testpostw90_fe_kubo_Axy/benchmark…`: `0.00000000E+00  0.30466638E+03 -0.12328027E-13`
  (second line `0.50000000E+00  0.11920666E+03  0.15615522E+03`)
- `testpostw90_fe_kubo_Szz/benchmark…`: `0.00000000E+00  0.17355951E+05 -0.14529221E-11`
- `testpostw90_fe_kubo_jdos/benchmark…`: `0.00000000E+00  0.71619828E+00`

Parsers: `parse_kubo_dat.py` needs exactly 3 fields → `energy, recond, imcond`;
`parse_jdos_dat.py` needs exactly 2 → `energy, jdos`. Tolerances (`tests/userconfig`):

```
[POSTW90_KUBO_OK]  tolerance = ( (1.0e-6, 5.0e-6, 'energy'),
                                 (1.0e-4, 1.0e+2, 'recond'),
                                 (1.0e-4, 1.0e+2, 'imcond'))
[POSTW90_JDOS_OK]  tolerance = ( (1.0e-6, 5.0e-6, 'energy'),
                                 (1.0e-4, 1.0e-4, 'jdos'))
```

(abs, rel; strict → for kubo the binding constraint is the absolute 1e-4 — rel 1e+2 is a
"disabled" relative check except that it fails when benchmark==0 and test/=0.)

---

## 3. Orbital magnetisation (berry_task=morb)

### 3.1 Operators required

`berry.F90:310-338`: morb calls get_HH_R, get_AA_R, then **get_BB_R** and **get_CC_R**.
Definitions (`get_oper.F90:800-801, 1105-1106`):

```
BB_a(R) = <0n|H(r-R)|Rm>   = FT of BB_a(k) = i<u|H|del_a u>       (a = x,y,z)
CC_ab(R) = <0|r_a.H.(r-R)_b|R> = FT of CC_ab(k) = <del_a u|H|del_b u>  (a,b = x,y,z)
```

`transl_inv` default `.false.`, `transl_inv_full` default `.false.`
(`postw90_types.F90:175-176`); both mutually exclusive
(`postw90_readwrite.F90:857-860`). **`transl_inv = T` is a hard error for morb**
(`berry.F90:556-559`: `'transl_inv=T disabled for morb'`). So by default get_BB_R/get_CC_R
apply **no** translational-invariance correction; the M-V band-diagonal correction
(transl_inv) exists only for AA_R (berry-ahc.md §A.3). The `transl_inv_full` scheme (r0
phases + real-space corrections below) IS supported for morb and is what the
fe_morb_transl_inv test exercises.

### 3.2 get_BB_R (get_oper.F90:795-1097)

Root re-reads `seedname.mmn` (formatted; header line + `nb nkp nntot`, checks each,
`get_oper.F90:908-929`). For each of the `num_kpts*nntot` blocks: reads
`ik ik2 nnl nnm nnn` then the `num_bands×num_bands` overlap `S_o(m,n)` (m inner; one
`(re,im)` pair per line, 938-944), matches the neighbour index nn (947-967), then

```fortran
call get_gauge_overlap_matrix(num_bands, num_wann, eigval, v_matrix, dis_manifold, &
                              ik, num_states(ik), kmesh_info%nnlist(ik, nn), &
                              num_states(kmesh_info%nnlist(ik, nn)), S_o, &
                              have_disentangled, H=H_q_qb)
```

(973-976). `get_gauge_overlap_matrix` (get_oper.F90:3235-3272) windows `S_o` with
`get_win_min` (first `lwindow=.true.` band per k, 3199-3232) and calls
`utility_zgemmm(v_matrix(1:ns_a,1:num_wann,ik_a),'C', S_o(win_a,win_b),'N',
v_matrix(1:ns_b,1:num_wann,ik_b),'N', S, eigval(win_a,ik_a), H)`; per
`src/utility.F90:189-249`, `S = V_a†·S_o·V_b` and `H = V_a†·diag(eigval_a)·S_o·V_b`. So
`H_q_qb = V(q)† diag(ε(q)) M(q,q+b) V(q+b)` = <ũ_q|H_q|ũ_{q+b}> (Wannier gauge; H acts via
the **bra** k-point's eigenvalues).

Accumulation (987-993):

```fortran
do idir = 1, 3
  nno = kmesh_info%nninv(nn, ik)
  BB_q_b(:, :, ik, nno, idir) = BB_q_b(:, :, ik, nno, idir) &
                                + cmplx_i*phase1(:, :)*kmesh_info%wb(nn)*kmesh_info%bk(idir, nn, ik) &
                                *H_q_qb(:, :)
end do
```

i.e. `BB_a(q) = i Σ_b w_b b_a <ũ_q|H_q|ũ_{q+b}>` (the finite-difference of i<u|H|∂_a u>).
`phase1 = 1` unless transl_inv_full, in which case
`phase1 = exp(+i r0·b)` with `r0(i,j,:) = (⟨r⟩_i + ⟨r⟩_j)/2` from
`wannier_centres_from_AA_R` (978-985, 890-898). `nninv` reorders shell indices so
`bk(:,nninv(nn,ik),1) = bk(:,nn,ik)` (`src/types.F90:198`); harmless in the default path
because `BB_q = sum(BB_q_b, 4)` (1056). Default path then Fourier-transforms q→R
(`fourier_loc_q_to_R`, phase `exp(-i R·q)/num_kpts`, get_oper.F90:3190-3194) and applies
degeneracy/ws reordering via `operator_wigner_setup` (3275-3327). transl_inv_full extra
steps (1016-1051): per-shell R-space phase `phase2 = exp(-i(R_c·b)/2)` (1027-1032) and the
final centre correction

```fortran
BB_R(:, :, ir, idir) = BB_R(:, :, ir, idir) + &
                       (r0(:, :, idir) - 0.5_dp*wigner_seitz%crvec_pw90(idir, ir))*HH_R(:, :, ir)
```

(1045-1050). Scissors shift is fatal for BB_R (879-882). Result broadcast (1085).

### 3.3 get_CC_R and the .uHu file (get_oper.F90:1100-1456)

.uHu open/read (1214-1242), controlled by `pw90_oper_read%uHu_formatted` (default
`.false.`, `postw90_types.F90:70-71`; keyword `uHu_formatted`,
`postw90_readwrite.F90:483`):

- **Formatted**: `open(file=trim(seedname)//".uHu", form='formatted')`;
  `read (uHu_in, *) header` (line 1219); `read (uHu_in, *) nb_tmp, nkp_tmp, nntot_tmp`
  (1221).
- **Unformatted**: `read (uHu_in) header` with `character(len=60) :: header` (1227,
  1161); `read (uHu_in) nb_tmp, nkp_tmp, nntot_tmp` (1229).
- Checks: nb=num_bands, nkp=num_kpts, nntot=kmesh_info%nntot (1231-1242, errors
  `'…has not the right number of bands/k-points/nearest neighbours'`).

Block loop order (1245-1309): `ik = 1..num_kpts` (outer), `nn2 = 1..nntot`,
`nn1 = 1..nntot` (inner); `qb1 = nnlist(ik,nn1)`, `qb2 = nnlist(ik,nn2)`. Per block —
comment at 1254: "Read from .uHu file the matrices <u_{q+b1}|H_q|u_{q+b2}> between the
original ab initio eigenstates":

```fortran
if (pw90_oper_read%uHu_formatted) then
  do m = 1, num_bands
    do n = 1, num_bands
      read (uHu_in, *, err=106, end=106) c_real, c_img
      Ho_qb1_q_qb2(n, m) = cmplx(c_real, c_img, dp)
    end do
  end do
else
  read (uHu_in, err=106, end=106) &
    ((Ho_qb1_q_qb2(n, m), n=1, num_bands), m=1, num_bands)
end if
! pw2wannier90 is coded a bit strangely, so here we take the transpose
Ho_qb1_q_qb2 = transpose(Ho_qb1_q_qb2)
```

(1257-1269). So after the transpose, `Ho_qb1_q_qb2(m,n) = <u_{m,q+b1}|H_q|u_{n,q+b2}>` — the
**first (row/bra) index belongs to q+b1 and is the conjugated side**. Gauge transform
(1280-1282) passes `H_qb1_q_qb2` as the `S` (prod1) argument — no eigval multiplication,
since H is already inside uHu: `H_qb1_q_qb2 = V(qb1)†·Ho·V(qb2)`.

Accumulation (1297-1305), upper triangle a<=b only:

```fortran
do b = 1, 3
  do a = 1, b
    nn1o = kmesh_info%nninv(nn1, ik)
    nn2o = kmesh_info%nninv(nn2, ik)
    CC_q_b(:, :, ik, nn1o, nn2o, a, b) = CC_q_b(:, :, ik, nn1o, nn2o, a, b) &
                                         + phase1(:, :)*kmesh_info%wb(nn1)*kmesh_info%bk(a, nn1, ik) &
                                         *kmesh_info%wb(nn2)*kmesh_info%bk(b, nn2, ik)*H_qb1_q_qb2(:, :)
```

i.e. `CC_ab(q) = Σ_{b1,b2} w_{b1} b1_a w_{b2} b2_b <ũ_{q+b1}|H_q|ũ_{q+b2}>` — **no i
factors** ((−i)(+i)=1 from bra/ket derivatives); **b1 ↔ row/bra/conjugated ↔ first
Cartesian index a**; b2 ↔ column/ket ↔ b. `phase1 = 1` by default; transl_inv_full:
`phase1 = exp(i r0·(b2 − b1))` (1284-1295). Default path completes the lower triangle in q
space by Hermiticity (1403-1415):

```fortran
CC_q = sum(sum(CC_q_b, 5), 4)
do b = 1, 3
  do a = 1, b
    do ik = 1, num_kpts
      CC_q(:, :, ik, b, a) = conjg(transpose(CC_q(:, :, ik, a, b)))
```

then q→R FT for all 9 (a,b) (1425-1437). transl_inv_full path (1317-1399) instead FTs each
(nn1,nn2,a,b) with R-space phase `phase2 = exp(-i R_c·(b1+b2)/2)` (1353-1360) and adds three
correction terms per R (1377-1398):

```fortran
CC_R(:, :, ir, a, b) = CC_R(:, :, ir, a, b) + (r0(:, :, a) + 0.5_dp*crvec(a, ir))*BB_R(:, :, ir, b)
! + for the R' = -R partner index ir2:
CC_R(:, :, ir, a, b) = CC_R(:, :, ir, a, b) + conjg(transpose(BB_R(:, :, ir2, a)))* &
                                              (r0(:, :, b) - 0.5_dp*crvec(b, ir))
CC_R(:, :, ir, a, b) = CC_R(:, :, ir, a, b) + (r0(:, :, a) + 0.5_dp*crvec(a, ir))* &
                                              crvec(b, ir)*HH_R(:, :, ir)
```

(requires allocated HH_R and BB_R, 1317-1324). Scissors fatal (1186-1189). Both paths
broadcast CC_R (1444).

The .uHu is **not** needed by fe_kubo/* at runtime, but those test dirs ship `Fe.uHu.bz2`
anyway; morb tests need it. All Fe tests set `uHu_formatted = .true.`.

### 3.4 k-space assembly: berry_get_imfgh_klist (berry.F90:1873-2139)

Docstring (1882-1890): computes `-2Im[f(k)]` [Eq.33 CTVR06 / Eq.6 LVTS12], `-2Im[g(k)]`
[Eq.34/7], `-2Im[h(k)]` [Eq.35/8], stored axial-vector form. imf J0/J1/J2 (2025-2049) is
in berry-ahc.md §C. For img/imh (only when BOTH `img_k_list` and `imh_k_list` are present,
2054): BB(k), CC(k) via Fourier —

```fortran
call pw90common_fourier_R_to_k_vec(..., BB_R, ..., OO_true=BB)
do j = 1, 3
  do i = 1, j
    call pw90common_fourier_R_to_k(..., CC(:, :, i, j), CC_R(:, :, :, i, j), ...)
    CC(:, :, j, i) = conjg(transpose(CC(:, :, i, j)))
```

(2065-2079). Workspaces (2058-2063): `tmp(:,:,1) = HH·AA(:,:,alpha_A(i))`,
`tmp(:,:,2) = LLambda` (Eq.(37) LVTS12 as pseudovector), `tmp(:,:,3) = HH·OOmega(:,:,i)`.
Then per pseudovector component i (α=alpha_A(i), β=beta_A(i); (α,β) = (y,z),(z,x),(x,y)):

```fortran
call utility_zgemm_new(HH, AA(:, :, alpha_A(i)), tmp(:, :, 1))
call utility_zgemm_new(HH, OOmega(:, :, i), tmp(:, :, 3))
tmp(:, :, 2) = cmplx_i*(CC(:, :, alpha_A(i), beta_A(i)) &
                        - conjg(transpose(CC(:, :, alpha_A(i), beta_A(i)))))
do ife = 1, nfermi_loc
  ! J0 terms for -2Im[g] and -2Im[h]
  ! tmp(:,:,5) = HH . AA(:,:,alpha_A(i)) . f_list(:,:,ife) . AA(:,:,beta_A(i))
  call utility_zgemm_new(tmp(:, :, 1), f_list(:, :, ife), tmp(:, :, 4))
  call utility_zgemm_new(tmp(:, :, 4), AA(:, :, beta_A(i)), tmp(:, :, 5))
  s = 2.0_dp*utility_im_tr_prod(f_list(:, :, ife), tmp(:, :, 5))
  img_k_list(1, i, ife) = utility_re_tr_prod(f_list(:, :, ife), tmp(:, :, 2)) - s
  imh_k_list(1, i, ife) = utility_re_tr_prod(f_list(:, :, ife), tmp(:, :, 3)) + s
  ! J1 terms
  call utility_zgemm_new(HH, JJm_list(:, :, ife, alpha_A(i)), tmp(:, :, 4))
  img_k_list(2, i, ife) = -2.0_dp* &
      ( utility_im_tr_prod(JJm_list(:, :, ife, alpha_A(i)), BB(:, :, beta_A(i))) &
      - utility_im_tr_prod(JJm_list(:, :, ife, beta_A(i)), BB(:, :, alpha_A(i))) )
  imh_k_list(2, i, ife) = -2.0_dp* &
      ( utility_im_tr_prod(tmp(:, :, 1), JJp_list(:, :, ife, beta_A(i))) &
      + utility_im_tr_prod(tmp(:, :, 4), AA(:, :, beta_A(i))) )
  ! J2 terms
  call utility_zgemm_new(JJm_list(:, :, ife, alpha_A(i)), HH, tmp(:, :, 4))
  call utility_zgemm_new(HH, JJm_list(:, :, ife, alpha_A(i)), tmp(:, :, 5))
  img_k_list(3, i, ife) = -2.0_dp*utility_im_tr_prod(tmp(:, :, 4), JJp_list(:, :, ife, beta_A(i)))
  imh_k_list(3, i, ife) = -2.0_dp*utility_im_tr_prod(tmp(:, :, 5), JJp_list(:, :, ife, beta_A(i)))
```

(`berry.F90:2084-2135`; comments: "Trace formula for -2Im[g], Eq.(66) LVTS12", "Trace
formula for -2Im[h], Eq.(56) LVTS12"). In words (H, A, Ω, f, J± all Wannier-gauge; tr
products from berry-ahc.md §C):

- img J0 = Re tr[f·Λ_i] − 2 Im tr[f·H·A_α·f·A_β], with Λ_i = i(CC_αβ − CC_αβ†)
- imh J0 = Re tr[f·H·Ω_i] + 2 Im tr[f·H·A_α·f·A_β]
- img J1 = −2( Im tr[J⁻_α·B_β] − Im tr[J⁻_β·B_α] )
- imh J1 = −2( Im tr[(H·A_α)·J⁺_β] + Im tr[(H·J⁻_α)·A_β] )
- img J2 = −2 Im tr[(J⁻_α·H)·J⁺_β]
- imh J2 = −2 Im tr[(H·J⁻_α)·J⁺_β]

morb has NO adaptive refinement: plain accumulation `imf_list2/img_list/imh_list +=
(imf/img/imh)_k_list*kweight` (`berry.F90:905-918`), reduced at 1197-1204.

### 3.5 M_orb assembly, units, output (berry.F90:1419-1495)

Unit derivation comment (1420-1445): "At this point X=img_ab(:)-fermi_energy*imf_ab(:) …
X(k)=-2*Im[g(k)-E_F.f(k)] … \tilde{M}^LC=-(e/2.hbar) int dk/(2.pi)^3 X(k) dk … (i) The
summand is an energy in eV times a Berry curvature in Ang^2. To convert to a.u., divide by
27.2 and by 0.529^2 (ii) Multiply by -(e/2.hbar)=-1/2 in atomic units (iii) … 1 Bohr
magneton = 1/2 atomic unit, so need to multiply by 2".

```fortran
fac = -physics%eV_au/physics%bohr**2
...
do if = 1, fermi_n
  LCtil_list(:, :, if) = (img_list(:, :, if) - fermi_energy_list(if)*imf_list2(:, :, if))*fac
  ICtil_list(:, :, if) = (imh_list(:, :, if) - fermi_energy_list(if)*imf_list2(:, :, if))*fac
  Morb_list(:, :, if) = LCtil_list(:, :, if) + ICtil_list(:, :, if)
```

(1447, 1459-1464). So **M = −(eV_au/bohr²)·[img + imh − 2·E_F·imf]** per J-term and
component, in **Bohr magnetons per cell** (NOT A/m). CODATA2006 defaults:
`eV_au = 3.674932540e-2` (`src/constants.F90:178`), `bohr = 0.52917720859`
(`bohr_angstrom_internal`, `src/constants.F90:182`; CODATA2006 selected by the `#define`
fallback at constants.F90:92-100; `bohr = bohr_angstrom_internal` unless
USE_WANNIER90_V1_BOHR, constants.F90:217-224). Note fac
combines the −1/2·(1/eV_au·bohr²)⁻¹·2 chain: −(1/2)·eV_au/bohr²·2 = −eV_au/bohr².

.wpout block (1465-1494): if fermi_n>1 a file `seedname-morb-fermiscan.dat` is written with
`'(4(F12.6,1x))') fermi_energy_list(if), sum(Morb_list(1:3,1,if)), sum(…,2,if)),
sum(…,3,if))`. stdout per Fermi energy:

```fortran
write (stdout, '(/,/,1x,a,F12.6)') 'Fermi energy (ev) =', fermi_energy_list(if)
write (stdout, '(/,/,1x,a)') 'M_orb (bohr magn/cell)        x          y          z'
! iprint <= 1:
write (stdout, '(1x,a22,2x,3(f10.4,1x),/)') '======================', &
  sum(Morb_list(1:3, 1, if)), sum(Morb_list(1:3, 2, if)), sum(Morb_list(1:3, 3, if))
```

(iprint>1 instead prints `'Local circulation :'` = ΣLCtil, `'Itinerant circulation:'` =
ΣICtil, then `'Total   :'`, 1472-1486). The sums run over the J0+J1+J2 index.

### 3.6 fe_morb tests

`tests/jobconfig:292-307`: all three use `program = POSTW90_WPOUT_OK`,
`inputs_args = ('Fe.win','')`, `output = Fe.wpout`.
`tests/userconfig [POSTW90_WPOUT_OK]` (parser `parse_wpout.py`, regex
`^\s*======================\s+(...)\s+(...)\s+(...)` after a line containing `M_orb`):

```
(1.0e-3, 2.0e-3, 'morb_x'), (1.0e-3, 2.0e-3, 'morb_y'), (1.0e-3, 2.0e-3, 'morb_z')
```

(abs, rel; the same program also checks ahc_x/y/z (1e-3, 2e-3) and spin fields).

- `testpostw90_fe_morb/Fe.win`: same Fe system as the ahc test; `use_ws_distance = .false.`,
  `fermi_energy = 12.6279`, `uHu_formatted = .true.`, `berry = true`, `berry_task = morb`,
  `berry_kmesh = 10 10 10`. Benchmark (`benchmark…:241-245`):

  ```
   Fermi energy (ev) =   12.627900
   M_orb (bohr magn/cell)        x          y          z
   ======================      0.0000     0.0000     0.0431
  ```

- `testpostw90_fe_morb_transl_inv/Fe.win` differs from fe_morb ONLY by (verified diff):
  `use_ws_distance = .true.` and added `transl_inv_full = .true.`. Benchmark
  (`benchmark…:224-228`): `Fermi energy (ev) =   12.627900` …
  `======================      0.0000    -0.0000     0.0415`.
- `testpostw90_fe_morb_transl_inv_higher/Fe.win` is a rewritten input (cell 2.71 bohr,
  projections `Fe:dxz;dyz;dxy` + `Fe:sp3d2`, `kmesh = 10 10 10`, `fermi_energy = 12.6631`,
  `transl_inv = F`, `transl_inv_full=T`, `higher_order_n=2`, `use_ws_distance = T`).
  Benchmark: `======================      0.0000    -0.0000    -0.0617` at E_F 12.6631.

---

## 4. geninterp (src/postw90/geninterp.F90)

### 4.1 Keywords

`postw90_readwrite.F90:436-437`: `geninterp` (logical, activates the module);
`postw90_readwrite.F90:1354-1379`: `geninterp_alsofirstder`, `geninterp_single_file`.
Defaults (`postw90_types.F90:230-236`): `alsofirstder = .false.`, `single_file = .true.`.

### 4.2 Input file seedname_geninterp.kpt (geninterp.F90:179-201, 302-321)

Read on root only (formatted, status='old'):
1. Line 1: free comment, `read (kpt_unit, '(A500)') commentline` (183) — echoed to output.
2. Line 2: coordinate-type token `cdum` (184-198): contains `'crystal'` or `'frac'` →
   `absoluteCoords = .false.`; contains `'cart'` or `'abs'` → `.true.`; else error
   `'Error on second line of file …_geninterp.kpt: unable to recognize keyword'`.
3. Line 3: `read (kpt_unit, *) nkinterp` — number of k-points (200-201).
4. Then nkinterp lines: `read (kpt_unit, *) kpointidx(i), kpt` — one **integer identifier**
   plus 3 reals (305-306). Fractional coords used as-is; cartesian (units of 2π/Å, i.e.
   Å⁻¹) are converted via
   `kpoints(j,i) = Σ_l real_lattice(j,l)*kpt(l); kpoints(:,i) = kpoints(:,i)/(2π)`
   (308-318 — "I use the real_lattice (transposed) and a factor of 2pi instead of inverting
   again recip_lattice").

### 4.3 Computation

`get_HH_R` once (226-230). Work is scattered over MPI ranks (`comms_array_split`,
`comms_scatterv`, 283-335). Per k-point (360-382):
- `alsofirstder = .true.`: `wham_get_eig_deleig` (wan_ham.F90:442-543) → HH(k) by Fourier
  interpolation, `utility_diagonalize`, plus `delHH(:,:,a)` (alpha=1,2,3 Fourier with iR_a
  factor) and `wham_get_deleig_a` per direction (wan_ham.F90:342-439):
  non-degenerate `deleig_a(i) = real((U†·∂H_a·U)(i,i))` [Eq.(27) YWVS07]; with
  `use_degen_pert` (default .false., degen_thr 1e-4, postw90_types.F90:97-104) degenerate
  groups (ΔE < degen_thr) diagonalize the sub-block of U†·∂H_a·U [Eq.(31) YWVS07]. Units:
  eV·Å (R in Å).
- else: `pw90common_fourier_R_to_k` (alpha=0) + `utility_diagonalize` only (374-379).

### 4.4 Output seedname_geninterp.dat

`single_file = .true.` (default): root gathers (`comms_gatherv`) and writes
`trim(seedname)//'_geninterp.dat'` (338-344). Else each rank writes
`trim(seedname)//'_geninterp_'//I5.5 rank//'.dat'` (`'(a,a,I5.5,a)'`, or `I0` if >99999
nodes; 346-351).

Header (`internal_write_header`, geninterp.F90:53-82):

```fortran
write (outdat_unit, '(A)') "# Written on "//cdate//" at "//ctime ! Date and time
write (outdat_unit, '(A)') "# Input file comment: "//trim(commentline)
! alsofirstder:
write (outdat_unit, '(A)') "#  Kpt_idx  K_x (1/ang)       K_y (1/ang)        K_z (1/ang)       Energy (eV)"// &
  "      EnergyDer_x       EnergyDer_y       EnergyDer_z"
! else the same line without the EnergyDer_* fields
```

Data lines (398-419, and 424-445 for multi-file): k printed in **Cartesian Å⁻¹**,
`frac(j) = Σ_l recip_lattice(l,j)*kpt(l)` (403-406; b-vectors are rows of recip_lattice),
one line per band (enidx = 1..num_wann), bands innermost:

```fortran
write (outdat_unit, '(I10,7G18.10)') kpointidx(i), frac, globaleig(enidx, i), globaldeleig(enidx, :, i)
! or without derivatives:
write (outdat_unit, '(I10,4G18.10)') kpointidx(i), frac, globaleig(enidx, i)
```

Note the k-point identifier written is the user-supplied integer from the .kpt file, not a
sequence number. Error handlers: 105/106/107 for open/read .kpt and open .dat
(471-476).

### 4.5 Test testpostw90_si_geninterp

`tests/jobconfig:262-265`:

```
[testpostw90_si_geninterp/]
program = POSTW90_GENINTERPDAT_OK
inputs_args = ('silicon.win', '')
output = silicon_geninterp.dat
```

(`testpostw90_si_geninterp_wsdistance` is identical except silicon.win has
`use_ws_distance = true` instead of `.false.`.)

`silicon.win` keys: `num_bands = 12`, `num_wann = 8`, `use_ws_distance = .false.`,
`search_shells=12`, **`geninterp = true`**, **`geninterp_alsofirstder = true`**,
`dis_win_max = 17.0`, `dis_froz_max = 6.4`, Si fcc cell (a/2 = 2.6988 Å entries),
`mp_grid = 4 4 4`, projections `Si : sp3`. No `geninterp_single_file` → single file.

`silicon_geninterp.kpt` (complete):

```
Sample points for testing implementation
crystal
3
1 0. 0.2 0.
2 0.15 0.15 0.15
3 0. 0.88 0.88
```

Benchmark (`benchmark.out.default.inp=silicon.win`, 27 lines = 3 header + 3 kpts × 8 bands);
first data line:

```
         1  0.2328140398      0.2328140398      0.2328140398      -5.137728841      0.7721365899      0.7720555811      0.7719609989
```

Parser `parse_geninterp_dat.py`: skips `#` lines; accepts 5 or 8 fields →
`bandidx, bandkptx/y/z, bandenergy[, bandderivx/y/z]`. Tolerances
(`tests/userconfig [POSTW90_GENINTERPDAT_OK]`):

```
tolerance = ( (1.0e-3, 5.0e-3, 'bandenergy'),
              (1.0e-6, None, 'bandkptx'),
              (1.0e-6, None, 'bandkpty'),
              (1.0e-6, None, 'bandkptz'),
              (1.0e-2, 1.0e-2, 'bandderivx'),
              (1.0e-2, 1.0e-2, 'bandderivy'),
              (1.0e-2, 1.0e-2, 'bandderivz'),
              (1.0e-6, None, 'bandidx'))
```

(abs, rel; `None` = no relative check — comment in userconfig: "ignore relative error
comparison of values ~zero").

---

## Implementation checklist (condensed)

1. **Fermi list**: min/max/step → n = nint(|Δ|/step)+1, step recomputed exactly; kubo
   requires n==1; ahc/morb scan all n.
2. **Adaptive (ahc)**: trigger |Σ_J imf_k(J,·,ife)|₂ > thresh (ang² default; /bohr² if
   curv_unit=bohr2); per-ife discard-and-replace with n³ cell-centred sub-mesh of weight
   kweight/n³; count triggers per ife; defaults n=1, thresh=100, unit ang2.
3. **Kubo**: A_band = U†A_WU + iD_h; σ^H += −π(f_m−f_n)(ε_m−ε_n)δ_η(ε_mn−ω)A_nm,iA_mn,j;
   σ^AH += i(f_m−f_n)(ε_m−ε_n)/(ε_mn−ω−iη)A_nm,iA_mn,j; jdos += f_n(1−f_m)δ_η;
   η adaptive = min(√2·|∇ε_m−∇ε_n|·Δk, 1 eV) by default (Gaussian w0gauss);
   ω grid 0..max step 0.01 default; cutoff kubo_eigval_max; ×1e8·e²/(ħV_c) → S/cm;
   files -kubo_S_ab.dat (Re→H sym, Im→AH sym), -kubo_A_ab.dat (Re→AH asym, Im→H asym),
   -jdos.dat, all `E16.8`, no headers.
4. **morb**: BB_R from .mmn (i·w_b·b_a·V†diag(ε)MV), CC_R from .uHu
   (w_b1 b1_a w_b2 b2_b V(qb1)†·uHu^T·V(qb2), a<=b + Hermitian completion; bra=b1
   conjugated); img/imh J0/J1/J2 traces above; M = −(eV_au/bohr²)(im{g,h}−E_F·imf) summed →
   μ_B/cell; transl_inv forbidden, transl_inv_full optional (r0 phases + 3 R-space
   corrections).
5. **geninterp**: .kpt = comment / crystal|frac|cart|abs / N / (idx kx ky kz)×N; output
   header 3 `#` lines; data `(I10,4G18.10)` or `(I10,7G18.10)` with k in Cartesian Å⁻¹;
   E_n eV; dE/dk eV·Å via U†∂HU diagonal (degenerate-perturbation optional).
