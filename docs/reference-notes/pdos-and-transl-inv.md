# postw90 reference notes: Projected DOS (dos module) and transl_inv / transl_inv_full

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Key paper: YWVS07 = Yates, Wang, Vanderbilt, Souza, PRB 75, 195121 (2007) — adaptive smearing;
MV97 = Marzari & Vanderbilt, PRB 56, 12847 (1997) — Eq.(31) Im-log band-diagonal position element.

---

# Part 1 — DOS / projected DOS (`dos = true`, src/postw90/dos.F90)

## 1.A Control flow (`dos_main`, dos.F90:54-368)

Called from `src/postw90/postw90.F90:438-439` only when
`pw90_calcs%dos .and. index(dos_data%task,'dos_plot') > 0`. The `find_fermi_energy` task is
parsed but the implementation is **commented out** (dos.F90:377-557, postw90.F90:447-448) — it
does nothing; setting it together with `fermi_energy` is an input error
(postw90_readwrite.F90:1250-1253).

Inputs needed: `HH_R` from `get_HH_R` (dos.F90:176-179) — built from chk (`u_matrix`,
`v_matrix`) + `seedname.eig`; if `spin_decomp=T` also `SS_R` from `get_SS_R`
(dos.F90:183-193), which reads **`seedname.spn`**. **No `.mmn` and no `.amn` are ever read by
the dos module** (get_AA_R is never called here) — this is why the pdos test ships only
`.chk.fmt + .eig + .win`.

Energy grid (dos.F90:146-158):

```
num_freq = nint((energy_max - energy_min)/energy_step) + 1;  if (num_freq==1) num_freq=2
d_omega  = (energy_max - energy_min)/(num_freq - 1)          ! step is RECOMPUTED
dos_energyarray(ifreq) = energy_min + (ifreq-1)*d_omega      ! ifreq = 1..num_freq
```

`ndim = 3` if `spin_decomp` else `1` (dos.F90:183-193); `dos_k(num_freq, ndim)`,
`dos_all(num_freq, ndim)`.

k-sampling (two branches):
- `wanint_kpoint_file=T` (dos.F90:237-281): irreducible-BZ points from `kpoint.dat`, each with
  its own weight; `dos_all += dos_k * kpoint_dist%weight(loop_tot)`.
- default full BZ (dos.F90:283-332): flat MPI-strided loop
  `loop_tot = my_node_id, product(mesh)-1, num_nodes`, decomposed with **z fastest**
  (lines 289-296), `kpt(i) = loop_{x,y,z}/mesh(i)` (unshifted Γ-centred grid),
  `kweight = 1/product(pw90_dos%kmesh%mesh)` (line 287), `dos_all += dos_k*kweight` (line 329).

Per k-point: if `smearing%use_adaptive` → `wham_get_eig_deleig` (eig + band velocities;
wan_ham.F90:442-543) + `dos_get_levelspacing` + `dos_get_k(..., levelspacing_k, UU)`;
else → `pw90common_fourier_R_to_k` (alpha=0) + `utility_diagonalize` + `dos_get_k(..., UU)`
(dos.F90:248-279 / 297-328). Finally `comms_reduce(dos_all,'SUM')` (line 336) and root writes
the file (lines 339-350).

`UU` from `utility_diagonalize` (src/utility.F90:652-696, ZHPEVX 'V','A','U'):
`H(k) = UU · diag(eig) · UU†`, i.e. **column** `UU(:,i)` is the i-th eigenvector expressed in
the Wannier basis: `UU(iw, ib) = <W_iw | psi_ib^(H)>`, eigenvalues ascending.

## 1.B Per-k contribution (`dos_get_k`, dos.F90:600-765)

For each interpolated band `i = 1..num_wann` (eigenvalue `eig_k(i)` in eV):

Smearing width (dos.F90:690-696):

```
fixed:    eta_smr = smearing%fixed_width                                       ! line 691
adaptive: eta_smr = min(levelspacing_k(i)*smearing%adaptive_prefactor,
                        smearing%adaptive_max_width)                           ! Eq.(35) YWVS07, line 694
```

with (dos_get_levelspacing, dos.F90:768-795 + postw90_common.F90:1013-1029):

```
Delta_k        = max_i( |b_i| / mesh(i) )        ! b_i = recip_lattice row i, Ang^-1
levelspacing_k(band) = sqrt(dot_product(del_eig(band,:), del_eig(band,:))) * Delta_k
```

`del_eig(band, 1:3) = dE/dk_a` in eV·Å from `wham_get_deleig_a` (wan_ham.F90:342-439):
`Re diag(UU† dH_a UU)`; with `use_degen_pert=T` (default **F**, `degen_thr=1e-4`,
postw90_types.F90:102-103) degenerate groups (gap < degen_thr) re-diagonalised.

Histogram-vs-smearing switch and bin range (dos.F90:698-715), with
`binwidth = EnergyArray(2)-EnergyArray(1)`, `smearing_cutoff = 10._dp` and
`min_smearing_binwidth_ratio = 2._dp` (src/constants.F90:80,82):

```
if (eta_smr/binwidth < 2.0) then              ! histogram: only nearest bin
  min_f = max(nint((eig_k(i)-E(1))/(E(N)-E(1))*(N-1)) + 1, 1)
  max_f = min(same, N);   DoSmearing = .false.
else                                          ! smear over eig_k(i) +/- 10*eta_smr
  min_f = max(nint((eig_k(i)-10*eta_smr-E(1))/(E(N)-E(1))*(N-1)) + 1, 1)
  max_f = min(nint((eig_k(i)+10*eta_smr-E(1))/(E(N)-E(1))*(N-1)) + 1, N)
  DoSmearing = .true.
end if
```

(If eig is far outside the window, max_f < min_f and the band contributes nothing.)

Delta-function value per bin (dos.F90:717-725):

```
DoSmearing: rdum = utility_w0gauss((E(loop_f)-eig_k(i))/eta_smr, type_index)/eta_smr
else:       rdum = 1/(E(2)-E(1))              ! whole state dumped in one bin
```

`utility_w0gauss(x, n)` (src/utility.F90:1008-1091): n=0 Gaussian `exp(-x^2)/sqrt(pi)`
(arg clamped at 200); n>0 Methfessel-Paxton order n; n=-1 cold/M-V
`exp(-(x-1/sqrt(2))^2)*(2-sqrt(2)x)/sqrt(pi)`; n=-99 Fermi-Dirac `1/(2+e^x+e^-x)` (0 for
|x|>36). String→index map (src/readwrite.F90:1864-1914): 'gauss'→0, 'm-p'/'m-pN'→N,
'm-v'/'cold'→-1, 'f-d'→-99.

**Accumulation — this is the projection formula** (dos.F90:730-761). With
`r_num_elec_per_state = real(num_elec_per_state)`:

```
if (num_project == num_wann) then                             ! total DOS fast path
  dos_k(loop_f,1) += rdum * r_num_elec_per_state
  if (spin_decomp) dos_k(loop_f,2) += rdum*alpha_sq           ! NO num_elec_per_state factor
  if (spin_decomp) dos_k(loop_f,3) += rdum*beta_sq
else                                                          ! 0 < num_project < num_wann
  do j = 1, num_project
    dos_k(loop_f,1) += rdum * r_num_elec_per_state * abs(UU(project(j), i))**2
    if (spin_decomp) dos_k(loop_f,2) += rdum*alpha_sq*abs(UU(project(j), i))**2
    if (spin_decomp) dos_k(loop_f,3) += rdum*beta_sq *abs(UU(project(j), i))**2
  end do
end if
```

So the projection weight of band i onto the selected WFs is

```
w_i(k) = sum_{j=1..num_project} |UU(project(j), i)|^2       ! FIRST index = WF, SECOND = band
```

i.e. |<W_p|psi_i^(H)(k)>|² with `UU` the eigenvector matrix of the interpolated H(k)
(column = band). Since Σ_p over all num_wann WFs gives 1 (unitarity), the fast path is the
num_project==num_wann limit. Note this holds even if the user explicitly lists all WFs.

Spin decomposition (dos.F90:671-688 + spin.F90:223-305): `spin_get_nk` interpolates
`SS_n(k) = Σ_a alpha_a S_a(k)` from SS_R (a=x,y,z Pauli),
`alpha = (sin θ cos φ, sin θ sin φ, cos θ)` from `spin_axis_polar/azimuth` (degrees),
`spn_nk(i) = Re[(UU† SS_n UU)_{ii}]`; then

```
alpha_sq = (1 + spn_nk(i))/2      ! |up|^2
beta_sq  = 1 - alpha_sq           ! |down|^2
```

`spin_decomp=T` requires `num_elec_per_state==1` (postw90_readwrite.F90:690-693).

## 1.C Output file: `seedname-dos.dat`

Root only (`iprint>0`), dos.F90:339-350. stdout announces
`'Output data files:'` then `'   '//trim(seedname)//'-dos.dat'`. File opened
STATUS='UNKNOWN', FORM='FORMATTED'. **No header lines.** One line per energy point:

```
write (dos_unit, '(4E16.8)') omega, dos_all(ifreq, :)
```

- spin_decomp=F: 2 columns → `E [eV], dos_total [states/eV/cell]`
- spin_decomp=T: 4 columns → `E, dos_total, dos_spin_up, dos_spin_down`
- Each field E16.8 (e.g. `  0.80000000E+01`); the `4E16.8` descriptor simply truncates to the
  actual list length.

Normalisation: ∫dos dE = num_elec_per_state × num_wann per cell (full projection); DOS is per
eV per unit cell (k-weights sum to 1, no volume factor).

stdout block (dos.F90:206-231), formats verbatim:

```
'(/,/,1x,a)'  'Properties calculated in module  d o s'
'(1x,a)'      '--------------------------------------'
num_project==num_wann: '(/,3x,a)' '* Total density of states (_dos)'
else: '(/,3x,a)' '* Density of states projected onto selected WFs (_dos)'
      '(3x,a)'   'Selected WFs |Rn> are:'   + per WF: '(5x,a,2x,i3)' 'n =', project(i)
'(/,5x,a,f9.4,a,f9.4,a)' 'Energy range: [', e_min, ',', e_max, '] eV'
'(/,5x,a,(f6.3,1x))'     'Adaptive smearing width prefactor: ', adaptive_prefactor
'(/,/,1x,a20,3(i0,1x))'  'Interpolation grid: ', mesh(1:3)
'(/,1x,a)'  'Sampling the full BZ'  (or 'Sampling the irreducible BZ only')
```

## 1.D Keywords (src/postw90/postw90_readwrite.F90, defaults from postw90_types.F90)

| keyword | default | units | notes |
|---|---|---|---|
| `dos` | `.false.` (postw90_types.F90:55) | — | activates module (readwrite 416-417) |
| `dos_task` | `'dos_plot'` (postw90_types.F90:147) | — | only 'dos_plot' / 'find_fermi_energy' accepted (1241-1254); latter is dead code |
| `dos_energy_min` | `minval(eigval) − 0.6667` (or `win_min−0.6667` w/o eig) (1675-1682) | eV | |
| `dos_energy_max` | `froz_max+0.6667` if frozen_states, else `maxval(eigval)+0.6667`, else `win_max+0.6667` (1664-1673) | eV | |
| `dos_energy_step` | `0.01` (postw90_types.F90:151) | eV | grid step, recomputed to fit range exactly |
| `dos_project` | all WFs `1..num_wann` (1327-1339) | — | range vector (`1:5`, `1,3,5-7` syntax); out-of-range → error (1322-1326) |
| `dos_kmesh` / `dos_kmesh_spacing` | falls back to global `kmesh`/`kmesh_spacing` (get_module_kmesh, 1894-1897/1902+) | — / Å⁻¹ | one int → n×n×n; three ints; error if module active and none set |
| `dos_adpt_smr` | global `adpt_smr` = `.true.` (pw90_smearing_type, postw90_types.F90:133) | — | 1268-1270 |
| `dos_adpt_smr_fac` | global `adpt_smr_fac` = `sqrt(2)` (postw90_types.F90:134) | — | must be >0 (1273-1279) |
| `dos_adpt_smr_max` | global `adpt_smr_max` = `1.0` (postw90_types.F90:137) | eV | must be >0 (1282-1288) |
| `dos_smr_fixed_en_width` | global `smr_fixed_en_width` = `0.0` (postw90_types.F90:136) | eV | 0 → pure histogram (1290-1297) |
| `dos_smr_type` | global `smr_type` → index 0 Gaussian (postw90_types.F90:135) | — | 1341-1349 |
| `smr_type`,`adpt_smr`,`adpt_smr_fac`,`adpt_smr_max`,`smr_fixed_en_width` | as above | | global versions, read at 580-636 |
| `spin_decomp` | `.false.` | — | needs `.spn` + num_elec_per_state=1 (687-693) |
| `spin_axis_polar`, `spin_axis_azimuth` | `0.0`, `0.0` (postw90_types.F90:84-85) | degrees | |
| `num_elec_per_state` | 2; forced 1 if `spinors=T` (src/readwrite.F90:448-465, types.F90:68) | — | multiplies total DOS only |
| `wanint_kpoint_file` | `.false.` (postw90_types.F90:174) | — | irreducible-BZ mode |
| `scissors_shift`, `num_valence_bands` | 0.0 / unset | eV / — | applied inside get_HH_R, affects DOS |

## 1.E Traps

- The energy step is **recomputed** (`d_omega`), so the last point is exactly `energy_max`;
  `num_freq` forced ≥ 2.
- The histogram fallback triggers whenever `eta_smr < 2*binwidth` — including the adaptive
  path at k-points with tiny band velocity. It puts the whole state (weight `1/binwidth`)
  into the single nearest bin (nint rounding).
- Smearing evaluated only within ±10·eta_smr (`smearing_cutoff=10`), bin indices via nint —
  reproduce the nint+clamp arithmetic exactly for file-precision matching.
- Adaptive smearing needs band derivatives → `use_ws_distance` and the R→k conventions must
  match get_HH_R (postw90 default `use_ws_distance=.true.`, src/types.F90:87-94; the copper
  pdos test sets it **false**).
- Spin-decomp channels (columns 3,4) are **not** multiplied by num_elec_per_state (comment
  dos.F90:735-737); with spinors num_elec_per_state=1 anyway.
- `dos_project` weight uses `UU(project(j), i)` — WF index is the **row**, band the column;
  no conjugation subtlety since only |·|² enters.
- No eigenvalue cutoff and no Fermi level enters dos_plot at all; `fermi_energy` is
  irrelevant to the .dat contents.
- CODATA2006 vs 2010/2018 constants affect nothing here except through recip_lattice
  (Ang) — the .eig energies are used as-is in eV.

## 1.F Reference test (`test-suite/tests/testpostw90_example04_pdos/`)

Files: `copper.win`, `copper.eig`, `copper.chk.fmt.bz2` (symlink →
`../../checkpoints/cu_postw90/`), `Makefile` (bunzip2 + `w90chk2chk.x -f2u copper`),
benchmark = `benchmark.out.default.inp=copper.win` which is the **copper-dos.dat content**.
Confirmed: no .mmn/.amn/.spn shipped or needed.

`copper.win` keywords (quoted): `num_bands = 12`, `num_wann = 7`,
**`use_ws_distance = .false.`**, `search_shells=12`, `dis_win_max = 38.0`,
`dis_froz_max = 13.0`, **`dos = true`**, **`kmesh = 10`** (global → dos mesh 10×10×10),
**`dos_energy_max = 10`**, **`dos_energy_min = 8`**, **`dos_energy_step = 0.25`**,
**`dos_project 1:5`** (the five Cu:d WFs; projections block `Cu:d` + two `f=…:s`),
`mp_grid : 4 4 4`. NOT set: `dos_adpt_smr` (→ default T, fac √2, max 1 eV, Gaussian),
`spin_decomp` (F), `spinors` (F → num_elec_per_state=2).

num_freq = nint((10−8)/0.25)+1 = 9 lines; benchmark file verbatim:

```
  0.80000000E+01  0.16146066E+01
  0.82500000E+01  0.22074413E+01
  0.85000000E+01  0.23050405E+01
  0.87500000E+01  0.19915528E+01
  0.90000000E+01  0.27906176E+01
  0.92500000E+01  0.44862488E+01
  0.95000000E+01  0.27430925E+01
  0.97500000E+01  0.22709949E+01
  0.10000000E+02  0.28830387E+01
```

Harness: `tests/jobconfig:431-434` `[testpostw90_example04_pdos/]`,
`program = POSTW90_DOS_OK`, `output = copper-dos.dat`. Tolerances
(`tests/userconfig:139-145`, parser `tools/parsers/parse_dos_dat.py`: 2-col → keys
energy/dos, 4-col adds dos_spin1/dos_spin2): `energy (1.0e-6, 5.0e-6)`,
`dos (1.0e-4, 1.0e-4)`, `dos_spin1/2 (1.0e-4, 1.0e-4)`.

---

# Part 2 — transl_inv and transl_inv_full (src/postw90/get_oper.F90)

## 2.A Keywords

| keyword | default | where read | notes |
|---|---|---|---|
| `transl_inv` | `.false.` (postw90_types.F90:175) | postw90_readwrite.F90:851-852 → `pw90_berry%transl_inv` | MV97 Eq.(31) band-diagonal Im-log |
| `transl_inv_full` | `.false.` (postw90_types.F90:176) | postw90_readwrite.F90:854-855 | full phase-corrected scheme (all elements) |
| both T | **error** `'Error: If transl_inv_full=T, transl_inv=T is not recommended'` (857-860) | | mutually exclusive |
| `guiding_centres` | `.false.` (postw90_types.F90:177) | 862-863 | relaxes the Wannier-centre consistency check in get_AA_R |
| `higher_order_n` | `1` (src/types.F90:141) | src/readwrite.F90:911-919 (kmesh_input) | higher-order finite-difference b-shells (kmesh.F90); used by the _higher test |

(There is no `berry_transl_inv` keyword — the input token is plain `transl_inv`; it lives in
`pw90_berry_mod_type`.)

## 2.B `transl_inv = T`: exactly what changes

Only inside `get_AA_R` (get_oper.F90:403-792), i.e. only the position operator A(R). Baseline
(both flags F) for every k-point ik and neighbour b (shell nn), after projecting the raw .mmn
block to the Wannier gauge `S = V(k)† S_o V(k+b)` (get_gauge_overlap_matrix, 3235-3272):

```
AA_q_b(:,:,ik,nn,idir) += cmplx_i * wb(nn) * bk(idir,nn,ik) * S(:,:)        ! lines 611-613
```

With `transl_inv=T` the **band-diagonal elements only** are rewritten a la MV97 Eq.(31)
(code verbatim, lines 606, 614-632; `nno==nn` here since transl_inv_full=F):

```
if (pw90_berry%transl_inv .and. ik .ne. ik_prev) AA_q_b_diag(:, :, :) = cmplx_0   ! line 606
...
if (pw90_berry%transl_inv) then
  ! Rewrite band-diagonal elements a la Eq.(31) of MV97
  do i = 1, num_wann
    AA_q_b_diag(i, nno, idir) = AA_q_b_diag(i, nno, idir) &
                                - kmesh_info%wb(nn)*kmesh_info%bk(idir, nn, ik) &
                                *aimag(log(S(i, i)))                              ! lines 618-622
  end do
end if
...
if (pw90_berry%transl_inv) then
  do n = 1, num_wann
    AA_q_b(n, n, ik, nno, idir) = AA_q_b_diag(n, nno, idir)                       ! lines 627-631
  end do
end if
```

Net effect per k after summing shells (line 734, `AA_q = sum(AA_q_b, 4)`):

```
off-diagonal (n≠m):  A_α,nm(k) =  i Σ_b w_b b_α S_nm(k,b)             (unchanged)
diagonal:            A_α,nn(k) = −Σ_b w_b b_α Im ln S_nn(k,b)         (replaced, real)
```

No factor i on the diagonal; the minus sign comes with Im ln; `aimag(log(z))` is the principal
branch arg(z) ∈ (−π, π]. Everything downstream is identical to the default path:
hermitization `0.5*(A + A†)` per (ik, idir) (lines 743-749), Fourier q→R with
`e^{−i2πq·R}/N_q` (fourier_loc_q_to_R), `operator_wigner_setup` degeneracy division
(see berry-ahc.md §A). Note the diagonal is already real, so hermitization preserves it.

Independent of the flag, get_AA_R always accumulates the Im-log centres
`wannier_centres_from_AA_R(:,i) −= wb·bk(:)·Im ln S_ii /num_kpts` (lines 598-602) and fails
with `'Computed and read Wannier centres different.'` if
`sum((computed − chk_centres)**2) > 1e-8`, unless `guiding_centres=T` (lines 639-647).

### Which modules honour / refuse `transl_inv`

- Honoured wherever `get_AA_R` supplies AA_R: berry module tasks ahc/kubo/sc/shc
  (berry.F90:301-424), kpath (kpath.F90:218/243), kslice (kslice.F90:217/243), gyrotropic
  (gyrotropic.F90:234). dos/boltzwann never touch AA_R.
- **Refused for morb**: `if (eval_morb) → error 'transl_inv=T disabled for morb'`
  (berry.F90:556-560); otherwise berry prints
  `'Using a translationally-invariant discretization for the'` /
  `'band-diagonal Wannier matrix elements of r, etc.'` (561-562).
- **Refused for the gyrotropic K-tensor**: `'transl_inv=T disabled for K-tensor'`
  (gyrotropic.F90:333-337), same message printed otherwise (338-342).
- Comparison: wannier90.x's `hamiltonian_write_tb` (src/hamiltonian.F90:962-979) ALWAYS uses
  the Im-log diagonal — `seedname_tb.dat` corresponds to postw90 with transl_inv=T (and no
  hermitization).

## 2.C `transl_inv_full = T`: the phase-corrected scheme (all elements)

Affects `get_AA_R`, `get_BB_R`, `get_CC_R`, and the SHC-Ryoo operators `get_SBB_R`
(get_oper.F90:2523-2819) / `get_SAA_R` (2822-3117). Explicitly **not implemented for
shc_method=qiao**: `'Error: transl_inv_full=T not implemented for shc_method=qiao'`
(postw90_readwrite.F90:1117-1124).

Common ingredients:
- `r0(i,j,:) = (rbar_i + rbar_j)/2` with `rbar = wigner_seitz%wannier_centres_from_AA_R`
  (Im-log centres accumulated in get_AA_R; e.g. get_oper.F90:650-657, 890-898, 1195-1203).
  Hence get_AA_R must run first (berry calls it before BB/CC).
- b-vectors are reordered to a k-independent ordering:
  `nno = kmesh_info%nninv(nn, ik)` with `bk(:, nninv(nn,ik), 1) = bk(:, nn, ik)`
  (src/types.F90:198, kmesh.F90:911-922) so that all k share the b-list of ik=1 (line 609).
- q-space phase per (i,j) and shell: `phase1 = exp(+i b·r0_ij)` (AA: lines 659-673 applied
  after the loop; BB: 978-985 applied inside; CC: `phase1 = exp(+i (b2−b1)·r0_ij)`,
  1284-1295).
- The b-sum is deferred: each shell is Fourier-transformed q→R separately, then multiplied by
  the R-dependent phase `phase2 = exp(−i R_cart·b/2)` (AA 701-706, BB 1027-1032; CC uses
  `exp(−i R·(b1+b2)/2)`, 1353-1362) and accumulated into the R-space operator.
- Final real-space corrections on root:
  - AA_R: overwrite the R=0 diagonal with the Im-log centres:
    `AA_R(i,i,ir0,:) = wannier_centres_from_AA_R(:,i)` (lines 719-728). No hermitization step
    in this branch (contrast lines 743-749 of the default branch).
  - BB_R: `BB_R(:,:,ir,idir) += (r0(:,:,idir) − 0.5*R_idir) * HH_R(:,:,ir)` (1045-1050);
    requires HH_R (error `'transl_inv_full=T for CC_R needs HH_R'` variant checks at
    1317-1324 for CC).
  - CC_R (1377-1398): for all a,b and every R (ir0(−R) located via utility_compar):
    `CC_R(:,:,R,a,b) += (r0_a + R_a/2)·BB_R(:,:,R,b)
                      + BB_R(:,:,−R,a)†·(r0_b − R_b/2)
                      + (r0_a + R_a/2)·R_b·HH_R(:,:,R)`.
    In this branch only a≤b slots of the uHu part are filled (1297-1305) and downstream only
    a≤b is used (berry.F90:2070-2079 reconstructs `CC(j,i) = CC(i,j)†` in k-space); in the
    default branch CC_q is hermitized in q-space instead (1407-1413).
- SBB_R requires SH_R (`'transl_inv_full=T for SBB_R needs SH_R'`, 2710-2712); SAA_R requires
  SS_R (`'transl_inv_full=T for SAA_R needs SS_R'`, 3009-3011).

Consequence: with transl_inv_full the operators (hence morb, ahc, shc-ryoo results) become
independent of a rigid shift of all Wannier centres; results match TB-convention
interpolation around the actual centres.

## 2.D The two Fe morb tests

Both use `program = POSTW90_WPOUT_OK` (tests/jobconfig:298-307), `output = Fe.wpout`;
tolerances (tests/userconfig:110-122, parser `tools/parsers/parse_wpout.py`, regex on the
`======================` row printed by berry.F90:1488-1491, format
`'(1x,a22,2x,3(f10.4,1x),/)'`): `morb_x/y/z (1.0e-3, 2.0e-3)`.

Runtime inputs: Fe.win, Fe.chk (from .chk.fmt.bz2 via w90chk2chk), Fe.eig, Fe.mmn (get_AA_R +
get_BB_R), Fe.uHu (get_CC_R). Fe.amn present but unused.

### testpostw90_fe_morb_transl_inv/Fe.win (quoted)

`num_bands = 28`, `num_wann = 18`, `use_ws_distance = .true.`, `search_shells=12`,
`spinors = true`, `fermi_energy = 12.6279`, **`uHu_formatted = .true.`** (pw90_oper_read
default F, postw90_types.F90:70), **`berry = true`**, **`berry_task = morb`**,
**`berry_kmesh = 10 10 10`**, **`transl_inv_full = .true.`**, cell 2.71175 bohr bcc,
`mp_grid = 2 2 2` (8 explicit kpoints). Note: `transl_inv` itself is NOT set — it stays F
(it would be an input error with morb, berry.F90:556-560). The b-shell list: 1 shell,
12 neighbours, wb = 0.104321 Ang².

Benchmark (`benchmark.out.default.inp=Fe.win:224-228`):

```
 Fermi energy (ev) =   12.627900
 M_orb (bohr magn/cell)        x          y          z
 ======================      0.0000    -0.0000     0.0415
```

Checked values: morb_x=0.0000, morb_y=-0.0000, morb_z=0.0415.

### testpostw90_fe_morb_transl_inv_higher/Fe.win (quoted)

`num_wann = 18`, `num_bands = 28`, `num_iter = 0`, `dis_froz_max = 30.0000`,
**`guiding_centres = T`** (relaxes the centre-consistency check — needed because the
higher-order Im-log centres differ from the chk centres), `use_ws_distance = T`,
`spinors = .true.`, `mp_grid = 2 2 2`, **`kmesh = 10 10 10`** (global keyword this time),
**`berry = T`**, **`berry_task = morb`**, **`fermi_energy = 12.6631`**, **`transl_inv = F`**,
**`transl_inv_full=T`**, **`higher_order_n=2`** (doubles the b-shell set with higher-order FD
weights, kmesh.F90:453-620; `nntot` becomes 2× the first-order count, weights from the
order-2 stencil). uHu_formatted not set → binary .uHu (the shipped Fe.uHu.bz2 is unformatted).

Benchmark (`benchmark.out.default.inp=Fe.win:268-272`):

```
 Fermi energy (ev) =   12.663100
 M_orb (bohr magn/cell)        x          y          z
 ======================      0.0000    -0.0000    -0.0617
```

Checked values: morb_x=0.0000, morb_y=-0.0000, morb_z=-0.0617.

## 2.E Traps

- `transl_inv` is a **diagonal-only** modification of A(k) before the q→R transform; do not
  touch off-diagonal elements, do not multiply by i, and use `Im ln` (principal branch) — a
  wrong branch cut flips centres by lattice vectors.
- The two flags are exclusive; the pair (morb, transl_inv=T) and (gyrotropic K-tensor,
  transl_inv=T) and (shc qiao, transl_inv_full=T) are hard errors.
- transl_inv_full changes the order of operations: per-shell FT then R-dependent phase — you
  cannot sum over b in q-space first. It also requires the `nninv` reordering so that the
  phase2 factors (built from `bk(:,nn,1)`) match the shell.
- The wannier-centre consistency check (1e-8 on the summed squared difference, Å²) is active
  in every get_AA_R call; with guiding_centres=T the computed centres are REPLACED by the chk
  centres (line 643) — this feeds into r0 for BB/CC in the _higher test.
- morb unit factor (berry.F90:1447): `fac = −eV_au/bohr**2` with CODATA2006 values
  (`eV_au = 27.21138386`-family, `bohr = 0.52917720859`; postw90 header prints
  'Using CODATA 2006 constant values'), then M = (g + h − 2·E_F·f)·fac summed as
  LC + IC (lines 1459-1467); printed with `3(f10.4,1x)`.
- Both tests interpolate on a 10×10×10 berry mesh from a 2×2×2 ab-initio grid — results are
  extremely sensitive to the WS R-vector set: the first test sets use_ws_distance **true**
  (postw90 default), unlike the fe_ahc test.
