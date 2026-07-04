# postw90 reference notes: spin module (spin_moment, spin_get_nk, dos spin_decomp)

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Modules covered:
- `src/postw90/spin.F90` (`w90_spin`: `spin_get_moment`, `spin_get_nk`, `spin_get_S`,
  private `spin_get_moment_k`)
- `src/postw90/get_oper.F90` (`get_SS_R`, lines 1652-1836)
- `src/postw90/dos.F90` spin_decomp branch (`dos_main` 54-368, `dos_get_k` 600-765)
- keyword parsing `src/postw90/postw90_readwrite.F90`, defaults `src/postw90/postw90_types.F90`

Convention used throughout: the "spin" operator is the **Pauli matrix vector sigma**
(eigenvalues in [-1,1]), NOT S = sigma/2. This is why no g=2 factor appears anywhere
(explicit comment spin.F90:198-199).

---

## A. SS(R): Wannier-gauge Pauli matrices from .spn (get_SS_R, get_oper.F90:1652-1836)

### A.1 .spn file reading (root node only)

The `.spn` file stores `<psi_nk | sigma_i | psi_mk>` (i = x,y,z) between the ORIGINAL
ab-initio eigenstates over the full `num_bands` set (upper triangle only, hermiticity
implied). Selected by `spn_formatted` (default `.false.`, postw90_types.F90:68; parsed at
postw90_readwrite.F90:477-478).

**Formatted layout** (get_oper.F90:1730-1737, 1755-1771):
```
line 1: header string (echoed to stdout; NOTE list-directed read into character(len=60)
        only captures the FIRST whitespace-delimited token, e.g. "Created")
line 2: num_bands  num_kpts        ! both validated against the run, fatal on mismatch
then for ik = 1..num_kpts:          ! outer loops: do m = 1, num_bands; do n = 1, m
  3 lines per (n,m) pair, each "re im" (list-directed):
    line A -> spn_o(n,m,ik,1)  (sigma_x)
    line B -> spn_o(n,m,ik,2)  (sigma_y)
    line C -> spn_o(n,m,ik,3)  (sigma_z)
  and the lower triangle is filled as spn_o(m,n,ik,i) = conjg(spn_o(n,m,ik,i))
```
Pair order = packed upper triangle, column-major: (1,1),(1,2),(2,2),(1,3),(2,3),(3,3),...
The stored element with n<=m is the (row=n, col=m) entry.

**Unformatted layout** (lines 1739-1745, 1772-1798): record 1 = header (60 chars), record 2 =
`nb, nkp`, then ONE record per ik containing `((spn_temp(s,m),s=1,3), m=1,(nb*(nb+1))/2)` —
i.e. 3*nb(nb+1)/2 complex(dp), sigma index fastest, same packed-triangle order.

Errors: `'...spn has wrong number of bands'` / `'...k-points'` (1747-1754); open/read failures
map to `'Error: Problem opening/reading input file <seed>.spn'` (1831-1834).

Guard: if `SS_R` is already allocated the routine returns immediately ("been here before",
lines 1707-1711) — relevant when spin_moment + spin_decomp/kpath run in one postw90 invocation.

### A.2 Windowing and gauge rotation (lines 1802-1812)

`num_states(ik) = dis_manifold%ndimwin(ik)` if `have_disentangled` else `num_wann`
(lines 1719-1725). For each ik and each sigma component is=1..3,
`get_gauge_overlap_matrix` (get_oper.F90:3235-3272) computes

```
SS_q(:,:,ik,is) = V(ik)^dagger . spn_o(wm:wm+ns-1, wm:wm+ns-1, ik, is) . V(ik)
```

- `wm` = index of the FIRST band with `lwindow(j,ik)=.true.` (`get_win_min`,
  get_oper.F90:3199-3232; wm=1 when not disentangled), `ns = num_states(ik)`.
  **TRAP: the window is taken as the contiguous block wm..wm+ns-1, not the lwindow mask.**
- `v_matrix` as in the AHC note: `v = u_matrix` (no disentanglement) or
  `v(m,j,k) = sum_i u_matrix_opt(m,i,k) u_matrix(i,j,k)` (postw90_common.F90:811-938).
- utility_zgemmm call with 'C','N','N' → prod = op(a).op(b).op(c) = V† S_o V, a num_wann×num_wann
  complex matrix per component.

### A.3 Fourier q→R and degeneracy weights (lines 1814-1823)

Identical machinery to H(R) (see berry-ahc.md §A.1):

```
SS_R_temp(:,:,ir,is) = (1/N_kpts) sum_q e^{-i 2pi q·R_ir} SS_q(:,:,q,is)   ! fourier_q_to_R, 3124-3157
```
then `operator_wigner_setup` (3275-3327) divides by `ndegen(ir)` (and, if
`use_ws_distance=.true.`, additionally by `ws_distance%ndeg(i,j,ir)` while rescattering onto
the pw90 R-grid `ir_ind_ws_to_pw90`). Result `SS_R(num_wann, num_wann, nrpts_pw90, 3)`,
broadcast to all nodes (line 1823). Index 4 = sigma component (1=x, 2=y, 3=z).

k-interpolation everywhere below uses `pw90common_fourier_R_to_k` with alpha=0
(postw90_common.F90:1032-1093): `S_i(k) = sum_R e^{+i 2pi k·R} SS_R(:,:,ir,i)` — no further
degeneracy weights.

---

## B. spin_get_nk — per-band spin projection (spin.F90:223-305)

Computes `spn_nk(m) = <psi_mk^(H) | sigma·n | psi_mk^(H)>`, m = 1..num_wann:

1. `H(k) = sum_R e^{i2pi k·R} HH_R`; `utility_diagonalize` → eig, UU with
   H(k) = UU·diag(eig)·UU† (lines 279-283).
2. `S_i(k) = sum_R e^{i2pi k·R} SS_R(:,:,:,i)`, i=1..3 (lines 285-290).
3. Quantization axis from the two input angles, **degrees** (lines 294-297):
```
conv = 180/pi
n_x = sin(axis_polar/conv) * cos(axis_azimuth/conv)
n_y = sin(axis_polar/conv) * sin(axis_azimuth/conv)
n_z = cos(axis_polar/conv)
```
4. `SS_n = n_x S_x + n_y S_y + n_z S_z` (line 301), then (line 303)
```
spn_nk(m) = real( (UU^dagger . SS_n . UU)(m,m) )      ! utility_rotate_diag, utility.F90:782-802
```

No clamping inside spin_get_nk. Clamping to ±(1−eps8), eps8 = 1.0e-8 (constants.F90:76),
is applied ONLY by the callers kpath (kpath.F90:348-354) and kslice (kslice.F90:367-373) —
NOT by the DOS (see §D).

Companion `spin_get_S` (spin.F90:387-452): same, but returns all three Cartesian components
`S(m,i) = real((UU† S_i UU)(m,m))` (used by berry/gyrotropic spin tasks).

Callers for colouring:
- kpath (`kpath_bands_colour = 'spin'`, default `'none'`, postw90_types.F90:113):
  `seedname-bands.dat` rows `(3E16.8)` = xval, eig(i,kpt), colour(i,kpt), band-major with a
  blank line between bands (kpath.F90:505-519); gnuplot palette `(-1 "blue", 0 "green", 1 "red")`,
  zrange [-1:1] (kpath.F90:550-555).
- kslice (`kslice_fermi_lines_colour = 'spin'`): `seedname-kslice-fermi-spn.dat`,
  `(3E16.8)` = kpt_x, kpt_y, spn masked to bands with |eig−E_F| < Delta_E (kslice.F90:397-405,
  546-551).

---

## C. spin_get_moment — total spin magnetic moment (spin.F90:54-220)

### C.1 Formula

Requires nfermi ≤ 1: error `'Routine spin_get_moment requires nfermi=1'` if
`size(fermi_energy_list) > 1` (lines 120-125). Uses `fermi_energy_list(1)`.
**TRAP: `fermi_energy_list` is ALWAYS allocated with n=1, value 0.0 eV, when no fermi_energy
keyword is given** (readwrite.F90:609-692) — an unset Fermi level silently means E_F = 0.

Calls `get_HH_R` + `get_SS_R` (lines 127-136). Per k-point (`spin_get_moment_k`,
spin.F90:311-384):

```
diagonalize H(k) -> eig, UU
occ(n) = 1 if eig(n) < ef else 0                  ! strict <, T=0 step; pw90common_get_occ,
                                                  ! postw90_common.F90:942-985
spn_nk(n,is) = aimag( i * (UU† S_is(k) UU)(n,n) ) ! == real part; lines 372-378
spn_k(is)    = sum_n occ(n) * spn_nk(n,is)        ! lines 379-381
```

BZ sum (full-BZ branch, lines 170-191; header `'Sampling the full BZ (not using symmetry)'`):

```
kweight = 1 / product(spin_kmesh)
loop_tot = my_node_id, product(mesh)-1, num_nodes      ! flat MPI-strided loop
loop_x = loop_tot/(m2*m3); loop_y = (loop_tot - loop_x*m2*m3)/m3; loop_z = remainder
kpt = (loop_x/m1, loop_y/m2, loop_z/m3)                ! unshifted Gamma-centred grid
spn_all += spn_k * kweight
```
(With `wanint_kpoint_file=T` the irreducible-wedge kpoint.dat weights are used instead,
lines 146-168, with an IBZ warning.) Then `comms_reduce(spn_all,3,'SUM')` and

```
spn_mom(1:3) = − spn_all(1:3)          ! line 201; sign flip: m = −<sigma>·mu_B
```

No factor of 2 (sigma not S), no division by degeneracy; result is in **Bohr magnetons per
cell** directly. Angles (defined "as in pwscf", lines 210-217):

```
magnitude = sqrt(mx² + my² + mz²)
theta = acos(mz/magnitude) * 180/pi
phi   = atan(my/mx)        * 180/pi     ! plain atan, NOT atan2 -> phi in (−90°,90°);
                                        ! NaN/garbage if mx == 0 exactly
```

### C.2 Exact .wpout output lines (root, iprint>0; lines 138-218)

```
write (stdout, '(/,/,1x,a)') '------------'
write (stdout, '(1x,a)')     'Calculating:'
write (stdout, '(1x,a)')     '------------'
write (stdout, '(/,3x,a)')   '* Spin magnetic moment'
write (stdout, '(/,1x,a)')   'Sampling the full BZ (not using symmetry)'
write (stdout, '(/,1x,a)')   'Spin magnetic moment (Bohr magn./cell)'
write (stdout, '(1x,a,/)')   '===================='
write (stdout, '(1x,a18,f11.6)')   'x component:', spn_mom(1)
write (stdout, '(1x,a18,f11.6)')   'y component:', spn_mom(2)
write (stdout, '(1x,a18,f11.6)')   'z component:', spn_mom(3)
write (stdout, '(/,1x,a18,f11.6)') 'Polar theta (deg):', theta
write (stdout, '(1x,a18,f11.6)')   'Azim. phi (deg):', phi
```
`a18` right-justifies the label in an 18-char field. Rendered example (benchmark):

```
 Spin magnetic moment (Bohr magn./cell)
 ====================

       x component:  -0.000003
       y component:  -0.000000
       z component:   3.090787

 Polar theta (deg):   0.000059
   Azim. phi (deg):   5.831720
```

There is NO data file for spin_moment; everything goes to stdout/.wpout.

---

## D. dos spin_decomp (src/postw90/dos.F90)

### D.1 Setup (dos_main, 54-368)

Energy grid (lines 146-158):
```
num_freq = nint((dos_energy_max − dos_energy_min)/dos_energy_step) + 1;  if 1 -> 2
d_omega  = (max − min)/(num_freq − 1)                 ! actual step, may differ from input step
E(ifreq) = dos_energy_min + (ifreq−1)*d_omega
```

`spin_decomp=T` → `ndim = 3` and `get_SS_R` is called (lines 183-193); else ndim=1.
k-sum identical in structure to spin_get_moment but over `dos_kmesh`
(kweight = 1/product(mesh), flat strided loop, lines 283-332; `'Sampling the full BZ'`).
Adaptive smearing branch computes band gradients via `wham_get_eig_deleig` and
`dos_get_levelspacing` (levelspacing(n) = |∇E_n| * Delta_k, lines 768-795); fixed branch just
diagonalizes H(k).

### D.2 Per-k contribution (dos_get_k, 600-765)

Spin weights per band i (lines 671-688), from `spn_nk = spin_get_nk(...)` (i.e. sigma·n̂ along
the SAME spin_axis_polar/azimuth axis):

```
alpha_sq = (1 + spn_nk(i)) / 2       ! |alpha|^2, weight of the spin-UP channel
beta_sq  = 1 − alpha_sq              ! |beta|^2 = spin-DOWN weight
```

**NO clamping of |spn_nk| > 1 in the DOS path** — if interpolation slightly overshoots,
alpha_sq can exceed 1 and beta_sq go negative by the same margin (only kpath/kslice clamp).
Up + down always sums exactly to the (per-state) total.

Smearing width per band (lines 690-696):
```
fixed:    eta = dos_smr_fixed_en_width
adaptive: eta = min(levelspacing(i)*dos_adpt_smr_fac, dos_adpt_smr_max)   ! Eq.(35) YWVS07
```
Bin range optimisation (lines 699-715): if `eta/binwidth < min_smearing_binwidth_ratio` (= 2.0,
constants.F90:82) NO smearing — the whole state is dumped in the single nearest bin
`nint((eig−E1)/(E_last−E1)*(nbins−1))+1` (clamped to [1,nbins]) with weight `1/binwidth`;
otherwise bins within `± smearing_cutoff*eta` (= 10*eta, constants.F90:80) get
`rdum = utility_w0gauss((E_bin−eig)/eta, type_index)/eta` (w0gauss: utility.F90:1008-1091;
type_index 0 = Gaussian exp(−x²)/sqrt(pi), −1 = Marzari-Vanderbilt cold, −99 = Fermi-Dirac,
n>0 = Methfessel-Paxton order n).

Accumulation, full-projection case `num_project == num_wann` (lines 730-743):

```
dos_k(loop_f, 1) += rdum * num_elec_per_state        ! total DOS
if (spin_decomp):
  dos_k(loop_f, 2) += rdum * alpha_sq                ! spin-up   (NO num_elec_per_state factor)
  dos_k(loop_f, 3) += rdum * beta_sq                 ! spin-down (NO num_elec_per_state factor)
```

Comment at 735-737: the up/down columns deliberately omit num_elec_per_state because
spin_decomp implies spinor calculation ⇒ num_elec_per_state = 1 (enforced at parse time, see
§E). Projected case (`dos_project` set, lines 749-760): every term additionally multiplied by
`abs(UU(project(j), i))**2` summed over the selected WFs j.

### D.3 Output file (dos_main, lines 339-350)

Filename `seedname-dos.dat` (stdout announces `'Output data files:'` then
`'   <seed>-dos.dat'`). One row per energy bin, no header line:

```
write (dos_unit, '(4E16.8)') omega, dos_all(ifreq, :)
```
Column order: `E [eV] | total DOS | spin-up DOS | spin-down DOS` (columns 3-4 only when
spin_decomp; with ndim=1 the same format writes just 2 columns). Units: states/eV/cell
(the k-weights 1/N sum the normalized BZ average).

---

## E. Input keywords (parse: postw90_readwrite.F90; defaults: postw90_types.F90)

| keyword | default | where | notes |
|---|---|---|---|
| `spin_moment` | `.false.` | postw90_types.F90:60; parsed 676 | activates spin_get_moment |
| `spin_decomp` | `.false.` | postw90_types.F90:61; parsed 687 | used by dos (also kubo/boltzwann); **error `'spin_decomp can be true only if num_elec_per_state is 1'`** (690-693) |
| `spin_axis_polar` | `0.0` | postw90_types.F90:84; parsed 679-680 | degrees |
| `spin_axis_azimuth` | `0.0` | postw90_types.F90:85; parsed 683-684 | degrees |
| `spin_kmesh` / `spin_kmesh_spacing` | falls back to global `kmesh`/`kmesh_spacing` | get_module_kmesh, prefix `'spin'`, postw90_readwrite.F90:1889-1892 (routine 1902-1998) | 1 or 3 ints; mutually exclusive with _spacing; error if spin_moment=T and neither local nor global mesh given |
| `spn_formatted` | `.false.` (binary) | postw90_types.F90:68; parsed 477-478 | |
| `fermi_energy` | list = [0.0] when absent | readwrite.F90:609-692 | spin_get_moment errors only for nfermi>1 |
| `num_elec_per_state` | 2; forced 1 by `spinors=T` | types.F90:68; readwrite.F90:448-465 | only 1 or 2 allowed |
| `dos` | `.false.` | pw90_calculation | activates dos_main (dos_task default `'dos_plot'`) |
| `dos_kmesh` / `_spacing` | global kmesh | postw90_readwrite.F90:1894-1897 | |
| `dos_energy_step` | `0.01` eV | postw90_types.F90:151; parsed 1263 | grid is re-derived, see §D.1 |
| `dos_energy_min` | `minval(eigval) − 0.6667` | postw90_readwrite.F90:1675-1682 | eV |
| `dos_energy_max` | `froz_max + 0.6667` if frozen window else `maxval(eigval) + 0.6667` | postw90_readwrite.F90:1664-1673 | eV |
| `dos_adpt_smr` | global `adpt_smr` (default `.true.`) | postw90_types.F90:133; parsed 1267-1270 | |
| `dos_adpt_smr_fac` | global (default `sqrt(2)` ≈ 1.414) | postw90_types.F90:134; parsed 1272-1279 | |
| `dos_adpt_smr_max` | global (default `1.0` eV) | postw90_types.F90:137; parsed 1281-1288 | |
| `dos_smr_fixed_en_width` | global (default `0.0` eV) | postw90_types.F90:136; parsed 1290-1297 | with the ratio-2 rule, width 0 ⇒ histogram bins |
| `dos_smr_type` | global `smr_type` (Gaussian, type_index 0) | parsed 1341-1349 | |
| `dos_project` | all WFs 1..num_wann | parsed 1305-1339 | range vector |
| `kpath_bands_colour` | `'none'` | postw90_types.F90:113 | `'spin'` colours bands by spn_nk |
| `kslice_fermi_lines_colour` | `'none'` | postw90_types.F90:126 | `'spin'` |
| `use_ws_distance` | `.true.` | types.F90:87-94 | both Fe tests set `.false.` |

Parameter-block echo in .wpout (postw90_readwrite.F90:2123-2131):
```
'|  Spin decomposition                        :', L8   ! '(1x,a46,10x,L8,13x,a1)'
'|  Compute Spin moment                       :', L8
'|  Polar angle of spin quantisation axis     :', f8.3 ! '(1x,a46,10x,f8.3,13x,a1)'
'|  Azimuthal angle of spin quantisation axis :', f8.3
'          |  Spn file-type                   :', 'formatted'/'unformatted'
```

---

## F. Reference tests

### F.1 testpostw90_fe_spin (spin moment)

`tests/testpostw90_fe_spin/Fe.win` — relevant keywords:
`num_bands = 28`, `num_wann = 18`, **`use_ws_distance = .false.`** (l.3), `spinors = true`
(l.15), **`fermi_energy = 12.6279`** (l.20), **`spin_moment=true`** (l.22),
**`spn_formatted=true`** (l.23), **`kmesh = 4`** (l.24 — GLOBAL kmesh, inherited as
spin_kmesh 4×4×4), bcc cell 2.71175 bohr, `mp_grid = 2 2 2`, 8 explicit kpoints.
spin_axis_polar/azimuth not set → 0/0 (ẑ axis). Files present: Fe.win, Fe.spn (formatted),
Fe.eig, Fe.chk.fmt.bz2 (→ w90chk2chk.x -f2u in Makefile), Fe.mmn.bz2 (not needed by spin),
Fe.amn (unused).

Benchmark (`benchmark.out.default.inp=Fe.win`), quantities checked:

```
       x component:  -0.000003
       y component:  -0.000000
       z component:   3.090787

 Polar theta (deg):   0.000059
   Azim. phi (deg):   5.831720
```

Harness: `tests/jobconfig:365-368` — program `POSTW90_WPOUT_OK`, output `Fe.wpout`.
`tests/userconfig:110-123` parses with `tools/parsers/parse_wpout.py` regexes
(`x\ component:\s*(...)`, `Polar\ theta\ \(deg\):`, `Azim.\ phi\ \(deg\):`); tolerances
(rel, abs): spin_x/y/z (1.0e-3, 2.0e-3), spin_p/spin_a (1.0e-2, 2.0e-2). Note x,y and the
angles are numerical noise around the z axis (hence loose angle tolerances; phi = atan of a
0/0-ish ratio).

Benchmark also fixes the incidental stdout: `'Reading spin matrices from Fe.spn in
get_SS_R : Created'` (header first token) and CODATA2006 constants banner
(bohr = 0.52917720859 Å converts the 2.71175-bohr cell to 1.434996 Å).

### F.2 testpostw90_fe_dos_spin (spin-decomposed DOS)

`tests/testpostw90_fe_dos_spin/Fe.win` = the fe_spin input (including `spin_moment=true`,
`kmesh = 4`) plus:

```
dos = true
spin_decomp = true
dos_energy_max = 13.0
dos_energy_min = 10.0
dos_energy_step = 0.2
dos_adpt_smr = false
dos_smr_fixed_en_width = 0.5
```

→ 16 energy bins (nint(3.0/0.2)+1), fixed Gaussian smearing eta = 0.5 eV, binwidth 0.2,
eta/binwidth = 2.5 ≥ 2 ⇒ smeared, bins within ±5 eV of each state. num_elec_per_state = 1
(spinors). Axis = ẑ (defaults).

Benchmark (`benchmark.out.default.inp=Fe.win`) IS the `Fe-dos.dat` file, 16 rows × 4 cols in
`(4E16.8)`; first and last rows:

```
  0.10000000E+02  0.83154571E+00  0.11547299E+00  0.71607272E+00
  ...
  0.13000000E+02  0.16904725E+01  0.11508848E+01  0.53958769E+00
```
(columns: E, total, up = (1+⟨sigma_z⟩)/2-weighted, down; up+down = total since
num_elec_per_state=1). Harness: `tests/jobconfig:371-374` — program `POSTW90_DOS_OK`, output
`Fe-dos.dat`; tolerances `tests/userconfig:139-145`: energy (1.0e-6, 5.0e-6), dos / dos_spin1 /
dos_spin2 each (1.0e-4, 1.0e-4).

---

## Implementation checklist (condensed)

1. Read `.spn` (formatted or record-based binary; packed upper triangle, sigma fastest in the
   binary, per-(n,m) x/y/z lines in the formatted file); hermitize.
2. Per k, per component: SS_q = V† S_o[wm:wm+ns−1, wm:wm+ns−1] V (contiguous window block);
   FT (1/N_q) Σ e^{−i2πq·R}; divide ndegen (+ ws_distance rescatter).
3. spn_nk(m) = Re[(U† (n̂·S)(k) U)_mm] with n̂ from polar/azimuth in DEGREES; clamp to
   ±(1−1e-8) only for kpath/kslice plotting, never for DOS.
4. Moment: m = −Σ_k w_k Σ_n θ(E_F−E_nk) ⟨sigma⟩_nk (each component; w_k = 1/N; strict <,
   E_F defaults to 0 when unset); print in μ_B with f11.6, angles via acos/atan (not atan2).
5. DOS: per band, up-weight (1+spn)/2, down = 1−up; total column ×num_elec_per_state, spin
   columns ×1; smear with eta (fixed or min(|∇E|Δk·fac, max)), single-bin dump when
   eta < 2·binwidth, else ±10·eta window of w0gauss/eta; write `E total up down` in 4E16.8.
