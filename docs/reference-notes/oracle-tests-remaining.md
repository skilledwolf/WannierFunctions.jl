# postw90 reference notes: remaining oracle tests — validation-target inventory

Inventory of the postw90 test-suite cases not yet covered by a dedicated milestone target.
All paths relative to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`.
Test dirs live under `test-suite/tests/`. This note is deliberately **inventory-only**
(what is checked, to which tolerance, in which file format, driven by which keywords);
the physics/formula specs live in the sibling notes:

- AHC / A(R),H(R),B(R),C(R) operators → `berry-ahc.md`
- morb, Kubo, adaptive kmesh, Fermi scan → `kubo-morb-geninterp.md`
- SHC (qiao/ryoo), DOS, kslice machinery → `shc-dos-boltz-kslice.md`
- harness mechanics + `wannier90.x` tests → `test-suite-and-targets.md`

---

## 1. Harness recap (what "checked" means)

- `test-suite/tests/jobconfig` maps each `[testpostw90_*]` section to a `program`
  (= parser + tolerance profile in `test-suite/tests/userconfig`), the input
  (`inputs_args = ('<seed>.win', '')` → run `postw90.x <seed>`), and ONE primary `output` file.
- The golden reference in each dir is the file `benchmark.out.default.inp=<seed>.win` —
  it is a byte-copy of the primary output file (for `POSTW90_WPOUT_OK` it is the full `.wpout`).
  The parser is run on both benchmark and fresh output; extracted dicts are compared per key.
- Comparison rule (`test-suite/testcode/lib/testcode2/validation.py`): pass iff
  `abs(test-bench) < abs_tol` AND (`bench==0` or `abs((test-bench)/bench) < rel_tol`);
  a `None` entry skips that check. Tuple order in userconfig is `(absolute, relative, 'label')`.
- **Default tolerance trap**: any extracted key whose label is NOT in the `tolerance` list is
  compared with the default `(1.e-10, None)` (testcode2 `config.py:84`) — i.e. the printed value
  must match the benchmark essentially exactly. This bites e.g. `efermi` in the gyrotropic
  parser (the userconfig lists `energy`, but `parse_gyro_dat.py` emits key `efermi`) and the
  `near_neigh_*` / `b_k_*` / `w_b` keys of `parse_wpout.py`.
- Every value is round-tripped through the printed ASCII, so matching "to file precision"
  means matching after Fortran formatting (e.g. `E16.8` → 8 significant digits).

### 1.1 Parser → key map (only those used below; `test-suite/tools/parsers/`)

| program (userconfig) | parser | file layout expected | keys extracted |
|---|---|---|---|
| `POSTW90_WPOUT_OK` | `parse_wpout.py` | full `.wpout`; regex scrape | `ahc_x/y/z` (rows after a line containing `AHC`, matching `^\s*==========\s+f f f`), `morb_x/y/z` (rows after `M_orb`, matching `^\s*======================\s+f f f`), `spin_x/y/z` (`x component:` etc.), `spin_p` (`Polar theta (deg):`), `spin_a` (`Azim. phi (deg):`), plus `near_neigh_dist/mult`, `b_k_x/y/z`, `w_b`, `final_*` if present |
| `POSTW90_DOS_OK` | `parse_dos_dat.py` | `#`-comments skipped; rows of 2 cols → `energy`,`dos`; rows of 4 cols → + `dos_spin1`,`dos_spin2`; any other width = error |
| `POSTW90_CURVDAT_OK` | `parse_curv_dat.py` | 4 cols → `bandpath`,`bandcurvx/y/z`; 3 cols → no path |
| `POSTW90_MORBDAT_OK` | `parse_morb_dat.py` | 4 cols → `bandpath`,`bandmorbx/y/z`; 3 cols → `bandmorbx/y/z` (kslice case) |
| `POSTW90_SHCFREQDAT_OK` | `parse_shc_dat.py` | line 1 header, `split()[1]` MUST be `Frequency(eV)`; rows of 4 tokens → `frequency`(col2),`shc_re`(col3),`shc_im`(col4) |
| `POSTW90_SHCFERMIDAT_OK` | `parse_shc_dat.py` | header `split()[1]` MUST be `Fermi`; rows of 3 tokens → `energy`(col2),`shc`(col3) |
| `POSTW90_SHCKPATHDAT_OK` | `parse_shc_kpath_dat.py` | rows of 2 → `path`,`shc` |
| `POSTW90_SHCKPATHBANDSDAT_OK` | `parse_shc_kpath_bandsdat.py` | rows of 3 → `path`,`energy`,`shc` |
| `POSTW90_SHCKSLICEDAT_OK` | `parse_shc_kslice_dat.py` | rows of 1 → `shc` |
| `POSTW90_GYRO_OK` | `parse_gyro_dat.py` | `#`-comments skipped; rows of exactly 11 → `efermi`,`omega`,`gyro_xx/yy/zz/xy/xz/yz`,`gyro_x/y/z` |

### 1.2 Tolerance profiles (from `test-suite/tests/userconfig`, `(abs, rel)`)

| profile | tolerances |
|---|---|
| `POSTW90_WPOUT_OK` | ahc_x/y/z (1e-3, 2e-3); morb_x/y/z (1e-3, 2e-3); spin_x/y/z (1e-3, 2e-3); spin_p, spin_a (1e-2, 2e-2); everything else default (1e-10, None) |
| `POSTW90_DOS_OK` | energy (1e-6, 5e-6); dos, dos_spin1, dos_spin2 (1e-4, 1e-4) |
| `POSTW90_CURVDAT_OK` | bandpath (1e-6, 5e-6); bandcurvx/y/z (1e-6, **1e+3**) — rel effectively off (curvature spans ~1e-19…1e+1) |
| `POSTW90_MORBDAT_OK` | bandpath (1e-6, 5e-6); bandmorbx/y/z (1e-4, 1e-4) |
| `POSTW90_SHCFREQDAT_OK` | frequency (1e-6, 5e-6); shc_re (**1e+1, 1e+1**); shc_im (1e-1, 1e-1) |
| `POSTW90_SHCFERMIDAT_OK` | energy (1e-6, 5e-6); shc (1e-1, 1e-1) |
| `POSTW90_SHCKPATHDAT_OK` | path (1e-6, 5e-6); shc (1e-1, None) |
| `POSTW90_SHCKPATHBANDSDAT_OK` | path (1e-6, 5e-6); energy (1e-3, 2e-3); shc (1e-1, None) |
| `POSTW90_SHCKSLICEDAT_OK` | shc (1e-1, None) |
| `POSTW90_GYRO_OK` | omega (1e-4, 1e-2); gyro_xx/yy/zz/xy/xz/yz (1e-4, 1e-2); gyro_x/y/z (1e-4, None); efermi → **default (1e-10, None)** (label mismatch, see §1) |

Input decompression (per-dir `Makefile`): `*.mmn.bz2/*.amn.bz2/*.uHu.bz2/*.spn.bz2/*.sHu.bz2/*.sIu.bz2`
→ `bunzip2` in place; `*.chk.fmt.bz2` → `bunzip2` to `<seed>.chk.fmt` then
`w90chk2chk.x -f2u <seed>` → binary `<seed>.chk`.

---

## 2. Per-test inventory

### 2.1 Common Fe block (bcc Fe, shared by all 7 `testpostw90_fe_*` below)

`num_bands=28`, `num_wann=18`, `spinors=true` (proj `Fe: sp3d2;dxy;dxz;dyz`),
`fermi_energy=12.6279` eV, `search_shells=12`, disentanglement `dis_win_min=-8.0`,
`dis_win_max=70.0`, `dis_froz_min=-8.0`, `dis_froz_max=30.0`; `mp_grid = 2 2 2`
(8 explicit kpoints); cell (bohr units in `.win`): a1=(2.71175,2.71175,2.71175),
a2=(−2.71175,2.71175,2.71175), a3=(−2.71175,−2.71175,2.71175); `kpoint_path`
G(0,0,0)→H(0.5,−0.5,−0.5)→P(0.75,0.25,−0.25). Differences are listed per test.

### 2.2 testpostw90_fe_spin

| item | value |
|---|---|
| postw90 keys | `spin_moment=true`; `spn_formatted=true`; `kmesh = 4` (global → spin kmesh 4×4×4); `use_ws_distance = .false.`; `fermi_energy=12.6279` |
| checked output | `Fe.wpout` (program `POSTW90_WPOUT_OK`) |
| quantities | `spin_x = -0.000003`, `spin_y = -0.000000`, `spin_z = 3.090787` (tol 1e-3/2e-3); `spin_p = 0.000059`, `spin_a = 5.831720` (tol 1e-2/2e-2); plus b_k/near-neigh tables at default tol |
| inputs | `Fe.win`, `Fe.eig`, `Fe.chk.fmt.bz2` (→ .chk), `Fe.spn` (formatted, uncompressed), `Fe.mmn.bz2` (present, unused by spin), `Fe.amn` (unused) |

`.wpout` block format (`src/postw90/spin.F90:204-217`):
```
 Spin magnetic moment (Bohr magn./cell)
 ====================

       x component:  -0.000003        ! '(1x,a18,f11.6)'
 ...
 Polar theta (deg):   0.000059        ! '(/,1x,a18,f11.6)'
   Azim. phi (deg):   5.831720        ! '(1x,a18,f11.6)'
```

### 2.3 testpostw90_fe_dos_spin

| item | value |
|---|---|
| postw90 keys | `dos = true`; `spin_decomp = true`; `dos_energy_min=10.0`, `dos_energy_max=13.0`, `dos_energy_step=0.2` (→ 16 grid points); `dos_adpt_smr = false`; `dos_smr_fixed_en_width = 0.5` eV; `kmesh = 4` (→ dos kmesh 4×4×4); also `spin_moment=true`, `spn_formatted=true`; `use_ws_distance=.false.` |
| checked output | `Fe-dos.dat` (program `POSTW90_DOS_OK`) |
| quantities | 16 rows × 4 cols: `energy` (1e-6/5e-6), `dos`, `dos_spin1` (up), `dos_spin2` (down) (1e-4/1e-4) |
| inputs | `Fe.win`, `Fe.eig`, `Fe.chk.fmt.bz2`, `Fe.spn` (formatted; required by spin_decomp), `Fe.mmn.bz2` (unused), `Fe.amn` (unused) |

### 2.4 testpostw90_fe_kpathcurv

| item | value |
|---|---|
| postw90 keys | `kpath = true`; `kpath_task = bands+morb+curv`; `kpath_num_points = 10`; `uHu_formatted = .true.`; `use_ws_distance = .false.`; `fermi_energy=12.6279`; (`berry_task=ahc+morb`, `berry_kmesh=10 10 10` present but `berry` is commented out — inert) |
| checked output | `Fe-curv.dat` (program `POSTW90_CURVDAT_OK`); the run also writes `Fe-bands.dat`, `Fe-morb.dat`, `Fe-path.kpt` (not checked) |
| quantities | 19 data rows × 4 cols: `bandpath` (Å⁻¹, 1e-6/5e-6), `bandcurvx/y/z` = **−Ω(k)** summed over occupied bands, Å² (1e-6 abs) |
| inputs | `Fe.win`, `Fe.eig`, `Fe.chk.fmt.bz2`, `Fe.mmn.bz2` (→ A(R), B(R)), `Fe.uHu.bz2` (formatted text after bunzip2 → C(R)), `Fe.amn` (unused) |

### 2.5 testpostw90_fe_kpathmorbcurv

Same `.win` as 2.4 (identical file). Checked output `Fe-morb.dat` (program
`POSTW90_MORBDAT_OK`): 19 rows × 4 cols `bandpath` (1e-6/5e-6), `bandmorbx/y/z`
(1e-4/1e-4). Values are the k-resolved morb integrand in **eV·Å²** (see §3 kpath —
no Bohr-magneton conversion on the kpath). Inputs identical to 2.4.

### 2.6 testpostw90_fe_kpathmorbcurv_ws

Identical to 2.5 except `use_ws_distance = true`. Checked output `Fe-morb.dat`,
same parser/tolerances. This is the (only) kpath-morb oracle for the ws_distance
branch of the R→k Fourier interpolation.

### 2.7 testpostw90_fe_kslicemorb

| item | value |
|---|---|
| postw90 keys | `kslice = true`; `kslice_task = morb+fermi_lines`; `kslice_2dkmesh = 5 5`; `kslice_corner = 0 0 0`; `kslice_b1 = 0.5 -0.5 -0.5`; `kslice_b2 = 0.5 0.5 0.5`; `uHu_formatted=.true.`; `use_ws_distance=.false.`; `fermi_energy=12.6279` |
| checked output | `Fe-kslice-morb.dat` (program `POSTW90_MORBDAT_OK`, 3-column branch); run also writes `Fe-kslice-coord.dat`, `Fe-kslice-bands.dat`, gnu/py scripts (unchecked) |
| quantities | (5+1)·(5+1)=36 rows × 3 cols: `bandmorbx/y/z` (1e-4/1e-4), eV·Å², one row per slice point, i1 fast / i2 slow, both endpoints included (`k1=i1/n1`, `kslice.F90:341-346`) |
| inputs | as 2.4 |

### 2.8 testpostw90_fe_morb_transl_inv

| item | value |
|---|---|
| postw90 keys | `berry = true`; `berry_task = morb`; `berry_kmesh = 10 10 10`; `transl_inv_full = .true.`; `use_ws_distance = .true.`; `uHu_formatted=.true.`; `fermi_energy=12.6279` |
| checked output | `Fe.wpout` (program `POSTW90_WPOUT_OK`) |
| quantities | `morb_x= 0.0000, morb_y=-0.0000, morb_z= 0.0415` Bohr magneton/cell (tol 1e-3/2e-3) scraped from the `M_orb` block; plus b_k/near-neigh at default tol |
| inputs | as 2.4 (mmn for A/B, uHu for C, eig, chk) |

`.wpout` M_orb block (`src/postw90/berry.F90:1468-1492`, single Fermi level branch):
```
 Fermi energy (ev) =   12.627900                                ! '(/,/,1x,a,F12.6)'

 M_orb (bohr magn/cell)        x          y          z
 ======================      0.0000    -0.0000     0.0415       ! '(1x,a22,2x,3(f10.4,1x),/)'
```
(with `fermi_n > 1` a `<seed>-morb-fermiscan.dat` with rows `'(4(F12.6,1x))'` would be
written instead — not exercised here.)

### 2.9 Common Pt "qiao-kpath" block (fcc Pt; tests 2.10–2.12)

`num_bands=40`, `num_wann=18`, `spinors=true` (proj `Pt: d;s;p`), `fermi_energy=17.9919`,
`dis_win_min=0`, `dis_win_max=60`, `dis_froz_min=0`, `dis_froz_max=30`; `mp_grid=4 4 4`
(64 kpoints); cell (bohr): a=3.703863220455861 fcc set; `spn_formatted=true`;
`shc_method=qiao`; `shc_freq_scan=false`; `shc_alpha=1`, `shc_beta=2`, `shc_gamma=3`;
`kubo_adpt_smr=false`; `kubo_smr_fixed_en_width=1` (eV); `berry_curv_unit=ang2`;
`berry_task=eval_shc` + `berry_kmesh=10` present with `berry` commented out (inert, but
`berry_task` containing `shc` FORCES `shc_method` to be set — readwrite check at
`src/postw90/postw90_readwrite.F90:1108-1117`). kpoint_path: W(0.75,0.5,0.25)→L(0.5,0,0)→
G→X(0.5,0.5,0)→W→G (5 segments).

### 2.10 testpostw90_pt_kpathshc

| item | value |
|---|---|
| postw90 keys | §2.9 + `kpath = true`; `kpath_task = shc`; `kpath_num_points = 10` |
| checked output | `Pt-shc.dat` (program `POSTW90_SHCKPATHDAT_OK`) |
| quantities | 59 data rows × 2 cols: `path` (Å⁻¹, 1e-6/5e-6), `shc` (1e-1 abs, rel skipped) = k-resolved SHC Berry-curvature term at E_F, Å² (curv_unit=ang2) |
| inputs | `Pt.win`, `Pt.eig`, `Pt.chk.fmt.bz2`, `Pt.mmn.bz2`, `Pt.spn.bz2` (**formatted** text after bunzip2, since `spn_formatted=true`), `Pt.amn.bz2` (unused) |

### 2.11 testpostw90_pt_kpathbandsshc

| item | value |
|---|---|
| postw90 keys | §2.9 + `kpath = true`; `kpath_task = bands`; `kpath_bands_colour = shc`; `kpath_num_points = 10` |
| checked output | `Pt-bands.dat` (program `POSTW90_SHCKPATHBANDSDAT_OK`) |
| quantities | num_wann(18) blocks × 60 rows of 3 cols (blank line between bands): `path` (1e-6/5e-6), `energy` = interpolated eigenvalue (1e-3/2e-3), `shc` = band-resolved SHC colour (1e-1 abs, rel skipped). Band-loop outer, kpt inner (`kpath.F90:508-518`) |
| inputs | as 2.10 |

### 2.12 testpostw90_pt_ksliceshc

| item | value |
|---|---|
| postw90 keys | §2.9 + `kslice = true`; `kslice_task = shc`; `kslice_2dkmesh = 10` (one integer → 10×10); `kslice_corner = 0 0 0`; `kslice_b1 = 1 0 0`; `kslice_b2 = 0.3535533905932738 1.0606601717798214 0` |
| checked output | `Pt-kslice-shc.dat` (program `POSTW90_SHCKSLICEDAT_OK`) |
| quantities | 11·11 = 121 rows × 1 col: `shc` (1e-1 abs, rel skipped), Å² |
| inputs | as 2.10 |

### 2.13 testpostw90_pt_shc_ryoo

| item | value |
|---|---|
| postw90 keys | `berry = true`; `berry_task = shc`; `shc_method = ryoo`; `shc_freq_scan = .true.`; `berry_kmesh = 9 9 9`; `kubo_adpt_smr=.false.`; `kubo_smr_fixed_en_width = 0.1`; `kubo_eigval_max = 1000`; `kubo_freq_min=0.00`, `kubo_freq_max=7.00`, `kubo_freq_step=0.1` (→ 71 freqs); **`shc_alpha=3`, `shc_beta=2`, `shc_gamma=1`** (σ^{spin-x}_{zy} — note the permuted indices!); `use_ws_distance=.false.`; `spn_formatted=true`; `fermi_energy=18.3823`; `num_bands=24`, `num_wann=18` (disentanglement WITHOUT explicit windows → win_min/max default to eigval range), `spinors=T` (proj `Pt:l=0;l=1;l=2`); cell bohr a=3.6963 fcc; `mp_grid=4 4 4` |
| checked output | `Pt-shc-freqscan.dat` (program `POSTW90_SHCFREQDAT_OK`) |
| quantities | 71 rows: `frequency` (1e-6/5e-6), `shc_re` (**10/10**), `shc_im` (0.1/0.1), units (ħ/e)·S/cm |
| inputs | `Pt.win`, `Pt.eig`, `Pt.chk.fmt.bz2`, `Pt.mmn.bz2`, `Pt.spn.bz2` (formatted), `Pt.sHu.bz2`, `Pt.sIu.bz2` (**Fortran unformatted** after bunzip2 — no `sHu_formatted` keyword exists; readers hardcode `form='unformatted'`, `src/postw90/get_oper.F90:2633,2932`), `Pt.amn.bz2` (unused) |

### 2.14 testpostw90_pt_shc_ryoo_transl_inv

Identical to 2.13 except `transl_inv_full = .true.` and `use_ws_distance=.true.`.
Same checked output `Pt-shc-freqscan.dat`, same parser/tolerances, same inputs.
(`transl_inv_full=T` together with `transl_inv=T` is an input error; `transl_inv_full`
with `shc_method=qiao` is an input error — `postw90_readwrite.F90:854-861,1118-1125`.)

### 2.15 testpostw90_gaas_shc

| item | value |
|---|---|
| postw90 keys | `berry = true`; `berry_task = shc`; `shc_method = qiao`; `shc_freq_scan = true`; `shc_alpha=1`, `shc_beta=2`, `shc_gamma=3`; `berry_kmesh = 10` (→ 10×10×10); `fermi_energy = 7.9366`; `kubo_freq_min=0.0`, `kubo_freq_max=8.0`, `kubo_freq_step=0.01` (→ 801 freqs); `kubo_adpt_smr=false`; `kubo_smr_fixed_en_width=0.05`; **`scissors_shift = 1.117`** eV with `num_valence_bands = 8`; `exclude_bands = 1-10`; `spn_formatted=true`; `spinors=true` (proj `As:sp3`, `Ga:sp3`); `num_bands=16`, `num_wann=16` (no disentanglement); `mp_grid=4 4 4`; cell bohr fcc a=5.342256 |
| checked output | `GaAs-shc-freqscan.dat` (program `POSTW90_SHCFREQDAT_OK`) |
| quantities | 801 rows: `frequency` (1e-6/5e-6), `shc_re` (10/10), `shc_im` (0.1/0.1) |
| inputs | `GaAs.win`, `GaAs.eig`, `GaAs.chk.fmt.bz2`, `GaAs.mmn.bz2`, `GaAs.spn.bz2` (formatted), `GaAs.amn.bz2` (unused) |

### 2.16 Common Te gyrotropic block (trigonal Te; the 7 `testpostw90_te_gyrotropic*` dirs)

All seven dirs share a byte-identical `.win` **except the `gyrotropic_task` line**:

`gyrotropic=true`; `fermi_energy_min=2`, `fermi_energy_max=10`, `fermi_energy_step=2`
(→ fermi list [2,4,6,8,10] eV, 5 entries); `gyrotropic_freq_min=0.0`, `_max=0.1`,
`_step=0.05` (→ freq list [0.0,0.05,0.1] eV; each entry gets `+ i*gyrotropic_smr_fixed_en_width`,
`postw90_readwrite.F90:1752-1757`); `gyrotropic_smr_fixed_en_width=0.1` eV;
`gyrotropic_smr_max_arg=5`; `gyrotropic_degen_thresh=0.001` eV;
`gyrotropic_box_b1=0.2 0 0`, `_b2=0 0.2 0`, `_b3=0 0 0.2`,
`gyrotropic_box_center=0.33333 0.33333 0.5` (→ `box_corner = center − 0.5·(b1+b2+b3)`,
reduced coords, `postw90_readwrite.F90:736-742`); `gyrotropic_kmesh=5 5 5`;
`use_ws_distance=.false.`; `uHu_formatted=.true.`; `num_bands=12`, `num_wann=9`,
`spinors=.false.` (3 Te p-triplet projections with explicit z/x axes);
`dis_win_min=-0.5`, `dis_win_max=10`, `dis_froz_min=0.0`, `dis_froz_max=8`;
`mp_grid=2 2 2` (8 kpoints with explicit weights column); hexagonal cell in Å
(4.457, γ=120°, c=5.9581176). Also `bands_plot=true`+`kpoint_path` (wannier90-side, inert
for postw90). NB: no spin variants — `-spin` task would abort since `spinors=false`.

Per-variant task and checked file (jobconfig):

| test | `gyrotropic_task` | checked output | program |
|---|---|---|---|
| `testpostw90_te_gyrotropic` | `-C-dos-D0-Dw-K-NOA` (computes everything) | `Te-gyrotropic-C.dat` | `POSTW90_GYRO_OK` |
| `testpostw90_te_gyrotropic_C` | `-C` | `Te-gyrotropic-C.dat` | `POSTW90_GYRO_OK` |
| `testpostw90_te_gyrotropic_D0` | `-D0` | `Te-gyrotropic-D.dat` | `POSTW90_GYRO_OK` |
| `testpostw90_te_gyrotropic_Dw` | `-Dw` | `Te-gyrotropic-tildeD.dat` | `POSTW90_GYRO_OK` |
| `testpostw90_te_gyrotropic_K` | `-K` | `Te-gyrotropic-K_orb.dat` | `POSTW90_GYRO_OK` |
| `testpostw90_te_gyrotropic_NOA` | `-NOA` | `Te-gyrotropic-NOA_orb.dat` | `POSTW90_GYRO_OK` |
| `testpostw90_te_gyrotropic_dos` | `-dos` | `Te-gyrotropic-DOS.dat` | **`POSTW90_DOS_OK`** (2-col branch: `energy`,`dos`) |

Task-string matching is on the lower-cased task (`gyrotropic.F90:196-202`):
`-k`→K, `-c`→C, `-d0`→D0, `-dw`→Dw, `-noa`→NOA, `-dos`→DOS, `-spin`→spin part, `all`→all.
Inputs (all 7 dirs): `Te.win`, `Te.eig`, `Te.chk.fmt.bz2`, `Te.mmn` (uncompressed),
`Te.uHu.bz2` (formatted text), `Te.amn` (unused). `testpostw90_te_gyrotropic` additionally
carries provenance dirs `input/` (QE scf/nscf/pw2wan + pseudo + run.sh) and `reference/`
(all 6 gyrotropic dat files + `Te.wpout` + band files).

Quantities per file: 5 Fermi rows × 11 cols per (freq-)block. For `C`, `D`, `K_orb`
(static): one block at `omega=0.000000E+00`. For `tildeD`, `NOA_orb` (ω-dependent):
3 consecutive blocks (ω = 0, 0.05, 0.1), each with its own header lines. Tolerances §1.2.
For `DOS`: 5 rows × 2 cols.

### 2.17 testpostw90_example04_pdos (Cu projected DOS)

| item | value |
|---|---|
| postw90 keys | `dos = true`; `kmesh = 10` (→ dos kmesh 10×10×10); `dos_energy_min = 8`, `dos_energy_max = 10`, `dos_energy_step = 0.25` (→ 9 grid points); `dos_project 1:5` (range-vector syntax, no `=`; → project onto WFs 1–5 = the Cu d + 2 interstitial s... WFs 1-5); `use_ws_distance = .false.`; `search_shells=12`; NO smearing keys → **adaptive smearing defaults on** (Gaussian, prefactor √2, max width 1 eV); `num_bands=12`, `num_wann=7`, no spinors; `dis_win_max=38.0`, `dis_froz_max=13.0` (win_min/froz_min default); `mp_grid : 4 4 4` |
| checked output | `copper-dos.dat` (program `POSTW90_DOS_OK`, 2-col branch) |
| quantities | 9 rows × 2 cols: `energy` (1e-6/5e-6), `dos` = projected DOS in eV⁻¹ per cell (1e-4/1e-4) |
| inputs | `copper.win`, `copper.eig`, `copper.chk.fmt.bz2` — nothing else (DOS needs only H(R) + U(k); no mmn/amn/spn shipped) |

Note: NO `fermi_energy` in this `.win` — allowed because plain dos_plot never references it.

---

## 3. Exact output-file formats (Fortran writers)

All formats verified against both source and the shipped benchmarks.

### 3.1 kpath files (`src/postw90/kpath.F90`)

| file | rows | format | content, column order |
|---|---|---|---|
| `<seed>-bands.dat` (colour ≠ none) | per band, per kpt; blank line (`write(*,*) ' '`) after each band | `'(3E16.8)'` (line 513) | `xval(kpt)`, `eig(band,kpt)`, `color(band,kpt)`; colour=none → `'(2E16.8)'` (line 511) |
| `<seed>-curv.dat` | one per kpt + trailing blank line | `'(4E16.8)'` (line 646) | `xval`, `−Ω_x`, `−Ω_y`, `−Ω_z` (sign flip `curv = -curv` at line 641; Å², ÷bohr² only if `berry_curv_unit=bohr2`) |
| `<seed>-morb.dat` | one per kpt + trailing blank line | `'(4E16.8)'` (line 735) | `xval`, `Morb_k(x,y,z)` where `Morb_k = −(img + imh − 2·E_F·imf)/2` summed over first index (lines 380-386); **eV·Å², no curv_unit and no Bohr-magneton conversion** |
| `<seed>-shc.dat` | one per kpt + trailing blank line | `'(2E16.8)'` (line 817) | `xval`, `shc_k_fermi(1)` (Å²; ÷bohr² iff bohr2) |
| `<seed>-path.kpt` | header = total_pts (list-directed), then rows | `'(3f12.6,3x,f4.1)'` (line 469) | k in reduced coords + weight 1.0 |

kpath grid (`kpath.F90:1263-1315`): segment lengths in Å⁻¹ via the reciprocal metric;
`kpath_pts(1) = kpath_num_points`, `kpath_pts(i>1) = nint(num_points·len_i/len_1)`
(**relative to segment 1, not the longest**); `total_pts = sum + 1`;
`xval(1)=0`, increments `len_i/pts_i`, last point pinned to `sum(kpath_len)`.
Bands colour=`spin` is clamped to ±(1−eps8); colour=`shc` is NOT clamped
(`kpath.F90:344-370`). kpath demands exactly one Fermi energy (lines 189-196).

### 3.2 kslice files (`src/postw90/kslice.F90`)

| file | rows | format | content |
|---|---|---|---|
| `<seed>-kslice-coord.dat` | (n1+1)(n2+1) | `'(2E16.8)'` (line 521) | in-plane Cartesian coords `kpt_x, kpt_y` (Å⁻¹) |
| `<seed>-kslice-bands.dat` | one eigenvalue per line, all kpts of band 1, then band 2… | `'(E16.8)'` (line 527) | interpolated eigenvalues (fermi_lines w/o colour) |
| `<seed>-kslice-morb.dat` | (n1+1)(n2+1) + trailing blank | `'(4E16.8)'` writing 3 values (line 571) | `Morb_k(x,y,z)` per point, eV·Å² (same combination as kpath; kslice.F90:437-443) |
| `<seed>-kslice-shc.dat` | (n1+1)(n2+1) + trailing blank | `'(1E16.8)'` (line 567) | `shc_k_fermi(1)`; ÷bohr² iff `berry_curv_unit=bohr2` (line 565) |

Grid: `itot = 0 … (n1+1)(n2+1)−1`, `i2 = itot/(n1+1)` slow, `i1` fast;
`k1 = i1/n1`, `k2 = i2/n2` → **k = corner + k1·b1 + k2·b2 covers both endpoints**
(duplicated boundary for periodic slices). One Fermi level required (lines 189-196).

### 3.3 dos file (`src/postw90/dos.F90:340-348`)

`<seed>-dos.dat`, no header. Rows `'(4E16.8)'`: `omega, dos_all(ifreq, 1:ndim)` —
ndim = 1 (plain or projected DOS) or 3 (with `spin_decomp`: total, spin-up, spin-down).
Energy grid: `num_freq = nint((emax−emin)/step) + 1`,
`d_omega = (emax−emin)/(num_freq−1)`, `omega_i = emin + (i−1)·d_omega`
(dos.F90:146-157) — the step is re-fitted, not taken literally.

### 3.4 SHC scan files (`src/postw90/berry.F90:1704-1727`)

`<seed>-shc-fermiscan.dat`:
```
#No.   Fermi energy(eV)   SHC((hbar/e)*S/cm)          ! '(a,3x,a,3x,a)'
   1    17.991900     0.12345678E+04                  ! '(I4,1x,F12.6,1x,E17.8)'
```
`<seed>-shc-freqscan.dat`:
```
#No.   Frequency(eV)   Re(sigma)((hbar/e)*S/cm)   Im(sigma)((hbar/e)*S/cm)   ! '(a,3x,a,3x,a,3x,a)'
   1     0.000000    -0.25760024E+04    0.00000000E+00                       ! '(I4,1x,F12.6,1x,1x,2(E17.8,1x))'
```
The parser keys off `split()[1]` of line 1 (`Fermi` vs `Frequency(eV)`) — reproduce the
header token-for-token. Unit conversion applied immediately before writing:
`fac = 1.0e8·e²/(ħ·V_cell)/2` (berry.F90:1690; the /2 is the ħ/2 of the spin operator),
result in (ħ/e)·S/cm.

### 3.5 gyrotropic files (`src/postw90/gyrotropic.F90:1146-1290`)

Filename pattern: `<seed>-gyrotropic-<NAME>.dat` with NAME ∈
{`C`, `D`, `tildeD`, `K_orb`, `K_spin`, `NOA_orb`, `NOA_spin`, `DOS`}.
Header (per file): `write(*,*) "#"//comment` then `write(*,*) "# in units of [ <units> ] "`
(list-directed → leading blank). Then per frequency block:

- symmetrized tensors (C, D, tildeD, K_*): banner line
  `'(a1,29x,a1,38x,a14,37x,a2,14x,a15,14x,a1)'` →
  `#                             |                                      symmetric part                                     ||              asymmetric part              |`
  then `'(11a15)'` → `   # EFERMI(eV)      omega(eV)             xx …  yz              x              y              z`;
  rows `'(11E15.6)'`: `E_f, ω, xx, yy, zz, (xy+yx)/2, (xz+zx)/2, (yz+zy)/2, (T_23−T_32)/2, (T_31−T_13)/2, (T_12−T_21)/2`.
- NOA (symmetrize=.false.): single header `'(11a15)'` with labels
  `yzx zxy xyz yzy yzz zxz xyy xyx zxx`; rows `'(11E15.6)'`:
  `E_f, ω, T11, T22, T33, T12, T13, T23, T32, T31, T21` (raw tensor, no sym/antisym split).
- DOS: header `'(2a15)'` → `   # EFERMI(eV) `; rows `'(11E15.6)'` with 2 values: `E_f, dos`.
- After each block: two empty list-directed writes (two blank-ish lines).

Comments/units written per quantity (gyrotropic.F90:425-570): C `"the C tensor -- Eq. B6 of TAS17"` /
`Ampere/cm`; D `"the D tensor -- Eq. 2 of TAS17"` / `dimensionless`; tildeD
`"the tildeD tensor -- Eq. 12 of TAS17"` / `dimensionless`; K_orb
`"orbital part of the K tensor -- Eq. 3 of TAS17"` / `Ampere`; NOA_orb
`"the tensor $gamma_{abc}^{orb}$ (Eq. C12,C14 of TAS17)"` / `Ang`; DOS
`"density of states"` / `eV^{-1}.Ang^{-3}`. (Shipped Te benchmarks carry an older NOA
comment "Eq. C10" — comments are parser-skipped, only column data are compared.)
Post-processing factors before writing (SI constants §5): K_orb `e²/(2ħV_c)`;
K_spin `−10²⁰·eħ/(2m_e·V_c)`; D, tildeD `1/V_c`; C `10⁸·e²/(2πħ·V_c)`;
NOA_orb `10¹⁰·e/(V_c·ε₀)`; NOA_spin `10³⁰·ħ²/(V_c·ε₀·m_e)`; DOS `1/V_c` (V_c in Å³).

### 3.6 `.wpout` scraped blocks

See 2.2 (spin) and 2.8 (M_orb). Both are fixed-point (`f11.6` / `f10.4`) — the regexes in
`parse_wpout.py` do **not** match exponent notation, so any reimplementation must print
these blocks in `F`-format for the harness to extract them.

---

## 4. Keyword defaults (source of truth)

From `src/postw90/postw90_types.F90` (compile-time defaults) and
`src/postw90/postw90_readwrite.F90` / `src/readwrite.F90` (read-time defaults & derived values):

| keyword | default | units / notes |
|---|---|---|
| `berry`, `gyrotropic`, `dos`, `kpath`, `kslice`, `spin_moment`, `spin_decomp` | `.false.` | postw90_types.F90:53-61 |
| `use_ws_distance` | `.true.` | src/types.F90 (**postw90 default is T; most of these tests force F**) |
| `spn_formatted`, `uhu_formatted` | `.false.` | binary Fortran unformatted by default; `.sHu`/`.sIu` are ALWAYS unformatted (no keyword) |
| `fermi_energy` | 0.0 (list [0.0] if nothing set) | eV; `fermi_energy_min/max/step` → list with `n = nint(|max−min|/step)+1`, step re-fitted (readwrite.F90:630-693); `fermi_energy_max` default `min+1.0`, `fermi_energy_step` default 0.01 |
| `kmesh` / `kmesh_spacing` (global) | unset | 1 int → n×n×n; 3 ints → n1×n2×n3; module override `berry_kmesh`, `dos_kmesh`, `spin_kmesh`, `gyrotropic_kmesh` (same syntax, postw90_readwrite.F90:1778-1998); required for the module used |
| `kpath_task` | `'bands'` | any of bands/curv/morb/shc (+ combos with `+`) |
| `kpath_num_points` | 100 | points on FIRST path segment |
| `kpath_bands_colour` | `'none'` | none/spin/shc |
| `kslice_task` | `'fermi_lines'` | fermi_lines/curv/morb/shc combos; curv+morb, shc+morb, shc+curv forbidden |
| `kslice_2dkmesh` | 50 50 | 1 or 2 ints |
| `kslice_corner` | 0 0 0 | reduced |
| `kslice_b1` / `kslice_b2` | (1,0,0) / (0,1,0) | reduced |
| `berry_task` | `' '` (required when berry=T) | ahc/morb/kubo/sc/shc/kdotp |
| `berry_curv_unit` | `'ang2'` | `'bohr2'` divides curv/shc outputs by bohr² (bohr=0.52917720859 Å) |
| `berry_curv_adpt_kmesh` | 1 | adaptive refinement OFF by default |
| `berry_curv_adpt_kmesh_thresh` | 100.0 | Å² |
| `transl_inv` | `.false.` | MV97 diagonal A(R) |
| `transl_inv_full` | `.false.` | full translation-invariant scheme; incompatible with `transl_inv=T` and with `shc_method=qiao` |
| `smr_type` / module `*_smr_type` | Gaussian (`type_index=0`) | |
| `adpt_smr` (and dos/kubo copies) | `.true.` | adaptive smearing ON by default |
| `adpt_smr_fac` | √2 | |
| `adpt_smr_max` | 1.0 | eV |
| `smr_fixed_en_width` | 0.0 | eV; with adaptive off and width 0 → pure histogram |
| `kubo_smr_fixed_en_width` | inherits `smr_fixed_en_width` | eV |
| `kubo_adpt_smr` etc. | inherit global | |
| `kubo_freq_min` | 0.0 | eV |
| `kubo_freq_max` | frozen: `froz_max − E_F(1) + 0.6667`; else `max(eig)−min(eig)+0.6667` | eV (postw90_readwrite.F90:1687-1697) |
| `kubo_freq_step` | 0.01 | eV; `nfreq = nint((max−min)/step)+1` (min 2), step re-fitted |
| `kubo_eigval_max` | frozen: `froz_max + 0.6667`; else `max(eig)+0.6667` | eV — band cutoff for sums (postw90_readwrite.F90:1760-1767) |
| `shc_freq_scan` | `.false.` | fermiscan by default |
| `shc_alpha`, `shc_beta`, `shc_gamma` | 1, 2, 3 | Cartesian indices of σ^{spin-γ}_{αβ} |
| `shc_method` | `' '` — REQUIRED (qiao or ryoo) whenever `berry_task` contains `shc` even if `berry=F` | |
| `shc_bandshift`(_firstband/_energyshift) | F / 0 / 0.0 | rigid conduction shift; mutually exclusive with scissors_shift |
| `scissors_shift` | 0.0 | eV; active only with `num_valence_bands > 0` and abs > 1e-7 |
| `dos_task` | `'dos_plot'` | |
| `dos_energy_min` | `min(eig) − 0.6667` | eV |
| `dos_energy_max` | frozen: `froz_max + 0.6667`; else `max(eig)+0.6667` | eV |
| `dos_energy_step` | 0.01 | eV |
| `dos_project` | all WFs (`num_project = num_wann`) | range-vector, e.g. `dos_project 1:5` |
| `spin_axis_polar`, `spin_axis_azimuth` | 0.0, 0.0 | degrees |
| `gyrotropic_task` | `'all'` | matched lower-case: `-k -c -d0 -dw -spin -noa -dos` |
| `gyrotropic_kmesh` | falls back to global `kmesh` | |
| `gyrotropic_freq_*` | min 0.0 / max = kubo_freq_max default / step 0.01 | eV; each freq gets `+ i·gyrotropic_smr_fixed_en_width` |
| `gyrotropic_smr_fixed_en_width` | inherits `smr_fixed_en_width` | eV; gyrotropic adaptive smearing is hard-disabled |
| `gyrotropic_smr_max_arg` | 5.0 (via `smr_max_arg` default 5.0) | dimensionless cutoff of smearing argument |
| `gyrotropic_degen_thresh` | 0.0 | eV |
| `gyrotropic_box_b1/2/3` | identity vectors | reduced; `gyrotropic_box_center` default → corner at 0 |
| `gyrotropic_band_list` | all `1…num_wann` | |
| `gyrotropic_eigval_max` | = kubo_eigval_max default | eV |
| `exclude_bands` | none | GaAs test excludes 1-10 (already baked into eig/mmn/spn band counts) |

---

## 5. Traps for the Julia reimplementation

1. **CODATA2006** constants (`src/constants.F90:161-187`, the build default):
   `e=1.602176487e-19`, `m_e=9.10938215e-31`, `ħ=1.054571628e-34`,
   `ε₀=8.854187817e-12`, `bohr=0.52917720859 Å`, `eV_au=3.674932540e-2`.
   All SHC/gyrotropic/morb unit factors use these; a different CODATA set shifts
   results by ~1e-8 relative — visible at the 1e-10 default tolerance keys.
2. **use_ws_distance default flip**: postw90 defaults to `.true.`, but every test here
   except `fe_kpathmorbcurv_ws`, `fe_morb_transl_inv`, `pt_shc_ryoo_transl_inv`
   explicitly sets `.false.`. Read it from the `.win`, never assume.
3. **kpath/kslice morb has no unit conversion**: `-morb.dat` / `-kslice-morb.dat` carry
   `−(g+h−2E_F·f)/2` in eV·Å² directly; only the `.wpout` M_orb block applies
   `fac = −eV_au/bohr²` to reach Bohr magnetons (berry.F90:1447). `berry_curv_unit`
   affects curv and shc outputs but NOT morb.
4. **Sign of curv output**: kpath `-curv.dat` and kslice `-kslice-curv.dat` print
   **minus** the Berry curvature; kpath/kslice shc values are printed as-is.
5. **kpath point counts scale to segment 1** (`nint(num_points·len_i/len_1)`), not to the
   longest segment; total = Σ+1; last point pinned to end. Fe path (10) → 19 pts;
   Pt path (10, 5 segments) → 59 pts (bands file 60 rows/band incl. duplicate? count from
   benchmark: 60 rows/band, 59+1).
6. **kslice covers both endpoints** (`k1=i1/n1` for `i1=0…n1`) — (n+1)² points with
   duplicated zone boundary; row order i1-fast.
7. **SHC ryoo index permutation**: the Pt ryoo tests compute σ^{spin-x}_{zy}
   (`shc_alpha=3, shc_beta=2, shc_gamma=1`), NOT the usual σ^{spin-z}_{xy}.
8. **shc_method is mandatory** whenever `berry_task` contains `shc` — even when the
   `berry` flag itself is false (Pt kpath tests keep `berry_task=eval_shc` around and
   must therefore set `shc_method=qiao`).
9. **SHC overall factor** `1.0e8·e²/(ħ·V_cell)/2` → (ħ/e)S/cm; the ½ compensates using
   σ (Pauli) instead of ħσ/2 as spin operator. Same factor for fermiscan and freqscan.
10. **kubo_eigval_max default** couples to the disentanglement window
    (`froz_max+0.6667` when frozen states exist). Pt-ryoo overrides with 1000 (i.e. no
    cutoff); GaAs (no disentanglement) defaults to `max(eig)+0.6667`. Reproduce exactly —
    it gates which bands enter the SHC sums.
11. **Frequency/energy grids are re-fitted**: `nfreq = nint((max−min)/step)+1` and the
    step is recomputed from the endpoints (kubo, gyrotropic, dos, fermi list alike);
    a literal `arange` with the input step drifts for non-commensurate ranges.
12. **Gyrotropic freq_list is complex**: `ω_i + i·gyrotropic_smr_fixed_en_width`
    (Te tests: 0.1i eV). The imaginary part is the NOA/tildeD broadening.
13. **Gyrotropic smearing** is always non-adaptive; `gyrotropic_smr_max_arg` (Te: 5)
    truncates the delta-smearing argument — values beyond give exactly 0 contribution,
    which is why several Te benchmark rows are exact zeros.
14. **gyrotropic box**: integration restricted to the reduced-coordinate box
    `corner + m1·b1 + m2·b2 + m3·b3` with corner = center − ½Σb_i; the 5×5×5
    `gyrotropic_kmesh` samples the box, not the full BZ.
15. **spin_decomp requires spinor .spn** and `num_elec_per_state==1`; dos ndim jumps
    1→3 and the dat file gains 2 columns — same file name.
16. **`dos_project 1:5`** uses the colon range syntax without `=`; projection uses
    |U(k)|² weights of the listed WFs (dos_get_k), still one `dos` column.
17. **File formats of matrix inputs differ per test**: `.spn` formatted (Pt/GaAs/Fe:
    `spn_formatted=true`), `.uHu` formatted (Fe/Te), but `.sHu`/`.sIu` are ALWAYS
    Fortran unformatted (sequential records with 4-byte markers) — `w90chk2chk.x -f2u`
    is needed for the `.chk.fmt.bz2` → `.chk` conversion.
18. **parse_shc_dat is strict about the header**: second whitespace token must be
    `Fermi` or `Frequency(eV)`; row token counts must be exactly 3 / 4.
19. **The wpout parser's F-format-only regexes** (no `E`-notation) for spin/morb/AHC
    blocks; and unlisted keys (b_k vectors, shell distances) fall to the 1e-10 default
    tolerance — the kmesh/b-vector machinery must be bit-compatible in print.
20. **Scissors shift in GaAs SHC**: applied inside get_oper H(R) construction
    (needs `num_valence_bands=8`); `exclude_bands=1-10` means the shipped eig/mmn/spn
    already contain only 16 bands — do not re-exclude.
21. **kpath colour clamps**: spin colour clamped to ±(1−1e-8); shc colour unclamped.
22. **Degeneracy**: `gyrotropic_degen_thresh=0.001` groups near-degenerate bands
    (band-derivative regularisation); kpath/kslice SHC uses
    `pw90_band_deriv_degen` defaults (`use_degen_pert=.false.`, `degen_thr=1e-4`).

---

## 6. Summary: feature → tests → checked files

| feature (module / branch) | tests | checked file(s) | profile |
|---|---|---|---|
| spin magnetic moment (spin.F90) | fe_spin | `Fe.wpout` | WPOUT (spin_* 1e-3) |
| DOS + spin_decomp (dos.F90) | fe_dos_spin | `Fe-dos.dat` (4 col) | DOS 1e-4 |
| DOS projected (dos_project, adaptive smr) | example04_pdos | `copper-dos.dat` (2 col) | DOS 1e-4 |
| kpath Berry curvature (−Ω) | fe_kpathcurv | `Fe-curv.dat` | CURV 1e-6 abs |
| kpath morb integrand | fe_kpathmorbcurv | `Fe-morb.dat` | MORB 1e-4 |
| kpath morb + use_ws_distance | fe_kpathmorbcurv_ws | `Fe-morb.dat` | MORB 1e-4 |
| kslice morb (+fermi_lines) | fe_kslicemorb | `Fe-kslice-morb.dat` (3 col) | MORB 1e-4 |
| berry_task=morb + transl_inv_full + ws | fe_morb_transl_inv | `Fe.wpout` (M_orb) | WPOUT morb_* 1e-3 |
| kpath SHC (qiao) | pt_kpathshc | `Pt-shc.dat` | SHCKPATH 1e-1 |
| kpath bands coloured by SHC | pt_kpathbandsshc | `Pt-bands.dat` (3 col) | SHCKPATHBANDS |
| kslice SHC | pt_ksliceshc | `Pt-kslice-shc.dat` (1 col) | SHCKSLICE 1e-1 |
| SHC freq-scan, ryoo (sHu/sIu) | pt_shc_ryoo | `Pt-shc-freqscan.dat` | SHCFREQ (10/0.1) |
| SHC ryoo + transl_inv_full + ws | pt_shc_ryoo_transl_inv | `Pt-shc-freqscan.dat` | SHCFREQ |
| SHC freq-scan, qiao + scissors + exclude_bands | gaas_shc | `GaAs-shc-freqscan.dat` | SHCFREQ |
| gyrotropic all tasks | te_gyrotropic | `Te-gyrotropic-C.dat` | GYRO 1e-4 |
| gyrotropic C (conductivity tensor) | te_gyrotropic_C | `Te-gyrotropic-C.dat` | GYRO |
| gyrotropic D0 (static D) | te_gyrotropic_D0 | `Te-gyrotropic-D.dat` | GYRO |
| gyrotropic Dw (ω-dependent D̃, 3 blocks) | te_gyrotropic_Dw | `Te-gyrotropic-tildeD.dat` | GYRO |
| gyrotropic K (orbital gme) | te_gyrotropic_K | `Te-gyrotropic-K_orb.dat` | GYRO |
| gyrotropic NOA (orbital, non-symmetrized) | te_gyrotropic_NOA | `Te-gyrotropic-NOA_orb.dat` | GYRO |
| gyrotropic DOS-at-E_F | te_gyrotropic_dos | `Te-gyrotropic-DOS.dat` (2 col) | DOS |

Related already-specced oracle tests (context): `testpostw90_fe_ahc*`, `testpostw90_fe_morb`,
`testpostw90_pt_shc` → `berry-ahc.md` / `kubo-morb-geninterp.md` / `shc-dos-boltz-kslice.md`.
