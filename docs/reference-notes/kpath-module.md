# postw90 reference notes: kpath module (kpath = true)

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Main source: `src/postw90/kpath.F90` (module `w90_kpath`, entry point `k_path`, lines 62-1144;
private helpers `k_path_print_info` 1149-1222 and `k_path_get_points` 1225-1315).

Cross-references to existing notes:
- `berry-ahc.md` §C — `berry_get_imf_klist` / `berry_get_imfgh_klist` (J0/J1/J2 trace formulas)
- `shc-dos-boltz-kslice.md` §1.3 — `spin_get_nk`; §1.5-1.6 — `get_SHC_R`, `berry_get_shc_klist`
- `interpolation.md` — `pw90common_fourier_R_to_k`, `utility_diagonalize`
- `berry-ahc.md` §A — `get_HH_R`, `get_AA_R`, Wigner-Seitz setup, `use_ws_distance`

---

## A. Input keywords and defaults

| keyword | default | type / units | where |
|---|---|---|---|
| `kpath` | `.false.` | logical | `postw90_types.F90:53` (`pw90_calculation_type%kpath`), parsed `postw90_readwrite.F90:424-425` |
| `kpath_task` | `'bands'` | string (character(len=20)) | `postw90_types.F90:111`, parsed `postw90_readwrite.F90:1170-1171` |
| `kpath_num_points` | `100` | integer (points in FIRST segment) | `postw90_types.F90:112`, parsed `postw90_readwrite.F90:1185-1186` |
| `kpath_bands_colour` | `'none'` | string: none/spin/shc | `postw90_types.F90:113`, parsed `postw90_readwrite.F90:1193-1194` |
| `berry_curv_unit` | `'ang2'` | string: ang2/bohr2 | `postw90_types.F90:168` (`pw90_berry_mod_type%curv_unit`), parsed `postw90_readwrite.F90:893-898` |
| `fermi_energy` | unset | eV; builds `fermi_energy_list` of length `fermi_n` | `src/readwrite.F90:609-692` (see berry-ahc.md §D) |
| `kubo_adpt_smr` | inherits `adpt_smr` (default `.true.`) | logical | `postw90_readwrite.F90:908-911`, `postw90_types.F90:133` |
| `kubo_smr_fixed_en_width` | inherits `smr_fixed_en_width` (default `0.0`) | eV | `postw90_readwrite.F90:931-936`, `postw90_types.F90:136` |
| `use_ws_distance` | `.true.` | logical | `src/types.F90:87-94` |

Task/colour matching is by SUBSTRING (`index(...)>0`), so `kpath_task = bands+curv`,
`bands+morb+curv`, etc. all work; flags set in `kpath.F90:175-178`:

```
plot_bands = index(task,'bands')>0 ; plot_curv = index(task,'curv')>0
plot_morb  = index(task,'morb')>0  ; plot_shc  = index(task,'shc')>0
```

Validation (`postw90_readwrite.F90:1152-1208`, subroutine `w90_wannier90_readwrite_read_pw90_kpath`):
- if `kpath=T`, task must contain at least one of bands/curv/morb/shc, else
  `'Error: value of kpath_task not recognised in w90_wannier90_readwrite_read'` (1173-1178)
- if `kpath=T` and no `kpoint_path` block:
  `'Error: a kpath plot has been requested but there is no kpoint_path block'` (1180-1183)
- `kpath_num_points < 0` → `'Error: kpath_num_points must be positive'` (1188-1191; note 0 passes)
- bands_colour must contain none/spin/shc (1196-1200)
- `kpath_task` may not contain BOTH 'shc' and 'spin' (1202-1206):
  `"Error: kpath_task cannot include both 'shc' and 'spin'"`

Runtime checks in `kpath.F90`:
- shc plotting (`plot_shc` OR `bands_colour=='shc'`) requires FIXED smearing:
  `if (pw90_berry%kubo_smearing%use_adaptive)` →
  `'Error: Must use fixed smearing when plotting spin Hall conductivity'` (181-188).
  Since `kubo_adpt_smr` defaults true, the user MUST set `kubo_adpt_smr = false` (Pt tests do).
- `plot_shc` requires exactly one Fermi level: `fermi_n==0` → `'Error: must specify Fermi
  energy'`; `fermi_n/=1` → `'Error: kpath plot only accept one Fermi energy, use fermi_energy
  instead of fermi_energy_min'` (189-198)
- `k_path_print_info` (1149-1222): curv → `'Must specify one Fermi level when kpath_task=curv'`
  if `fermi_n /= 1` (1194-1197); same for morb (1202-1205) and shc (1215-1218).

### A.1 kpoint_path block parsing

`w90_readwrite_read_kpath` (`src/readwrite.F90:479-542`) →
`w90_readwrite_get_keyword_kpath` (`src/readwrite.F90:4662-4757`).

- Block delimited by lines starting with `begin kpoint_path` / `end kpoint_path`.
- Each interior line defines ONE segment = TWO special points:
  `read(dummy,*) labels(2i-1), (points(:,2i-1)), labels(2i), (points(:,2i))`
  (readwrite.F90:4733-4738), i.e. `label1 kx ky kz label2 kx ky kz`, fractional coords of the
  reciprocal lattice.
- `kpoint_path%labels(:)` is `character(len=20)` (`src/types.F90:260`); the second index of
  `points(3, :)` is 2×(number of segments).
- Labels are UPPERCASED character-by-character (`x` → `X`) (readwrite.F90:4741-4749).
- The separate keyword `bands_num_points` fills `kpoint_path%num_points_first_segment`
  (default 100, types.F90:259) — that is wannier90.x's band-plot density and is NOT used by
  postw90 kpath, which uses `kpath_num_points` (`pw90_kpath%num_points`).

---

## B. Path construction: k_path_get_points (kpath.F90:1225-1315)

Inputs: `kpoint_path%points(3, 2*num_paths)` (fractional), `recip_lattice(3,3)` from
`utility_recip_lattice_base` (rows = b_i in Å⁻¹, includes 2π), `pw90_kpath%num_points`.

```
num_paths = size(kpoint_path%labels)/2                              ! 1255-1260
recip_metric = recip_lattice . recip_lattice^T                      ! utility_metric, utility.F90:416-442
                                                                    ! metric(i,j) = sum_l lattice(i,l)*lattice(j,l)
vec = points(:,2p) - points(:,2p-1)                                 ! 1264-1265
kpath_len(p) = sqrt( vec . recip_metric . vec )                     ! 1266-1267   [Å^-1, includes 2π]

kpath_pts(1) = num_points                                           ! 1272-1273
kpath_pts(p) = nint( real(num_points,dp) * kpath_len(p)/kpath_len(1) )   ! p>1, 1275-1276
               ! Fortran nint = round half away from zero
total_pts    = sum(kpath_pts) + 1                                   ! 1279
```

Point placement (1292-1313), counter runs 1..sum(kpath_pts), then one extra final point:

```
do p = 1, num_paths
  do i = 1, kpath_pts(p)
    counter += 1
    xval(1) = 0.0
    xval(counter) = xval(counter-1) + kpath_len(p)/real(kpath_pts(p),dp)     ! 1296-1301
    plot_kpoint(:,counter) = points(:,2p-1) + (points(:,2p)-points(:,2p-1)) *
                             real(i-1,dp)/real(kpath_pts(p),dp)              ! 1302-1306
  end do
end do
xval(total_pts)          = sum(kpath_len)                                    ! 1312  (FORCED, not accumulated)
plot_kpoint(:,total_pts) = points(:, 2*num_paths)                            ! 1313
```

Consequences (verified against Fe benchmark, see §G):
- Segment p is sampled at fractions 0, 1/n_p, ..., (n_p−1)/n_p — the segment ENDPOINT is never
  sampled within the segment; it appears once as the first point (fraction 0) of the NEXT
  segment. **No duplicated vertices.** Only the global last point closes the path.
- **xval trap**: crossing from segment p to p+1, the increment used is the NEW segment's step
  `kpath_len(p+1)/kpath_pts(p+1)`, so the x-coordinate of an interior vertex is
  `Σ_{q≤p} kpath_len(q) − kpath_len(p)/n_p + kpath_len(p+1)/n_{p+1}` — NOT the cumulative
  length `Σ kpath_len(1:p)`. Example (Fe test, n₁=10, n₂=9, len₁=2.1892688, len₂=1.8959624):
  vertex H is row 11 at x = 9·len₁/10 + len₂/9 = 2.1810044, while the gnuplot tick for H sits
  at len₁ = 2.1892688. Only the final row is forced to x = len₁+len₂ = 4.0852312.
- The gnuplot/python tick locations use exact cumulative `sum(kpath_len(1:j))` (kpath.F90:
  538-540, 606-610), inconsistent with the .dat x values by design.
- xval units: Å⁻¹ including the 2π factor (distance in reciprocal Cartesian space).

MPI: points are scattered with `comms_array_split`/`comms_scatterv` and gathered back in
original order (289-303, 419-454) — serial replication just evaluates points in order.

---

## C. Per-k quantities for each task

All tasks first need `get_HH_R` (207-211). Additional operators (all built once, before the
k loop):

| condition | operators (files read) | kpath.F90 lines |
|---|---|---|
| `plot_curv .or. plot_morb` | `get_AA_R` (.mmn) | 213-224 |
| `plot_morb` | `get_BB_R` (.mmn+eig), `get_CC_R` (.uHu; `uHu_formatted` flag) | 225-238 |
| `plot_shc .or. bands_colour=='shc'` | `get_AA_R` (.mmn), `get_SS_R` (.spn), `get_SHC_R` (.spn+.mmn → SH_R, SHR_R, SR_R) | 240-263 |
| `plot_bands .and. bands_colour=='spin'` | `get_SS_R` (.spn) | 265-271 |

**Trap**: kpath only ever calls `get_SHC_R` (the Qiao-method operators). The Ryoo operators
(`get_SH_R`/`get_SAA_R`/`get_SBB_R`, arrays SAA_R/SBB_R) are passed to `berry_get_shc_klist`
but never populated here → kpath SHC effectively requires `shc_method = qiao` (which is also
what the Pt tests set). `shc_method` itself is only mandatory when `berry_task` contains shc.

Main loop (322-416), one k-point at a time; `kpt(:)` in fractional coords.

### C.1 bands (plot_bands, 325-368)

```
H(k) = fourier_R_to_k(HH_R, kpt, deriv=0)                 ! 326-328
eig(:,k), UU = utility_diagonalize(H, num_wann)           ! 330  (ascending eigenvalues, eV)
```

Colour column (only if `bands_colour /= 'none'`):

- **spin** (337-354): `spin_get_nk` → `spn_k(n) = diag(UU† (α·σ_x + β·σ_y + γ·σ_z interpolated
  from SS_R) UU)_n` = ⟨ψ_nk|σ·n̂|ψ_nk⟩ ∈ [−1,1] with quantization axis from
  `spin_axis_polar/azimuth` (defaults 0,0 → ẑ). See shc-dos-boltz-kslice.md §1.3.
  Then clamped to prevent bands vanishing in the gnuplot z-range (348-354):
  `color = min(max(color, −1+eps8), 1−eps8)` with `eps8 = 1.0e-8_dp` (`src/constants.F90:76`).
- **shc** (355-367): `berry_get_shc_klist(..., shc_k_band=shc_k_band)` — band-resolved
  Berry-curvature-like term of QZYZ18 (PRB 98, 214402):

```
Ω^{spin γ}_{n,αβ}(k) = Σ_{m≠n} (−2) Im[ js_k(n,m) · i(e_m−e_n) · A_{β}(m,n) ] /
                       ((e_m−e_n)² + η²)                          ! berry.F90:2869-2890
shc_k_band(n) = Ω^{spin γ}_{n,αβ}(k)     (no occupation factor)   ! berry.F90:2900-2901
```

  with `js_k(n,m) = ⟨ψ_n|½{σ_γ, v_α}|ψ_m⟩` (Qiao branch, see shc-dos-boltz-kslice.md §1.6),
  `η = kubo_smr_fixed_en_width` (adaptive forbidden), α/β/γ = `shc_alpha/beta/gamma`
  (defaults 1/2/3 = σ^{spin z}_{xy}), and pairs skipped when `eig(m) > kubo_eigval_max .or.
  eig(n) > kubo_eigval_max` (berry.F90:2867). `kubo_eigval_max` default = `dis_froz_max +
  0.6667` if frozen window set, else `maxval(eigval)+0.6667` (postw90_readwrite.F90:1758-1768).
  Units: Å². No clamping.

### C.2 curv (plot_curv, 388-402)

If morb is NOT also requested, calls `berry_get_imf_klist` (T=0 step occupations at
`fermi_energy_list`; occupancy `occ(i)=1` iff `eig(i) < ef`, `pw90common_get_occ`,
postw90_common.F90:942-985); otherwise reuses `imf_k_list` from the imfgh call (389-398):

```
my_curv(k, i) = sum(imf_k_list(:, i, 1))     ! i = axial component: 1=Ω_yz(=Ω_x), 2=Ω_zx, 3=Ω_xy
                                             ! sum over first index = J0+J1+J2 terms   399-401
```

i.e. curv(k,i) = Ω_i(k) = Σ_{n occ} Ω_{n,i}(k) in Å² (WYSV06/LVTS12 traces, berry-ahc.md §C).
Post-processing on root:

```
if (berry_curv_unit == 'bohr2') curv = curv/bohr**2       ! 474   bohr = 0.52917720859 Å (CODATA2006 default, constants.F90:92-99,182,224)
...
curv = -curv    ! 641  "It is conventional to plot the negative curvature"
```

**The .dat file contains −Ω** (sign flip AFTER unit conversion; both linear so order moot).

### C.3 morb (plot_morb, 370-386)

`berry_get_imfgh_klist` returns the three LVTS12 trace matrices at `fermi_energy_list(1)`
(each `(3 terms J0/J1/J2, 3 axial comps, fermi_n)`): F = imf (Ω), G = img (LC-tilde), H = imh
(IC-tilde). Then:

```
Morb_k = img_k_list(:,:,1) + imh_k_list(:,:,1) − 2.0_dp*fermi_energy_list(1)*imf_k_list(:,:,1)   ! 380-381
Morb_k = −Morb_k/2.0_dp        ! 382  "differs by −1/2 from Eq.97 LVTS12"
my_morb(k, i) = sum(Morb_k(:, i))    ! i = axial component, sum over J0+J1+J2    383-385
```

So morb(k,i) = −½ [ G_i + H_i − 2 E_F F_i ](k), the k-space integrand of the orbital
magnetization (LVTS12 Eq. 97 times −½). Units: **eV·Å² always** — `berry_curv_unit` is NOT
applied to morb, and there is NO sign flip beyond the −½. (The python-script y-label
"Ry·Å²" at kpath.F90:802-803 is a known cosmetic lie; stdout correctly says
`'* Orbital magnetization k-space integrand in eV.Ang^2'`, 1200-1201.)

### C.4 shc (plot_shc, 404-415)

```
berry_get_shc_klist(..., shc_k_fermi=shc_k_fermi)     ! fermi mode
my_shc(k) = shc_k_fermi(1)                            ! 414
```

`shc_k_fermi(1) = Σ_n occ_n(E_F) · Ω^{spin γ}_{n,αβ}(k)` (berry.F90:2894-2897) with the same
band formula as C.1-shc and step occupations at `fermi_energy_list(1)`. Units Å²; then
`if (berry_curv_unit=='bohr2') shc = shc/bohr**2` (479-481). **No sign flip.** This is the
raw k-resolved curvature-like term — none of the −e²/ħ/V, ħ/2/e, or 1e8 S/cm factors that
the berry module applies to `-shc-fermiscan.dat` (berry.F90:1690-1695).

Same conversion for the bands colour column: `if (bands_colour=='shc' .and.
curv_unit=='bohr2') color = color/bohr**2` (476-478).

### C.5 Allowed combinations observed

Any substring combination works; the tests use `bands+morb+curv` and `bands`+colour shc.
When both morb and curv are requested, `berry_get_imfgh_klist` is called once and imf reused
for curv (389 guard). `plot_bands .and. plot_shc` produces the extra `-bands+shc.py` script;
`plot_bands .and. (plot_curv .or. plot_morb)` produces `-bands+curv_*.py` / `-bands+morb_*.py`
(curv wins the filename if both, 1027-1031).

---

## D. Output files (root only, 456-1131)

stdout lists them under `'Output files:'` (497). Data formats below are the complete records —
none of the .dat files has a header line.

### D.1 seedname-path.kpt (only when plot_bands; 466-472)

```
write(dataunit,*) total_pts                                  ! list-directed integer, first line
write(dataunit,'(3f12.6,3x,f4.1)') plot_kpoint(1:3,k), 1.0   ! one line per k; 1.0 = pwscf weight
```

### D.2 seedname-bands.dat (plot_bands; 505-519)

Loop **band-major**: `do i = 1, num_wann; do k = 1, total_pts`:

```
colour none:   write(dataunit,'(2E16.8)') xval(k), eig(i,k)
colour spin/shc: write(dataunit,'(3E16.8)') xval(k), eig(i,k), color(i,k)
```

After each band: `write(dataunit,*) ' '` (517) — a whitespace-only separator line (list-directed
blank + `' '`). Energies in eV; colour = ⟨σ·n̂⟩ (clamped) or Ω^{spinγ}_{n,αβ}(k) (Å²/bohr²).

### D.3 seedname-curv.dat (plot_curv; 642-650)

```
write(dataunit,'(4E16.8)') xval(k), curv(k,1), curv(k,2), curv(k,3)   ! per k
write(dataunit,*) ' '                                                 ! ONE trailing blank line
```

Columns: x, −Ω_x, −Ω_y, −Ω_z (axial: Ω_x=Ω_yz, Ω_y=Ω_zx, Ω_z=Ω_xy), occupied-sum at E_F,
in Å² (or bohr² if `berry_curv_unit=bohr2`).

### D.4 seedname-morb.dat (plot_morb; 731-738)

```
write(dataunit,'(4E16.8)') xval(k), morb(k,1), morb(k,2), morb(k,3)
write(dataunit,*) ' '
```

Columns: x, M^orb integrand components (eV·Å², see C.3).

### D.5 seedname-shc.dat (plot_shc; 813-821)

```
write(dataunit,'(2E16.8)') xval(k), shc(k)
write(dataunit,*) ' '
```

### D.6 Plot scripts (safe to SKIP for test parity — never benchmark-checked)

Written only in the stated combinations:
- `-bands.gnu` + `-bands.py` — only when bands alone (`.not.(curv|morb|shc)`, 522-637)
- `-curv_{x,y,z}.gnu` / `-curv_{x,y,z}.py` — curv without bands (653-728); `achar(119+i)` = x/y/z
- `-morb_{x,y,z}.gnu` / `-morb_{x,y,z}.py` — morb without bands (741-810)
- `-shc.gnu` / `-shc.py` — shc without bands (824-904)
- `-bands+shc.py` — bands AND shc (906-1020)
- `-bands+curv_{x,y,z}.py` or `-bands+morb_{x,y,z}.py` — bands AND (curv|morb) (1022-1129)

Gnuplot format statements (1133-1142), quoted for completeness:

```
701 format('set style data dots',/,'unset key',/,'set xrange [0:',F8.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
702 format('set xtics (',:20('"',A3,'" ',F8.5,','))
703 format(A3,'" ',F8.5,')')
705 format('set arrow from ',F16.8,',',F16.8,' to ',F16.8,',',F16.8,' nohead')
706 format('unset key',/,'set xrange [0:',F9.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
707 format('set style data lines',/,'set nokey',/,'set xrange [0:',F8.5,']',/,'set yrange [',F16.8,' :',F16.8,']')
```

Tick labels `glabel` are `character(len=3)` (159): `glabel(1)=' '//labels(1)//' '` etc.
(485-493) — since labels are len=20, the len-3 truncation keeps a leading space + first two
label chars, and the `'label_end/label_start'` discontinuity form (488) degenerates to the
first 3 chars of the end label. ` G ` is special-cased to `$\Gamma$` in the python scripts.
Tick positions = exact `sum(kpath_len(1:j))`; y-ranges = data min/max ∓ 1 eV (bands) or
∓0.02·range (others); the shc scripts apply a signed-log10 transform for |z|>10.

---

## E. Degeneracies, Fermi energy, smearing

- **Fermi level enters through T=0 step occupations only** (`pw90common_get_occ`,
  postw90_common.F90:964-972: `occ=1` iff `eig < ef`, no smearing, strict `<`). Used by
  imf/imfgh (curv, morb: f_list/g_list projectors, berry-ahc.md §B-C) and by the shc fermi
  sum. Exactly one Fermi energy required for curv/morb/shc (§A errors); `fermi_energy_min/max`
  scans are rejected for shc.
- **shc degeneracy regularization**: the energy denominator is
  `rfac = −2/((e_m−e_n)² + η²)` with η = `kubo_smr_fixed_en_width` (eV). Default η is 0.0 —
  degenerate/near-degenerate pairs then blow up; the Pt tests set 1 eV. Adaptive smearing is
  a hard error in kpath (kpath.F90:181-188) because Δk comes from `berry_kmesh` (berry.F90
  comment 2844-2845).
- **shc band cutoff**: pairs with either eigenvalue above `kubo_eigval_max` are skipped
  (berry.F90:2867) — with defaults this is `dis_froz_max + 0.6667` (Pt: 30+0.6667).
- curv/morb have no eigenvalue cutoff and no smearing; exact degeneracies at E_F simply
  produce large J1/J2 values (adaptive-kmesh refinement exists only in the berry module,
  not in kpath).
- `scissors_shift`/`num_valence_bands` pass through to `get_HH_R` and `get_SHC_R` as usual.
- `use_ws_distance` affects every interpolation through `operator_wigner_setup` /
  `pw90common_fourier_R_to_k` (berry-ahc.md §A.1-A.2). postw90 default is `.true.`; the two
  Fe kpath tests set it `.false.` explicitly (the `_ws` variant sets `true`); the Pt tests
  leave the default (**true**).

---

## F. Reference tests (test-suite/tests/)

All four run `postw90.x` on a pre-converged formatted checkpoint
(`Makefile`: `bunzip2 *.chk.fmt.bz2` → `w90chk2chk.x -f2u`, plus bunzip of .mmn/.uHu/.spn).
Only the single `output=` file is compared (benchmark file `benchmark.out.default.inp=*.win`
IS that .dat file, no header).

### F.1 testpostw90_fe_kpathcurv — checks `Fe-curv.dat`

`Fe.win` kpath-relevant settings (identical file to fe_kpathmorbcurv):

```
num_bands = 28, num_wann = 18, use_ws_distance = .false., search_shells=12
dis_win_min=-8.0 dis_win_max=70.0 dis_froz_min=-8.0 dis_froz_max=30.0
spinors = true, fermi_energy = 12.6279, uHu_formatted = .true.
kpath = true
kpath_task = bands+morb+curv
kpath_num_points=10
begin kpoint_path
G 0.0000 0.0000 0.0000   H 0.500 -0.5000 -0.5000
H 0.500 -0.5000 -0.5000  P 0.7500 0.2500 -0.2500
end kpoint_path
mp_grid = 2 2 2 (8 explicit kpoints); bcc cell 2.71175 bohr
(berry = true is COMMENTED OUT; berry_task/berry_kmesh present but inert)
```

n₁=10, n₂=nint(10·len₂/len₁)=9, total_pts=20; len₁=2.1892688, len₁+len₂=4.0852312 (see §B).
jobconfig (line 328-332): `program = POSTW90_CURVDAT_OK, output = Fe-curv.dat`. userconfig
(101-108): parser `parse_curv_dat` (fields from 4-token rows: bandpath, bandcurvx/y/z),

```
tolerance = ( (1.0e-6, 5.0e-6, 'bandpath'), (1.0e-6, 1.0e+3, 'bandcurvx'),
              (1.0e-6, 1.0e+3, 'bandcurvy'), (1.0e-6, 1.0e+3, 'bandcurvz'))
```

Benchmark first rows (`0.21892688E+00` spacing; curv_z hits −5.6 at k₃):

```
  0.00000000E+00 -0.47753935E-19  0.52889654E-20  0.47482879E-21
  0.21892688E+00 -0.60208601E-06 -0.17632139E-06  0.34409405E-03
  0.43785377E+00 -0.10052614E-04  0.32238014E-04 -0.56272413E+01
```

### F.2 testpostw90_fe_kpathmorbcurv — checks `Fe-morb.dat`

Same Fe.win as F.1 (diff confirms identical); jobconfig 316-320:
`program = POSTW90_MORBDAT_OK, output = Fe-morb.dat`. `testpostw90_fe_kpathmorbcurv_ws/`
is the same but with `use_ws_distance = true` (jobconfig 322-326). Tolerances (userconfig
93-99): `((1.0e-6,5.0e-6,'bandpath'),(1.0e-4,1.0e-4,'bandmorbx'),(...y),(...z))`.
Benchmark first rows:

```
  0.00000000E+00  0.71495070E-07  0.23833843E-06  0.32522033E+00
  0.21892688E+00  0.50309754E-06  0.97327112E-06  0.35041958E+00
```

Note both Fe tests also write Fe-bands.dat, Fe-path.kpt, Fe-curv.dat, Fe-morb.dat and the
bands+curv/bands+morb python scripts — only the named .dat is compared.

### F.3 testpostw90_pt_kpathshc — checks `Pt-shc.dat`

`Pt.win` kpath-relevant settings:

```
shc_freq_scan = false, shc_alpha = 1, shc_beta = 2, shc_gamma = 3
spn_formatted = true, shc_method = qiao
kpath = true
kpath_task = shc                      (#kpath_bands_colour = shc commented)
kpath_num_points = 10
kubo_adpt_smr = false
kubo_smr_fixed_en_width = 1
fermi_energy = 17.9919
begin kpoint_path
W  0.75  0.50  0.25    L  0.50  0.00  0.00
L  0.50  0.00  0.00    G  0.00  0.00  0.00
G  0.00  0.00  0.00    X  0.50  0.50  0.00
X  0.50  0.50  0.00    W  0.75  0.50  0.25
W  0.75  0.50  0.25    G  0.00  0.00  0.00
end kpoint_path
berry_curv_unit = ang2, num_bands = 40, num_wann = 18, spinors = true
dis_win 0..60, dis_froz 0..30 (→ kubo_eigval_max = 30.6667), mp_grid = 4 4 4
use_ws_distance NOT set → default .true.
(kpath needs Pt.chk, Pt.eig, Pt.mmn, Pt.spn — spn formatted)
```

5 segments, n₁=10, total_pts=60; first-segment spacing 0.11333885, final x=6.7178249.
jobconfig 503-506: `program = POSTW90_SHCKPATHDAT_OK, output = Pt-shc.dat`. userconfig
221-225: parser `parse_shc_kpath_dat` (2-token rows → path, shc),
`tolerance = ((1.0e-6, 5.0e-6, 'path'), (1.0e-1, None, 'shc'))`.
Benchmark first rows:

```
  0.00000000E+00  0.73385296E-01
  0.11333885E+00  0.10001571E+00
```

### F.4 testpostw90_pt_kpathbandsshc — checks `Pt-bands.dat`

Same Pt.win as F.3 except:

```
kpath_task = bands
kpath_bands_colour = shc
```

jobconfig 497-500: `program = POSTW90_SHCKPATHBANDSDAT_OK, output = Pt-bands.dat`.
userconfig 214-219: parser `parse_shc_kpath_bandsdat` (3-token rows → path, energy, shc;
blank separator lines skipped),
`tolerance = ((1.0e-6,5.0e-6,'path'),(1.0e-3,2.0e-3,'energy'),(1.0e-1,None,'shc'))`.
Benchmark first rows (x, E [eV], Ω^{spin z}_{xy,band} [Å²]):

```
  0.00000000E+00  0.12074410E+02 -0.38012782E+00
  0.11333885E+00  0.12061117E+02 -0.59918484E+00
```

18 bands × 60 points + blank separators ≈ 1098 lines.

---

## G. Trap checklist for the Julia port

1. **xval accumulation quirk** (§B): per-point step uses the CURRENT segment's spacing even
   at segment boundaries; final point forced to `sum(kpath_len)`. Do not "fix" this — the
   benchmark x column encodes it (path tolerance is 5e-6 relative).
2. `nint` rounding for per-segment point counts (round half away from zero).
3. Segment lengths via the metric `recip_lattice·recip_latticeᵀ` in Å⁻¹ with 2π included.
4. curv .dat stores **−Ω**; shc and morb are NOT sign-flipped; morb has the extra −½ and
   −2·E_F·imf term; morb ignores `berry_curv_unit`.
5. bohr² conversion divides by `bohr**2` with bohr = 0.52917720859 Å (CODATA2006 — the
   compile-time DEFAULT, constants.F90:92-99; CODATA2018/2022 differ in the 10th digit).
6. Step occupations with strict `eig < ef`; exactly one `fermi_energy`.
7. shc: fixed smearing mandatory (`kubo_adpt_smr=false`), η in eV, `kubo_eigval_max =
   dis_froz_max + 0.6667` default cutoff, qiao operators only (SR/SHR/SH from .spn+.mmn,
   no .uHu/.sHu/.sIu in kpath), `shc_alpha/beta/gamma` defaults 1/2/3.
8. spin colour clamped to ±(1−1e-8); shc colour not clamped.
9. `use_ws_distance` defaults **true** in this code base; Fe kpath tests override to false,
   Pt tests do not.
10. -bands.dat is band-major with a whitespace separator line after each band; curv/morb/shc
    .dat have exactly one trailing separator line; -path.kpt starts with a bare list-directed
    integer count and uses `(3f12.6,3x,f4.1)` rows with weight 1.0.
11. All .dat numeric records use `E16.8` (e.g. ` 0.21892688E+00`, `-0.56272413E+01`).
12. .gnu/.py scripts are cosmetic — not compared by the test suite; skip or stub them.
13. Fermi-energy list parsing: `fermi_energy` gives fermi_n=1; a min/max scan gives fermi_n>1
    which kpath rejects for curv/morb/shc.
