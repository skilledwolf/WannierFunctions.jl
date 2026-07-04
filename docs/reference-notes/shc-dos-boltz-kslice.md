# postw90 reference notes: spin/SHC, DOS, BoltzWann, kslice

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Shared machinery is **already documented** and only referenced here:
- `docs/reference-notes/berry-ahc.md`: get_HH_R, get_AA_R (mmn reader, `v_matrix` windowing,
  `get_gauge_overlap_matrix`, `fourier_q_to_R`, `operator_wigner_setup`), wan_ham k-space
  machinery (`wham_get_eig_deleig`, `wham_get_D_h`, degenerate perturbation), JJ± matrices,
  `berry_get_imf(gh)_klist` (J0/J1/J2), AHC assembly, `adkpt` adaptive sub-mesh construction.
- `docs/reference-notes/kubo-morb-geninterp.md`: adaptive smearing / `utility_w0gauss`
  (smearing type indices), `pw90common_kmesh_spacing`, `pw90common_get_occ` occupations,
  module kmesh resolution (`<mod>_kmesh` / `<mod>_kmesh_spacing` falling back to global
  `kmesh`), Fermi-scan basics, testcode tolerance semantics: tuples are
  `(abs_tol, rel_tol, name)`, strict ⇒ **both** must pass
  (`test-suite/testcode/lib/testcode2/validation.py:114-158`).

Key papers cited (`src/postw90/berry.F90:43-44`): QZYZ18 = Qiao-Zhou-Yuan-Zhao, PRB 98,
214402 (2018); RPS19 = Ryoo-Park-Souza, PRB 99, 235113 (2019).

Physical constants (`src/constants.F90:108-124`, 2018 CODATA branch):
`elem_charge_SI = 1.602176634e-19` C, `hbar_SI = 1.054571817e-34` J·s,
`k_B_SI = 1.380649e-23` J/K, `eV_seconds = 6.582119569e-16`,
`smearing_cutoff = 10._dp` (line 80), `min_smearing_binwidth_ratio = 2._dp` (line 82).

---

## 1. SPIN operator and Spin Hall conductivity

### 1.1 The .spn file (reader in get_SS_R, get_oper.F90:1727-1798)

Selected by `spn_formatted` (default `.false.`, `postw90_types.F90:68`; keyword read at
`postw90_readwrite.F90:477-478`). Content = ⟨ψ_nk|σ_i|ψ_mk⟩ between ab-initio eigenstates
(σ = Pauli matrices), **upper triangle only** (n ≤ m), 3 spin components x,y,z per pair.

Formatted variant (`get_oper.F90:1730-1737, 1755-1771`):
```
line 1:  header                          (read (spn_in,*) header; character(len=60))
line 2:  nb_tmp, nkp_tmp                 (must equal num_bands, num_kpts; 1747-1754)
then, loop ik=1,num_kpts { loop m=1,num_bands { loop n=1,m {
   line: s_real s_img    -> spn_o(n,m,ik,1)   (sigma_x)
   line: s_real s_img    -> spn_o(n,m,ik,2)   (sigma_y)
   line: s_real s_img    -> spn_o(n,m,ik,3)   (sigma_z)
}}}
```
i.e. per (m,n) pair three separate "re im" lines; loop order: n fastest (1..m), then m,
then ik. Lower triangle rebuilt as `spn_o(m,n,ik,i) = conjg(spn_o(n,m,ik,i))` (1766-1768).

Unformatted variant (`get_oper.F90:1739-1745, 1772-1798`): record 1 = `header`
(character(len=60)); record 2 = `nb_tmp, nkp_tmp` (integers); then **one record per
k-point**:
```
read (spn_in) ((spn_temp(s, m), s=1, 3), m=1, (num_bands*(num_bands + 1))/2)   ! line 1779
```
with `spn_temp` `complex(kind=dp)` of shape `(3, nb(nb+1)/2)` — i.e. for each packed
upper-triangular pair (counter runs n=1..m inner, m outer; 1780-1791) the 3 spin components
are contiguous (spin index fastest). Unpacked identically with conjugate for (m,n).

### 1.2 get_SS_R (get_oper.F90:1652-1836)

`SS_R(num_wann, num_wann, nrpts_pw90, 3)` = ⟨0n|σ_x,y,z|Rm⟩. Steps (root node only, then
broadcast at 1823):
1. If `SS_R` already allocated, return ("been here before", 1707-1711).
2. `num_states(ik) = dis_manifold%ndimwin(ik)` if disentangled else `num_wann` (1719-1725).
3. Read spn_o from `.spn` (above).
4. Gauge transform to Wannier gauge with the same windowed rotation used for HH_R:
   `get_gauge_overlap_matrix(..., spn_o(:,:,ik,is), have_disentangled, SS_q(:,:,ik,is))`
   for each ik, is (1804-1812) — effectively SS_q = V† spn_o(window) V (see berry-ahc.md §A).
5. `fourier_q_to_R` per spin component (1814-1816), then `operator_wigner_setup` per
   component (1818-1820) (ws_distance/pw90 R-set mapping, see berry-ahc.md).
No `transl_inv`/phase correction is applied to SS_R (unlike AA_R).

### 1.3 spin_get_nk (spin.F90:223-305)

Computes ⟨ψ_mk^(H)|σ·n̂|ψ_mk^(H)⟩ for m=1..num_wann:
- Interpolate H(k) (`pw90common_fourier_R_to_k` with deriv flag 0), diagonalize → UU, eig
  (279-283).
- Interpolate each SS(:,:,is) from SS_R (285-290).
- Quantization axis from `spin_axis_polar`, `spin_axis_azimuth` (degrees; defaults 0.0, 0.0,
  `postw90_types.F90:84-85`; keywords at `postw90_readwrite.F90:679-686`):
```
alpha(1) = sin(polar/conv)*cos(azimuth/conv); alpha(2) = sin(polar/conv)*sin(azimuth/conv)
alpha(3) = cos(polar/conv)                      ! conv = 180/pi   (294-297)
SS_n = alpha(1)*SS_x + alpha(2)*SS_y + alpha(3)*SS_z              ! (301)
spn_nk(:) = real(utility_rotate_diag(SS_n, UU, num_wann), dp)     ! (303)  = diag(UU^† SS_n UU)
```
Default axis = ẑ, so `spn_nk` = ⟨σ_z⟩_n. Values span [-1,1] (Pauli σ, not S; spin.F90:199-201).
`spin_get_moment` sums `spn_mom = -Σ_k Σ_n occ_n(E_F) spn_nk / N_k` (spin.F90:174-201);
requires nfermi=1 (122-125). `spin_decomp` keyword requires `num_elec_per_state == 1`
(`postw90_readwrite.F90:690-693`).

### 1.4 SHC input parameters

`pw90_spin_hall_type` defaults (`postw90_types.F90:198-210`): `freq_scan = .false.`
(→ Fermi-scan mode is the default), `alpha = 1`, `beta = 2`, `gamma = 3` (i.e. σ^{spin z}_{xy}),
`bandshift = .false.`, `method = ' '`. Keyword parsing `postw90_readwrite.F90:1047-1124`:
`shc_freq_scan`, `shc_alpha/beta/gamma` (each must be 1,2,3), `shc_bandshift(+_firstband,
+_energyshift)` (mutually exclusive with `scissors_shift`, 1080-1083), and **`shc_method` is
mandatory when `berry_task` contains `shc`** ("Error: berry_task=shc and shc_method is not
set", 1108-1111); value must contain 'qiao' or 'ryoo' (1112-1116). So there is no default
method — the pt_shc test sets `shc_method = qiao` explicitly. `transl_inv_full=T` is rejected
with qiao (1117-1124). Activation: `if (index(pw90_berry%task, 'shc') > 0) eval_shc=.true.`
(berry.F90:286). Fermi scan allowed (multiple fermi energies); freq_scan requires nfermi=1:
`not_scannable = eval_kubo .or. (eval_shc .and. pw90_spin_hall%freq_scan)` (berry.F90:342-347).

Fermi list construction (`src/readwrite.F90:648-691`): `fermi_energy_min/max/step`;
`n = nint(abs((max-min)/step)) + 1`; `step` recomputed as `(max-min)/(n-1)`;
`fermi_energy_list(i) = min + (i-1)*step`. Default step 0.01 if omitted.

### 1.5 Operators required (berry_main, berry.F90:413-454)

For eval_shc: `get_HH_R`, `get_AA_R` (reads `.mmn`), `get_SS_R` (reads `.spn`); then
- **qiao**: `get_SHC_R` (berry.F90:435-438) → SR_R, SHR_R, SH_R.
- **ryoo**: `get_SH_R` (reads `.spn`), `get_SAA_R` (reads **`.sIu`**, get_oper.F90:2932-2935),
  `get_SBB_R` (reads **`.sHu`**, get_oper.F90:2633-2636).

`get_SHC_R` (get_oper.F90:1839-2259) builds, from `.spn` + `.mmn` only (no uHu!):
```
SR_R  = <0n|sigma_{x,y,z}.(r-R)_alpha|Rm>       (1846)
SHR_R = <0n|sigma_{x,y,z}.H.(r-R)_alpha|Rm>     (1847)
SH_R  = <0n|sigma_{x,y,z}.H|Rm>                 (1848)
```
- H_o = diag(eigval) with scissors_shift on bands > num_valence_bands, or shc_bandshift on
  bands >= bandshift_firstband (2045-2063).
- SH_o(:,:,ik,is) = spn_o·H_o (QZYZ18 Eq.48, 2109); gauge-projected → SH_q (2111-2113).
- From `.mmn` overlaps S_o=⟨u_ik|u_ik2⟩: SM_o = spn_o·S_o (Eq.50), SHM_o = SH_o·S_o (Eq.51)
  (2170-2172); gauge-projected with bra at ik and ket at neighbour
  (`get_gauge_overlap_matrix(..., ik, ..., nnlist(ik,nn), ...)`, 2181-2189); finite-difference
  sums with `wb·bk`:
```
SR_q(:,:,ik,is,idir)  += wb(nn)*bk(idir,nn,ik)*(SM_q(:,:,is) - SS_q(:,:,is))     (2195-2197)
SHR_q(:,:,ik,is,idir) += wb(nn)*bk(idir,nn,ik)*(SHM_q(:,:,is) - SH_q(:,:,ik,is)) (2199-2201)
```
- fourier_q_to_R + operator_wigner_setup per (is, idir) (2210-2236), then
  `SR_R = cmplx_i*SR_R; SHR_R = cmplx_i*SHR_R` (2238-2239).

### 1.6 berry_get_shc_klist (berry.F90:2684-3099)

Returns the "Berry curvature-like term" of QZYZ18 Eqs.(3),(4): σ^{γ}_{αβ}(k) contribution
per k, **in Å²** (2698-2699). ħ factors: js_k lacks ħ/2 (spin op) and 1/ħ (velocity), the 2nd
velocity lacks 1/ħ — all cancelled by ħ² in Eq.(3) except an overall ħ/2/(-e) handled in the
final unit factor (2701-2705, 3046-3052).

Per k-point:
1. `wham_get_eig_deleig` → eig, del_eig, delHH, UU (2814-2820); `wham_get_D_h` → D_h (2822).
2. Optional `shc_bandshift`: `eig(firstband:) += energyshift` (2825-2827).
3. AA = rotate(FT AA_R) + i·D_h ("Eq.(25) WYSV06", 2829-2837).
4. `berry_get_js_k` (2917-3097) computes js_k(n,m) = ⟨ψ|½(σ_γ v_α + v_α σ_γ)|ψ⟩ (QZYZ18 Eq.23):
```
S_k  = rotate(FT SS_R(:,:, :,gamma))                                   ! Eqs.(25),(36),(30)  2994-3000
(qiao branch, 3002-3052):
SR_alpha_k  = -i * rotate(FT SR_R(:,:,:,gamma,:))(alpha)               ! Eq.(31)  3006-3013
K_k  = SR_alpha_k + matmul(S_k, D_alpha_h)                             ! Eq.(26)  3014
SHR_alpha_k = -i * rotate(FT SHR_R(:,:,:,gamma,:))(alpha)              ! Eq.(32)  3019-3025
SH_k = rotate(FT SH_R(:,:,:,gamma))                                    ! Eq.(32)  3027-3033
L_k  = SHR_alpha_k + matmul(SH_k, D_alpha_h)                           ! Eq.(27)  3034
B_k  = del_eig_mat*S_k + eig_mat*K_k - L_k     ! element-wise *, eig_mat(i,:)=eig(:)  3038-3044
js_k = 1/2*(B_k + conjg(transpose(B_k)))                               ! Eq.(23)  3052
(ryoo branch, 3054-3095): spinVel0 = VV0·S_k + S_k·VV0 with VV0 = rotate(delHH_alpha);
js_k(n,m) = spinVel0(n,m) - i(eig(m)·SAA(n,m) - SBB(n,m)) + i(eig(n)·conjg(SAA(m,n)) -
conjg(SBB(m,n))); js_k = js_k/2                       ! RPS19 Eqs.(21),(26)
```
5. Band sums (2858-2903), for each n: Ω^γ_{n,αβ}(k) accumulated over m ≠ n, skipping pairs
   with `eig(m) > kubo_eigval_max .or. eig(n) > kubo_eigval_max` (2867):
```
prod = js_k(n, m) * cmplx_i * (eig(m)-eig(n)) * AA(m, n, shc_beta)     ! 2869-2872
fermi/band mode: rfac = -2/((eig(m)-eig(n))**2 + eta_smr**2);  omega += rfac*aimag(prod)  (2888-2890)
freq mode:  cdum = freq(ifreq) + i*eta_smr; cfac = -2/(rfac**2 - cdum**2)
            omega_list(ifreq) += cfac*aimag(prod)                       (2882-2887)
```
   (the commented line 2871 shows the AHC analogue `prod = -rfac*i*AA(n,m,α) * rfac*i*AA(m,n,β)`.)
6. Smearing η (2873-2881): if `kubo_smearing%use_adaptive` (default: inherits global
   `adpt_smr`, default true; `postw90_readwrite.F90:908-911`):
   `eta_smr = min(|del_eig(m,:)-del_eig(n,:)|·Δk·adaptive_prefactor, adaptive_max_width)`
   ("Eq.(35) YWVS07"), with `Δk = pw90common_kmesh_spacing(berry_kmesh, recip_lattice)`
   (2846-2849); else `eta_smr = kubo_smearing%fixed_width`. Note comment 2844-2845: adaptive
   only works with `berry_kmesh` — do **not** use adaptive in kpath/kslice.
7. Occupations: fermi mode weights `shc_k_fermi(i) += occ_fermi(n,i)*omega` with
   `pw90common_get_occ(fermi_energy_list(i), eig, ...)` per Fermi energy (2852-2856,
   2894-2897); freq mode uses occ at `fermi_energy_list(1)` (2850-2851, 2898-2899); band mode
   `shc_k_band(n) = omega` (2900-2901). `kubo_eigval_max` default: `froz_max + 0.6667` if
   frozen states, else `maxval(eigval)+0.6667` (postw90_readwrite.F90:1758-1768).

### 1.7 berry_main SHC accumulation, adaptive kmesh, units, output

Fermi-scan loop (regular grid branch berry.F90:968-1046; kpoint-file branch 760-836):
`kweight = db1*db2*db3` (=1/N_k), `kweight_adpt = kweight/curv_adpt_kmesh**3` (843-844).
After the first `berry_get_shc_klist` call, adaptive refinement is triggered when
`berry_curv_adpt_kmesh > 1` and `abs(shc_k_fermi(if))` (converted `/bohr**2` if
`berry_curv_unit='bohr2'`) exceeds `berry_curv_adpt_kmesh_thresh` **at any Fermi energy**
(788-803); then the k-point is replaced by `curv_adpt_kmesh**3` points `kpt + adkpt(:,i)`
accumulated with `kweight_adpt` (804-820); else `shc_fermi += kweight*shc_k_fermi` (822).
Defaults `curv_adpt_kmesh = 1`, thresh `= 100.0` (postw90_types.F90:166-167) — so **no
refinement by default**. Freq-scan branch has no adaptive kmesh (824-836).
After `comms_reduce` (1228-1234):

```
! (i) multiply -e^2/hbar/(V*N_k) as in QZYZ18 Eq.(5) (1/N_k already in kweight)
! (ii) spin current: overall -hbar/2/e ;  (iii) 1e8 -> S/cm
fac = 1.0e8_dp*physics%elem_charge_SI**2/(physics%hbar_SI*cell_volume)/2.0_dp   ! berry.F90:1690
shc_fermi = shc_fermi*fac   (or shc_freq)                                       ! 1691-1695
```
Final unit: **(ħ/e)·S/cm** (1677-1688). `cell_volume` in Å³.

Output files (berry.F90:1704-1727):
- Fermi scan → `seedname-shc-fermiscan.dat`; header
  `write(...,'(a,3x,a,3x,a)') '#No.', 'Fermi energy(eV)', 'SHC((hbar/e)*S/cm)'` (1713-1714);
  rows `'(I4,1x,F12.6,1x,E17.8)') n, fermi_energy_list(n), shc_fermi(n)` (1716-1717).
- Freq scan → `seedname-shc-freqscan.dat`; header
  `'#No.', 'Frequency(eV)', 'Re(sigma)...', 'Im(sigma)...'`; rows
  `'(I4,1x,F12.6,1x,1x,2(E17.8,1x))'` (1720-1725).
stdout also prints "Qiao's SHC (Phys.Rev.B 98.214402)" / "Ryoo's ..." (541-546).

### 1.8 testpostw90_pt_shc

jobconfig (`test-suite/tests/jobconfig:472-476`):
```
[testpostw90_pt_shc/]
program = POSTW90_SHCFERMIDAT_OK
inputs_args = ('Pt.win', '')
output = Pt-shc-fermiscan.dat
```
So the harness compares the **`.dat` file**, not the `.wpout`. Parser
`test-suite/tools/parsers/parse_shc_dat.py`: line 0 must have `pieces[1] == 'Fermi'`
(fermiscan) or `'Frequency(eV)'` (freqscan); fermiscan rows must have exactly 3 tokens →
`energy = float(pieces[1])`, `shc = float(pieces[2])`; blank lines skipped. Tolerances
(`test-suite/tests/userconfig:201-205`):
```
[POSTW90_SHCFERMIDAT_OK]
extract_fn = tools parsers.parse_shc_dat.parse
tolerance = ( (1.0e-6, 5.0e-6, 'energy'),
              (1.0e-1, 1.0e-1, 'shc'))        # (abs, rel), strict: both must pass
```
Input files shipped: `Pt.win`, `Pt.spn.bz2`, `Pt.chk.fmt.bz2`, `Pt.mmn.bz2`, `Pt.amn.bz2`,
`Pt.eig`. **No uHu** — the default-method test uses qiao, which needs only .spn + .mmn
(+ chk/eig). Pt.win key settings: `shc_freq_scan = false`, `shc_alpha/beta/gamma = 1/2/3`,
`spn_formatted = true`, `shc_method = qiao`, `berry = true`, `berry_task = eval_shc`,
`berry_kmesh = 15` (→15×15×15), `berry_curv_unit = ang2`, `fermi_energy_min = 6`,
`fermi_energy_max = 26`, `fermi_energy_step = 0.1` (→ 201 Fermi energies),
`num_bands = 40`, `num_wann = 18`, spinors, FCC Pt, mp_grid 4 4 4. Adaptive-kmesh keywords
commented out (default 1 = off); kubo smearing keywords commented out → adaptive smearing
with defaults (prefactor √2, max 1.0 eV, Gaussian), kubo_eigval_max default = froz_max
(30.0)+0.6667. Benchmark `benchmark.out.default.inp=Pt.win` = the full 201-row fermiscan
file; e.g. row 121 `18.000000 0.14194102E+04`, row 1-17 exactly 0.0 (bands below window).

---

## 2. DOS module (dos.F90)

### 2.1 Parameters and defaults

- `dos_task` default `'dos_plot'` (postw90_types.F90:147; other value: `find_fermi_energy`,
  postw90_readwrite.F90:1241-1254; find_fermi_energy is dead code, dos.F90:372-557 commented).
- `dos_energy_step` default `0.01` eV (postw90_types.F90:151).
- `dos_energy_max` default: `froz_max + 0.6667` if frozen states, else `maxval(eigval)+0.6667`,
  else `win_max + 0.6667` (postw90_readwrite.F90:1664-1673). `dos_energy_min` default:
  `minval(eigval) - 0.6667` else `win_min - 0.6667` (1675-1682).
- Smearing: `pw90_dos%smearing` inherits the global smearing then per-key overrides
  `dos_adpt_smr`, `dos_adpt_smr_fac`, `dos_adpt_smr_max`, `dos_smr_fixed_en_width`,
  `dos_smr_type` (postw90_readwrite.F90:1267-1349). Global defaults
  (postw90_types.F90:129-141): `use_adaptive = .true.`, `adaptive_prefactor = sqrt(2)`,
  `type_index = 0` (Gaussian), `fixed_width = 0`, `adaptive_max_width = 1.0` eV.
- `dos_project` range vector; default = all WFs 1..num_wann (1305-1339).
- k-mesh: `dos_kmesh`/`dos_kmesh_spacing` else global `kmesh` (module kmesh mechanism, see
  kubo-morb-geninterp.md).

### 2.2 dos_main (dos.F90:60-368)

Energy grid (146-158):
```
num_freq = nint((energy_max - energy_min)/energy_step) + 1;  if (num_freq==1) num_freq=2
d_omega  = (energy_max - energy_min)/(num_freq - 1)
dos_energyarray(ifreq) = energy_min + (ifreq-1)*d_omega
```
Operators: `get_HH_R` always; `get_SS_R` iff `spin_decomp` (then ndim=3 else ndim=1)
(176-193). k loop (full-BZ branch 285-330): same flattened loop as spin/berry;
`kweight = 1/PRODUCT(dos_kmesh)`; per point:
- adaptive: `wham_get_eig_deleig` → del_eig; `dos_get_levelspacing` (768-795):
  `levelspacing(band) = sqrt(dot_product(del_eig(band,:), del_eig(band,:)))*Delta_k` with
  `Delta_k = pw90common_kmesh_spacing(kmesh, recip_lattice)`; then `dos_get_k(...,
  levelspacing_k=..., UU=UU)`.
- fixed: plain Fourier + diagonalize, `dos_get_k(..., UU=UU)` without levelspacing.
`dos_all += dos_k*kweight` (329); `comms_reduce` (336); wanint_kpoint_file branch (237-281)
uses IBZ points/weights from kpoint.dat.

Output (339-348): file `seedname-dos.dat`:
```
write (dos_unit, '(4E16.8)') omega, dos_all(ifreq, :)      ! dos.F90:346
```
→ 2 columns (E, DOS) or 4 columns (E, total, spin-up, spin-down) when spin_decomp. Units:
states/eV/cell, normalized so ∫DOS dE = num electrons (comment 563-566).

### 2.3 dos_get_k (dos.F90:600-765)

Mutually-exclusive checks: levelspacing_k present ⟺ smearing%use_adaptive (653-665).
`spin_decomp` → `spin_get_nk` → per band `alpha_sq = (1+spn_nk(i))/2`, `beta_sq = 1-alpha_sq`
(671-688). Per band i:
```
eta_smr = smearing%fixed_width                                    ! fixed (691)
eta_smr = min(levelspacing_k(i)*adaptive_prefactor, adaptive_max_width)  ! adaptive, Eq.(35) YWVS07 (694)
```
Bin-range optimization (699-715): if `eta_smr/binwidth < min_smearing_binwidth_ratio` (=2):
**no smearing** — only the single bin nearest eig_k(i), weight `rdum = 1/binwidth` (724);
else bins within `± smearing_cutoff*eta_smr` (=10η) and
`rdum = utility_w0gauss((E_f - eig_k(i))/eta_smr, type_index)/eta_smr` (720-721).
Accumulation (730-761): total DOS `dos_k(loop_f,1) += rdum*num_elec_per_state`; spin channels
`+= rdum*alpha_sq` / `rdum*beta_sq` (no num_elec_per_state — spinor calc has 1 e/state,
comment 735-737). Projected DOS (num_project < num_wann): weight `*abs(UU(project(j), i))**2`
per selected WF (749-760).

### 2.4 DOS tests

`testpostw90_example04_dos` (jobconfig:424-428): program POSTW90_DOS_OK,
`inputs_args = ('copper.win','')`, `output = copper-dos.dat`. Ships only
`copper.chk.fmt.bz2`, `copper.eig`, `copper.win` (+Makefile) — DOS needs no mmn/amn/spn.
copper.win: `num_bands=12, num_wann=7, use_ws_distance=.false., search_shells=12,
dis_win_max=38.0, dis_froz_max=13.0, dos=true, kmesh=10` (global → 10×10×10),
`dos_energy_max=10, dos_energy_min=8, dos_energy_step=0.25`, mp_grid 4 4 4. Adaptive smearing
by default. Benchmark: 9 rows × 2 cols, first row `0.80000000E+01 0.19268028E+01`.

`testpostw90_fe_dos_spin` (jobconfig:370-374): program POSTW90_DOS_OK, input Fe.win,
output `Fe-dos.dat`. Ships `Fe.spn` (uncompressed, formatted), `Fe.chk.fmt.bz2`, `Fe.eig`,
`Fe.amn`, `Fe.mmn.bz2`. Fe.win: `num_bands=28, num_wann=18, spinors=true,
fermi_energy=12.6279, spin_moment=true, spn_formatted=true, kmesh=4, dos=true,
spin_decomp=true, dos_energy_max=13.0, dos_energy_min=10.0, dos_energy_step=0.2,
dos_adpt_smr=false, dos_smr_fixed_en_width=0.5`, use_ws_distance=.false., mp_grid 2 2 2.
Benchmark 16 rows × 4 cols, first row
`0.10000000E+02 0.83154571E+00 0.11547299E+00 0.71607272E+00`.

Tolerances (`userconfig:139-145`):
```
[POSTW90_DOS_OK]
extract_fn = tools parsers.parse_dos_dat.parse
tolerance = (  (1.0e-6, 5.0e-6, 'energy'),
               (1.0e-4, 1.0e-4, 'dos'),
               (1.0e-4, 1.0e-4, 'dos_spin1'),
               (1.0e-4, 1.0e-4, 'dos_spin2'))
```
parse_dos_dat.py: skips `#`/blank lines; 2 tokens → energy,dos; 4 tokens → energy, dos,
dos_spin1, dos_spin2; other lengths → error.

---

## 3. BoltzWann (boltzwann.F90)

Reference paper: Pizzi et al., Comp. Phys. Comm. 185, 422 (2014) (boltzwann.F90:50-52).
Packed tensor indices `XX=1, XY=2, YY=3, XZ=4, YZ=5, ZZ=6` — mapping (i,j)→`i+((j-1)*j)/2`
for i≤j (boltzwann.F90:70-75, 850-855).

### 3.1 Parameters and defaults (postw90_types.F90:241-265; postw90_readwrite.F90:1382-1634)

- `boltz_relax_time` default **10.0 fs** (types:261; "By default: 10 fs relaxation time",
  readwrite:1610-1613). Constant, band- and k-independent (boltzwann.F90:1491-1493).
- `boltz_tdf_energy_step` default **0.001 eV** (types:258; "the energy step for the TDF is
  1 meV", readwrite:1571-1578; must be > 0).
- `boltz_tdf_smr_fixed_en_width`: default = global `smr_fixed_en_width` (default 0 → **no
  smearing / plain histogram binning**; readwrite:1580-1589); `boltz_tdf_smr_type` default =
  global. TDF smearing is **never adaptive**: `pw90_boltzwann%TDF_smearing%use_adaptive =
  .false.` (readwrite:1412) and an error if set (boltzwann.F90:268-271).
- (μ,T) grid — all **required** when boltzwann=true (errors otherwise, readwrite:1505-1567):
  `boltz_mu_min`, `boltz_mu_max`, `boltz_mu_step` (eV), `boltz_temp_min/max/step` (K;
  temp_min > 0). Grids: `TempNumPoints = int(floor((temp_max-temp_min)/temp_step)) + 1`,
  `TempArray(i) = temp_min + (i-1)*temp_step` (239-247); same for Mu (258-266).
  `KTArray = TempArray*k_B_SI/elem_charge_SI` (eV; 249-256).
- `boltz_kmesh`/`boltz_kmesh_spacing` else global `kmesh`.
- `boltz_2d_dir` ∈ {no,x,y,z} → dir_num_2d 0..3 (readwrite:1420-1438); 2×2 submatrix inversion
  for Seebeck in 2D (boltzwann.F90:414-453).
- `boltz_calc_also_dos` default `.false.` (types:245); DOS keywords `boltz_dos_energy_step`
  (default 0.001), `boltz_dos_energy_min` default `minval(eigval) - 0.6667`,
  `boltz_dos_energy_max` default `maxval(eigval) + 0.6667` (readwrite:1448-1471); DOS smearing
  `boltz_dos_adpt_smr(_fac,_max)`, `boltz_dos_smr_fixed_en_width`, `boltz_dos_smr_type`
  default to the global smearing (1473-1503, 1601-1608).
- `boltz_bandshift(_firstband,_energyshift)`: rigid shift `eig(firstband:) += energyshift`
  before TDF accumulation (boltzwann.F90:1076-1079).

### 3.2 TDF energy grid (boltzwann.F90:272-290)

```
TDF_exceeding_energy = max(TDF_exceeding_energy_times_smr*tdf_smearing%fixed_width, 0.2_dp)   ! (278); constant = 3._dp (182)
TDFEnergyNumPoints = int(floor((dis_win_max - dis_win_min + 2*TDF_exceeding_energy)/tdf_energy_step)) + 1   ! (279-280)
TDFEnergyArray(i)  = dis_win_min - TDF_exceeding_energy + (i-1)*tdf_energy_step               ! (287-290)
```
(dis window defaults to eigenvalue range when not specified.)

### 3.3 TDF_kpt (boltzwann.F90:1300-1495) and assembly

Velocities: `deleig_k = del_eig` from `wham_get_eig_deleig` (1067-1071) — band derivatives
∇_k ε in eV·Å (degeneracy-corrected via `pw90_band_deriv_degen`; see berry-ahc.md). Binning
identical to dos_get_k (fixed smearing only, `smear = tdf_smearing%fixed_width`; unsmeared
single bin with `rdum = 1/binwidth` when `smear/binwidth < min_smearing_binwidth_ratio`,
1413-1441):
```
TDF_k(XX,loop_f,1) += rdum * num_elec_per_state * deleig_k(BandIdx,1)*deleig_k(BandIdx,1)   ! 1445-1446
... (XY: v1*v2, YY: v2*v2, XZ: v1*v3, YZ: v2*v3, ZZ: v3*v3)                                 ! 1447-1456
(spin_decomp: channels 2/3 weighted alpha_sq/beta_sq via spin_get_nk, no num_elec_per_state) ! 1460-1487
TDF_k = TDF_k * pw90_boltzwann%relax_time                                                   ! 1493
```
Assembly: `TDF = Σ_k TDF_k * kweight / cell_volume` with `kweight = 1/PRODUCT(boltz_kmesh)`
(1039, 1086-1090), then `comms_allreduce` (1176). So
**TDF_ij(ε) = (1/(V_cell·N_k)) Σ_{n,k} v_i v_j τ δ(ε−ε_nk) · (spin degeneracy)**, stored as
ħ²·TDF in units **eV·fs/Å** ("The TDF array contains now the TDF, or more precisely hbar^2 *
TDF in units of eV * fs / angstrom", 313-314).

### 3.4 Transport integrals (boltzwann.F90:381-567)

Fermi derivative (1268-1297):
```
MyExp = (E - mu)/KT; if (abs(MyExp) > 36) 0 else 1/KT*exp(MyExp)/((exp(MyExp)+1)**2)   ! MaxExp=36 (1287-1295)
```
Per (μ,T) pair (383-520), using only the total (spin-unresolved) TDF (392):
```
IntegrandArray = TDF(:,:,1) * MinusFermiDerivative(E,mu,KT)                       (393-396)
LocalElCond = sum(IntegrandArray, DIM=2) * tdf_energy_step                        (399)   ! = K0; sigma/e^2
IntegrandArray *= (E - mu)                                                        (472-474)
LocalSigmaS = sum(IntegrandArray, DIM=2)*tdf_energy_step/TempArray(TempIdx)       (477)   ! = K1/T = (sigma*S)*T/e /T
ThisSeebeck = -matmul(ElCondInverse, SigmaS_FP)     ! sign: electron charge < 0   (494)
IntegrandArray *= (E - mu)   ! again -> (E-mu)^2                                  (514-516)
LocalKappa = sum(IntegrandArray, DIM=2)*tdf_energy_step/TempArray(TempIdx)        (517)   ! = K2/T
```
i.e. σ = e²∫dε(−∂f/∂ε)Σ(ε); σS = (e/T)∫dε(−∂f/∂ε)Σ(ε)(ε−μ); S = σ⁻¹·(σS);
K = (1/T)∫dε(−∂f/∂ε)Σ(ε)(ε−μ)² (the "K coefficient", ingredient of κ). ElCond/SigmaS/Kappa
are symmetric 6-component packed; **Seebeck is a full 3×3, 9 components, row-major
(xx,xy,xz,yx,yy,yz,zx,zy,zz)** (495-501). Zero-determinant σ ⇒ Seebeck=0 + warning (457-467,
524-533).

Unit conversions (535-567):
```
LocalElCond = LocalElCond*physics%elem_charge_SI**3/(physics%hbar_SI**2)*1.e-5_dp   ! (546) -> 1/Ohm/m
LocalSigmaS = LocalSigmaS*  (same factor)                                           ! (553) -> Ampere/m/K
! Seebeck already in volt/kelvin (556; derivation 503-510)
LocalKappa  = LocalKappa *  (same factor)                                           ! (566) -> W/m/K
```
(derivation of `e³/ħ²·1e-5`: eV·fs/Å × conversion; comments 538-545, 558-565).
TDF is written **without** any conversion (in 1/ħ² · eV·fs/Å).

### 3.5 Output files (formats 101-104 at boltzwann.F90:774-777)

- `seedname_tdf.dat` (318-333): header lines
  `# Energy TDF_xx TDF_xy TDF_yy TDF_xz TDF_yz TDF_zz` (+12 extra columns if spin_decomp:
  6 up then 6 down); rows `101 FORMAT(7G18.10)` (or `102 FORMAT(19G18.10)`).
- `seedname_elcond.dat` (633-642): `# Mu(eV) Temp(K) ElCond_xx ElCond_xy ElCond_yy ElCond_xz
  ElCond_yz ElCond_zz`; rows `103 FORMAT(8G18.10)`; loop MuIdx outer, TempIdx inner.
- `seedname_sigmas.dat` (646-655): `# Mu(eV) Temp(K) (Sigma*S)_xx ... (Sigma*S)_zz`
  (Ampere/m/K), format 103.
- `seedname_seebeck.dat` (659-669): `# Mu(eV) Temp(K) Seebeck_xx Seebeck_xy Seebeck_xz
  Seebeck_yx Seebeck_yy Seebeck_yz Seebeck_zx Seebeck_zy Seebeck_zz`; rows
  `104 FORMAT(11G18.10)`.
- `seedname_kappa.dat` (673-686): `# Mu(eV) Temp(K) Kappa_xx Kappa_xy Kappa_yy Kappa_xz
  Kappa_yz Kappa_zz` (W/m/K), format 103.
- `seedname_boltzdos.dat` (994, 1179-1212): header comments incl. smearing info and
  `# Cell volume (ang^3):`, `# Energy(eV) DOS [DOS DOS ...]`; rows
  `'(1X,<1+ndim>G18.10)') DOS_EnergyArray(EnIdx), dos_all(EnIdx,:)`. DOS from `dos_get_k`
  with `pw90_boltzwann%dos_smearing`; **adaptive-smearing zero-velocity refinement**: if any
  `|levelspacing_k| < SPACING_THRESHOLD = 1.e-3` (899), the point is replaced by its 8
  neighbours at ±1/4 grid spacing in each direction, each weighted kweight/8
  (1094-1143). DOS in states/eV (not divided by cell volume), like dos.F90.

### 3.6 testpostw90_boltzwann

jobconfig:273-277: program `POSTW90_BOLTZWANN_ELCOND_OK`, `inputs_args = ('silicon.win','')`,
`output = silicon_elcond.dat` — only elcond is compared. Parser
`parse_boltzwann.parse_elcond`: skips `#`, requires exactly 8 tokens → mu, temp, elcond_xx,
xy, yy, xz, yz, zz. Tolerances (userconfig:232-241):
`tolerance = ( (10., 1.0e-4, 'elcond_xx'), ... same for xy, yy, xz, yz, zz)` — abs 10 (S/m),
rel 1e-4 ("Abs tolerances are big (10) because the absolute values are also big").
Ships: `silicon.win`, `silicon.chk.fmt.bz2`, `silicon.eig`, `silicon.amn`, `silicon.mmn.bz2`.
silicon.win sets: `num_bands=12, num_wann=8, use_ws_distance=true, boltzwann=true,
boltz_calc_also_dos=true, boltz_dos_energy_step=0.1, smr_type=gauss,
boltz_dos_adpt_smr=false, boltz_dos_smr_fixed_en_width=0.03, kmesh=20 (→20³),
boltz_mu_min=5., boltz_mu_max=5., boltz_mu_step=0.01, boltz_temp_min=300.,
boltz_temp_max=300., boltz_temp_step=50, boltz_relax_time=10.`, `dis_win_max=17.0,
dis_froz_max=6.4`, mp_grid 4 4 4. Benchmark single data row:
`5.000000000 300.0000000 6504317.803 -207173.9111 6669730.422 200281.9835 209188.4251 6478066.135`.

---

## 4. kslice (kslice.F90)

### 4.1 Parameters

Defaults set at postw90_readwrite.F90:161-164 (and 335-338):
`kslice_corner = 0.0 0.0 0.0`, `kslice_b1 = 1 0 0`, `kslice_b2 = 0 1 0`,
`kslice_2dkmesh = 50 50` (a single value is duplicated to both; must be > 0;
readwrite:530-552). `kslice_task` default `'fermi_lines'` (types:121); allowed tokens
fermi_lines/curv/morb/shc; pairwise exclusions curv+morb, shc+morb, shc+curv
(readwrite:507-528). `kslice_fermi_lines_colour` ∈ {none, spin}, default none (566-574).
Task flags (kslice.F90:174-179): `plot_fermi_lines = index(task,'fermi_lines')>0` etc.;
`heatmap = plot_curv .or. plot_morb .or. plot_shc`; spin-coloured fermi_lines forbidden with
heatmap (180-184). SHC extra constraints (185-198): **fixed smearing required**
("Error: Must use fixed smearing when plotting spin Hall conductivity") and exactly one
Fermi energy.

### 4.2 Geometry and grid (kslice.F90:269-357)

b1, b2 are **fractional** (reciprocal-lattice) coordinates; Cartesian
`bvec(1,:) = matmul(b1, recip_lattice)` etc. (272-273); zvec = b1×b2, area = |zvec| (must be
≠ 0); yvec = zvec×b1; dual vectors avec_2d from the triad (275-295). `square='True'` iff
b1⊥b2 and |b1|=|b2| (307-311).

```
nkpts = (kmesh2d(1) + 1)*(kmesh2d(2) + 1)                 ! (313) — INCLUSIVE endpoints, N+1 points per side
i2 = itot/(kmesh2d(1) + 1)   ! slow index (341)
i1 = itot - i2*(kmesh2d(1) + 1)   ! fast index (342)
k1 = i1/real(kmesh2d(1), dp);  k2 = i2/real(kmesh2d(2), dp)      ! in [0,1] (345-346)
kpt = kslice_corner + k1*kslice_b1 + k2*kslice_b2                 ! fractional (347)
```
2D plot coordinates: k1,k2 shifted by the corner's in-plane projection (351-353), then
`kpt_x = k1*b1mod + k2*b2mod*cosb1b2; kpt_y = k2*b2mod*cosyb2` (356-357).

### 4.3 Per-point evaluation (kslice.F90:361-455)

Operators loaded up front: HH_R always (206); AA_R for curv/morb (212-223); BB_R+CC_R for
morb (224-236); AA_R+SS_R+SHC_R for shc (238-260); SS_R for spin-coloured fermi lines
(262-268).
- fermi_lines (plain): Fourier + diagonalize → `my_bandsdata(:,iloc) = eig(:)` (385-395).
- fermi_lines (colour=spin): `spin_get_nk` → spn_k clamped to ±(1−eps8); `wham_get_eig_deleig`;
  mask points with `|eig(n) − E_F| < Delta_E` where `Delta_E = |v_in-plane|·Delta_k`,
  `Delta_k = max(b1mod/N1, b2mod/N2)` (362-406).
- curv: `berry_get_imf_klist` at `fermi_energy_list` → `curv(i) = sum(imf_k_list(:, i, 1))`;
  `/bohr**2` if `berry_curv_unit='bohr2'`; stored as **−curv** ("Print _minus_ the Berry
  curvature", 409-424). Same −2Im⟨∂u|∂u⟩ AHC integrand as berry (see berry-ahc.md §C).
- morb: `berry_get_imfgh_klist` →
  `Morb_k = img_k_list(:,:,1) + imh_k_list(:,:,1) - 2*fermi_energy_list(1)*imf_k_list(:,:,1);
   Morb_k = -Morb_k/2` ("differs by -1/2 from Eq.97 LVTS12");
  `morb(i) = sum(Morb_k(:, i))` (425-441). Requires exactly one Fermi energy (1029-1032).
- shc: `berry_get_shc_klist(..., shc_k_fermi=shc_k_fermi)`;
  `my_zdata(1, iloc) = shc_k_fermi(1)` — the raw Å² curvature-like term, no S/cm factor;
  `/bohr**2` at write time if curv_unit=bohr2 (442-453, 564-568).

### 4.4 Output files (kslice.F90:511-577, helpers 1049-1115)

Point ordering in all files = grid order above (i1 fast along b1, i2 slow along b2).
- `seedname-kslice-coord.dat` (unless colour): rows `'(2E16.8)') kpt_x, kpt_y`, written by
  `write_data_file` (519-521, 1049-1068); ends with one blank line (`write (fileunit,*) ''`).
- `seedname-kslice-bands.dat` (fermi_lines, no colour): **one value per line**,
  `'(E16.8)'`, reshaped from bandsdata(num_wann, nkpts) column-major → all num_wann energies
  of point 1, then point 2, … (524-528).
- `seedname-bnd_NNN.dat` per band (fermi_lines without heatmap, gnuplot 'grid data'):
  `'(3E16.8)') kpt_x kpt_y E_n`, with a blank line after every `kmesh2d(1)+1` rows
  (blocklen; 531-543, 1097-1112) and a final blank line.
- `seedname-kslice-fermi-spn.dat` (colour=spin): `'(3E16.8)') kpt_x kpt_y spn` only for
  masked points near E_F (546-551, 1089-1096).
- `seedname-kslice-curv.dat` / `-morb.dat`: **3 values per line, no coordinates**,
  `write (dataunit, '(4E16.8)') zdata(:, loop_kpt)` (569-572), then one blank line
  (`write (dataunit, *) ' '`, 574). curv rows = (−Ω_x, −Ω_y, −Ω_z) [Å²];
  morb rows = (M_x, M_y, M_z) local morb integrand [eV·Å²].
- `seedname-kslice-shc.dat`: 1 value per line `'(1E16.8)') zdata(1, loop_kpt)` (564-568).
- Plot scripts always written on root for the relevant combos (a greenfield port can skip
  these; nothing in the test suite parses them): `-kslice-fermi_lines.gnu/.py` (579-659,
  661-719), heatmap python `-kslice-{curv,morb}_{x,y,z}[+fermi_lines].py` (721-832),
  `-kslice-shc[+fermi_lines].py` (834-965); helpers `script_common` (1117-1163),
  `script_fermi_lines` (1165-1200).

### 4.5 kslice tests

jobconfig:333-343:
```
[testpostw90_fe_kslicecurv/]  program = POSTW90_CURVDAT_OK  inputs_args = ('Fe.win','')  output = Fe-kslice-curv.dat
[testpostw90_fe_kslicemorb/]  program = POSTW90_MORBDAT_OK  inputs_args = ('Fe.win','')  output = Fe-kslice-morb.dat
```
Tolerances (userconfig:93-108); parsers parse_curv_dat / parse_morb_dat treat 3-token rows
as (bandcurvx, bandcurvy, bandcurvz) / (bandmorbx,...) with no path column:
```
[POSTW90_MORBDAT_OK] tolerance = ( (1.0e-6, 5.0e-6, 'bandpath'), (1.0e-4, 1.0e-4, 'bandmorbx'),
                                   (1.0e-4, 1.0e-4, 'bandmorby'), (1.0e-4, 1.0e-4, 'bandmorbz'))
[POSTW90_CURVDAT_OK] # relative error not helpful for e-20ish values expected
                     tolerance = ( (1.0e-6, 5.0e-6, 'bandpath'), (1.0e-6, 1.0e+3, 'bandcurvx'),
                                   (1.0e-6, 1.0e+3, 'bandcurvy'), (1.0e-6, 1.0e+3, 'bandcurvz'))
```
Both Fe dirs ship `Fe.win, Fe.chk.fmt.bz2, Fe.eig, Fe.amn, Fe.mmn.bz2, Fe.uHu.bz2` (uHu is
needed for morb's CC_R; curv only needs mmn but the dir ships uHu anyway). Fe.win (identical
except `kslice_task`): `num_bands=28, num_wann=18, spinors=true, use_ws_distance=.false.,
fermi_energy=12.6279, uHu_formatted=.true.`, `kslice=true`,
`kslice_task = curv+fermi_lines` (curv test) / `morb+fermi_lines` (morb test),
`kslice_2dkmesh = 5 5` → 6×6 = 36 points, `kslice_corner = 0 0 0`,
`kslice_b1 = 0.5 -0.5 -0.5`, `kslice_b2 = 0.5 0.5 0.5`, mp_grid 2 2 2. Benchmarks: 36 rows ×
3 cols (benchmark files stored without the trailing blank line); e.g. curv row 2
`-0.10052614E-04 0.32238014E-04 -0.56272413E+01`, morb row 1
`0.71495070E-07 0.23833843E-06 0.32522033E+00`.
