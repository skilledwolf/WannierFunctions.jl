# postw90 reference notes: SHC via Ryoo–Park–Souza method (shc_method = ryoo)

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Companion notes (shared machinery documented there, only referenced here):
- `docs/reference-notes/berry-ahc.md`: `v_matrix` construction, `get_HH_R`, `get_AA_R`
  (mmn reader), `fourier_q_to_R`, `operator_wigner_setup` (ndegen / use_ws_distance),
  `pw90common_fourier_R_to_k_new/_vec` conventions, `wham_get_eig_deleig`, `wham_get_D_h`.
- `docs/reference-notes/shc-dos-boltz-kslice.md` §1: `.spn` format, `get_SS_R`, qiao-side
  `get_SHC_R`, generic `berry_get_shc_klist` band-sum, SHC units/output, testpostw90_pt_shc.

Papers (berry.F90:43-44): QZYZ18 = Qiao, Zhou, Yuan, Zhao, PRB 98, 214402 (2018) [qiao];
RPS19 = Ryoo, Park, Souza, PRB 99, 235113 (2019) [ryoo]; GP22 = PRB 106, 075126 (2022)
(tetrahedron method, berry.F90:46 — separate `berry_get_shc_tetrahedron`, berry.F90:3101ff,
only when `tetrahedron_method=T`; not covered here).

---

## A. Method dispatch and operator inventory

`eval_shc` when `index(berry_task,'shc') > 0` (berry.F90:286). Operator setup
(berry.F90:413-454) — for BOTH methods first:
`get_HH_R`, `get_AA_R` (reads `.mmn`), `get_SS_R` (reads `.spn`). Then:

- **qiao** (`index(pw90_spin_hall%method,'qiao') > 0`, berry.F90:434-439):
  `get_SHC_R` → `SR_R`, `SHR_R`, `SH_R` (from `.spn` + `.mmn`, no extra files).
- **ryoo** (anything else that passed the qiao/ryoo validation, berry.F90:440-453):
  1. `get_SH_R` (get_oper.F90:2272-2520) → `SH_R(num_wann,num_wann,nrpts_pw90,3)`
     = ⟨0n|σ_{x,y,z}·H|Rm⟩, from `.spn` only (same spn reader as get_SS_R;
     `SH_o(:,:,ik,is) = matmul(spn_o(:,:,ik,is), H_o(:,:,ik))` with
     `H_o = diag(eigval(:,ik))` + scissors_shift on bands > num_valence_bands OR
     shc_bandshift on bands ≥ bandshift_firstband, get_oper.F90:2443-2461, 2476-2484;
     then `get_gauge_overlap_matrix` V†·SH_o(window)·V at same k, `fourier_q_to_R`,
     `operator_wigner_setup` per is, 2486-2496; NO factor of i).
  2. `get_SAA_R` (get_oper.F90:2822-3117) → `SAA_R(num_wann,num_wann,nrpts_pw90,3,3)`
     from **`.sIu`**; comment (2827-2828): `SAA_ab(R) = <0|s_a.(r-R)_b|R>`, FT of
     `SAA_ab(k) = <u|s_a|del_b u>` (a = spin index, b = Cartesian derivative).
  3. `get_SBB_R` (get_oper.F90:2523-2819) → `SBB_R(num_wann,num_wann,nrpts_pw90,3,3)`
     from **`.sHu`**; comment (2528-2529): `SBB_ab(R) = <0|s_a.H.(r-R)_b|R>`, FT of
     `SBB_ab(k) = <u|s_a.H|del_b u>`.

ryoo does NOT use SR_R/SHR_R; the interpolation-time formula uses HH_R, AA_R, SS_R,
SAA_R, SBB_R (SH_R is needed only inside get_SBB_R's transl_inv_full correction).
Runtime files for ryoo: `seedname.win`, `.chk`, `.eig`, `.mmn` (for AA_R), `.spn`
(for SS_R and SH_R), `.sHu`, `.sIu`.

Validation (postw90_readwrite.F90:1105-1125): `shc_method` is MANDATORY with
berry_task=shc ("Error: berry_task=shc and shc_method is not set", 1108-1111); the value
must contain 'qiao' or 'ryoo' (1112-1116); `transl_inv_full=T` is fatal with qiao
("Error: transl_inv_full=T not implemented for shc_method=qiao", 1117-1124).

---

## B. .sHu and .sIu file formats (readers: get_oper.F90:2633-2705 / 2932-3005)

Both files are **Fortran unformatted sequential only** in this codebase — the `open` is
unconditionally `form='unformatted'` (get_oper.F90:2633-2634 for `.sHu`,
2932-2933 for `.sIu`). There is NO formatted variant and NO `sHu_formatted`/`sIu_formatted`
keyword (`pw90_oper_read_type` has only `spn_formatted` and `uHu_formatted`,
postw90_types.F90:64-72). (pw2wannier90's `write_sHu/write_sIu` also has formatted writers,
but this postw90 cannot read them.)

Record layout (identical for both files; `.sHu` shown, lines 2637-2678):

```
record 1: header                          ! character(len=60)
record 2: nb_tmp, nkp_tmp, nntot_tmp      ! 3 default integers; must equal
                                          ! num_bands, num_kpts, kmesh_info%nntot (2640-2651)
then, loop ik = 1, num_kpts
       loop nn2 = 1, nntot                ! neighbour index IN nnlist ORDER
         loop ipol = 1, 3                 ! sigma_x, sigma_y, sigma_z
  record: ((X(n, m, ipol), n=1,num_bands), m=1,num_bands)     ! complex(kind=dp), n fastest
```

so there are `num_kpts*nntot*3` matrix records, one full `num_bands×num_bands`
complex(dp) matrix each; loop order ipol fastest, then nn2, then ik (2654-2703).

Immediately after each read the code transposes (2676-2677, comment "pw2wannier90 is coded
a bit strangely, so here we take the transpose"):

```
Ho_q_qb2(:,:,ipol) = transpose(Ho_q_qb2(:,:,ipol))
```

After the transpose, `Ho_q_qb2(ii, jj, ipol)` with ii a band at q (bra) and jj a band at
q+b2 (ket) is (comments 2670-2672 / 2969-2971):

```
.sHu:  Ho_q_qb2(ii,jj,ipol) = <u_{ii,q} | sigma_ipol * H_q | u_{jj,q+b2}>
.sIu:  Ho_q_qb2(ii,jj,ipol) = <u_{ii,q} | sigma_ipol         | u_{jj,q+b2}>
```

⇒ **on disk (before transpose) the fast index n runs over the KET band at q+b2 and the
slow index m over the BRA band at q**: file element (n,m) = ⟨u_{m,q}|σ H_q|u_{n,q+b2}⟩.
`H_q` is the ab-initio Hamiltonian at the BRA k-point q. σ are Pauli matrices
(dimensionless; the ħ/2 of the spin operator is restored in the final unit factor).
q+b2 = `kmesh_info%nnlist(ik, nn2)` — the reader does NOT verify neighbour identity per
block (unlike the `.mmn` reader); records are trusted to be in postw90's own
nnlist order, which must coincide with the `.nnkp` order used by pw2wannier90.

Error handling: open failure → `'Error: Problem opening input file seedname.sHu'`
(get_oper.F90:2814, label 111); read failure/EOF → `'Error: Problem reading input
file seedname.sHu'` (2816, label 112); dims mismatch → `'seedname.sHu has not the right
number of bands / k-points / nearest neighbours'` (2640-2651). Same texts with `.sIu`
(3112-3115, labels 113/114; dims 2939-2950).

---

## C. get_SAA_R / get_SBB_R — R-space construction (default transl_inv_full = F)

Both routines are structurally identical; only the input file and the correction operator
differ (SBB↔`.sHu`/SH_R, SAA↔`.sIu`/SS_R). Early-exit if the target array is already
allocated ("been here before" pattern, 2587-2593 / 2886-2892). Fatal if
`|scissors_shift| > 1e-7`: `'Error: scissors correction not yet implemented for SBB_R'`
(2604-2607) / `'... for SAA_R'` (2903-2906) — **ryoo + scissors_shift is unusable**.

Per (ik, nn2, ipol), after reading + transposing `Ho_q_qb2`:

1. **Windowed gauge rotation** V†·X·V (hand-coded, get_oper.F90:2680-2695 / 2979-2994) —
   the exact parallel of `get_gauge_overlap_matrix` with bra at ik, ket at qb2:

```
H_q_qb2(n,m) = sum_{i=1..num_states(ik)} sum_{j=1..num_states(qb2)}
               conjg(v_matrix(i,n,ik)) * Ho_q_qb2(winmin_q+i-1, winmin_qb2+j-1, ipol)
                                       * v_matrix(j,m,qb2)
```

   `num_states(ik) = ndimwin(ik)` if disentangled else num_wann (2624-2631);
   `winmin` = first band with lwindow=T (`get_win_min`, 3199-3232; assumes a contiguous
   window). NOTE: `H_q_qb2 = cmplx_0` sits BEFORE `do ipol = 1,3` (2680-2681 / 2979-2980)
   and the sum is `H_q_qb2(n,m) = H_q_qb2(n,m) + ...` — it is NOT reset between ipols,
   so successive spin components pile up (σ_x, then σ_x+σ_y, then σ_x+σ_y+σ_z). See
   trap T8.

2. **b-sum accumulation** (2696-2700 / 2995-2999), with `w = kmesh_info%wb(nn2)` (Å²),
   `b = kmesh_info%bk(1:3, nn2, ik)` (Cartesian Å⁻¹), canonical slot
   `nn2o = kmesh_info%nninv(nn2, ik)` (index such that `bk(:,nn2o,1) = bk(:,nn2,ik)`;
   types.F90:198, built by `kmesh_bvectors_perm`, kmesh.F90:2280-2329, tol 1e-7):

```
do b_dir = 1, 3
  SBB_q_b(:,:,ik,nn2o,ipol,b_dir) += cmplx_i * phase1(:,:) * wb(nn2) * bk(b_dir,nn2,ik) * H_q_qb2(:,:)
```

   `phase1 = 1` when transl_inv_full=F (2663-2665 / 2962-2964). Note the **+i factor is
   applied here at accumulation** (contrast qiao's get_SHC_R, which multiplies SR_R/SHR_R
   by `cmplx_i` only at the very end, get_oper.F90:2238-2239). There is NO subtraction of
   the on-site term (qiao subtracts `SS_q`/`SH_q` inside the b-sum; here the analogous
   term vanishes by the Σ_b w_b·b = 0 sum rule and is simply omitted). There is NO
   hermitization (these operators are not Hermitian).

3. **q→R Fourier + WS setup** (default branch 2770-2803 / 3069-3102):

```
SBB_q = sum(SBB_q_b, dim=4)                          ! sum over neighbour slots (2774/3073)
per (ipol, b_dir):  fourier_loc_q_to_R:  O(R) = (1/N_q) sum_q e^{-i 2π q·R} O(q)
                    operator_wigner_setup: /ndegen (+ ws_distance remap)      (2786-2798)
comms_bcast SBB_R (2805/3104)
```

Resulting index convention: `SAA_R(n, m, ir, ipol, b_dir)`, `SBB_R(...)` — 4th index =
spin direction γ∈{x,y,z}, 5th = Cartesian derivative/position direction.

---

## D. What transl_inv_full = T changes (test "pt_shc_ryoo_transl_inv")

Keyword `transl_inv_full` (logical, default `.false.`, postw90_types.F90:176; read at
postw90_readwrite.F90:854-856; fatal together with `transl_inv=T`, 857-860; fatal with
shc_method=qiao, 1117-1124). **The test directory name says "transl_inv" but the keyword
in its Pt.win is `transl_inv_full = .true.`** (`transl_inv` itself — the MV97 Im-log
diagonal for AA_R, default F, postw90_types.F90:175 — is untouched by these tests).

With transl_inv_full = T, three operators change:

**(1) AA_R** (get_AA_R, get_oper.F90:608-609, 649-729): b-contributions are stored in
canonical slots `nno = nninv(nn,ik)` (609); after the mmn loop each slot is multiplied by
`phase1 = exp(+i b·r0)` with `r0(i,j,:) = (τ_i + τ_j)/2` (Cartesian Å,
`wannier_centres_from_AA_R`; 649-673); the q→R transform is then done PER SLOT and each
R-matrix gets `phase2 = exp(-(i/2) B·R_cart)` with `B = bk(:,nn,1)` and
`R_cart = crvec_pw90(:,ir)` (690-709); **no hermitization** (the 0.5(A+A†) of the default
branch, 743-749, is skipped); finally the R=0 diagonal is overwritten:
`AA_R(i,i,ir0,:) = wannier_centres_from_AA_R(:,i)` (719-728).

**(2) SAA_R / SBB_R** (get_oper.F90:2658-2665+2710-2769 / 2957-2964+3009-3068):
- accumulation phase: `phase1(i,j) = exp(+i bk(:,nn2,ik)·r0(i,j))` multiplies the
  gauge-rotated matrix inside step C.2 (2658-2662/2957-2961); r0 as above
  (2615-2623/2914-2922 — requires get_AA_R to have run first, which berry_main guarantees).
- per-slot FT: instead of summing slots in q-space, each canonical slot nn2 is FT'd
  separately and multiplied by
  `phase2(ir) = exp(-(i/2) * dot(bk(:,nn2,1), crvec_pw90(:,ir)))` before being added into
  SBB_R/SAA_R (2728-2751/3027-3050).
- **diagonal-correction term** (2761-2768 for SBB with SH_R; 3060-3067 for SAA with SS_R):

```
SBB_R(:,:,ir,ipol,b_dir) += (r0(:,:,b_dir) - 0.5*crvec_pw90(b_dir,ir)) * SH_R(:,:,ir,ipol)
SAA_R(:,:,ir,ipol,b_dir) += (r0(:,:,b_dir) - 0.5*crvec_pw90(b_dir,ir)) * SS_R(:,:,ir,ipol)
```

  (element-wise product over (n,m); r0 in Å). Fatal if the partner operator is missing:
  `'transl_inv_full=T for SBB_R needs SH_R'` (2711-2713), `'...for SAA_R needs SS_R'`
  (3010-3012). Note: berry_main's ryoo branch always builds SH_R and SS_R first, so in
  practice this never fires.

Everything downstream (interpolation, js_k, unit factor, output) is IDENTICAL — the
transl_inv_full flag acts only at R-matrix construction time.

---

## E. berry_get_shc_klist (berry.F90:2684-3099) — per-k Ω-like term

Common part (both methods), per interpolation k:
1. `wham_get_eig_deleig` → eig(:), del_eig(:,1:3), delHH(:,:,1:3), HH, UU (2814-2820);
   `wham_get_D_h` → `D_h(n,m,i) = (U†∂_iH U)(n,m)/(eig(m)-eig(n))`, zero when n==m or
   |Δε|<1e-7 (2822; wan_ham.F90:102-142).
2. `shc_bandshift`: `eig(bandshift_firstband:) += bandshift_energyshift` (2825-2827).
3. `AA(:,:,i) = U†[Σ_R e^{i2πk·R} AA_R(:,:,:,i)]U + i*D_h(:,:,i)` (2829-2837, "Eq.(25) WYSV06").
4. `berry_get_js_k` fills `js_k(num_wann,num_wann)` = matrix of
   ⟨ψ_nk|½(σ_γ v_α + v_α σ_γ)|ψ_mk⟩ **without** the ħ/2 of the spin operator and without
   the 1/ħ of each velocity (comments 2700-2705; restored in the final `fac`).

### E.1 js_k, qiao branch (berry.F90:3002-3052) — for reference

```
S_k         = U† [FT SS_R(:,:,:,gamma)] U                            ! QZYZ18 (25)(36)(30)
SR_alpha_k  = -i * U† [FT_vec SR_R(:,:,:,gamma,:)](alpha) U          ! (31)
K_k         = SR_alpha_k + matmul(S_k, D_h(:,:,alpha))               ! (26)
SHR_alpha_k = -i * U† [FT_vec SHR_R(:,:,:,gamma,:)](alpha) U         ! (32)
SH_k        = U† [FT_vec SH_R](gamma) U                              ! (32)
L_k         = SHR_alpha_k + matmul(SH_k, D_h(:,:,alpha))             ! (27)
B_k(i,j)    = del_eig(j,alpha)*S_k(i,j) + eig(j)*K_k(i,j) - L_k(i,j) ! elementwise; eig_mat(i,:)=eig(:)
js_k        = 1/2 * (B_k + conjg(transpose(B_k)))                    ! (23)
```

### E.2 js_k, ryoo branch (berry.F90:3054-3095) — RPS19 Eqs.(21),(26),(37)-(40)

Only the (γ=shc_gamma, α=shc_alpha) component of SAA/SBB is interpolated:

```
SAA = U† [FT SAA_R(:,:,:,gamma,alpha)] U          ! scalar fourier_R_to_k_new, 3056-3060, 3076-3077
SBB = U† [FT SBB_R(:,:,:,gamma,alpha)] U          ! 3062-3066, 3078-3079
(HH, delHH re-interpolated from HH_R at 3069-3073 — redundant recomputation, same values)
VV0 = U† delHH(:,:,alpha) U                       ! bar-velocity  = U†(∂_α H^W)U, 3075
spinVel0 = matmul(VV0, S_k) + matmul(S_k, VV0)    ! 3081-3082   (S_k as in E.1)

do n; do m:                                       ! 3084-3093
  js_k(n,m) = spinVel0(n,m)
              - i*( eig(m)*SAA(n,m) - SBB(n,m) )
              + i*( eig(n)*conjg(SAA(m,n)) - conjg(SBB(m,n)) )
js_k = js_k / 2                                   ! 3094
```

Index/conjugation detail: the second line uses the (n,m) elements, the third the complex
conjugates of the TRANSPOSED elements (m,n); eig(m) multiplies SAA(n,m), eig(n) multiplies
conjg(SAA(m,n)). Since SAA/SBB already carry the +i from construction (C.2), the −i/+i here
produce the Hermitian combination of ⟨u|σ(H−ε)|∂u⟩-type terms; no further hermitization is
applied. vs. qiao: qiao builds B_k from band-diagonal `del_eig`·S plus ε·K − L then takes
½(B+B†); ryoo uses the full matrix VV0 (not just diagonal velocities) and the ε-weighted
SAA/SBB difference directly.

### E.3 Band sum (both methods, berry.F90:2846-2903)

For n = 1..num_wann, m ≠ n, skipping any pair with
`eig(m) > kubo_eigval_max .or. eig(n) > kubo_eigval_max` (2867):

```
rfac = eig(m) - eig(n)
prod = js_k(n,m) * cmplx_i * rfac * AA(m,n,shc_beta)                 ! 2869-2872
eta_smr: adaptive (kubo_adpt_smr=T): min(|del_eig(m,:)-del_eig(n,:)|*Δk*kubo_adpt_smr_fac,
                                         kubo_adpt_smr_max)          ! Eq.(35) YWVS07, 2873-2878
         fixed:    eta_smr = kubo_smr_fixed_en_width                 ! 2880
Fermi scan / band:  omega += (-2/(rfac**2 + eta_smr**2)) * aimag(prod)          ! 2888-2891
Freq scan:  do ifreq: cdum = Re(kubo_freq_list(ifreq)) + i*eta_smr
            omega_list(ifreq) += (-2/(rfac**2 - cdum**2)) * aimag(prod)         ! 2882-2887
```

then occupation weighting (T=0 step, `pw90common_get_occ`: occ=1 iff eig<E_F):

```
Fermi scan: shc_k_fermi(i) += occ_fermi(n,i) * omega          ! per fermi_energy_list(i), 2894-2897
Freq scan:  shc_k_freq(:)  += occ_freq(n) * omega_list(:)     ! occ at fermi_energy_list(1), 2898-2899
Band mode:  shc_k_band(n) = omega                             ! kpath colouring, 2900-2901
```

Δk = `pw90common_kmesh_spacing(berry_kmesh, recip_lattice)` (2846-2849) — adaptive smearing
is only meaningful on the berry_kmesh grid (comment 2844-2845; kpath/kslice must use fixed).
The per-k result is in Å² (comment 2698-2699). shc_k_freq is COMPLEX; shc_k_fermi real.

---

## F. Assembly, frequency scan, units, output

- Fermi scan is the default (`shc_freq_scan = .false.`); freq scan requires exactly ONE
  Fermi energy (`not_scannable = eval_kubo .or. (eval_shc .and. freq_scan)`,
  berry.F90:342-347).
- Regular-grid loop (berry.F90:968-1046; kpoint-file branch 760-836): `kweight = 1/N_k`.
  Fermi-scan mode supports adaptive k-mesh refinement: trigger when
  `berry_curv_adpt_kmesh > 1` and `|shc_k_fermi(if)|` (converted /bohr² if
  `berry_curv_unit='bohr2'`, using `physics%bohr`) `> berry_curv_adpt_kmesh_thresh` at ANY
  Fermi level (788-803/996-1011); the point is then recomputed on the
  `curv_adpt_kmesh**3` sub-cluster with `kweight_adpt = kweight/curv_adpt_kmesh**3`,
  REPLACING the coarse contribution (804-823). **Freq-scan mode has NO adaptive kmesh**
  (824-836). MPI reduce at 1228-1237.
- Unit conversion (berry.F90:1675-1695, comments quoted):

```
! (i)   multiply -e^2/hbar/(V*N_k) as in the QZYZ18 Eq.(5)  (1/N_k already in kweight)
! (ii)  convert charge current to spin current: overall -hbar/2/e
! (iii) multiply 1e8 to convert to S/cm
fac = 1.0e8_dp * physics%elem_charge_SI**2 / (physics%hbar_SI * cell_volume) / 2.0_dp   ! 1690
shc_freq = shc_freq*fac      (freq scan)    /    shc_fermi = shc_fermi*fac              ! 1691-1695
```

  **POSITIVE sign** (AHC uses −1e8·e²/(ħV), berry.F90:1360). cell_volume in Å³
  (det(real_lattice), berry.F90:262-267). Final unit: (ħ/e)·S/cm.
  Constants default to the **CODATA2006** branch (preprocessor default, constants.F90:92-100):
  `elem_charge_SI = 1.602176487e-19` C, `hbar_SI = 1.054571628e-34` J·s,
  `bohr = 0.52917720859` Å (constants.F90:161-187, 224). Build flags -DCODATA2010/2018/2022
  select other sets (2018/2022: e=1.602176634e-19, ħ=1.054571817e-34) — sub-1e-7 relative
  effect on fac, but flag it for bit-reproducibility.

Output files (berry.F90:1704-1727), written by root:

```
Fermi scan → seedname-shc-fermiscan.dat
  header: write(unit,'(a,3x,a,3x,a)') '#No.', 'Fermi energy(eV)', 'SHC((hbar/e)*S/cm)'
  rows:   write(unit,'(I4,1x,F12.6,1x,E17.8)') n, fermi_energy_list(n), shc_fermi(n)

Freq scan → seedname-shc-freqscan.dat
  header: write(unit,'(a,3x,a,3x,a,3x,a)') '#No.', 'Frequency(eV)',
          'Re(sigma)((hbar/e)*S/cm)', 'Im(sigma)((hbar/e)*S/cm)'
  rows:   write(unit,'(I4,1x,F12.6,1x,1x,2(E17.8,1x))') n,
          real(kubo_freq_list(n),dp), real(shc_freq(n),dp), aimag(shc_freq(n))
```

stdout (berry.F90:540-551): `'* Spin Hall Conductivity'`, then
`"  Qiao's SHC (Phys.Rev.B 98.214402)"` or `"  Ryoo's SHC (Phys.Rev.B 99.235113)"`, then
`'  Frequency scan'` / `'  Fermi energy scan'`. get_SAA_R/get_SBB_R announce
`' Reading sIu overlaps from seedname.sIu in get_SAA_R: '` + header (2934-2937), similarly
for sHu (2635-2638).

kpath/kslice reuse: `berry_get_shc_klist(..., shc_k_band=...)` for band colouring and
`shc_k_fermi` (fermi_n must be 1) for the shc curve (kpath.F90:355-366, 404-414).

---

## G. Keywords (defaults from postw90_types.F90 / parse in postw90_readwrite.F90)

| keyword | default | notes |
|---|---|---|
| `shc_method` | `' '` — **mandatory** for berry_task=shc | must contain 'qiao' or 'ryoo' (readwrite 1105-1116) |
| `shc_freq_scan` | `.false.` (types 202) | T = ac-SHC ω scan; F = Fermi scan (readwrite 1047-1049) |
| `shc_alpha` | 1 (types 203) | spin-current flow direction α ∈ {1,2,3} (readwrite 1051-1057) |
| `shc_beta` | 2 (types 204) | E-field direction β (1059-1065) |
| `shc_gamma` | 3 (types 205) | spin polarisation γ (1067-1073); default = σ^{z}_{xy} |
| `shc_bandshift` | `.false.` (types 206) | rigid shift of upper bands; only active with berry & freq_scan (1075-1079); exclusive with scissors_shift (1080-1083) |
| `shc_bandshift_firstband` | 0 (types 207) | 1-based ab-initio band index, required if bandshift (1085-1095) |
| `shc_bandshift_energyshift` | 0.0 (types 208) | eV, required if bandshift (1097-1103) |
| `transl_inv_full` | `.false.` (types 176) | ryoo-only one-shell scheme, §D (851-860) |
| `transl_inv` | `.false.` (types 175) | MV97 diag of AA_R; NOT the keyword of the ryoo_transl_inv test |
| `spn_formatted` | `.false.` (types 68) | .spn text vs stream (readwrite 477-478); no analogue for sHu/sIu |
| `kubo_freq_min` | 0.0 eV (readwrite pw90_extra_io, line 56) | (1684-1686) |
| `kubo_freq_max` | froz_max−E_F(1)+0.6667 if frozen states else maxval(eigval)−minval(eigval)+0.6667 (1688-1694) | eV |
| `kubo_freq_step` | 0.01 eV (line 58) | `kubo_nfreq = nint((max−min)/step)+1` (≥2); step recomputed = (max−min)/(nfreq−1); list linear, real (1708-1724) |
| `kubo_eigval_max` | froz_max+0.6667 if frozen states, else maxval(eigval)+0.6667, else win_max+0.6667 (1758-1768) | eV; pair-skip cutoff E.3 |
| `kubo_adpt_smr` | inherits `adpt_smr` = `.true.` (types 133; readwrite 908-911) | |
| `kubo_adpt_smr_fac` | inherits `adpt_smr_fac` = √2 (types 134) | |
| `kubo_adpt_smr_max` | inherits `adpt_smr_max` = 1.0 eV (types 137) | |
| `kubo_smr_fixed_en_width` | inherits `smr_fixed_en_width` = 0.0 eV (types 136) | η enters ANALYTICALLY (iη / η²) — `kubo_smr_type` is irrelevant to SHC |
| `berry_kmesh` | falls back to global `kmesh` | 1 value → n×n×n, or 3 values (get_module_kmesh, readwrite 1902-1998) |
| `berry_curv_adpt_kmesh` | 1 = off (types 166) | Fermi-scan only for SHC |
| `berry_curv_adpt_kmesh_thresh` | 100.0 (types 167) | Å² unless curv_unit=bohr2 |
| `berry_curv_unit` | 'ang2' (types 168) | affects threshold comparison only |
| `fermi_energy(_min/_max/_step)` | — (src/readwrite.F90:609-692) | list linear, step default 0.01; freq scan needs exactly 1 |
| `use_ws_distance` | `.true.` (src/types.F90:87-94) | via operator_wigner_setup on ALL R-operators |
| `scissors_shift` | 0.0 | **fatal with ryoo** (get_SAA_R/get_SBB_R §C) |

---

## H. Test cases

Parser for both dat flavours: `test-suite/tools/parsers/parse_shc_dat.py` (header token
`'Fermi'` vs `'Frequency(eV)'`); tolerances `test-suite/tests/userconfig:201-212`:

```
[POSTW90_SHCFERMIDAT_OK]  tolerance = ((1.0e-6, 5.0e-6, 'energy'), (1.0e-1, 1.0e-1, 'shc'))
[POSTW90_SHCFREQDAT_OK]   tolerance = ((1.0e-6, 5.0e-6, 'frequency'),
                                       (1.0e+1, 1.0e+1, 'shc_re'), (1.0e-1, 1.0e-1, 'shc_im'))
```

### H.1 testpostw90_pt_shc_ryoo (jobconfig:485-488, POSTW90_SHCFREQDAT_OK, output Pt-shc-freqscan.dat)

Ships `Pt.win, Pt.eig, Pt.amn.bz2, Pt.chk.fmt.bz2, Pt.mmn.bz2, Pt.spn.bz2, Pt.sHu.bz2,
Pt.sIu.bz2`. Pt.win keywords (full list of non-structural ones):

```
num_bands = 24, num_wann = 18, spinors = T, projections Pt:l=0;l=1;l=2
fermi_energy = 18.3823
berry = true, berry_task = shc
shc_freq_scan = .true., shc_method = ryoo
berry_kmesh = 9 9 9
kubo_adpt_smr = .false., kubo_smr_fixed_en_width = 0.1
kubo_eigval_max = 1000
kubo_freq_min = 0.00, kubo_freq_max = 7.00, kubo_freq_step = 0.1      ! kubo_nfreq = 71
use_ws_distance = .false.
shc_gamma = 1, shc_alpha = 3, shc_beta = 2        ! sigma^{spin-x}_{z y}
spn_formatted = true
fcc Pt, a/2 = 3.6963 bohr; mp_grid = 4 4 4 (64 explicit kpoints)
```

num_bands(24) > num_wann(18) with no dis_win keywords → disentangled with default (full)
window. Benchmark `benchmark.out.default.inp=Pt.win` = header + 71 rows; first/last rows:

```
   1     0.000000    -0.25760024E+04    0.00000000E+00
   2     0.100000    -0.25950607E+04   -0.40023211E+02
  71     7.000000    -0.67676662E+02   -0.37714785E+03
```

### H.2 testpostw90_pt_shc_ryoo_transl_inv (jobconfig:491-494, POSTW90_SHCFREQDAT_OK)

Same Pt.win as H.1 **except** (diff of the two files):

```
transl_inv_full = .true.
use_ws_distance = .true.        ! (ryoo test has .false.)
```

Benchmark first rows:

```
   1     0.000000    -0.21304230E+04    0.00000000E+00
   2     0.100000    -0.21322109E+04   -0.35314936E+02
```

(≈17% shift of the ω→0 value vs H.1 — transl_inv_full + ws_distance together; useful as a
sensitivity check that both switches are actually wired in.)

### H.3 testpostw90_gaas_shc (jobconfig:479-482, POSTW90_SHCFREQDAT_OK, output GaAs-shc-freqscan.dat)

GaAs uses **shc_method = qiao** with **shc_freq_scan = true** (the ac-SHC test for qiao;
NO .sHu/.sIu shipped). GaAs.win keywords:

```
shc_freq_scan = true, shc_alpha = 1, shc_beta = 2, shc_gamma = 3     ! sigma^{z}_{xy}
spn_formatted = true
berry = true, berry_task = shc, berry_kmesh = 10                      ! → 10 10 10
fermi_energy = 7.9366
kubo_freq_min = 0.0, kubo_freq_max = 8.0, kubo_freq_step = 0.01       ! kubo_nfreq = 801
scissors_shift = 1.117, num_valence_bands = 8                         ! (shc_bandshift commented out)
shc_method = qiao
kubo_adpt_smr = false, kubo_smr_fixed_en_width = 0.05
exclude_bands = 1-10, num_bands = 16, num_wann = 16                   ! num_bands == num_wann → NO disentanglement
spinors = true, As:sp3 + Ga:sp3, GaAs fcc a/2 = 5.342256 bohr, mp_grid 4 4 4
```

kubo_eigval_max not set → default maxval(eigval)+0.6667. Benchmark = header + 801 rows:

```
   1     0.000000    -0.42820457E+03    0.00000000E+00
   2     0.010000    -0.42820226E+03    0.23134881E-01
 801     8.000000     0.40466629E+03    0.38906067E+02
```

scissors_shift works here because qiao's get_SHC_R applies it to H_o
(get_oper.F90:2053-2056) — the same input would be FATAL with ryoo (§C).

---

## I. Traps checklist

T1. **.sHu/.sIu are unformatted-only** in this reference; complex(kind=dp), one full
    nb×nb matrix per record, record order (ipol fastest, nn2, ik), n (ket at k+b) fastest
    inside a record; every matrix must be TRANSPOSED after reading.
T2. Neighbour records are consumed in postw90's `nnlist` order without verification —
    the Julia kmesh must reproduce wannier90's shell search + ordering bit-for-bit,
    including `nninv` (same-b-vector index at ik=1, tol 1e-7).
T3. Factor placement: SAA/SBB get `+i` at accumulation (C.2); qiao's SR/SHR get `+i`
    only after the q→R FT (get_oper.F90:2238-2239) and then `−i` at interpolation
    (berry.F90:3013, 3025). Net objects agree; do not double-count.
T4. js_k lacks ħ/2 (spin) and 1/ħ per velocity; all restored by the single factor
    `fac = +1.0e8·e²/(ħ·V_c)/2` — note the sign is + (AHC's is −) and the /2.
T5. `kubo_smr_type` has NO effect on SHC (η enters as +iη / η² analytically). Defaults:
    adaptive smearing ON (√2 prefactor, 1.0 eV cap) unless kubo_adpt_smr=F.
T6. Pair skip: `eig(m) > kubo_eigval_max .OR. eig(n) > kubo_eigval_max` — applied
    ONLY to the (n,m) pair loop, not to js_k construction; degeneracy guard only inside
    D_h (|Δε| < 1e-7 → element zeroed). Occupations are T=0 steps at eig<E_F.
T7. Freq scan: occupations from `fermi_energy_list(1)`; exactly one Fermi energy allowed;
    output real AND imaginary parts; Im ≡ 0 at ω=0. No adaptive kmesh in freq scan.
T8. `H_q_qb2` in get_SAA_R/get_SBB_R is zeroed BEFORE the ipol loop
    (get_oper.F90:2680/2979-2980) and the windowed rotation ACCUMULATES (`+=`) with no
    per-ipol reset — so the matrix used for ipol=2 contains σ_x+σ_y and for ipol=3
    contains σ_x+σ_y+σ_z. Physically this is a bug (spin components mixed), but it is
    what the reference computes. **Both Pt ryoo tests set shc_gamma=1**, which consumes
    only the clean ipol=1 slot, so their benchmarks are unaffected; matching the
    reference for shc_gamma=2 or 3 would require reproducing the accumulation verbatim.
T9. `transl_inv_full` (NOT `transl_inv`) is the keyword of the ryoo_transl_inv test;
    it also changes AA_R (§D.1: canonical slots, e^{ib·r0}, e^{−ib·R/2}, no hermitization,
    R=0 diagonal ← Wannier centres). AA enters through `prod` (E.3), so both AA and
    SAA/SBB differences contribute to the H.2 benchmark.
T10. Ryoo needs get_AA_R run before get_SAA_R/get_SBB_R because r0 uses
    `wannier_centres_from_AA_R` (berry_main order guarantees this); centres are checked
    against chk centres within 1e-8 sum of squares unless guiding_centres=T.
T11. scissors_shift: fatal with ryoo; with shc_bandshift instead, the shift is applied to
    interpolated eig (berry.F90:2825-2827) and to H_o inside get_SH_R
    (get_oper.F90:2455-2459), but the H inside the `.sHu` matrix elements is UNSHIFTED
    (frozen into the file) — the ryoo bandshift is therefore only partially consistent.
T12. Constants: preprocessor default is CODATA2006 (e=1.602176487e-19, ħ=1.054571628e-34,
    bohr=0.52917720859) — other notes in this repo quote 2018/2022 values from the same
    file; those are the non-default #ifdef branches (constants.F90:92-100 picks 2006).
T13. use_ws_distance defaults to TRUE in postw90; H.1 turns it off, H.2 turns it on —
    every R-space operator (HH, AA, SS, SH, SAA, SBB) goes through the same
    `operator_wigner_setup` remap, so the flag must be honoured uniformly.
T14. E17.8 output fields print `-0.25760024E+04`-style (0.x mantissa) — match the exact
    Fortran descriptors `(I4,1x,F12.6,1x,E17.8)` / `(I4,1x,F12.6,1x,1x,2(E17.8,1x))`
    for file-precision comparison.
