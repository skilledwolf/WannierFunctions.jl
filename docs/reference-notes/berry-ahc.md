# postw90 reference notes: Anomalous Hall Conductivity (berry_task = ahc)

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Key papers referenced by the code (`src/postw90/berry.F90:38-45`):
- WYSV06 = Wang, Yates, Souza, Vanderbilt, PRB 74, 195118 (2006) — AHC by Wannier interpolation
- LVTS12 = Lopez, Vanderbilt, Thonhauser, Souza, PRB 85, 014435 (2012) — J0/J1/J2 trace formulas
- MV97 = Marzari & Vanderbilt, PRB 56, 12847 (1997) — Eq.(31) band-diagonal position element

---

## A. Real-space operators: get_HH_R and get_AA_R (src/postw90/get_oper.F90)

### A.0 Inputs from the checkpoint

`v_matrix` is precomputed once in `pw90common_wanint_data_dist`
(`src/postw90/postw90_common.F90:811-938`):

- no disentanglement: `v_matrix(1:num_wann,:,:) = u_matrix(1:num_wann,:,:)` (line 874)
- disentangled: `v_matrix(m,j,k) = sum_i u_matrix_opt(m,i,k) * u_matrix(i,j,k)`,
  `m = 1..ndimwin(k)` (lines 877-886)

Eigenvalues `eigval(:, :)` come from `seedname.eig`; `u_matrix`, `u_matrix_opt`, `lwindow`,
`ndimwin`, `wannier_centres` from `seedname.chk`. The overlaps M are re-read from `seedname.mmn`
by `get_AA_R` (NOT from the chk `m_matrix`).

### A.1 get_HH_R (get_oper.F90:64-291)

Wannier-gauge H(q) on each ab-initio k-point (get_oper.F90:237-247):

```
HH_q(n,m,ik) = sum_{i=1..num_states(ik)} conjg(v_matrix(i,n,ik)) * eigval(winmin_q+i-1,ik) * v_matrix(i,m,ik)
HH_q(m,n,ik) = conjg(HH_q(n,m,ik))          ! filled for n<=m, hermitized
```

`num_states(ik) = ndimwin(ik)` if disentangled else `num_wann` (lines 230-234); `winmin_q` is the
lowest band index with `lwindow(j,ik)=.true.` (`get_win_min`, get_oper.F90:3199-3232).

Then Fourier q→R (get_oper.F90:250) via `fourier_q_to_R` (get_oper.F90:3124-3157):

```
O_ij(R) = (1/N_kpts) sum_q e^{-i q·R} O_ij(q)      ! rdotq = twopi*dot_product(kpt_latt(:,ik), irvec(:,ir))
                                                    ! phase_fac = exp(-cmplx_i*rdotq); op_R = op_R/real(num_kpts,dp)
```

R is in **reduced** coordinates here (`irvec`), phase = exp(−i 2π k·R), normalisation 1/N_q.

Optional scissors shift (only if `num_valence_bands>0` and `|scissors_shift|>1e-7`,
get_oper.F90:255-278): builds P_valence in the Wannier gauge from `u_matrix` and adds
`scissors_shift * (P + delta at R=0)`.

Finally `operator_wigner_setup` (get_oper.F90:281, defined 3275-3327) divides by the R-vector
degeneracy and (if use_ws_distance) rescatters onto the per-(i,j)-shifted R grid:

```
use_ws_distance=T:  op_R_opt_ws(i,j,jr) += op_R(i,j,ir) / (ndegen(ir) * ws_distance%ndeg(i,j,ir))
                    with jr = ir_ind_ws_to_pw90(ideg,i,j,ir)              (lines 3307-3318)
use_ws_distance=F:  op_R_opt_ws(:,:,ir) = op_R(:,:,ir) / ndegen(ir)      (lines 3322-3324)
```

After this, all k-interpolation uses `nrpts_pw90 / irvec_pw90 / crvec_pw90` with NO further
degeneracy weights (they are already baked in).

### A.2 R-vector list and degeneracies

`wignerseitz` (postw90_common.F90:1622-1794): the R list is all lattice vectors inside the
Wigner-Seitz supercell of the BvK cell (mp_grid(i)*a_i), searched over
`(2*ws_search_size+1)^3` supercells with tolerance `ws_distance_tol**2` (line 1729);
`ndegen(ir)` counts ties on the WS boundary (lines 1732-1737); sum rule
`sum 1/ndegen = product(mp_grid)` enforced (lines 1779-1790).

Defaults (`src/types.F90:87-94`): `use_ws_distance = .true.`, `ws_distance_tol = 1e-5`,
`ws_search_size = 2`. **Note: postw90 default is use_ws_distance=T, but the Fe AHC test
explicitly sets `use_ws_distance = .false.`** (see E). Cartesian
`crvec(:,ir) = matmul(transpose(real_lattice), irvec(:,ir))` (postw90_common.F90:170).

### A.3 get_AA_R (get_oper.F90:403-792) — Berry connection A(R)

Reads `seedname.mmn` on root (lines 516-554): for each (ik, ik2=neighbour b) block the raw
`S_o(m,n) = <u_mk|u_n,k+b>` over num_bands. Matches (ik2,nnl,nnm,nnn) against
`kmesh_info%nnlist/nncell` to identify the shell index nn (lines 563-586).

Projects into the Wannier gauge via `get_gauge_overlap_matrix` (lines 591-594, defined
3235-3272): `S = V(k)^dagger . S_o(win_min windows) . V(k+b)` (a num_wann×num_wann matrix).

**A(k) finite-difference formula** (lines 611-632). For every direction idir and neighbour b:

```
AA_q_b(:,:,ik,nn,idir) += cmplx_i * wb(nn) * bk(idir,nn,ik) * S(:,:)        ! line 612-613
```

i.e. A_α(k) = i Σ_b w_b b_α S(k,b) — this is WYSV06 Eq.(44) **without** the "−1": the code uses
i·w_b·b·M, not i·w_b·b·(M−1) (the −1 term cancels exactly because Σ_b w_b b = 0). This same
linear form is applied to the **diagonal** as well when `transl_inv = .false.` (the DEFAULT,
`src/postw90/postw90_types.F90:175`).

Only if `transl_inv = .true.` are the band-diagonal elements rewritten à la MV97 Eq.(31)
(lines 614-631):

```
AA_q_b_diag(i,nn,idir) -= wb(nn) * bk(idir,nn,ik) * aimag(log(S(i,i)))      ! lines 619-621
...then AA_q_b(n,n,ik,nn,idir) = AA_q_b_diag(n,nn,idir)                     ! line 629
```

**Difference vs. wannier90.x `_tb.dat` writer** (`src/hamiltonian.F90:962-979`,
`hamiltonian_write_tb`): the tb.dat writer ALWAYS uses the Im-log form for the diagonal
(`pos_r -= wb*bk*aimag(log(m_matrix(i,i,nn,ik)))`, lines 966-975) and `+i*wb*bk*m_matrix(j,i,...)`
for off-diagonal (lines 977-978) — i.e. tb.dat matches postw90 with `transl_inv=T`, while
postw90's own default (`transl_inv=F`) uses the linear form for ALL elements including the
diagonal. Also the tb.dat writer does NOT hermitize; postw90 does (next paragraph).

Wannier centres consistency check: the code accumulates
`wannier_centres_from_AA_R(:,i) -= wb*bk*aimag(log(S(i,i)))/num_kpts` (lines 598-602) and
errors out if `sum((centres_from_AA_R − wann_data%centres)**2) > 1e-8` unless
`guiding_centres=T` (lines 639-647).

After the k/b loop (default `transl_inv_full=F` branch, lines 730-775):

```
AA_q = sum(AA_q_b, 4)                                            ! sum over neighbours, line 734
AA_q(:,:,ik,idir) = 0.5*(AA_q + conjg(transpose(AA_q)))          ! hermitization, lines 743-749
```

(comment at 737-741: Eq.(44) WYSV06 does not preserve Hermiticity, take Hermitean part).
Then the same `fourier_q_to_R` transform ((1/N_q) Σ_q e^{−i2πq·R}, via `fourier_loc_q_to_R`,
get_oper.F90:3160-3196) per Cartesian direction, then `operator_wigner_setup` (line 768) applies
1/ndegen (and ws_distance reordering) exactly as for H(R).

`transl_inv_full` (default F, postw90_types.F90:176) is a separate one-shell-per-b Marzari
scheme with extra e^{i b·r0} phases (lines 649-729); mutually exclusive with transl_inv
(postw90_readwrite.F90:857-860). Not needed for baseline AHC.

---

## B. k-space machinery (src/postw90/wan_ham.F90, postw90_common.F90)

### B.1 Fourier R→k and velocities

`pw90common_fourier_R_to_k_new` (postw90_common.F90:1096-1158), used by the AHC path:

```
rdotk = twopi * dot_product(kpt, irvec_pw90(:,ir));  phase_fac = exp(+i rdotk)
H(k)      = sum_R  phase_fac * HH_R(:,:,ir)                                   ! line 1149
dH_a(k)   = sum_R  cmplx_i * crvec_pw90(a,ir) * phase_fac * HH_R(:,:,ir)      ! lines 1150-1155
```

i.e. ∂_α H(k) = Σ_R e^{ik·R} (i R_α^cart) H(R) — R_α in **Cartesian Å** (crvec), phase from
**reduced** k·R. No 1/ndegen here (already folded into HH_R). Same convention in
`pw90common_fourier_R_to_k` (postw90_common.F90:1032-1093, alpha=0 → O(k); alpha=1,2,3 → i R_α O).

### B.2 A_W(k) and Omega_bar_W(k)

`pw90common_fourier_R_to_k_vec` (postw90_common.F90:1339-1403) with `OO_true=AA, OO_pseudo=OOmega`
(called from berry.F90:2020-2022):

```
A_a(k)          = sum_R e^{ik·R} AA_R(:,:,ir,a)                               ! lines 1386-1388
Omega_pseudo_1  = sum_R e^{ik·R} [ i R_2 A_3(R) − i R_3 A_2(R) ]              ! lines 1391-1393
Omega_pseudo_2  = sum_R e^{ik·R} [ i R_3 A_1(R) − i R_1 A_3(R) ]              ! lines 1394-1396
Omega_pseudo_3  = sum_R e^{ik·R} [ i R_1 A_2(R) − i R_2 A_1(R) ]              ! lines 1397-1399
```

i.e. Ω̄^W(k) = ∇×A_W(k) = Σ_R e^{ikR} iR × A(R), packed as pseudovector with the module-level
convention (berry.F90:66-73): `alpha_A=(2,3,1)`, `beta_A=(3,1,2)` → component 1=(y,z), 2=(z,x),
3=(x,y).

### B.3 Diagonalisation and gauge rotation; J± matrices

`wham_get_eig_UU_HH_JJlist` (wan_ham.F90:588-678): one Fourier call gives HH, delHH(x,y,z)
(lines 659-662), then `utility_diagonalize(HH, num_wann, eig, UU, ...)` (line 665) —
H(k) = UU · diag(eig) · UU†. Then per direction `wham_get_JJp_JJm_list` (wan_ham.F90:196-264):

```
call utility_rotate_new(delHH, UU, num_wann)            ! in-place: delHH ← UU† delHH UU (H gauge), line 236
for each Fermi level fe:
  if (eig(n) > fe .and. eig(m) < fe):                   ! n = empty, m = occupied
     JJp_list(n,m) = cmplx_i * delHH(n,m) / (eig(m) − eig(n))     ! line 251
     JJm_list(m,n) = cmplx_i * delHH(m,n) / (eig(n) − eig(m))     ! line 252
  else: 0
call utility_rotate_new(JJ±, UU, num_wann, reverse=.true.)        ! back to Wannier gauge: UU J UU†, lines 260-261
```

So J⁺ has nonzero (empty,occupied) blocks with denominator (E_occ − E_empty) (i.e. **E_m − E_n**
with n=row=empty, m=col=occupied), J⁻ = (J⁺)† structure with (E_n − E_m). These are the
D-matrix analogues; there is **no degeneracy threshold** in the AHC path — separation is purely
occupied-vs-empty at each fe (T=0 step: `pw90common_get_occ`, postw90_common.F90:942-985:
`occ=1 if eig<ef else 0`; the Fermi-Dirac branch is commented out).

For reference, the generic D-matrix routines (used by kubo/morb/sc, not by imf):
- `wham_get_D_h` (wan_ham.F90:102-142): `D_h(n,m,i) = (UU† dH_i UU)(n,m)/(eig(m) − eig(n))`,
  skipped when `n==m .or. abs(eig(m)−eig(n)) < 1.0e-7` (line 136) — sign: **denominator E_m − E_n**.
- `wham_get_deleig_a` (wan_ham.F90:342-439): band velocities dE/dk_a = Re diag(UU† dH_a UU)
  (line 435); with `use_degen_pert=T` (default F, `degen_thr=1e-4`, postw90_types.F90:102-103)
  degenerate subspaces (gap < degen_thr) are re-diagonalised (lines 392-422).
- `wham_get_eig_deleig` (wan_ham.F90:442-543): eig + del_eig wrapper (calls get_HH_R,
  fourier_R_to_k alpha=0..3, diagonalize, wham_get_deleig_a).

Occupation matrices `wham_get_occ_mat_list` (wan_ham.F90:267-339), in the **Wannier gauge**:

```
f_list(n,m,if) = sum_i UU(n,i) * occ_i(fe_if) * conjg(UU(m,i))    ! lines 330-331
g_list = 1 − f_list                                                ! lines 333-334
```

---

## C. berry_get_imf_klist / berry_get_imfgh_klist (src/postw90/berry.F90)

`berry_get_imf_klist` (berry.F90:1777-1870) is a thin wrapper forwarding to
`berry_get_imfgh_klist` (berry.F90:1873-2139) with only `imf_k_list` requested; optional `occ`
(fixed occupations) and `ladpt` (per-Fermi-level adaptive mask) pass through.

In `berry_get_imfgh_klist`: gather HH, UU, eig, JJ± (lines 1992-2018), f_list/g_list, then
A_W and Ω̄_W (lines 2020-2022). The **J0/J1/J2 decomposition** — trace formula for −2Im[f],
LVTS12 Eq.(51) (lines 2025-2049):

```
do i = 1, 3      ! axial component: 1=(y,z), 2=(z,x), 3=(x,y)  [alpha_A/beta_A, berry.F90:72-73]
  ! J0 term (Omega_bar term of WYSV06)
  imf_k_list(1,i,ife) =  Re Tr[ f · Ω̄_i ]                                      ! lines 2033-2034
  ! J1 term (D·A term of WYSV06)
  imf_k_list(2,i,ife) = −2 * ( Im Tr[ A_{alpha_A(i)} · J⁺_{beta_A(i)} ]
                             + Im Tr[ J⁻_{alpha_A(i)} · A_{beta_A(i)} ] )       ! lines 2037-2041
  ! J2 term (D·D term of WYSV06)
  imf_k_list(3,i,ife) = −2 *   Im Tr[ J⁻_{alpha_A(i)} · J⁺_{beta_A(i)} ]        ! lines 2044-2045
```

with `utility_re_tr_prod(a,b) = Re Σ_ij a(i,j)b(j,i)`, `utility_im_tr_prod = Im Σ_ij a(i,j)b(j,i)`
(`src/utility.F90:825-877`). All factors are in the Wannier gauge (f_list = U f U†, JJ± rotated
back). Everything is evaluated at T=0 step occupations for each entry of `fermi_energy_list`.

First index of `imf_k_list(1:3, 1:3, ife)` = J0/J1/J2 term; second index = axial component
x↔(y,z), y↔(z,x), z↔(x,y). So imf(:,3,:) is Ω_xy → σ_xy → AHC z-component.

`use_degen_pert` plays no role in imf (only in band-derivative/velocity routines).

---

## D. berry_main AHC assembly (src/postw90/berry.F90:96-1774)

- Task selection: `if (index(pw90_berry%task,'ahc') > 0) eval_ahc = .true.` (line 282).
- Setup: `get_HH_R` + `get_AA_R` (lines 291-308; `get_AA_R_effective` if effective_model).
- Requires at least one Fermi level: error `'Must specify one or more Fermi levels when
  berry=true'` if fermi_n==0 (lines 254-257). `fermi_energy_list` built in
  `w90_readwrite_read_fermi_energy` (`src/readwrite.F90:609-692`): single `fermi_energy`, or scan
  `fermi_energy_min/max/step` (default step 0.01, list linearly spaced, n = nint(|max−min|/step)+1).
- Cell volume computed as det(real_lattice) explicitly (lines 262-267).
- k-mesh: `berry_kmesh` (or `berry_kmesh_spacing`), read via `get_module_kmesh` with fallback to
  the global `kmesh`/`kmesh_spacing` (postw90_readwrite.F90:1902-1960). Mesh spacings
  `db_i = 1/mesh(i)` (lines 271-273). Main loop (no kpoint file, no tetrahedron; lines 841-903):

```
kweight = db1*db2*db3                                   ! line 843
loop_xyz = my_node_id, product(mesh)-1, num_nodes       ! flat MPI-strided loop, line 846
kpt = (loop_x*db1, loop_y*db2, loop_z*db3)              ! lines 852-854 (unshifted Γ-centred grid)
```

- **Adaptive refinement**: `berry_curv_adpt_kmesh` (integer, DEFAULT 1 = **off**;
  postw90_types.F90:166) and `berry_curv_adpt_kmesh_thresh` (DEFAULT 100.0; line 167), unit
  controlled by `berry_curv_unit` (default 'ang2'; line 168). Refinement sub-mesh (lines 617-626):

```
adkpt(1,ikpt) = db1*((i+0.5)/curv_adpt_kmesh − 0.5)   ! i=0..n−1, centred cluster spanning one coarse cell
kweight_adpt = kweight/curv_adpt_kmesh**3             ! line 844
```

  Trigger test per Fermi level (lines 867-880): vdum(j) = Σ_terms imf_k_list(:,j,if) (sum of
  J0+J1+J2 for each axial component), converted to bohr² if `curv_unit=='bohr2'` (line 872),
  `rdum = |vdum|`; if `rdum > curv_adpt_kmesh_thresh` the coarse point's contribution is
  **discarded** and replaced by the n³ refined points (each with kweight_adpt, evaluated only for
  the triggering Fermi levels via `ladpt`; lines 881-902). Otherwise
  `imf_list += imf_k_list*kweight` (line 878).
- MPI reduce: `comms_reduce(imf_list, 3*3*fermi_n, 'SUM')` (line 1191).
- **Unit conversion** (lines 1329-1361). At this point imf_list = (1/N)Σ_k Ω_αβ(k) ≈
  V_c ∫ dk/(2π)³ Ω(k) [Ω in Å²]; want σ_αβ = −(e²/ħ) ∫ dk/(2π)³ Ω(k) in S/cm:

```
fac = −1.0e8_dp * elem_charge_SI**2 / (hbar_SI * cell_volume)      ! line 1360
ahc_list(:,:,:) = imf_list(:,:,:) * fac                            ! line 1361
```

  with cell_volume in Å³, e = 1.602176634e-19 C, ħ = 1.054571817e-34 J·s
  (`src/constants.F90:108,112`). The 1e8 converts Å⁻¹ → cm⁻¹; sign is the **explicit (−1)** of
  σ = −(e²/ħ)∫Ω (comment lines 1342-1357).
- Output (lines 1373-1416): per Fermi level prints `'AHC (S/cm)       x          y          z'`
  (line 1395) then, at iprint≤1, one row `'=========='` with
  `sum(ahc_list(:,1,if)), sum(ahc_list(:,2,if)), sum(ahc_list(:,3,if))` — the J0+J1+J2 totals
  of σ for axial components x=σ_yz, y=σ_zx, z=σ_xy (lines 1411-1413). With `iprint>1` the J0/J1/J2
  rows are printed separately (lines 1396-1409). If fermi_n>1, also writes
  `seedname-ahc-fermiscan.dat` with `E_F, ahc_x, ahc_y, ahc_z` per line
  (lines 1369-1376, format `4(F12.6,1x)`).

---

## E. Fe AHC test case (test-suite/tests/testpostw90_fe_ahc/)

Directory contents:
- `Fe.win` (input), `Fe.eig`, `Fe.amn` (present but NOT read by postw90 for ahc)
- `Fe.chk.fmt.bz2` → symlink to `../../checkpoints/fe_postw90/Fe.chk.fmt.bz2`
- `Fe.mmn.bz2` → symlink to `../../checkpoints/fe_postw90/Fe.mmn.bz2`
- `Makefile`: `bunzip2 Fe.chk.fmt.bz2; w90chk2chk.x -f2u Fe` (formatted→unformatted chk) and
  `bunzip2 Fe.mmn.bz2 > Fe.mmn`
- `benchmark.out.default.inp=Fe.win` (reference .wpout)

Files actually required at runtime for `postw90.x Fe`: **Fe.win, Fe.chk (binary), Fe.eig,
Fe.mmn** (mmn read by get_AA_R, get_oper.F90:516; eig read via postw90 readwrite; chk gives
u_matrix/u_matrix_opt/lwindow/ndimwin/centres/kmesh info).

`Fe.win` key content (full file, lines cited from tests/testpostw90_fe_ahc/Fe.win):
- `num_bands = 28`, `num_wann = 18` (l.1-2), `spinors = true` (l.15), projections
  `Fe: sp3d2;dxy;dxz;dyz` (l.16-18)
- **`use_ws_distance = .false.`** (l.3) — overrides the postw90 default (.true.)
- `search_shells = 12` (l.4)
- disentanglement: `dis_win_min=-8.0 dis_win_max=70.0 dis_froz_min=-8.0 dis_froz_max=30.0` (l.6-9)
- **`fermi_energy = 12.6279`** (l.20)
- **`berry = true`, `berry_task = ahc`, `berry_kmesh = 10 10 10`** (l.24-26)
- bcc cell 2.71175 bohr half-lattice (l.43-48), `mp_grid = 2 2 2`, explicit 8 kpoints (l.54-65)
- NOT set (so defaults): `transl_inv` (F), `berry_curv_adpt_kmesh` (1 → benchmark prints
  "Adaptive refinement : none"), `berry_curv_adpt_kmesh_thresh` (100), `berry_curv_unit` (ang2)

Benchmark output (`benchmark.out.default.inp=Fe.win:234-239`):

```
 Interpolation grid: 10 10 10
 Fermi energy (ev):   12.6279
 AHC (S/cm)       x          y          z
 ==========     0.0334     0.0572  1222.1510
```

Harness: `tests/jobconfig:280-283` — `[testpostw90_fe_ahc/]`, `program = POSTW90_WPOUT_OK`,
`inputs_args = ('Fe.win','')`, `output = Fe.wpout`. Tolerances `tests/userconfig:110-115`:
`POSTW90_WPOUT_OK` parses the wpout and checks `ahc_x`, `ahc_y`, `ahc_z` each with
`(1.0e-3, 2.0e-3)` (relative, absolute).

Companion test `testpostw90_fe_ahc_adaptandfermi/` (jobconfig:285-289, program
POSTW90_FERMISCAN_OK, output `Fe-ahc-fermiscan.dat`, tolerances userconfig:125-131:
fermienergy (1e-6,5e-6), ahc_x/y/z (1e-4,2e-4)). Its Fe.win differs only by:
`fermi_energy_min = 11.6279`, `fermi_energy_max = 13.6279`, `fermi_energy_step = 0.2`,
`berry_curv_adpt_kmesh = 5`, `berry_curv_adpt_kmesh_thresh = 10`. First/benchmark row of
its fermiscan dat: `11.627900  50.845203  -50.619458  157.962756`; at E_F=12.6279:
`12.627900  41.926963  -41.913963  722.267035` (adaptive-refined values).

---

## Implementation checklist (condensed)

1. Read chk (+eig, mmn); build v = u_opt·u per k.
2. H(q)=v†·diag(eig,window)·v; A_α(q)= hermitize[ Σ_b i w_b b_α V†M_bV ] (diagonal Im-log only if
   transl_inv=T); FT with e^{−i2πq·R}/N_q onto WS irvec; divide by ndegen (ws_distance optional).
3. Per interpolation k: H(k), dH(k)=Σ iR_cart e^{i2πk·R}H(R); diagonalize → E,U;
   J± from (U†dHU)/(ΔE) with i factor, occupied/empty split at E_F, rotate back with U;
   f=U·occ·U†; A_W(k), Ω̄_W(k)=Σ e^{ikR} iR×A(R).
4. imf = [ReTr(fΩ̄_i), −2ImTr(A_a J⁺_b)−2ImTr(J⁻_a A_b), −2ImTr(J⁻_a J⁺_b)] for i↔(a,b) cyclic.
5. σ(S/cm) = −1e8 e²/(ħ V_c[Å³]) × (1/N_k)Σ_k Σ_terms imf; print x=σ_yz, y=σ_zx, z=σ_xy.
