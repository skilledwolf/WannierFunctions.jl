# postw90 reference notes: Gyrotropic module (gyrotropic = true)

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.

Key paper (quoted throughout `src/postw90/gyrotropic.F90:34-37`):
- TAS17 = S.S. Tsirkin, P. Aguado-Puente, I. Souza, arXiv:1710.03204 (2017) /
  PRB 97, 035158 (2018) — "Gyrotropic effects in trigonal tellurium studied from first principles".
- ZMS16 (comment at gyrotropic.F90:585) — Zhong, Moore, Souza, gyrotropic magnetic effect.
- WYSV06 Eq.(25) for the covariant Berry connection A + iD.

Module: `src/postw90/gyrotropic.F90` (module `w90_gyrotropic`, entry `gyrotropic_main`, called
from postw90.x when `gyrotropic = true`; flag default `.false.`, `src/postw90/postw90_types.F90:57`,
parsed `src/postw90/postw90_readwrite.F90:432-433`).

Related notes: berry-ahc.md (get_HH_R/get_AA_R, fourier R→k, JJ± machinery),
kubo-morb-geninterp.md §3 (get_BB_R/get_CC_R, imfgh trace formulas).

---

## A. Input keywords and defaults

Parsing in `w90_wannier90_readwrite_read_gyrotropic` (`src/postw90/postw90_readwrite.F90:698-818`)
and `w90_wannier90_readwrite_read_energy_range` (same file, 1638-1773). The whole .win is
lowercased on read (`src/readwrite.F90:2734`), so task strings are case-insensitive.

| keyword | default | units / notes |
|---|---|---|
| `gyrotropic` | `.false.` | activates module (postw90_types.F90:57) |
| `gyrotropic_task` | `'all'` | char(len=120), substring-matched (see B.1); type default postw90_types.F90:216 |
| `gyrotropic_kmesh` / `gyrotropic_kmesh_spacing` | fall back to global `kmesh`/`kmesh_spacing` | via `get_module_kmesh` (postw90_readwrite.F90:1883-1887, 1902-1998); **error** `'Error: gyrotropic module required, but no interpolation mesh given.'` if gyrotropic=T and neither local nor global mesh set. 1 int → replicated to 3. Mutually exclusive with `_kmesh_spacing`. |
| `gyrotropic_freq_min` | `0.0` | eV (pw90_extra_io default, postw90_readwrite.F90:53) |
| `gyrotropic_freq_max` | = default of kubo_freq_max: `froz_max − fermi_energy_list(1) + 0.6667` if frozen states, else `maxval(eigval) − minval(eigval) + 0.6667` | eV. **Copied from the kubo DEFAULT before the user's `kubo_freq_max` keyword is read** (lines 1688-1698), so setting kubo_freq_max does NOT change the gyro default. |
| `gyrotropic_freq_step` | `0.01` | eV (line 55). `nfreq = nint((max−min)/step)+1`; clamped to ≥2; step recomputed as `(max−min)/(nfreq−1)` (lines 1738-1744) |
| `gyrotropic_smr_fixed_en_width` | = global `smr_fixed_en_width` (default `0.0`) | eV; must be ≥ 0 (lines 797-804). **0 ⇒ division by zero in the delta; every real calculation must set it** |
| `gyrotropic_smr_type` | = global `smr_type` (default gauss, index 0) | index map (`src/readwrite.F90:1864-1914`): `m-v`/`cold` → −1, `m-pN` → N (bare `m-p` → 1), `f-d` → −99, `gauss` → 0 |
| `gyrotropic_smr_max_arg` | = `smr_max_arg` (default `5.0`) | dimensionless; must be > 0 (lines 779-795). Band/Fermi pair skipped when `|E−E_f|/η > max_arg` |
| `gyrotropic_degen_thresh` | `0.0` | eV (postw90_types.F90:221) |
| `gyrotropic_band_list` | `1..num_wann` | range vector (e.g. `1-4,7`); each entry must be in [1,num_wann] (lines 743-777) |
| `gyrotropic_box_b1/b2/b3` | rows of identity | reduced (fractional) coordinates of the box edge vectors; `box(i,:) = b_i` (lines 728-735; type default `box=0.0` then diag set to 1) |
| `gyrotropic_box_center` | not set → `box_corner=(0,0,0)` | reduced coords; if given, `box_corner = center − 0.5*(b1+b2+b3)` (lines 736-741) |
| `gyrotropic_eigval_max` | = default of kubo_eigval_max: `froz_max + 0.6667` if frozen states, else `maxval(eigval) + 0.6667` | eV; only used by NOA to cut off unoccupied bands (lines 1758-1772; note again: copied before the user's `kubo_eigval_max` is read) |
| `uhu_formatted` | `.false.` | read seedname.uHu as formatted text (postw90_types.F90:70, parsed postw90_readwrite.F90:482-483). Keyword spelled `uHu_formatted` in .win files (lowercased anyway) |
| `spn_formatted` | `.false.` | read seedname.spn formatted (postw90_types.F90:68) |
| `fermi_energy` or `fermi_energy_{min,max,step}` | REQUIRED | `'Must specify one or more Fermi levels when gyrotropic=true'` (gyrotropic.F90:171-174). List construction `src/readwrite.F90:609-692`: `n = nint(abs((max−min)/step))+1`, `step` recomputed `(max−min)/(n−1)`, `E_f(i) = min + (i−1)*step`; default step 0.01, default max = min+1.0 |
| `smearing use_adaptive` | forced `.false.` for gyrotropic (postw90_readwrite.F90:779) | belt-and-braces error `'Adaptive smearing not allowed in Gyrotropic'` in gyrotropic_get_k_list:690-693 |
| `transl_inv` | `.false.` | if `.true.` together with the K task: fatal `'transl_inv=T disabled for K-tensor'` (gyrotropic.F90:333-337; note the check sits inside the `iprint>0` root-only block) |
| `use_ws_distance` | `.true.` (postw90 default) | Te tests set `.false.`; affects all R→k sums as usual |

Frequency list is COMPLEX (postw90_readwrite.F90:1751-1756):

```
freq_list(i) = freq_min + (i-1)*(freq_max - freq_min)/(nfreq-1)  + cmplx_i * smearing%fixed_width
```

i.e. **Im(ω) = gyrotropic_smr_fixed_en_width** — the same η that smears the Fermi-surface delta
also broadens all frequency denominators of Dw and NOA (the read order guarantees this:
read_gyrotropic at postw90_readwrite.F90:209 runs before read_energy_range at 249).

---

## B. Task selection, k sampling, band selection

### B.1 Task flags (gyrotropic.F90:196-218)

`pw90_gyrotropic%task` is scanned with `index()` (already lowercase):

```
'-k'    → eval_K      (K tensor, Eq.3 TAS17: orbital + optional spin part)
'-c'    → eval_C      (C tensor, Eq.B6 TAS17)
'-d0'   → eval_D      (Berry-curvature dipole D, Eq.2 TAS17)
'-dw'   → eval_Dw     (frequency-dependent tildeD, Eq.12 TAS17)
'-spin' → eval_spn    (spin contributions to K and NOA)
'-noa'  → eval_NOA    (natural optical activity gamma, Eq.C12 TAS17)
'-dos'  → eval_DOS    (density of states)
'all'   → all of the above; eval_spn only if w90_system%spinors
```

`if (.not.(eval_K .or. eval_NOA)) eval_spn = .false.` (line 213). If `-spin` requested without
`spinors=true`: fatal `"spin contribution requested for gyrotropic, but the wavefunctions are
not spinors"` (215-218).

### B.2 Required matrix elements per task (gyrotropic.F90:220-263)

- always: `get_HH_R` (from chk + eig).
- `eval_D .or. eval_Dw .or. eval_K .or. eval_NOA`: `get_AA_R` (re-reads seedname.mmn), or
  `get_AA_R_effective` for effective_model.
- `eval_spn`: `get_SS_R` (reads seedname.spn; formatted if `spn_formatted`).
- `eval_K`: `get_BB_R` (mmn + eig; **fatal if |scissors_shift|>1e-7**, get_oper.F90:880-883) and
  `get_CC_R` (reads seedname.uHu; formatted if `uHu_formatted`). uHu contains
  `<u_{q+b1}|H_q|u_{q+b2}>` over num_bands; the code reads each (n,m) block then TRANSPOSES
  (`Ho_qb1_q_qb2 = transpose(...)`, get_oper.F90:1257-1274) and projects into the Wannier
  gauge WITHOUT eigenvalue weighting: `H_qb1_q_qb2 = V(qb1)† Ho V(qb2)`
  (`get_gauge_overlap_matrix`, positional `S` output, get_oper.F90:1286-1288, 3235-3272). Then
  (get_oper.F90:1296-1305, default transl_inv_full=F path)

  ```
  CC_q(:,:,ik,a,b) (a<=b) = sum_{b1,b2} wb(nn1)*bk(a,nn1,ik) * wb(nn2)*bk(b,nn2,ik) * H_qb1_q_qb2
  ```

  q→R Fourier + 1/ndegen as for H (kubo-morb-geninterp.md §3.2). BB per
  get_oper.F90:988-993: `BB_q_b(:,:,ik,nninv(nn,ik),idir) += i*wb(nn)*bk(idir,nn,ik)*H_q_qb`
  with `H_q_qb = V(k)† diag(eigval(k)) S_o V(k+b)` (i.e. B_a(k)=i<u|H|∂_a u>).
- Only `-c` and/or `-dos`: nothing beyond HH_R (no mmn/uHu/spn needed).

### B.3 k-point box sampling (gyrotropic.F90:355-381)

No `kpoint.dat` support (comment line 355) and no adaptive refinement. Flat MPI-strided loop
over `product(gyrotropic_kmesh)` points:

```
db_i    = 1/mesh(i)                                              ! lines 184-186
kweight = db1*db2*db3 * utility_det3(box)                        ! line 357
loop_x  = loop_xyz/(mesh(2)*mesh(3)); loop_y, loop_z analogous   ! lines 360-364
kpt     = (loop_x*db1, loop_y*db2, loop_z*db3)                   ! reduced coords in the box
kpt(:)  = box_corner(:) + matmul(kpt, box)                       ! line 368: k = corner + Σ_i f_i * box(i,:)
```

`utility_det3` (src/utility.F90:111-124) is a plain 3×3 determinant — a left-handed box gives a
NEGATIVE kweight (not guarded). Default box = identity, corner = 0 → unshifted Γ-based full-BZ
grid identical to the berry module. The k integral is normalised so that
`sum_k kweight = det(box)` (fraction of the BZ covered).

### B.4 Band list and degeneracy threshold (gyrotropic.F90:730-741)

Band loop runs over `band_list` (`n = band_list(n1)`). A band is SKIPPED entirely (cycle) when
it is quasi-degenerate with an adjacent band in the FULL band index space:

```
if (n > 1        .and. eig(n) - eig(n-1) <= degen_thresh) cycle
if (n < num_wann .and. eig(n+1) - eig(n) <= degen_thresh) cycle
```

Note `<=`: with the default `degen_thresh = 0.0` exact degeneracies are still skipped.
This threshold applies to the Fermi-surface quantities (K, D, Dw, C, DOS); the NOA part is
computed OUTSIDE this loop and does NOT use degen_thresh (it uses band_list + eigval_max only).

---

## C. Per-k evaluation (gyrotropic_get_k_list, gyrotropic.F90:574-849)

Common machinery per k-point:

1. `wham_get_eig_deleig` (wan_ham.F90:442-543) → `eig(:)`, `del_eig(:,3)` [eV·Å], `UU`, `HH`,
   `delHH` (see berry-ahc.md §B.3; `pw90_band_deriv_degen` defaults: use_degen_pert=F,
   degen_thr=1e-4).
2. If `eval_Dw .or. eval_NOA` (lines 710-722): covariant Berry connection, Eq.(25) WYSV06

   ```
   call wham_get_D_h(delHH, D_h, UU, eig, num_wann)   ! D_h(n,m,i) = (UU† dH_i UU)(n,m)/(eig(m)-eig(n)),
                                                      ! zero if n==m or |eig(m)-eig(n)| < 1.0e-7   (wan_ham.F90:102-142)
   AA(:,:,i) = utility_rotate(AA_W(:,:,i), UU, num_wann)   ! = UU† A_W(k) UU;  A_W from fourier_R_to_k_vec(AA_R)
   AA        = AA + cmplx_i * D_h
   ```

3. `eta_smr = gyrotropic_smr_fixed_en_width` (line 726). Per band n (band_list, degen filter) and
   per Fermi level `E_f(ifermi)`:

   ```
   arg = (eig(n) - E_f)/eta_smr ;   cycle if |arg| > smr_max_arg          ! lines 746-751
   delta = utility_w0gauss(arg, smr_type_index)/eta_smr * kweight         ! lines 805-806  (kweight folded in!)
   ```

   `utility_w0gauss` (src/utility.F90:1008-1091): Gaussian n=0 → `exp(-x²)/sqrt(pi)`
   (via Hermite recursion), n=−1 cold, n=−99 Fermi-Dirac `1/(2+e^x+e^-x)`, n≥1 Methfessel-Paxton.

4. Spin expectation (once per k, only when needed, lines 756-761): `spin_get_S`
   (src/postw90/spin.F90:387-452): diagonalises H(k) itself and returns
   `S(n,j) = Re[(UU† S_j(k) UU)(n,n)]` with `S_j(k) = Σ_R e^{ik·R} SS_R(:,:,:,j)`
   — dimensionless Pauli-matrix expectation values (factor ħ/2 NOT included).

5. Orbital/berry quantities per band n via FAKE OCCUPATIONS (`occ = 0; occ(n) = 1`, computed
   once per (k,n), lines 763-803):
   - `eval_K`: `berry_get_imfgh_klist(..., imf_k, img_k, imh_k, occ)` (berry.F90:1873-2139;
     trace formulas in berry-ahc.md §C and kubo-morb-geninterp.md §3.4; the occ branch uses
     `wham_get_JJp_JJm_list` occ logic, wan_ham.F90:241-249: nonzero only between occ(n)>0.5 and
     occ(m)<0.5 states, `JJm(n,m)=i·delHH̄(n,m)/(eig(m)−eig(n))`, `JJp(m,n)=i·delHH̄(m,n)/(eig(n)−eig(m))`,
     both rotated back to the W gauge; `f = UU·diag(occ)·UU†`). Then (lines 778-781)

     ```
     orb_nk(i)  = sum(imh_k(:,i,1)) - sum(img_k(:,i,1))      ! sum over the 3 J-terms (first index)
     curv_nk(i) = sum(imf_k(:,i,1))                          ! band-n Berry curvature, axial component i, Å²
     ```

   - `eval_D` only (no K): `berry_get_imf_klist(..., imf_k, ..., occ)` → `curv_nk` (lines 782-798).
   - `eval_Dw`: `gyrotropic_get_curv_w_k(eig, AA, curv_w_nk, ...)` (called at line 800; NOTE it
     fills curv_w_nk for ALL bands in band_list in one call, see C.1).

6. Accumulation (lines 812-828), all with the SAME broadened delta (i = row index runs over the
   velocity direction via array syntax `del_eig(n,:)`):

   ```
   gyro_K_spn(:,j,ifermi) += del_eig(n,:) * S(n,j)        * delta        ! [Å]        (only eval_K.and.eval_spn)
   gyro_K_orb(:,j,ifermi) += del_eig(n,:) * orb_nk(j)     * delta        ! [eV·Å³]
   gyro_D  (:,j,ifermi)   += del_eig(n,:) * curv_nk(j)    * delta        ! [Å³]
   gyro_Dw (i,j,ifermi,:) += del_eig(n,i) * delta * curv_w_nk(n,:,j)     ! [Å³]  per frequency
   gyro_C  (:,j,ifermi)   += del_eig(n,:) * del_eig(n,j)  * delta        ! [eV·Å²]  (comment says eV·Å³; actual eV·Å²... see F)
   gyro_DOS(ifermi)       += delta                                       ! [1/eV]
   ```

   First tensor index = velocity direction (dE/dk_i), second = curvature/spin/moment direction j.

### C.1 gyrotropic_get_curv_w_k (gyrotropic.F90:851-895) — tildeΩ(ω)

Band-resolved frequency-dependent Berry curvature (module-local; NOT shared with berry):

```
curv_w_k(n, iw, i) = - Σ_{m ∈ band_list, m≠n}  2 * Im[ AA(n,m,alpha_A(i)) * AA(m,n,beta_A(i)) ]
                                               * Re[ wmn² / (wmn² − freq_list(iw)²) ]
wmn = eig(m) − eig(n)
```

with `alpha_A=(2,3,1)`, `beta_A=(3,1,2)` (module constants, gyrotropic.F90:56-57) and COMPLEX
`freq_list` (ω + iη), so the real part implements Lorentzian broadening. At ω=0 the ratio → 1 and
tildeΩ reduces to the interband decomposition of Ω. AA is the covariant connection of step C.2.

### C.2 NOA: gyrotropic_get_NOA_k (gyrotropic.F90:897-1055)

Called ONCE per k after the band loop, when `eval_NOA` (with `SS_R`+`gyro_NOA_spn` optional args
when eval_spn, lines 833-847). If spin: `SS(:,:,j) = Σ_R e^{ikR} SS_R(:,:,:,j)` via
`pw90common_fourier_R_to_k_new`, then `S_h(:,:,j) = UU† SS_j UU` (lines 971-982).

Per Fermi level ifermi (lines 984-1011): partition band_list by

```
eig(n) <  E_f              → occupied  (occ_list, num_occ)
E_f <= eig(n) < eigval_max → unoccupied (unocc_list, num_unocc)
eig(n) >= eigval_max       → dropped entirely
```

If num_occ==0 or num_unocc==0: warning to stdout only when `iprint >= 2`, and the Fermi level is
skipped (cycle).

**B-matrices.** `gyrotropic_get_NOA_Bnl_orb` (1057-1104), for n∈occ, l∈unocc, a,c=1..3
[units eV·Å²]:

```
Bnl_orb(n1,l1,a,c) = -i*(del_eig(n,a) + del_eig(l,a)) * AA(n,l,c)
                     + Σ_{m ∈ band_list} [ (eig(n)-eig(m)) * AA(n,m,a) * AA(m,l,c)
                                          - (eig(l)-eig(m)) * AA(n,m,c) * AA(m,l,a) ]
```

`gyrotropic_get_NOA_Bnl_spin` (1106-1144) [dimensionless]: initialised to 0; for each b=1..3 with
`c = alpha_A(b), a = beta_A(b)` (so (a,b,c) is an EVEN permutation, e.g. b=1 → (a,c)=(3,2)):

```
Bnl_spin(n1,l1,a,c) = S_h(n,l,b) ;   then Bnl_spin = Bnl_spin * (-cmplx_i)
```

All other (a,c) entries stay zero (the comment `-i eps_{abc} <u_n|sigma_b|u_l>` refers only to
those cyclic entries; the anticyclic ones are NOT filled).

**Accumulation** (lines 1022-1050). For each occupied n, unoccupied l:

```
wln       = eig(l) − eig(n)
multW1(:) = 1 / (wln² − freq_list(:)²)                     ! complex (freq has +iη)
multWm(:) = Re(multW1) * kweight
multWe(:) = Re( −multW1*(2*wln²*multW1 + 1) ) * kweight    ! = −(3wln²−ω²)/(wln²−ω²)²  (broadened)

do ab = 1,3 ;  a = alpha_A(ab) ; b = beta_A(ab)            ! ab: 1=(y,z), 2=(z,x), 3=(x,y)
  do c = 1,3
    gyro_NOA_orb(ab,c,ifermi,:) += multWm * Re( AA(l,n,b)*Bnl_orb(n1,l1,a,c) − AA(l,n,a)*Bnl_orb(n1,l1,b,c) )
                                 + multWe * (del_eig(n,c) + del_eig(l,c)) * Im( AA(n,l,a)*AA(l,n,b) )
    gyro_NOA_spn(ab,c,ifermi,:) += multWm * Re( AA(l,n,b)*Bnl_spin(n1,l1,a,c) − AA(l,n,a)*Bnl_spin(n1,l1,b,c) )
  end do
end do
```

So the stored tensor `gyro_NOA_*(ab,c)` is γ_{abc} with (a,b) the antisymmetric pair encoded as
the axial index ab. Units before conversion: orb eV⁻¹·Å³ (comment line 523), spin eV⁻²·Å
(comment line 539).

MPI: all accumulators `comms_reduce(...,'SUM')` (lines 384-422); conversion+output on root only.

---

## D. Reused vs module-local routines

Reused as-is:
- `w90_get_oper`: `get_HH_R`, `get_AA_R`, `get_AA_R_effective`, `get_BB_R`, `get_CC_R`, `get_SS_R`.
- `w90_berry`: `berry_get_imf_klist`, `berry_get_imfgh_klist` (with the `occ` optional argument —
  single "Fermi index" 1, fake occupations; identical trace formulas as ahc/morb).
- `w90_wan_ham`: `wham_get_eig_deleig`, `wham_get_D_h` (and, through imfgh,
  `wham_get_eig_UU_HH_JJlist`, `wham_get_occ_mat_list`).
- `w90_spin`: `spin_get_S`.
- `w90_postw90_common`: `pw90common_fourier_R_to_k_vec` (A_W), `pw90common_fourier_R_to_k_new`
  (SS for NOA). (`pw90common_fourier_R_to_k_new_second_d` is imported at line 618 but unused.)
- `w90_utility`: `utility_rotate` (U†OU), `utility_w0gauss`, `utility_det3`, `utility_rotate_diag`
  (inside spin_get_S).

Module-local (do NOT exist in berry):
- `gyrotropic_get_k_list` — per-k dispatcher.
- `gyrotropic_get_curv_w_k` — tildeΩ(ω) (different from kubo's Kubo-formula curvature; sum
  restricted to band_list, real-part-of-complex-frequency regularisation).
- `gyrotropic_get_NOA_k`, `gyrotropic_get_NOA_Bnl_orb`, `gyrotropic_get_NOA_Bnl_spin`.
- `gyrotropic_outprint_tensor`, `gyrotropic_outprint_tensor_w` — all file output.
- Local copies of `alpha_A=(2,3,1)`, `beta_A=(3,1,2)` (gyrotropic.F90:56-57), same values as berry.

---

## E. Unit conversions (gyrotropic.F90:431-568, exact expressions)

`cell_volume` = explicit det(real_lattice) (lines 179-181), Å³. Constants from
`src/constants.F90`; the DEFAULT build is **CODATA2006** (constants.F90:92-96):
`elem_charge_SI = 1.602176487e-19` C, `elec_mass_SI = 9.10938215e-31` kg,
`hbar_SI = 1.054571628e-34` J·s, `eps0_SI = 8.854187817e-12` F/m. (CODATA2010/2018/2022 are
selectable by preprocessor define; values differ in the 7th digit — file-precision relevant.)

| quantity | Fortran factor (verbatim) | output units |
|---|---|---|
| K_spin | `fac = -1.0e20_dp*physics%elem_charge_SI*physics%hbar_SI/(2.*physics%elec_mass_SI*cell_volume)` (445-446) — note the LEADING MINUS (−g_s·e·ħ/(4m_e) ≈ −e·ħ/(2m_e)) | Ampere |
| K_orb | `fac = physics%elem_charge_SI**2/(2.*physics%hbar_SI*cell_volume)` (468) | Ampere |
| D | `fac = 1./cell_volume` (479) | dimensionless |
| tildeD (Dw) | `fac = 1./cell_volume` (490) | dimensionless |
| C | `fac = 1.0e+8_dp*physics%elem_charge_SI**2/(twopi*physics%hbar_SI*cell_volume)` (512) — i.e. (e/h)·e·1e8/V_c | Ampere/cm |
| NOA_orb | `fac = 1e+10_dp*physics%elem_charge_SI/(cell_volume*physics%eps0_SI)` (529) | Ang |
| NOA_spin | `fac = 1e+30_dp*physics%hbar_SI**2/(cell_volume*physics%eps0_SI*physics%elec_mass_SI)` (546) | Ang |
| DOS | `gyro_DOS = gyro_DOS/cell_volume` (562) | eV⁻¹·Å⁻³ |

(For C the code comment at 502-505 explains: raw accumulator eV·Å² per band → /V_c[Å³] → eV/Å,
×1e8·e → J/cm, ×e/h → A/cm; `twopi*hbar_SI = h`.)

---

## F. Output files (gyrotropic_outprint_tensor(_w), gyrotropic.F90:1146-1285)

Filename: `trim(seedname)//"-gyrotropic-"//trim(f_out_name)//".dat"` (line 1181). One file per
quantity; also `write(stdout,'(/,3x,a)') '* '//file_name` to the .wpout. `f_out_name` values:
`K_spin`, `K_orb`, `D`, `tildeD`, `C`, `NOA_orb`, `NOA_spin`, `DOS`.

Header lines 1 and 2 (list-directed `write(file_unit,*)`, hence ONE LEADING BLANK):

```
 #<comment>
 # in units of [ <units> ] 
```

comment strings (verbatim, lines 450/473/484/495/517/533/550/565):
- K_spin: `spin part of the K tensor -- Eq. 3 of TAS17`
- K_orb: `orbital part of the K tensor -- Eq. 3 of TAS17`
- D: `the D tensor -- Eq. 2 of TAS17`
- tildeD: `the tildeD tensor -- Eq. 12 of TAS17`
- C: `the C tensor -- Eq. B6 of TAS17`
- NOA_orb: `the tensor $gamma_{abc}^{orb}$ (Eq. C12,C14 of TAS17)`
- NOA_spin: `the tensor $gamma_{abc}^{spin}$ (Eq. C12,C15 of TAS17)`
- DOS: `density of states`

### F.1 Symmetrized 3×3 tensors — D, tildeD, C, K_orb, K_spin (symmetrize defaults to .true.)

Per frequency block (static quantities are a single block with omega=0.0; tildeD loops
`i=1..nfreq` passing `real(freq_list(i))`), two header lines:

```
write(file_unit,'(a1,29x,a1,38x,a14,37x,a2,14x,a15,14x,a1)') '#',"|","symmetric part","||","asymmetric part","|"
write(file_unit,'(11a15)') '# EFERMI(eV)', "omega(eV)", 'xx','yy','zz','xy','xz','yz','x','y','z'
```

then per Fermi level one row `write(file_unit,'(11E15.6)') efermi, omega, xx,yy,zz,xy,xz,yz,x,y,z`
with (T = gyro tensor, lines 1236-1247):

```
xx=T(1,1)  yy=T(2,2)  zz=T(3,3)
xy=(T(1,2)+T(2,1))/2  xz=(T(1,3)+T(3,1))/2  yz=(T(2,3)+T(3,2))/2
x =(T(2,3)-T(3,2))/2  y =(T(3,1)-T(1,3))/2  z =(T(1,2)-T(2,1))/2      ! antisymmetric part as polar vector
```

Each block is terminated by TWO blank lines (`write(file_unit,*)` twice, lines 1283-1284).

### F.2 NOA tensors — NOA_orb, NOA_spin (symmetrize=.false., lines 534-536, 551-553)

One header line per frequency block:

```
write(file_unit,'(11a15)') '# EFERMI(eV)', "omega(eV)", 'yzx','zxy','xyz','yzy','yzz','zxz','xyy','xyx','zxx'
```

then rows `'(11E15.6)'` with efermi, omega and, in order (lines 1249-1257, T(ab,c)=γ_{abc},
ab axial 1=(y,z),2=(z,x),3=(x,y)):

```
T(1,1)=γ_yzx  T(2,2)=γ_zxy  T(3,3)=γ_xyz  T(1,2)=γ_yzy  T(1,3)=γ_yzz
T(2,3)=γ_zxz  T(3,2)=γ_xyy  T(3,1)=γ_xyx  T(2,1)=γ_zxx
```

One block per frequency, each followed by two blank lines. NOTE: the shipped Te benchmarks were
generated by an older code revision whose header says `(Eq. C10 of TAS17)` and labels column 8
`yzz` (twice) instead of `xyx` — numbers are unchanged; the test parser skips `#` lines.

### F.3 DOS (arrEf1d path, lines 1276-1282)

```
write(file_unit,'(2a15)') '# EFERMI(eV) '            ! single item → "  # EFERMI(eV) "
write(file_unit,'(11E15.6)') fermi_energy_list(i), dos(i)     ! one row per Fermi level
```

followed by the same two blank lines. (Yes: the data format string is the same `'(11E15.6)'`
even though only 2 items are written.)

E15.6 renders like `   0.200000E+01` / `  -0.210217E+00`.

### F.4 stdout (.wpout) side output (gyrotropic.F90:302-353, 424-429)

`'Properties calculated in module  g y r o t r o p i c'`, per-task lines
(`'* D-tensor  --- Eq.2 of TAS17 '`, `'* K-tensor  --- Eq.3 of TAS17 '` +
`'    * including/excluding spin component '`, `'* Dw-tensor  --- Eq.12 of TAS17 '`,
`'* C-tensor  --- Eq.B6 of TAS17 '`, `'* gamma-tensor of NOA --- Eq.C12 of TAS17 '`,
`'* density of states '`), then `write(stdout,'(1x,a20,3(i0,1x))') 'Interpolation grid: ', mesh(1:3)`,
`'Calculation finished, writing results'`, and one `'* <file>'` line per written file.

---

## G. Te reference tests (test-suite/tests/testpostw90_te_gyrotropic*)

Seven directories; the `Te.win` files are IDENTICAL except for `gyrotropic_task`. Common
gyrotropic-relevant content (quoting testpostw90_te_gyrotropic/Te.win):

```
gyrotropic=true
gyrotropic_task=<per test, see below>
fermi_energy_step=2 ; fermi_energy_min=2 ; fermi_energy_max=10     → E_f = 2,4,6,8,10 eV (5 levels)
gyrotropic_freq_step=0.05 ; gyrotropic_freq_min=0.0 ; gyrotropic_freq_max=0.1   → ω = 0,0.05,0.1 (+0.1i)
gyrotropic_smr_fixed_en_width=0.1 ; gyrotropic_smr_max_arg=5
gyrotropic_degen_thresh=0.001
gyrotropic_box_b1=0.2 0.0 0.0 ; gyrotropic_box_b2=0.0 0.2 0.0 ; gyrotropic_box_b3=0.0 0.0 0.2
gyrotropic_box_center=0.33333 0.33333 0.5        → corner = center − (0.1,0.1,0.1); det(box)=0.008
gyrotropic_kmesh=5 5 5                           → kweight = (1/125)*0.008 = 6.4e-5
use_ws_distance = .false. ; search_shells=12
uHu_formatted=.true.
dis_win_min=-0.5 ; dis_win_max=10 ; dis_froz_min=0.0 ; dis_froz_max=8   → frozen → eigval_max default 8.6667
num_bands=12 ; num_wann=9 ; spinors = .false. ; mp_grid = 2 2 2 (8 explicit kpoints)
trigonal Te cell (a=4.457, c=5.9581176 Å), 3 Te atoms, p-projections
```

Smearing type not set → global default Gaussian (index 0). No `gyrotropic_band_list` → 1..9.
No `-spin` anywhere (spinors=false). Inputs: Te.win, Te.eig, Te.mmn, Te.uHu (bunzip2), Te.chk
(from Te.chk.fmt.bz2 via `w90chk2chk.x -f2u`); Te.amn unused by postw90.

Per-test task and benchmark-checked output (`tests/jobconfig:376-416`):

| test dir | gyrotropic_task | program | checked output |
|---|---|---|---|
| testpostw90_te_gyrotropic | `-C-dos-D0-Dw-K-NOA` | POSTW90_GYRO_OK | `Te-gyrotropic-C.dat` |
| testpostw90_te_gyrotropic_C | `-C` | POSTW90_GYRO_OK | `Te-gyrotropic-C.dat` |
| testpostw90_te_gyrotropic_D0 | `-D0` | POSTW90_GYRO_OK | `Te-gyrotropic-D.dat` |
| testpostw90_te_gyrotropic_Dw | `-Dw` | POSTW90_GYRO_OK | `Te-gyrotropic-tildeD.dat` |
| testpostw90_te_gyrotropic_K | `-K` | POSTW90_GYRO_OK | `Te-gyrotropic-K_orb.dat` |
| testpostw90_te_gyrotropic_NOA | `-NOA` | POSTW90_GYRO_OK | `Te-gyrotropic-NOA_orb.dat` |
| testpostw90_te_gyrotropic_dos | `-dos` | POSTW90_DOS_OK | `Te-gyrotropic-DOS.dat` |

(The combined test writes C, DOS, D, tildeD, K_orb, NOA_orb but only C.dat is compared.
K_spin/NOA_spin are never produced: no spinors.)

Tolerances (`tests/userconfig:155-169`, parser `tools/parsers/parse_gyro_dat.py` expects exactly
11 numeric columns, skips `#` lines and blanks): energy (1e-6, 5e-6); omega and gyro_xx..yz
(1e-4, 1e-2); gyro_x/y/z (1e-4, None). POSTW90_DOS_OK (userconfig:139-145): energy (1e-6, 5e-6),
dos (1e-4, 1e-4).

Benchmark first rows for sanity (E_f=2, ω=0):
- C:   `0.200000E+01 0.000000E+00 0.361959E+01 0.326754E+01 0.447568E+01 -0.210217E+00 -0.393533E+00 0.324436E+00 0. 0. 0.`
- D:   xx `0.472879E-02`; tildeD: xx `0.349534E-03`; K_orb: xx `-0.825797E-07`;
  NOA_orb: yzx `-0.443840E+02`; DOS: `0.277695E-03`.
  (The antisymmetric x,y,z columns of C, D, tildeD, K_orb are exactly 0 in these benchmarks —
  the accumulated tensors are symmetric… no: for D/tildeD/K they vanish because at these smearing
  settings the antisymmetric residue is below E15.6; C is exactly symmetric by construction.)

---

## H. Traps / gotchas

1. **kweight includes det(box)** and delta includes kweight — do not multiply twice.
2. **Complex frequencies**: ω_list = ω + i·gyrotropic_smr_fixed_en_width. Both Dw and NOA take
   `Re[...]` of complex-rational multipliers; a purely real-ω implementation will not match.
3. **η=0 default**: gyrotropic_smr_fixed_en_width inherits global smr_fixed_en_width default 0.0
   → `delta = w0gauss(arg)/0` = Inf. The reference gives no guard; tests always set it.
4. **degen_thresh uses `<=`** and compares against adjacent bands n±1 over ALL num_wann bands
   (not just band_list). Applies only to the Fermi-surface loop, NOT to NOA.
5. **Fake occupations** per band feed `berry_get_imfgh_klist(..., occ)`: occ-branch selects
   occ>0.5 vs occ<0.5 (no Fermi comparison), and only ONE fermi slot is computed.
6. **orb_nk = Σ(imh) − Σ(img)** summed over the three J-terms; curv_nk = Σ(imf). Signs and the
   −2Im[...] conventions live inside imfgh (berry-ahc.md §C, kubo-morb §3.4).
7. **K_spin factor has an explicit minus** (−1e20·e·ħ/(2m_e·V_c)); K_orb has +e²/(2ħV_c).
8. **CODATA2006** constants by default (e=1.602176487e-19, ħ=1.054571628e-34, m_e=9.10938215e-31,
   ε₀=8.854187817e-12) — different from the CODATA2022 numbers quoted in some other notes.
9. **uHu read order**: per (ik, nn2, nn1); each num_bands×num_bands block is TRANSPOSED after
   reading (`pw2wannier90 is coded a bit strangely`); the gauge projection uses NO eigenvalue
   weighting (the H is already inside uHu). BB uses eigval of the BRA k-point.
10. **transl_inv=T + K task is fatal**; transl_inv also changes AA_R diagonals for all other
    tasks (allowed, prints the usual message).
11. **eigval_max/freq_max defaults are captured from the kubo DEFAULTS before user kubo values
    are parsed** — user-set kubo_eigval_max/kubo_freq_max never leak into gyrotropic.
12. **NOA band cutoffs**: occupied = eig<E_f (regardless of eigval_max); unoccupied =
    E_f ≤ eig < eigval_max; empty lists → level silently skipped (warning only at iprint≥2).
13. **Bnl_spin fills only even-permutation entries** (a,c)=(beta_A(b),alpha_A(b)); everything
    else remains zero, and the antisymmetrized combination in the accumulation relies on that.
14. **Output rounding**: all data E15.6; header lines have a LEADING SPACE for the two
    list-directed comment/units lines but none for the format-written column headers.
15. **The task match is plain substring** on the lowercased string; `-d0` ≠ `-dos` ≠ `-dw`;
    `all` enables everything (+spin only if spinors).
16. **Adaptive smearing forbidden**, no adaptive kmesh refinement, no kpoint.dat path —
    unlike berry_main.
17. `use_ws_distance` default .true. in postw90, but all Te gyro tests set `.false.`.
18. Fermi levels REQUIRED; fermi list from min/max/step with `n=nint(abs((max−min)/step))+1`
    and step recomputed — reproduce exactly (fermi values are compared at 1e-6/5e-6).

---

## Implementation checklist (condensed)

1. Parse keywords (§A); build fermi list, freq list (complex, +iη), band_list, box.
2. Operators: HH_R always; AA_R if D/Dw/K/NOA; BB_R+CC_R (uHu) if K; SS_R (spn) if -spin.
3. Loop k over box grid, kweight = db1 db2 db3 det(box); per k: eig/del_eig/UU;
   AA = U†A_WU + iD_h if Dw/NOA.
4. Fermi-surface tasks: per band (degen filter), per E_f (max_arg cut),
   delta = w0gauss((E−E_f)/η)/η·kweight; accumulate D/Dw/C/DOS/K per §C.6 using per-band fake-occ
   imf/img/imh and spin_get_S; tildeΩ per §C.1.
5. NOA per §C.2 with occ/unocc split per E_f, Bnl_orb/Bnl_spin, multWm/multWe.
6. MPI-sum; scale by §E factors; write §F files (E15.6, exact headers, 2 blank lines per block).
