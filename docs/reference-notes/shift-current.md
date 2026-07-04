# postw90 reference notes: Shift current (berry_task = sc)

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.
Constants: the reference build defaults to **CODATA2006** (`src/constants.F90:92-103`).

Key papers referenced by the code (`src/postw90/berry.F90:38-46`, `berry.F90:2371-2386`,
`berry.F90:2626`):
- IATS18 = Ibañez-Azpiroz, Tsirkin, Souza, PRB 97, 245143 (2018) — shift current by Wannier
  interpolation (Eqs. 8, 30, 32, 34 implemented)
- SS00 = Sipe & Shkrebtii, PRB 61, 5337 (2000) — Eq. (57), the σ_abc definition
- WYSV06 = PRB 74, 195118 (2006) — Eq. (24), the D matrix
- Eta correction: Eq. (19) of PRB 103, 247101 (2021) (Comment on IATS18)

Notation correspondence quoted from the code header (`berry.F90:2375-2384`):

```
AA_da_bar              <-->   \mathbbm{b}
AA_bar                 <-->   \mathbbm{a}
HH_da_bar              <-->   \mathbbm{v}
HH_dadb_bar            <-->   \mathbbm{w}
D_h(n,m)               <-->   \mathbbm{v}_{nm} * Re[1/(E_{m}-E_{n}+i*sc_eta)]
D_h_no_eta(n,m)        <-->   \mathbbm{v}_{nm} / (E_{m}-E_{n})
sum_AD                 <-->   summatory of Eq. 32 IATS18
sum_HD                 <-->   summatory of Eq. 30 IATS18
eig_da(n)-eig_da(m)    <-->   \mathbbm{Delta}_{nm}
```

---

## A. Dispatch, setup, k-loop, reduction (src/postw90/berry.F90)

- Dispatch: `if (index(pw90_berry%task, 'sc') > 0) eval_sc = .true.` (berry.F90:285).
  Substring match on `berry_task` (parse at postw90_readwrite.F90:865-879; berry=T with an
  unrecognised task errors out).
- Requires ≥1 Fermi level: `'Must specify one or more Fermi levels when berry=true'` if
  fermi_n==0 (berry.F90:254-257). **sc is NOT in the `not_scannable` list** (berry.F90:342 covers
  only kubo and shc-freq-scan), so a Fermi *scan* is accepted — but `berry_get_sc_klist` uses
  **only `fermi_energy_list(1)`** (berry.F90:2531). Trap: extra Fermi levels silently ignored.
- Setup (berry.F90:392-411): `get_HH_R` + `get_AA_R` (or `get_AA_R_effective` for
  effective_model) — the sc task needs **HH_R and AA_R** (i.e. seedname.chk + .eig + .mmn), no
  BB/CC/spin operators. Arrays:

```
allocate (sc_k_list(3, 6, pw90_berry%kubo_nfreq));  sc_k_list = 0     ! lines 407-410
allocate (sc_list  (3, 6, pw90_berry%kubo_nfreq));  sc_list   = 0
```

  First index a = 1..3 (direction of the generalised derivative = DC current direction), second
  index jk = 1..6 packs the symmetric light-polarisation pair (b,c) via the module constants
  (berry.F90:76-89):

```
integer, dimension(6), parameter :: alpha_S = (/1, 2, 3, 1, 1, 2/)   ! b   (line 84)
integer, dimension(6), parameter :: beta_S  = (/1, 2, 3, 2, 3, 3/)   ! c   (line 85)
! jk: 1<->xx 2<->yy 3<->zz 4<->xy 5<->xz 6<->yz  (comment lines 77-82)
```

- Main k-loop, regular-grid branch (berry.F90:841-854, identical "CODE BLOCK 1" copy in the
  wanint_kpoint_file branch at 640-651/746-756): Γ-centred unshifted grid
  `kpt = (loop_x/N1, loop_y/N2, loop_z/N3)` over `berry_kmesh`, `kweight = 1/(N1*N2*N3)`,
  MPI-strided flat loop. Accumulation (berry.F90:953-964, kpoint-file copy 746-756):

```
call berry_get_sc_klist(...)                 ! per-k integrand
sc_list = sc_list + sc_k_list*kweight        ! line 963 (755)
```

  **No adaptive kmesh refinement for sc** (that machinery only wraps ahc/shc).
- MPI reduce: `comms_reduce(sc_list(1,1,1), 3*6*kubo_nfreq, 'SUM')` (berry.F90:1223-1226).
- Cell volume = det(real_lattice) expanded explicitly, **signed, no abs** (berry.F90:262-267).

### A.1 R→k machinery and use_ws_distance in this tree

In this reference tree ws_distance is pre-folded into the R-space operators:
`operator_wigner_setup` (get_oper.F90:3275-3327) rescatters `HH_R`/`AA_R` onto the union grid
`irvec_pw90/crvec_pw90` with weights `1/(ndegen*ndeg)` (see berry-ahc.md §A.1/A.2), and
`wigner_seitz_opt_setup` (postw90_common.F90:1797-1968) builds that union list. Consequently the
Fourier routines used by sc contain **no ws_distance branches** — with `use_ws_distance = .true.`
(the postw90 default, `src/types.F90:87-94`) the effect enters purely through the expanded
R-list; with `.false.` the plain WS list divided by ndegen is used.

---

## B. Keywords and defaults

| keyword | default | where |
|---|---|---|
| `sc_phase_conv` | **1** (TB convention) | default postw90_types.F90:170; parse postw90_readwrite.F90:940-946; must be 1 or 2 else error `'Error: sc_phase_conv must be either 1 or 2'` |
| `sc_eta` | **0.04** (eV) | default postw90_types.F90:171; parse postw90_readwrite.F90:962-964 |
| `sc_w_thr` | **5.0d0** (dimensionless, in units of eta_smr) | default postw90_types.F90:172; parse postw90_readwrite.F90:966-968 |
| `sc_use_eta_corr` | **.true.** | default postw90_types.F90:173; parse postw90_readwrite.F90:948-950 |

Reused kubo_* machinery (all shared with berry_task=kubo):

- Frequency list (postw90_readwrite.F90:1684-1724): `kubo_freq_min` (default 0.0),
  `kubo_freq_max` (default: `froz_max − fermi_energy_list(1) + 0.6667` if frozen window set,
  else `maxval(eigval) − minval(eigval) + 0.6667`, else `win_max − win_min + 0.6667`; lines
  1688-1694), `kubo_freq_step` (default 0.01, must be >0). Then:

```
kubo_nfreq = nint((max-min)/step) + 1;  if (kubo_nfreq <= 1) kubo_nfreq = 2   ! lines 1708-1710
step is recomputed = (max-min)/(nfreq-1)                                      ! lines 1711-1712
kubo_freq_list(i) = min + (i-1)*(max-min)/(nfreq-1)                           ! lines 1720-1724
```

  The list is stored complex but built purely real; sc uses `real(kubo_freq_list(:),dp)` only
  (berry.F90:2558). Trap: nint rounding — e.g. min=0, max=10, step=0.03 → nfreq=334 and actual
  step 10/333 = 0.030030….
- `kubo_eigval_max` (postw90_readwrite.F90:1758-1769): default `froz_max + 0.6667` if
  `dis_froz_max` was given (frozen_states flag, readwrite.F90:823-827), else
  `maxval(eigval) + 0.6667`, else `win_max + 0.6667`.
- Smearing (`kubo_smearing`): global defaults `pw90_smearing_type` (postw90_types.F90:129-141):
  `use_adaptive = .true.`, `adaptive_prefactor = sqrt(2)`, `type_index = 0` (Gaussian),
  `fixed_width = 0.0`, `adaptive_max_width = 1.0` eV. Global keywords `smr_type`, `adpt_smr`,
  `adpt_smr_fac`, `adpt_smr_max`, `smr_fixed_en_width` (postw90_readwrite.F90:580-636); kubo
  overrides `kubo_adpt_smr`, `kubo_adpt_smr_fac`, `kubo_adpt_smr_max`,
  `kubo_smr_fixed_en_width`, `kubo_smr_type` (postw90_readwrite.F90:908-960). Smearing-name →
  index map (readwrite.F90:1864-1914): 'gauss'→0, 'm-p'/'m-pN'→N, 'm-v'/'cold'→−1, 'f-d'→−99.
  **Trap: the vectorised delta used by sc (`utility_w0gauss_vec`, utility.F90:1093-1164) only
  implements Gaussian (n=0)** and raises an error for any other index at runtime.
- `berry_kmesh` / `berry_kmesh_spacing`: via `get_module_kmesh` with prefix 'berry'
  (postw90_readwrite.F90:1878-1881, routine 1902-1998); falls back to global
  `kmesh`/`kmesh_spacing`; setting both errors; 1-integer form expands to cubic.
- wpout parameter block (postw90_readwrite.F90:2361-2370), printed when task contains 'sc':

```
'|  Smearing factor for shift current         :'  sc_eta   (f8.3)
'|  Frequency theshold for shift current      :'  sc_w_thr (f8.3)   [sic: "theshold"]
'|  Bloch sums                                :'  'Tight-binding convention' | 'Wannier90 convention'
'|  Finite eta correction for shift current   :'  T/F (L8)
```

  Convention strings from `w90_readwrite_get_convention_type` (readwrite.F90:1847-1862).

---

## C. berry_get_sc_klist (berry.F90:2365-2681) — the per-k integrand

### C.1 Gather W-gauge matrices — phase conventions

**sc_phase_conv = 1 (TB convention, default)** (berry.F90:2481-2505):
- `wham_get_eig_UU_HH_AA_sc_TB_conv` (wan_ham.F90:681-768): calls get_HH_R + get_AA_R, then
  `pw90common_fourier_R_to_k_new_second_d_TB_conv` (postw90_common.F90:1235-1336) giving, with
  `τ_ij = τ_j − τ_i` (τ_i = Wannier centres, see below) and per-(i,j) phases:

```
r_sum   = irvec_pw90(:,ir) + τ_j^frac − τ_i^frac              ! fractional, line 1307-1308
phase   = exp(i·2π k·r_sum)                                    ! line 1309-1310
HH(i,j)          = Σ_R phase · HH_R(i,j,ir)                                     ! line 1311
HH_da(i,j,a)     = Σ_R  i·(crvec_pw90(a,ir)+τ_ja−τ_ia) · phase · HH_R(i,j,ir)   ! 1312-1318
HH_dadb(i,j,a,b) = Σ_R −(crvec+τ)_a·(crvec+τ)_b · phase · HH_R(i,j,ir)          ! 1319-1331
```

  (crvec/τ Cartesian Å; k·R phase from fractional coordinates.) Then
  `utility_diagonalize(HH,…,eig,UU)` (wan_ham.F90:765) — the TB-convention UU differs from
  convention 2 by k-dependent diagonal phases, eig identical.
- Position operator: `pw90common_fourier_R_to_k_vec_dadb_TB_conv` (postw90_common.F90:1473-1615)
  called with `OO_da=AA, OO_dadb=AA_da` (berry.F90:2496-2499):

```
AA(i,j,c)      = Σ_R phase · AA_R(i,j,ir,c)                                ! lines 1580-1582
AA_da(i,j,c,a) = Σ_R i·(crvec_pw90(a,ir)+τ_ja−τ_ia) · phase · AA_R(i,j,ir,c)  ! 1600-1608
```

  **Diagonal subtraction**: for the R=0, i=j term only, `AA_R(i,i,0,c)` is replaced by
  `AA_R(i,i,0,c) − τ_ic` in both AA and AA_da (lines 1566-1573 and 1587-1597) — the WF's own
  centre is removed ("matrix element is zero in this convention").
- Band derivatives: `wham_get_eig_deleig_TB_conv` (wan_ham.F90:546-585) = three calls to
  `wham_get_deleig_a` on the TB-convention HH_da with the TB UU.
- The centres τ are `wigner_seitz%wannier_centres_from_AA_R`, accumulated inside get_AA_R via
  the Im-log formula `τ_i −= w_b b · Im ln S(i,i) / N_q` (get_oper.F90:596-602) — independent of
  `transl_inv`. Consistency check (get_oper.F90:639-647): if
  `Σ(τ_from_AA_R − chk centres)² > 1e-8` → fatal `'Computed and read Wannier centres
  different.'` unless `guiding_centres=T` (then chk centres are used instead).

**sc_phase_conv = 2 (usual W90 convention)** (berry.F90:2506-2527):
- `wham_get_eig_UU_HH_AA_sc` (wan_ham.F90:770-839) →
  `pw90common_fourier_R_to_k_new_second_d` (postw90_common.F90:1161-1232):

```
rdotk = 2π k·irvec_pw90(:,ir);  phase = exp(+i rdotk)
HH        = Σ_R phase·HH_R;      HH_da(…,a) = Σ_R i·crvec_pw90(a,ir)·phase·HH_R
HH_dadb(…,a,b) = Σ_R −crvec(a)·crvec(b)·phase·HH_R                         ! lines 1212-1229
```

- `pw90common_fourier_R_to_k_vec_dadb` (postw90_common.F90:1406-1470):
  `AA(…,c) = Σ_R phase·AA_R(…,ir,c)`, `AA_da(…,c,a) = Σ_R i·crvec(a)·phase·AA_R(…,ir,c)`
  (lines 1451-1468). **No diagonal subtraction.**
- `wham_get_eig_deleig` (wan_ham.F90:442-543): recomputes HH/delHH via
  `pw90common_fourier_R_to_k` alpha=0..3 and diagonalises again, then `wham_get_deleig_a`.

Index convention in both cases: `AA_da(:,:,c,a)` = ∂_a A_c(k) (position component c, derivative
direction a; comment berry.F90:2495); `HH_dadb(:,:,c,a)` = ∂_c∂_a H (symmetric).

### C.2 Occupations, D matrices, gauge rotation

- `pw90common_get_occ(fermi_energy_list(1), eig, occ, num_wann)` (berry.F90:2531): T=0 step,
  `occ(i)=1 if eig(i) < ef else 0` (strict `<`; postw90_common.F90:969-972; Fermi-Dirac branch
  commented out).
- `wham_get_D_h_P_value` → **D_h** (wan_ham.F90:145-193): with `ΔE = eig(m) − eig(n)`,

```
D_h(n,m,a) = (U†·∂_aH·U)(n,m) · ΔE/(ΔE² + sc_eta²)      ! line 188, zero for n==m (2186)
```

  i.e. v_nm·Re[1/(E_m − E_n + i·sc_eta)] — principal-value regularisation with **sc_eta**, no
  degeneracy threshold.
- `wham_get_D_h` → **D_h_no_eta** (wan_ham.F90:102-142):

```
D_h_no_eta(n,m,a) = (U†·∂_aH·U)(n,m)/(eig(m) − eig(n))   ! line 137
skipped (left 0) when n==m .or. abs(eig(m)-eig(n)) < 1.0e-7   ! line 136
```

- Rotation to H (Hamiltonian) gauge (berry.F90:2544-2555), `utility_rotate(M,U) = U†·M·U`
  (utility.F90:699-716):

```
AA_bar(:,:,a)        = U†·AA(:,:,a)·U
HH_da_bar(:,:,a)     = U†·HH_da(:,:,a)·U
AA_da_bar(:,:,a,b)   = U†·AA_da(:,:,a,b)·U
HH_dadb_bar(:,:,a,b) = U†·HH_dadb(:,:,a,b)·U
```

- Band derivatives `eig_da(n,a)` from `wham_get_deleig_a` (wan_ham.F90:342-439):
  `dE_n/dk_a = Re[(U†∂_aH U)(n,n)]`; only with `use_degen_pert=T` (default **F**, degen_thr
  1e-4, postw90_types.F90:102-103) degenerate subspaces are re-diagonalised.

### C.3 Frequency grid and adaptive smearing

- `omega = real(kubo_freq_list)`, `wmin = omega(1)`, `wmax = omega(nfreq)`,
  `wstep = omega(2) − omega(1)` (berry.F90:2557-2561) — **uniform grid assumed**.
- If `kubo_smearing%use_adaptive` (berry.F90:2538-2541, 2574-2581):

```
Delta_k = max_i |b_i| / berry_kmesh(i)            ! kmesh_spacing_mesh, postw90_common.F90:1024-1027
vdum    = eig_da(m,:) − eig_da(n,:)
eta_smr = min( |vdum|·Delta_k · adaptive_prefactor , adaptive_max_width )   ! YWVS07 Eq.(34-35)
```

  else `eta_smr = kubo_smearing%fixed_width` (berry.F90:2580). No lower bound — parallel bands
  give eta_smr→0 (1/eta blowup guarded only by the occupation filter in practice).

### C.4 Band-pair loop, cutoffs

Double loop over **ordered** pairs (n,m), n,m = 1..num_wann (berry.F90:2564-2588):

```
if (n == m) cycle
if (eig(m) > kubo_eigval_max .or. eig(n) > kubo_eigval_max) cycle    ! line 2568
occ_fac = occ(n) − occ(m);  if (abs(occ_fac) < 1e-10) cycle          ! lines 2570-2571
! w_thr window: skip pair if BOTH E_nm and E_mn lie outside [wmin,wmax] padded by ±sc_w_thr*eta_smr
if (((eig(n)-eig(m)+sc_w_thr*eta_smr < wmin) .or. (eig(n)-eig(m)-sc_w_thr*eta_smr > wmax)) .and. &
    ((eig(m)-eig(n)+sc_w_thr*eta_smr < wmin) .or. (eig(m)-eig(n)-sc_w_thr*eta_smr > wmax))) cycle
```

Each unordered pair {v,c} is visited twice (as (v,c) and (c,v)) and each visit fills **two**
delta branches (±E_nm), see C.7.

### C.5 Intermediate-state sums (IATS18 Eqs. 30/32)

(berry.F90:2592-2604; `utility_zdotu(a,b) = Σ_p a(p)·b(p)`, **no conjugation**,
utility.F90:181-186.) For all (c,a) ∈ 3×3:

```
sum_AD(c,a) = [ Σ_p AA_bar(n,p,c)·D_h(p,m,a)  −  AA_bar(n,n,c)·D_h(n,m,a) ]
            − [ Σ_p D_h(n,p,a)·AA_bar(p,m,c)  −  D_h(n,m,a)·AA_bar(m,m,c) ]
sum_HD(c,a) = [ Σ_p HH_da_bar(n,p,c)·D_h(p,m,a) − HH_da_bar(n,n,c)·D_h(n,m,a) ]
            − [ Σ_p D_h(n,p,a)·HH_da_bar(p,m,c) − D_h(n,m,a)·HH_da_bar(m,m,c) ]
```

The explicit subtractions remove the p=n term of the first sum and the p=m term of the second
(diagonals of D_h are zero already), so effectively p ∉ {n,m}. Note D_h here is the
**η-regularised** one (P-value).

### C.6 Dipole and generalised derivative (IATS18 Eq. 34 + 30 + 32)

Dipole matrix element (berry.F90:2607) — note the **(m,n)** ordering:

```
r_mn(b) = AA_bar(m,n,b) + i·D_h_no_eta(m,n,b)
```

Generalised derivative, for derivative direction a and each Cartesian component c
(berry.F90:2610-2623; the division by (E_m − E_n) applies ONLY to the last i·(…) group —
Fortran left-to-right `*`/`/` precedence):

```
gen_r_nm(c) =  AA_da_bar(n,m,c,a)
             + (AA_bar(n,n,c) − AA_bar(m,m,c))·D_h_no_eta(n,m,a)
             + (AA_bar(n,n,a) − AA_bar(m,m,a))·D_h_no_eta(n,m,c)
             − i·AA_bar(n,m,c)·(AA_bar(n,n,a) − AA_bar(m,m,a))
             + sum_AD(c,a)
             + i·[ HH_dadb_bar(n,m,c,a) + sum_HD(c,a)
                   + D_h_no_eta(n,m,c)·(eig_da(n,a) − eig_da(m,a))
                   + D_h_no_eta(n,m,a)·(eig_da(n,c) − eig_da(m,c)) ] / (eig(m) − eig(n))
```

**Finite-eta correction** (berry.F90:2625-2641), applied iff `sc_use_eta_corr` (default T);
η = sc_eta; for every p with p ≠ n and p ≠ m:

```
gen_r_nm(c) −=  η²/((eig(p)−eig(m))² + η²) / (eig(n)−eig(m))
                · ( AA_bar(n,p,c)·HH_da_bar(p,m,a)
                  − (HH_da_bar(n,p,c) + i·(eig(n)−eig(p))·AA_bar(n,p,c))·AA_bar(p,m,a) )
gen_r_nm(c) +=  η²/((eig(n)−eig(p))² + η²) / (eig(n)−eig(m))
                · ( HH_da_bar(n,p,a)·AA_bar(p,m,c)
                  − AA_bar(n,p,a)·(HH_da_bar(p,m,c) + i·(eig(p)−eig(m))·AA_bar(p,m,c)) )
```

(Note the asymmetry: in the first (subtracted) term the free component c sits on the FIRST
factor and direction a on the second; in the second (added) term a sits on the first factor.)

Matrix element, symmetrised over (b,c) (berry.F90:2646-2650):

```
do bc = 1, 6;  b = alpha_S(bc); c = beta_S(bc)
  I_nm(a,bc) = aimag( r_mn(b)·gen_r_nm(c) + r_mn(c)·gen_r_nm(b) )   ! Im[r^b_mn r^c_nm;a + (b<->c)]
```

### C.7 Delta functions and accumulation

Two branches per ordered pair (berry.F90:2653-2676). Branch 1, ω ≈ E_nm = eig(n)−eig(m):

```
istart = max( int((eig(n)−eig(m) − sc_w_thr*eta_smr − wmin)/wstep + 1), 1 )        ! int() truncates
iend   = min( int((eig(n)−eig(m) + sc_w_thr*eta_smr − wmin)/wstep + 1), kubo_nfreq )
if (istart <= iend):
  delta(i) = utility_w0gauss( (eig(m)−eig(n)+omega(i))/eta_smr, type_index ) / eta_smr
           = exp(−((omega(i)−E_nm)/eta_smr)²) / (sqrt(pi)·eta_smr)        ! Gaussian, type 0
  DGER(18, iend−istart+1, occ_fac, I_nm, 1, delta(istart:iend), 1, sc_k_list(:,:,istart:iend), 18)
```

i.e. `sc_k_list(a,bc,i) += occ_fac · I_nm(a,bc) · delta(i)` (I_nm(3,6) flattened to 18). Branch
2 is identical with E_mn = eig(m)−eig(n) windows and argument `(eig(n)−eig(m)+omega(i))/eta_smr`
(berry.F90:2667-2676). The Gaussian being even, each branch is δ(E_nm − ω) resp. δ(E_mn − ω).
Two deltas per visit × two ordered visits per pair — this is why the prefactor carries **1/4
instead of SS00's 1/2** (comment berry.F90:1626-1627).

---

## D. Unit conversion and output (berry.F90:1608-1669)

Comment block (berry.F90:1613-1642): at this point
`sc_list = (1/N) Σ_k (r^b r^c_a + r^c r^b_a)(k) δ(w)` ≈ `V_c ∫ dk/(2π)³ …` and we want
`σ_abc = (π e³/(4 ħ²)) ∫ dk/(2π)³ Im[…] δ(w)` in A/V². Prefactor, quoted verbatim
(berry.F90:1638 and 1644):

```
! fac = eV_seconds.( pi.e^3/(4.hbar^2.V_c) )
fac = physics%eV_seconds*pi*physics%elem_charge_SI**3/(4*physics%hbar_SI**(2)*cell_volume)
```

with `V_c` in Å³ (signed det, berry.F90:262-267) and, for **CODATA2006**
(constants.F90:164-182):

```
elem_charge_SI = 1.602176487e-19   ! C
hbar_SI        = 1.054571628e-34   ! J*s
eV_seconds     = 6.582119e-16      ! (only 7 significant digits in 2006 set!)
bohr           = 0.52917720859     ! Å (bohr_angstrom_internal) — used for bohr cells
```

(CODATA2018/2022 have e=1.602176634e-19, ħ=1.054571817e-34, eV_seconds=6.582119569e-16 —
~1e-7-level relative differences that matter for file-precision matching.)

Output (root only, inside `if (print_output%iprint > 0)`), stdout banner (berry.F90:1645-1651):

```
 ----------------------------------------------------------
 Output data files related to shift current:               
 ----------------------------------------------------------
```

then per file `write (stdout,'(/,3x,a)') '* '//file_name`. File loop (berry.F90:1652-1667):

```
do i = 1, 3
  do jk = 1, 6
    j = alpha_S(jk);  k = beta_S(jk)
    file_name = trim(seedname)//'-sc_'//achar(119+i)//achar(119+j)//achar(119+k)//'.dat'
    do ifreq = 1, kubo_nfreq
      write (file_unit, '(2E18.8E3)') real(kubo_freq_list(ifreq),dp), fac*sc_list(i,jk,ifreq)
```

- **18 files always**, named `seedname-sc_<a><b><c>.dat` with a ∈ {x,y,z} the generalised
  derivative (current) direction and (b,c) the 6 symmetric pairs in order xx,yy,zz,xy,xz,yz:
  `-sc_xxx, -sc_xyy, -sc_xzz, -sc_xxy, -sc_xxz, -sc_xyz, -sc_yxx, …, -sc_zyz`.
- **No header lines**; kubo_nfreq rows; Fortran format `(2E18.8E3)` → each field width 18,
  8 mantissa digits, 3-digit exponent, e.g. `   0.00000000E+000   0.35531546E-006`.
- Column 1 = ħω in eV; column 2 = σ_abc(0;ω,−ω) in **A/V²** (SI). No spin-degeneracy factor 2
  anywhere (occupations are 0/1 per Wannier band).

---

## E. Fermi energy and occupation summary

- fermi_n ≥ 1 required (berry.F90:254-257); list built by `w90_readwrite_read_fermi_energy`
  (readwrite.F90:609-692; single `fermi_energy` or min/max/step scan).
- Only `fermi_energy_list(1)` is used by sc (berry.F90:2531).
- occ = step function at T=0, strict `eig < ef` (postw90_common.F90:969-972);
  `occ_fac = occ(n) − occ(m)`, pair skipped when |occ_fac| < 1e-10 (berry.F90:2570-2571).
- Also both eig(n), eig(m) ≤ kubo_eigval_max required (berry.F90:2568).
- scissors_shift (if num_valence_bands>0) enters through get_HH_R exactly as for AHC.

---

## F. GaAs test cases (test-suite/tests/testpostw90_gaas_sc_*)

All five share: `berry = true`, `berry_task=sc`, `fermi_energy=7.7414`, `search_shells=12`,
`kubo_freq_min=0.0`, `kubo_freq_max=10.0`, `kubo_adpt_smr=true`, `num_bands=12`, `num_wann=8`,
`exclude_bands: 1-5`, dis (`dis_win_max=24.0d0`, `dis_froz_max=14.0d0` → frozen_states=T →
default `kubo_eigval_max = 14.6667`), fcc GaAs cell `(-5.34 0 5.34 / 0 5.34 5.34 / -5.34 5.34 0)`
**bohr** (V_c = 2·(5.34·0.52917720859)³ ≈ 45.12 Å³), `mp_grid 4 4 4` with explicit 64 kpoints,
projections `As: s,p; Ga: p,s` at f=(¼,¼,¼)/(0,0,0). Files needed: gaas.win, gaas.chk
(from `gaas.chk.fmt.bz2` via `w90chk2chk.x -f2u`, per-test Makefile), gaas.eig, gaas.mmn
(bz2 in some dirs); gaas.amn present in some dirs but unused.

Differing keywords (verbatim from each gaas.win):

| test | berry_kmesh | kubo_freq_step | sc_phase_conv | sc_eta | sc_use_eta_corr | use_ws_distance | nfreq |
|---|---|---|---|---|---|---|---|
| `sc_xyz` | 25 25 25 | 0.03 | 1 | 0.040 | **.false.** (explicit) | .false. | 334 |
| `sc_eta_corr` | 25 25 25 | 0.03 | 1 | 0.040 | **.true.** (explicit) | .false. | 334 |
| `sc_xyz_ws` | 10 10 10 | 0.05 | 1 | 0.040 | (unset → default **T**) | true | 201 |
| `sc_xyz_scphase2` | 15 15 15 | 0.05 | 2 | 0.040 | (unset → default **T**) | false | 201 |
| `sc_xyz_scphase2_ws` | 11 11 11 | 0.05 | 2 | **0.10** | (unset → default **T**) | true | 201 |

Note: `sc_xyz` vs `sc_eta_corr` differ ONLY in sc_use_eta_corr F/T; the ws/scphase2 tests run
with the eta correction ON (default). For step 0.03 the actual grid step is 10/333 =
0.030030030… (nint readjustment, §B).

Harness (`tests/jobconfig:437-465`): all five use `program = POSTW90_SC_OK`,
`inputs_args = ('gaas.win','')`, **`output = gaas-sc_xyz.dat`** — only the a=x,(b,c)=(y,z)
component is benchmarked. `tests/userconfig:171-175`:

```
[POSTW90_SC_OK]
exe = ../../postw90.x
extract_fn = tools parsers.parse_sc_dat.parse
tolerance = (  (1.0e-6, 5.0e-6, 'energy'),
               (1.0e-6, 1.0e+2, 'shiftcurr'))
```

Tuple order is **(abs_tol, rel_tol, name)** (testcode/lib/testcode2/config.py:35-53, strict:
both must pass) — so shiftcurr is effectively an absolute 1e-6 A/V² check per row (rel 1e2 is
a no-op except benchmark==0 with diff≠0 → Inf; the abs check catches that first). The parser
(test-suite/tools/parsers/parse_sc_dat.py) reads 2-column rows into 'energy'/'shiftcurr'.
[Aside: the "(relative, absolute)" gloss in berry-ahc.md §E has the order swapped.]

Benchmark anchors (`benchmark.out.default.inp=gaas.win` = reference gaas-sc_xyz.dat rows;
row 1 / row 34 / row 67):

```
sc_xyz:              0.00000000E+000  0.36274114E-006
                     0.99099099E+000  0.19253485E-005
                     0.19819820E+001  0.66548988E-005
sc_eta_corr:         0.00000000E+000  0.35531546E-006
                     0.99099099E+000  0.18907956E-005
                     0.19819820E+001  0.66172886E-005
sc_xyz_ws:           0.00000000E+000  0.11553248E-005
                     0.16500000E+001  0.72455917E-005
                     0.33000000E+001  0.13550804E-004
sc_xyz_scphase2:     0.00000000E+000  0.34843095E-006
                     0.16500000E+001  0.43509545E-005
                     0.33000000E+001  0.10232568E-004
sc_xyz_scphase2_ws:  0.00000000E+000  0.99741350E-006
                     0.16500000E+001  0.57984762E-005
                     0.33000000E+001  0.11771852E-004
```

(row 100: sc_xyz 0.29729730E+001/0.80211781E-005, sc_eta_corr 0.29729730E+001/0.80159245E-005,
sc_xyz_ws 0.49500000E+001/0.13086964E-004, scphase2 0.49500000E+001/0.11094091E-004,
scphase2_ws 0.49500000E+001/0.11499692E-004.)

---

## Traps

1. **Two distinct broadenings**: `sc_eta` (fixed, default 0.04 eV) regularises D_h and the
   eta-correction; `eta_smr` (adaptive kubo smearing, per pair per k) is used for the delta
   functions AND the sc_w_thr windows. Do not conflate them.
2. **Gaussian only**: `utility_w0gauss_vec` errors for any smearing type except type_index 0
   (utility.F90:1138-1163). m-p/m-v/f-d smr types break berry_task=sc at the first k-point.
3. Only `fermi_energy_list(1)` is used; Fermi scans are not rejected for sc, just ignored.
4. `kubo_eigval_max` default depends on the frozen window: GaAs tests → 14.0 + 0.6667 =
   14.6667 eV; bands above it are dropped from the pair loop entirely.
5. Frequency-window indices use Fortran `int()` (truncation toward zero) with `+1` offset;
   nfreq is forced ≥ 2 and the step is recomputed from (max−min)/(nfreq−1) after nint rounding.
6. The `/(eig(m) − eig(n))` in gen_r_nm divides ONLY the `cmplx_i*(w + sum_HD + Δ-terms)`
   group (operator precedence), not the whole expression.
7. `r_mn` uses matrix element (m,n) while `gen_r_nm` uses (n,m); D_h_no_eta has a 1e-7
   degeneracy guard, D_h (P-value) has none.
8. `utility_zdotu` = plain Σ a·b — **no complex conjugation** in the intermediate-state sums.
9. TB convention (sc_phase_conv=1) subtracts each WF's own centre from the diagonal R=0
   elements of AA and dAA, and needs Im-log Wannier centres from get_AA_R; mismatch with chk
   centres > 1e-8 (sum of squares) is fatal unless guiding_centres=T.
10. Diagonalisation differs between conventions only by diagonal phases in UU; the full
    integrand is gauge-covariant, so conv 1 vs 2 differ only through finite-R interpolation
    error (compare only physical outputs, per the gauge-invariance note).
11. No spin factor, no adaptive-kmesh refinement, no symmetrisation over components — all 18
    (abc) files written with a,b,c independent, b↔c symmetric by construction.
12. `cell_volume` is the signed determinant (no abs).
13. eta_smr has no floor; with adaptive smearing and (near-)parallel bands the 1/eta_smr in
    delta can blow up (only |occ_fac|>1e-10 and eigval_max filter pairs).
14. CODATA2006 `eV_seconds = 6.582119e-16` has fewer digits than the 2018/2022 value —
    reproduce the constant set exactly for file-precision matches.
15. Each unordered band pair contributes 4 delta terms (2 ordered visits × 2 branches); the
    1/4 in the prefactor (vs 1/2 in SS00 Eq. 57) compensates the doubled delta count.

---

## Implementation checklist (condensed)

1. Read chk/eig/mmn; build HH_R, AA_R (+ Im-log centres) exactly as for AHC (berry-ahc.md §A);
   fold ndegen(+ws_distance) into R-space operators.
2. Per k of berry_kmesh: per sc_phase_conv build H(k), ∂H, ∂∂H, A, ∂A (with/without τ_ij in
   phases/prefactors and TB diagonal subtraction); diagonalise → eig, U; rotate all to H gauge;
   band velocities eig_da; occ at E_F(1).
3. Per ordered pair (n≠m, both ≤ eigval_max, occ_fac≠0): eta_smr (adaptive: min(√2·|Δv|Δk, 1eV));
   skip if outside padded [wmin,wmax]; sum_AD/sum_HD with η-regularised D_h; r_mn (m,n-element);
   gen_r_nm per direction a (8 terms + optional eta correction over p∉{n,m});
   I_nm(a,bc) = Im[r_mn(b)gen_r(c) + r_mn(c)gen_r(b)].
4. Accumulate occ_fac·I_nm·Gaussian((ω−E_nm)/eta_smr)/(√π eta_smr) over the truncated index
   window, both ±E_nm branches; weight 1/N_k; sum over k.
5. Multiply by fac = eV_seconds·π·e³/(4ħ²V_c) (CODATA2006) and write 18 files
   `seedname-sc_abc.dat`, rows `(2E18.8E3)`: ω [eV], σ_abc [A/V²]; no headers.
