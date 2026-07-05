# Tetrahedron-method Spin Hall Conductivity (Ghim–Park + Kawamura correction)

Reference: Wannier90 (greenfield ref tree), `src/postw90/tetrahedron.F90` (module
`w90_tetrahedron`, 620 lines, read in full) and `src/postw90/berry.F90`
subroutine `berry_get_shc_tetrahedron` (line 3101) plus its caller inside
`berry_main` (the `else ! tetrahedron_method` branch, lines 1049–1187) and the
output writer (lines 1690–1727). Physics reference: **Minsu Ghim & Cheol-Hwan
Park, PRB 106, 075126 (2022)** for the SHC tetrahedron integrals; **Takashi
Kawamura et al., PRB 89, 094515 (2014)** for the optimized-tetrahedron
higher-order correction (the P-matrix). All line numbers are for this ref tree.

This note is the tetrahedron **integration-weight** replacement for Gaussian
smearing in SHC. The band-resolved integrand `Im[j^{spin γ}_{α,nm} v_{β,mn}]`
(what our `shc_k_band` already produces per (n,m,k)) is fed in unchanged; only
the way it is summed over k and combined with the energy denominator changes.

Constants: **the ref build defaults to CODATA2006** (`constants.F90:94-96`
select `CODATA2006` when no macro is set). Relevant values:
`elem_charge_SI = 1.602176487e-19` C, `hbar_SI = 1.054571628e-34` J·s,
`bohr = 0.52917720859` Å (`constants.F90`, CODATA2006 block).

---

## 0. Units, conventions, index base (silent-mismatch axes)

- **`imjv(n,m) = Im[ spinVel_{γα}(n,m) · VV_β(m,n) ]`** is stored real, in
  **eV²·Å²** (velocity ~ eV·Å from `∂H/∂k` in eV with k in Å⁻¹). Band indices are
  1-based; `n` = row (spin-current band), `m` = column (velocity band). Only
  `n ≠ m` terms enter (line 1119 `if (n == m) cycle`).
- **Eigenvalues `eig(n)` are in eV.** The Fermi level `Ef` and frequencies `ω`
  are in eV. So `Ef < E1_s(1)` comparisons and `hw = ±ω` are all in eV.
- **The tetrahedron vectors `t(3,3)` (`ttet`) are edge vectors of a tetrahedron
  in reduced (fractional) k-coordinates**, not Cartesian. `ttet(i,k) =
  kptv(k+1,i) − kptv(1,i)` (berry.F90:1113) with `kptv` in units of `db1,db2,db3
  = 1/mesh(1..3)`. `Det_t` in `tetrahedron_integral` is the fractional-coordinate
  parallelepiped volume; the physical 1/N_k volume normalization is applied
  **once at the end** through `fac` (see §5). Do NOT also divide by N_k inside the
  tetrahedron — `Det_t` already carries the per-cube fractional volume (6 tets ×
  Det_t sum = 1/N_k fractional cube).
- Fortran arrays 1-based; loops `do i=1,4` over tetrahedron vertices, `do k=1,20`
  over correction points, `do itet=1,6` over the six tetrahedra per cube.
- **Sign of the final accumulation is negative**: `shc_fermi(ifreq) = shc_fermi(ifreq)
  - real(shc_k_tet)` (line 1175); `shc_freq(ifreq) = shc_freq(ifreq) - shc_k_tet`
  (line 1177). The overall `fac` (§5) is also negative for fermiscan is NOT — see §5;
  the two signs are independent and must both be reproduced.

---

## 1. k-cube → 6 tetrahedra decomposition and the P-matrix correction

### 1.1 The 6-tetrahedra split of one cube (`tetrahedron_array_small_init`, lines 78–96)

Each mesh cube has 8 corners labelled 1..8 (layer 1 = z-bottom: 1,2,3,4;
layer 2 = z-top: 5,6,7,8; local numbering `1=(0,0,0), 2=(1,0,0), 3=(0,1,0),
4=(1,1,0), 5=(0,0,1), 6=(1,0,1), 7=(0,1,1), 8=(1,1,1)`). The cube is cut into 6
tetrahedra **sharing the main diagonal 1–8**:

```
tet 1: 1 5 6 8
tet 2: 1 7 5 8
tet 3: 1 3 7 8
tet 4: 1 4 3 8
tet 5: 1 2 4 8
tet 6: 1 6 2 8
```

This "small" (uncorrected) array is **not used** — the ref errors out if
`tetrahedron_higher_correction = .false.` (see §3). It documents the diagonal
choice; the corrected version below embeds the same 6 diagonal-sharing tets in a
larger 20-point stencil.

### 1.2 The corrected 20-point stencil (`tetrahedron_array_init`, lines 98–122)

With the Kawamura correction each tetrahedron is described by **4 vertices + 16
surrounding points = 20 points**, indexed into a **4×4×4 = 64-point super-block**
around the cube. `tet_array(6,20)`: the first 4 columns of each row are the
vertices (still the 6 diagonal-sharing tets, now in 64-block numbering):

```
tet 1: 22 38 39 43    (+16 corr. pts: 6 37 35 64 5 33 56 48 1 54 40 47 59 23 42 18)
tet 2: 22 42 38 43    (+16: 2 46 33 64 6 41 54 44 1 62 34 48 63 18 47 17)
tet 3: 22 26 42 43    (+16: 18 10 41 64 2 9 62 60 1 30 58 44 47 38 27 21)
tet 4: 22 27 26 43    (+16: 17 28 9 64 18 11 30 59 1 32 25 60 48 21 44 5)
tet 5: 22 23 27 43    (+16: 21 19 11 64 17 3 32 63 1 24 31 59 44 26 39 6)
tet 6: 22 39 23 43    (+16: 5 55 3 64 21 35 24 47 1 56 7 63 60 6 59 2)
```

Main diagonal is **22–43**. The 8 cube vertices are the 64-block indices
`22,23,26,27,38,39,42,43`. Of the 64 points, 18 (`4,8,12,13,14,15,16,20,29,36,
45,49,50,51,52,53,57,61`) are never referenced.

**64-block layout** (berry.F90:1096–1106): the running index is
`idx = 16*l + 4*k + i + 1`, `i,k,l ∈ {0,1,2,3}`, mapping to physical k
```
kptc(1,idx) = (loop_x + i - 1)*db1 + 0.5*db1     ! i is the x-offset
kptc(2,idx) = (loop_y + k - 1)*db2 + 0.5*db2     ! k is the y-offset
kptc(3,idx) = (loop_z + l - 1)*db3 + 0.5*db3     ! l is the z-offset
```
so `i` fastest (x), then `k` (y), then `l` (z, slowest). The central cube
vertices `22,23,26,27,38,39,42,43` correspond to `i,k∈{1,2}`, `l∈{1,2}`
(`22 = 16·1+4·1+1+1`). The whole stencil is **shifted by half a cell**
(`mesh_shift = 0.5`, berry.F90:233) so the Γ point is avoided.

### 1.3 The P-matrix (`tetrahedron_P_matrix_init`, lines 40–76)

`P_matrix(4,20)` maps the 20 stencil values (vertices + 16 correction points) to
**4 corrected corner values** per tetrahedron (Kawamura Eq. (16)):

```
E1tet(i) = Σ_{k=1..20} P_matrix(i,k) · E1_opt(k)     i=1..4   (berry.F90:1131)
E2tet(i) = Σ_k P_matrix(i,k) · E2_opt(k)
Ftet(i)  = Σ_k P_matrix(i,k) · F_opt(k)
```
The raw literal matrix (before the final `/1260`) — 4 rows × 20 cols in 5 blocks:

```
cols 1:4    (vertices, near-identity)   /1260 → main weights 1440,30
 row1: 1440    0   30    0
 row2:    0 1440    0   30
 row3:   30    0 1440    0
 row4:    0   30    0 1440
cols 5:8
 row1: -38   7  17 -28
 row2: -28 -38   7  17
 row3:  17 -28 -38   7
 row4:   7  17 -28 -38
cols 9:12
 row1: -56   9 -46   9
 row2:   9 -56   9 -46
 row3: -46   9 -56   9
 row4:   9 -46   9 -56
cols 13:16
 row1: -38 -28  17   7
 row2:   7 -38 -28  17
 row3:  17   7 -38 -28
 row4: -28  17   7 -38
cols 17:20
 row1: -18 -18  12 -18
 row2: -18 -18 -18  12
 row3:  12 -18 -18 -18
 row4: -18  12 -18 -18
```
then **`P_matrix = P_matrix / 1260.0`** (line 74). Row sums are 1 (e.g. row1:
`(1440+30 −38+7+17−28 −56+9−46+9 −38−28+17+7 −18−18+12−18)/1260 = 1260/1260 = 1`),
so the correction preserves constants. The correction bumps the interpolation
from linear to a higher-order polynomial that removes the O(1/N²) tetrahedron
error on the eigenvalues **and** on the integrand `F`.

**Trap:** the 16 correction points are shared across neighbouring cubes, so the
build precomputes `imjv`/`eig` on a `-1 … mesh+1` padded grid over 4 z-layers
(`imjv(:,:,0:mesh1+2,0:mesh2+2,0:3)`, berry.F90:474) and slides the z-window one
layer per `loop_z` (lines 1075–1083: reuse layer `i+1→i`, only recompute the new
top layer). A from-scratch port can just evaluate `imjv`/`eig` on the full padded
grid once — the sliding is a memory optimization, not a numerical one.

---

## 2. Per-tetrahedron integration weight vs the four corner eigenvalues

The whole point: for each ordered band pair `(n,m)` and each tetrahedron, the
band-resolved integrand at the 4 corrected corners is `Ftet(1..4)` (this is our
`imjv`), the two bands' corrected corner energies are `E1tet(1..4)` (band n) and
`E2tet(1..4)` (band m). The Kubo SHC needs

  ∫_tet d³k  θ(Ef − E_n(k)) · F(k) / D(k)^p ,  D = E_n − E_m,

with `p=2` for ω=0 (fermiscan / static) and `p=1` (two shifted poles) for ω≠0.
The tetrahedron method does this **analytically per tetrahedron**; there is no
smearing width.

### 2.1 Top-level dispatch: `tetrahedron_spinhall` (lines 125–167)

```
tetrahedron_spinhall(F,E1,E2,t,hw,Ef,type,cutoff,avoid_deg)
   = tetrahedron_fermidirac(F,E1,E2,...)   [ occupy by band n = E1 ]
   - tetrahedron_fermidirac(F,E2,E1,...)   [ occupy by band m = E2 ]
```
i.e. the antisymmetrized `(f_{nk} − f_{mk})` factor of the Kubo formula. Both
terms integrate the **same** integrand `F/D^p` but with the θ-occupation stepping
on `E1` then on `E2`. A fast exit: if all four corners have `E1<Ef` AND `E2<Ef`
(both fully occupied, `flag1`) or all four have both unoccupied (`flag2`), the
tetrahedron returns 0 (lines 151–160) — `f_{nk}−f_{mk}=0`.

### 2.2 Fermi-surface occupation weights: `tetrahedron_fermidirac` (lines 169–267)

Sort the four vertices so `E1_s(1) ≤ E1_s(2) ≤ E1_s(3) ≤ E1_s(4)` (band n),
carrying `F`, `t`, and `D = E1_s − E2_s` along (`tetrahedron_sort`, §2.4). Then a
**5-case split on where `Ef` falls** (this is the tetrahedron analog of the
Fermi–Dirac step, replacing `occ_fermi(n,i)·Ω` of the Gaussian path):

- **Case 1** `Ef < E1_s(1)`: band n empty on this tet → 0.
- **Case 2** `E1_s(1) ≤ Ef < E1_s(2)`: one corner below Ef → **one small
  sub-tetrahedron** at the E1_s(1) corner. Linear intercepts
  `x(j) = (Ef − E1_s(1))/(E1_s(j+1) − E1_s(1))`, j=1..3. Interpolate `F`, `D`,
  and scale edges `t(:,j) → t(:,j)·x(j)`; integrate that one small tet
  (lines 201–213).
- **Case 3** `E1_s(2) ≤ Ef < E1_s(3)`: two corners below Ef → **large minus/plus
  three sub-tetrahedra** (lines 215–243): `+tet A − tet B + tet C` with the
  parametrization `x(1)=(Ef−E4)/(E2−E4)`, `x(2)=(Ef−E1)/(E3−E1)`,
  `x(3)=(Ef−E1)/(E4−E1)`, `y=(Ef−E3)/(E2−E3)` (all on sorted E1_s). Interpolated
  `F_small`, `D_small` and the explicit edge combinations in lines 223–242 must be
  reproduced verbatim — they are the exact prism-into-3-tets decomposition.
- **Case 4** `E1_s(3) ≤ Ef < E1_s(4)`: three corners below Ef → **full tet minus
  the empty small tet** at the E1_s(4) corner (lines 244–261):
  `x(1)=(Ef−E4)/(E2−E4)`, `x(2)=(Ef−E4)/(E3−E4)`, `x(3)=(Ef−E1)/(E4−E1)`;
  `+integral(full) − integral(small)`.
- **Case 5** `Ef ≥ E1_s(4)`: band n fully occupied → **the whole tetrahedron**
  (line 262–263), `+integral(F_s, D, t_s)`.

Each `integral(...)` call is `tetrahedron_integral` (§2.3), which supplies the
`1/D^p` energy-denominator weight over that (sub-)tetrahedron.

### 2.3 Energy-denominator integral: `tetrahedron_integral` (lines 270–486)

Inputs: corner integrand `F(1..4)`, corner denominators `D(1..4)=E_n−E_m`,
tetrahedron edges `t(3,3)`, shift `hw`, `type ∈ {1,2,3}`, `tet_cutoff`,
`avoid_deg`. First **sort by D** (`tetrahedron_sort`, carrying F and t). Three
regularizations then the analytic integral:

**(a) Degeneracy guard (type 3 only, lines 313–321):** if `|D(j)| < avoid_deg`,
replace `D(j) = avoid_deg·sign(D(j))` and **zero the whole `F`** (the near-degenerate
pair contributes nothing to the static SHC). `avoid_deg` default `3e-4` eV.

**(b) Cutoff regularization (lines 323–347):** three pairwise checks on the sorted
`D` to avoid `1/(D_i − D_j)` blow-ups when two denominators are close relative to
`(DAV + hw)`. If `|(D(2)−D(3))/(DAV+hw)| < tet_cutoff` (DAV=(D2+D3)/2), spread
D(2),D(3) symmetrically by `±0.5·|DAV+hw|·cutoff` and shift D(1),D(4) to keep
ordering. Then similar single-sided nudges for the (1,2) and (3,4) pairs.
`tet_cutoff` default `1e-4` (see §3 — the **Pt tests override to `1e-1`**).

**(c) Volume factor:** `Det_t = |det t|` (line 357–358).

**Type 1 — nondissipative, `1/((E_n−E_m)² − ω²)` linearized as `1/(D+hw)` (used
for ω≠0):** with `dd(i) = (D(4)−D(i))/(D(i)+hw)`, `ll(i)=log1p(dd(i))` (lines
350–353), and the `cc(4,3)`, `bb(4)` coefficient tables (lines 363–376):
```
Ans = ff/(D(4)+hw) · Σ_{i=1..4} F(i)·( cc(i,1) ll(1)+cc(i,2) ll(2)+cc(i,3) ll(3)+bb(i) )
tetrahedron_integral = Ans · Det_t              ! ff = -Π(1+dd)/(dd(dd−dd'))² / 6
```

**Type 3 — nondissipative, `1/(E_n−E_m)²` static (used for ω=0, fermiscan):**
same structure with `hw=0`, different `cc`/`bb`/`ff` tables (lines 461–478):
```
Ans = ff/(D(4)+hw)² · Σ_i F(i)·( Σ_a cc(i,a) ll(a) + bb(i) )
tetrahedron_integral = Ans · Det_t              ! ff = +Π(1+dd)/(dd(dd−dd'))² / 2
```

**Type 2 — dissipative (Im part, the δ-function on `E_n−E_m=ω`), used for ω≠0
only:** a **surface integral** over the plane `D(k)=hw` inside the tet, weighted
by `1/|∇D|` (lines 385–458). `∇D` from the inverse edge matrix `t_inverse` (line
389–400). Then a 4-way split on where `hw` sits between sorted `D(1..4)` (empty
below D(1)/above D(4); one triangle in [D1,D2] and [D3,D4]; two triangles in
[D2,D3]) with a Jacobian `tetrahedron_jacobian(t,x,type)` (lines 543–594) and the
linear surface-integral formula `Jac·(F_uv(0)/2 + (F_uv(1)+F_uv(2))/6)/|∇D|`.

`tetrahedron_log1p` (lines 596–618) is a numerically stable `log(1+x)`:
principal `log|1+x|` for `|x|>0.5`, else `x·log(y)/(y−1)` with `y=1+x`.

### 2.4 `tetrahedron_sort` (lines 489–540)

Bubble-sorts the size-4 key array `a` ascending, permutes the two payload arrays
`b1,b2` the same way, **and rebuilds the edge matrix `t`** so that after sorting
`t(:,i-1) = t_temp(:,ref(i)) − t_temp(:,ref(1))` (i=2..4) — i.e. edges are always
measured from the new first vertex. This is why the caller must pass edges, not
absolute vertex positions.

---

## 3. Keywords and defaults

Defaults are set in **`postw90_readwrite.F90`** (not in the type declaration; the
`postw90_types.F90:185-188` fields are declared without initializers). Read the
keyword-parse block (lines 974–996):

| keyword | default | set at | notes |
|---|---|---|---|
| `tetrahedron_method` | **`.false.`** | `postw90_readwrite.F90:974` | must be `.true.` to use tetrahedra; errors if the berry task is not `shc` (line 977–978) |
| `tetrahedron_higher_correction` | **`.true.`** | `postw90_readwrite.F90:980` | must stay `.true.` — the ref **errors out** (`set_error_input`, line 983–984) if `.false.`: "tetrahedron_method works only with correction". |
| `tetrahedron_cutoff` | **`1.e-4`** (dp) | `postw90_readwrite.F90:986` | must be `>0` (line 989–990). The `(2±cutoff)/(2∓cutoff)` and symmetric-spread regularizers of §2.3(b). |
| `tetrahedron_avoid_degeneracy` | **`3.e-4`** (dp) | `postw90_readwrite.F90:992` | must be `>0` (line 995–996). Type-3 degeneracy guard §2.3(a). |

**Known ref quirk (do NOT copy the bug, but be aware for byte-matching):**
`postw90_readwrite.F90:982` reads the `tetrahedron_higher_correction` keyword's
value into `pw90_berry%tetrahedron_method` (`l_value=pw90_berry%tetrahedron_method`),
apparently a copy-paste slip. In practice both are `.true.` in the tests so it is
harmless; a Julia port should bind `tetrahedron_higher_correction` to its own
field.

Note also `readwrite.F90:1610` registers a stray `tetrahedron_correction`
keyword name in the "known keywords" list (no effect; the active name is
`tetrahedron_higher_correction`).

---

## 4. Caller structure (`berry_main`, tetrahedron branch, lines 1049–1187)

1. Allocate padded `imjv(num_wann,num_wann,0:mesh1+2,0:mesh2+2,0:3)` and
   `eig(num_wann,0:mesh1+2,0:mesh2+2,0:3)` (lines 474–475); init `P_matrix`,
   `tet_array` (480–481). `nfreq = kubo_nfreq` (freqscan) or `fermi_n` (fermiscan)
   (482–486). MPI split over `mesh(3)` z-layers (497).
2. For each z-layer, evaluate `berry_get_shc_tetrahedron` on the padded
   `-1…mesh+1` xy-grid across 4 z-layers `i=0..3` (offset `+mesh_shift`), with BZ
   wrap-around on x,y,z (lines 1058–1084). Slide the z-window to reuse layers.
3. Summation loop over interior cubes `loop_x,loop_y = 0…mesh-1`: gather the 64
   super-block points into `imjv_tet`, `eig_tet` (1096–1106). For each of 6 tets
   (`itet=1..6`): build edges `ttet` (1111–1115), then for every ordered band pair
   `n≠m` (1117–1119) gather the 20-point `F_opt,E1_opt,E2_opt` (1120–1124) and
   apply the P-matrix to get `Ftet,E1tet,E2tet` (1129–1135).
4. For each scan point `ifreq` (1142–1179): set `(ω,Ef)`:
   - fermiscan: `ω = Re(kubo_freq_list(1))` (=0 in the tests), `Ef =
     fermi_energy_list(ifreq)`.
   - freqscan: `ω = Re(kubo_freq_list(ifreq))`, `Ef = fermi_energy_list(1)`.
   Then:
   - **ω = 0** (static): `shc_k_tet = tetrahedron_spinhall(Ftet,E1tet,E2tet,ttet,
     0, Ef, type=3, cutoff, avoid_deg)` (lines 1151–1155).
   - **ω ≠ 0**: (lines 1156–1172)
     ```
     shc_k_tet = [ spinhall(...,hw=-ω, type=1) - spinhall(...,hw=+ω, type=1) ] / (2ω)
               + iπ·[ spinhall(...,hw=-ω, type=2) + spinhall(...,hw=+ω, type=2) ] / (2ω)
     ```
     The real part is the two shifted principal-value poles `1/(D∓ω)`; the
     imaginary part is the `±ω` δ-functions (dissipative).
5. Accumulate: `shc_fermi(ifreq) -= real(shc_k_tet)` (1175) or
   `shc_freq(ifreq) -= shc_k_tet` (1177). MPI-reduce afterwards (1228–1238).

### 4.1 `berry_get_shc_tetrahedron` — the band-resolved integrand (lines 3101–3276)

Produces `imjv(n,m)` and `eig_out(n)` at one k. Identical operator content to the
Gaussian `berry_get_shc_klist` path, but returns the **raw** `Im[jv]` without any
energy denominator (the tetrahedron supplies that). Steps:
- Fourier `HH_R→HH` with `∂H/∂k` (`delHH`, 3177–3182); diagonalize → `eig`, `UU`
  (3189); Fourier `AA`, `SS` and (qiao branch) `SR_R→SAA`, `SHR_R→SBB`
  (3190–3217); rotate all to the eigenbasis by `UU` (3219–3227).
- **Velocity** (β component): `VV(n,m,β) = VV0(n,m,β) − i·AA(n,m,β)·(eig(m)−eig(n))`
  where `VV0 = U†·delHH·U` (3230–3235). This is the covariant velocity
  `v_β = ∂_β H − i[A_β, H]` in the Wannier gauge.
- **Spin velocity** (γ spin, α velocity), qiao/Qiao QZYZ18 Eq.(23):
  `spinVel0(:,:,γ,α) = VV0(:,:,α)·SS(:,:,γ) + SS(:,:,γ)·VV0(:,:,α)` (anticommutator
  ½{v_α, σ_γ}), then a gauge-covariant correction
  `spinVel(n,m,γ,α) = spinVel0 − i(eig(m)·SAA(n,m) − SBB(n,m)) + i(eig(n)·conj(SAA(m,n))
  − conj(SBB(m,n)))`, and finally **`spinVel = spinVel/2`** (3238–3254).
- **Integrand:** `imjv(n,m) = aimag( spinVel(n,m,γ,α) · VV(m,n,β) )` (3257–3261).
  `eig_out = eig` (3262).

Indices: `α = pw90_spin_hall%alpha` (default 1), `β = beta` (2), `γ = gamma` (3)
(`postw90_types.F90:203-205`). The `ryoo` method branch (3199–3207) uses
`SAA_R`/`SBB_R` directly; the default `qiao` branch (3208–3217) uses
`SR_R`/`SHR_R` (this note assumes `shc_method = qiao`, as in the Pt tests).

**Note on `shc_k_band`:** our Julia `shc_k_band` already computes the band-diagonal
`Ω^{spin γ}_{n} = Σ_{m≠n} (−2/((E_m−E_n)²+η²))·Im[jv]` for the Gaussian path. For
the tetrahedron path, feed the **per-pair** `Im[jv](n,m) = imjv(n,m)` (i.e. the
summand *before* the `−2/(…)` energy factor and *before* the m-sum) into `Ftet`,
and let §2 supply the energy denominator and occupation. Do not pre-multiply by
any `−2/(…)` smearing factor.

---

## 5. Output format and the overall unit factor (lines 1690–1727)

After the k-sum, multiply by
```
fac = 1.0e8 · e² / ( ħ · V_cell ) / 2.0          (berry.F90:1690)
    = 1.0e8_dp * physics%elem_charge_SI**2 / (physics%hbar_SI * cell_volume) / 2.0_dp
```
`cell_volume` is the real-space cell volume in **Å³** (triple product of
`real_lattice`, lines 262–267, Å). The `1e8` converts cm→(the Å/eV bookkeeping to
S/cm); the `/2` and `e²/ħ` give spin current in `(ħ/e)·S/cm`. Comment block
(1678–1688) spells it out: the k-weight already carries `1/N_k`; convert charge
current to spin current via `−ħ/2/e`; `×1e8` → S/cm. Final unit
**`(ħ/e)·S/cm`**. Apply to `shc_freq` (complex) or `shc_fermi` (real) (1691–1695).

Then write `seedname-shc-fermiscan.dat` (fermiscan) or `seedname-shc-freqscan.dat`
(freqscan) (1704–1707):

**Fermiscan** (`E17.8`, one signed exponent field):
```
#No.   Fermi energy(eV)   SHC((hbar/e)*S/cm)
   n    <F12.6 Ef>   <E17.8 shc_fermi(n)>
```
format `(I4,1x,F12.6,1x,E17.8)` (line 1716). `shc_fermi(n)` is real.

**Freqscan** (two `E17.8` fields, real & imag):
```
#No.   Frequency(eV)   Re(sigma)((hbar/e)*S/cm)   Im(sigma)((hbar/e)*S/cm)
   n    <F12.6 ω>   <E17.8 Re>   <E17.8 Im>
```
format `(I4,1x,F12.6,1x,1x,2(E17.8,1x))` (line 1723) — note the **double space**
after the frequency and the trailing space after each E-field.

Scan axes:
- `fermi_energy_list(i) = fermi_energy_min + (i−1)·step`, `n = nint(|(max−min)/step|)+1`,
  step rescaled to `(max−min)/(n−1)` (`readwrite.F90:679-691`).
- `kubo_freq_list(i) = freq_min + (i−1)·(freq_max−freq_min)/(nfreq−1)`,
  `kubo_nfreq = nint((max−min)/step)+1` (`postw90_readwrite.F90:1708-1724`).

---

## 6. Reference tests (Pt, spinors, `mp_grid = 4 4 4`, CODATA2006 build)

Common to both: `num_bands=40`, `num_wann=18`, `spinors=true`,
`shc_method=qiao`, `berry=true`, `berry_task=eval_shc`, `berry_kmesh=10`
(→ interpolation mesh `10×10×10`), `berry_curv_unit=ang2`, `guiding_centres=T`,
projections `Pt: d;s;p`. FCC Pt cell (unit_cell_cart in bohr, `a≈3.7039` bohr
edge vectors). `shc_alpha=1, shc_beta=2, shc_gamma=3` (σ_xy^z).

### 6.1 `testpostw90_pt_tetra_shcfermi/Pt.win`
```
shc_freq_scan = false
tetrahedron_method = true
tetrahedron_higher_correction = true
tetrahedron_cutoff = 1.e-1              # OVERRIDES the 1e-4 default
tetrahedron_avoid_degeneracy = 3.e-4   # equals the default
fermi_energy_min = 13
fermi_energy_max = 23
fermi_energy_step = 0.5                 # → 21 Fermi points
```
Benchmarked file: **`Pt-shc-fermiscan.dat`** (jobconfig `[testpostw90_pt_tetra_shcfermi/]`,
`program = POSTW90_SHCFERMIDAT_OK`). Parser `parsers.parse_shc_dat.parse` extracts
columns `energy` (col 2) and `shc` (col 3). **Tolerances (userconfig
`[POSTW90_SHCFERMIDAT_OK]`, lines 203–205):**
`energy` abs=`1.0e-6` rel=`5.0e-6`; **`shc` abs=`1.0e-1` rel=`1.0e-1`**.
Anchor rows (verbatim from `benchmark.out.default.inp=Pt.win`):
```
   1    13.000000   -0.10477881E+04
  11    18.000000    0.18635090E+04
  21    23.000000    0.52279635E+03
```
(i.e. row 1 `−1047.7881`, row 11 `+1863.5090`, row 21 `+522.79635`, all in
`(ħ/e)·S/cm`.)

### 6.2 `testpostw90_pt_tetra_shcfreq/Pt.win`
```
shc_freq_scan = true
tetrahedron_method = true
tetrahedron_higher_correction = true
tetrahedron_cutoff = 1.e-1              # (no tetrahedron_avoid_degeneracy → default 3e-4)
fermi_energy = 18.3823                  # single Fermi level
kubo_freq_min = 0.0
kubo_freq_max = 2.0
kubo_freq_step = 0.1                    # → 21 frequency points
```
Benchmarked file: **`Pt-shc-freqscan.dat`** (jobconfig `[testpostw90_pt_tetra_shcfreq/]`,
`program = POSTW90_SHCFREQDAT_OK`). Parser extracts `frequency` (col 2),
`shc_re` (col 3), `shc_im` (col 4). **Tolerances (userconfig
`[POSTW90_SHCFREQDAT_OK]`, lines 207–212):** `frequency` abs=`1.0e-6` rel=`5.0e-6`;
**`shc_re` abs=`1.0e+1` rel=`1.0e+1`** (loose — the real part is large/steep);
**`shc_im` abs=`1.0e-1` rel=`1.0e-1`**. Anchor rows:
```
   1     0.000000     0.91021053E+03    0.00000000E+00
   9     0.800000     0.24536549E+04    0.13659641E+04
  21     2.000000     0.13488542E+03    0.11288179E+04
```
(ω=0 imag part is exactly 0 — type-2 contributes nothing at ω=0 because the ω=0
branch uses only type 3.)

---

## 7. Traps

1. **`fac` sign vs accumulation sign are independent.** `fac` (§5) is **positive**;
   the negative sign in the ref output comes from `shc_fermi -= real(shc_k_tet)`
   (line 1175). If a port folds the two into one factor, verify the total sign
   against row 1 of §6.1 (`−1047.79`).
2. **Static path uses type 3 only** (`ω=0` branch, line 1151–1155); it does **not**
   call type 1 or type 2. Type 2 (surface δ) is only reached for `ω≠0`.
3. **`tetrahedron_cutoff` default is `1e-4`, but both Pt tests set `1e-1`.** The
   `1e-4` default gives sharper (less regularized) peaks; benchmark numbers were
   generated with `1e-1`. Match the value from the input, not the default.
4. **`tetrahedron_higher_correction=.false.` is a hard error** in the ref; the
   uncorrected 6-tet path (`tetrahedron_array_small_init`) is dead code. Only the
   20-point P-matrix path is implementable.
5. **Edges, not vertices.** `tetrahedron_sort` rebuilds `t` as edges from the first
   sorted vertex; passing absolute k positions instead of edge vectors silently
   corrupts `Det_t` and every case-split parametrization.
6. **`mesh_shift = 0.5`** shifts the whole stencil by half a cell to avoid Γ; the
   effective integration mesh is `berry_kmesh` (=10 here), and `db_i = 1/mesh(i)`.
   The P-matrix stencil reaches ±1 cell beyond each cube, hence the `-1…mesh+1`
   padded evaluation with BZ wrap-around.
7. **`E17.8` field widths and the extra spaces** in the freqscan format
   (`F12.6,1x,1x,2(E17.8,1x)`) must be reproduced for byte-exact `.dat` matching;
   the fermiscan uses a single `1x` and no trailing space (`I4,1x,F12.6,1x,E17.8`).
8. **CODATA build matters for the last digits.** `fac` uses `elem_charge_SI` and
   `hbar_SI`; CODATA2006 vs 2018/2022 differ in the 7th significant figure, which
   is within the loose SHC tolerances here but would matter for a tight
   file-precision comparison. The ref benchmarks assume CODATA2006.
9. **`n≠m` only, ordered pairs.** The double loop is over all ordered `(n,m)` with
   `n≠m` (not `m>n`); `tetrahedron_spinhall` antisymmetrizes `E1↔E2` internally, so
   each ordered pair is summed once — do not additionally halve.
