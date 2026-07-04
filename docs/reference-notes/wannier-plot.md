# Wannier90 reference spec: real-space WF plotting (`wannier_plot = .true.`)

Extracted from `reference/wannier90/src/plot.F90` (subroutine `plot_wannier`,
lines 1547–2485) plus `wannier90_types.F90`, `wannier90_readwrite.F90`,
`constants.F90`, and the docs at
`reference/wannier90/docs/docs/user_guide/wannier90/files.md`.

NOTE: this reference tree is a refactored/development Wannier90. Relative to
official v3.x, `plot_wannier` here hoists the band contraction out of the grid
loop (`c_wvfn`, plot.F90:1858–1876) and factors the k-phase per axis
(plot.F90:1878–1887). The mathematics and all file formats are identical to
upstream; only loop organization differs. Cited line numbers are for this tree.

All Wannier90 internal quantities are in **Angstrom** (`real_lattice`,
`wannier_data%centres`, `atom_data%pos_cart`, `wannier_plot_radius`).
`kpt_latt(3, num_kpts)` are k-points in **fractional** coordinates of the
reciprocal lattice. `real_lattice(i, :)` is lattice vector **a_i** (row =
vector, column = Cartesian component; see PRIMVEC writes, plot.F90:2450–2452).

---

## A. UNK input files

### A.1 Naming

plot.F90:1664–1665 (also used at 1757–1759 with `loop_kpt`):

```fortran
200   format('UNK', i5.5, '.', i1)      ! scalar case: UNK00001.1 / UNK00001.2
199   format('UNK', i5.5, '.', 'NC')    ! spinor case: UNK00001.NC
```

- The `i5.5` field is the k-point index `p` = 1..`num_kpts` (order must match
  the `kpoints` block / `kpt_latt` ordering).
- The trailing digit is the spin channel `wvfn_read%spin_channel`: 1 = up,
  2 = down, selected by input keyword `spin = up|down`
  (wannier90_readwrite.F90:1049–1060; default 1, wannier90_types.F90:134).
- If `spinors = .true.` (noncollinear), a single file `UNKppppp.NC` holds both
  spinor components (plot.F90:1647, 1759).
- Formatted vs unformatted is chosen by keyword `wvfn_formatted`
  (default `.false.` = unformatted, wannier90_types.F90:132,
  wannier90_readwrite.F90:1045–1046). One flag for all files.
- Existence of the k=1 file is checked up front; its header defines
  `ngx, ngy, ngz` (plot.F90:1649–1662).

### A.2 Header record

One record of 5 integers: `ngx ngy ngz ik nbnd`.

```fortran
! formatted (plot.F90:1657, re-read per k at 1763)
read (file_unit, *) ngx, ngy, ngz, nk, nbnd
! unformatted (plot.F90:1660, 1766)
read (file_unit) ngx, ngy, ngz, nk, nbnd
```

- `ngx/ngy/ngz`: real-space FFT grid dimensions of the periodic part
  u_nk(r).
- `ik`: the k-point index; **checked** against the loop index, along with the
  grid dims (plot.F90:1769–1775) — mismatch is a fatal error. `nbnd` is *not*
  cross-checked against `num_bands` (it is only printed in the warning).
- Formatted read is list-directed, so any whitespace separation works.
  Example from `test-suite/tests/testw90_example01/UNK00001.1`, line 1:
  `  20  20  20   1   4` (20^3 grid, k-point 1, 4 bands; the file has
  1 + 20^3*4 = 32001 lines).

### A.3 Per-band data, scalar case

`num_bands` blocks, in band order (ascending, same ordering as the `.mmn`/
`.amn`/eigenvalue bands). Each block is the periodic part u_{b,k}(r) on the
grid, **one complex value per grid point, x fastest, then y, then z**. The
linear index convention used everywhere is plot.F90:1908:

```fortran
npoint = nx + (ny - 1)*ngx + (nz - 1)*ngy*ngx     ! nx,ny,nz in 1..ngx etc.
```

- Formatted variant (plot.F90:1808–1816): one grid point per record,
  list-directed `read (file_unit, *) w_real, w_imag`; i.e. each line contains
  the real and imaginary part. Example data line from the test file:
  ` -0.566  0.407`. (Writers may use any float format; the reader is `*`.)
- Unformatted variant (plot.F90:1825): **one Fortran sequential record per
  band** containing the whole grid as complex double precision (16 bytes per
  point, plus the compiler's standard 4-byte record markers):

```fortran
read (file_unit) (r_wvfn(nx, loop_b), nx=1, ngx*ngy*ngz)
```

The header is likewise a single record of 5 default (4-byte) integers.
(Confirmed by docs/user_guide/wannier90/files.md:670–675.)

### A.4 Per-band data, spinor case (`UNKppppp.NC`)

Same header; each band contributes **two** consecutive blocks: full up-spinor
grid, then full down-spinor grid.

- Formatted (plot.F90:1810–1822): `ngx*ngy*ngz` lines of `(re, im)` for the
  up component, then `ngx*ngy*ngz` lines for the down component.
- Unformatted (plot.F90:1827–1828): two records per band:

```fortran
read (file_unit) (r_wvfn_nc(nx, loop_b, 1), nx=1, ngx*ngy*ngz) ! up-spinor
read (file_unit) (r_wvfn_nc(nx, loop_b, 2), nx=1, ngx*ngy*ngz) ! down-spinor
```

### A.5 Disentanglement windowing while reading

When `have_disentangled` (plot.F90:1777–1805): bands are read sequentially into
slot `counter` of `r_wvfn_tmp(:, counter)`; `counter` advances **only** when
`inc_band(loop_b) = dis_manifold%lwindow(loop_b, ik)` is true, i.e. a band
outside the outer energy window is overwritten by the next band (discarded).
The loop exits once `counter > num_inc = dis_manifold%ndimwin(ik)`, so trailing
bands after the last in-window band are never read. Result: `r_wvfn_tmp`
contains only the `ndimwin(ik)` window states, in ascending band order.

---

## B. WF assembly in `plot_wannier`

### B.1 Which U matrices

Two-step contraction; combined transform is `u_matrix_opt · u_matrix` when
disentangled, `u_matrix` alone otherwise.

Step 1 (only if `have_disentangled`), plot.F90:1839–1844 — project window
states onto the `num_wann`-dim optimal subspace:

```fortran
do loop_w = 1, num_wann
  do loop_b = 1, num_inc
    r_wvfn(:, loop_w) = r_wvfn(:, loop_w) + &
                        u_matrix_opt(loop_b, loop_w, loop_kpt)*r_wvfn_tmp(:, loop_b)
```

(spinor case identical via `zaxpy`, applied to both components with the same
coefficient, plot.F90:1846–1855).

Step 2 (always), plot.F90:1858–1876 — apply the gauge matrix, but only for the
WFs actually being plotted (`wannier_plot%list`):

```fortran
! Contract band index: c_wvfn(npoint, w) = sum_b u_matrix(b, list(w), k) * r_wvfn(npoint, b)
uw = u_matrix(loop_b, wannier_plot%list(loop_w), loop_kpt)
c_wvfn(npoint, loop_w) = c_wvfn(npoint, loop_w) + uw*r_wvfn(npoint, loop_b)
```

So the smooth-gauge periodic part is
u~_{n k}(r) = Σ_b [U_opt(k) U(k)]_{b n} u_{b k}(r), with U indexed
`u_matrix(band, wann, kpt)` and `u_matrix_opt(window_band, wann, kpt)`.

### B.2 Supercell grid and bounds

`ngs = wannier_plot%supercell` (default (2,2,2)). Supercell grid indices
(plot.F90:1674–1680):

```fortran
nxx_lo = -((ngs(1))/2)*ngx
nxx_hi = ((ngs(1) + 1)/2)*ngx - 1     ! integer division
```

(same for y, z). `wann_func(nxx_lo:nxx_hi, nyy_lo:nyy_hi, nzz_lo:nzz_hi, wann_plot_num)`
(plot.F90:1682). The **home cell is grid indices 1..ngx** (index n ↔ fractional
coordinate (n−1)/ngx, see phase below). Examples: ngs=1 → indices 0..ngx−1
(home cell, shifted down by one point); ngs=2 → −ngx..ngx−1 (cells −1, 0);
ngs=3 → −ngx..2·ngx−1 (cells −1, 0, +1). I.e. unit cells `-(ngs/2) ..
(ngs+1)/2 - 1`: the supercell is centred on the home cell, with the extra cell
on the negative side for even ngs and on the positive side for odd ngs.

### B.3 The Bloch sum and phase

Per k-point, the phase is factored per axis (plot.F90:1878–1887):

```fortran
phase_x(nxx) = exp(twopi*cmplx_i*kpt_latt(1, loop_kpt)*real(nxx - 1, dp)/real(ngx, dp))
```

so with fractional position f(nxx,nyy,nzz) = ((nxx−1)/ngx, (nyy−1)/ngy,
(nzz−1)/ngz) (which ranges over the whole supercell, negative included), the
total phase `catmp = phase_x(nxx)*phase_y(nyy)*phase_z(nzz)` equals
**e^{i k·r}** evaluated directly at the supercell point r — this *is* the
e^{i k·(r0+R)} construction, built as one exponential rather than a separate
cell-phase e^{i k·R}. Grid index 1 ↔ r = 0 (note the `nxx - 1`).

Grid accumulation (plot.F90:1898–1923): loop nzz (outer), nyy, nxx; the
home-cell image of each supercell point is

```fortran
nz = mod(nzz, ngz); if (nz .lt. 1) nz = nz + ngz     ! maps to 1..ngz
npoint = nx + (ny - 1)*ngx + (nz - 1)*ngy*ngx
wann_func(nxx, nyy, nzz, loop_w) = wann_func(nxx, nyy, nzz, loop_w) + c_wvfn(npoint, loop_w)*catmp
```

The innermost loop is over plotted WFs (a performance choice, see comment
plot.F90:1895–1897; upstream v3.x loops bands innermost with
`u_matrix` applied there — same result).

K-loop: `do loop_kpt = 1, num_kpts` (plot.F90:1746), MPI-distributed by
`dist_k` and summed with `comms_reduce` (plot.F90:1927–1933).

Full formula implemented:

w_n(r) = (1/num_kpts) Σ_{k} e^{i k·r} Σ_b [U_opt(k) U(k)]_{b n} u_{b k}(r mod cell)

### B.4 Normalization and phase fixing

- **1/num_kpts**: non-spinor at plot.F90:1983 (inside the max-modulus scan),
  spinor at plot.F90:1964 (inside the spinor-collapse loop). Exactly one
  division per element in either path. No other normalization (no cell-volume
  factor; values inherit the UNK normalization convention).
- **Global phase fixing (non-spinor only)**, plot.F90:1972–1997: find the grid
  point maximizing |w|², take `wmod = wann_func(...)` there, then

```fortran
wmod = wmod/sqrt(real(wmod)**2 + aimag(wmod)**2)
wann_func(:, :, :, loop_w) = wann_func(:, :, :, loop_w)/wmod
```

i.e. rotate the global phase so the max-|w| point becomes **real positive**.
Prints `Wannier Function Num: <n>  Phase Factor = <1/wmod>`
(plot.F90:1995–1996).
- **Reality check**, plot.F90:2002–2017: over points with |Re w| ≥ 0.01,
  report max |Im w|/|Re w| as `Maximum Im/Re Ratio = <ratmax>`.
- **Spinor collapse** (plot.F90:1937–1969): per point compute
  `upspinor = |w_up|²`, `dnspinor = |w_dn|²`; per `wannier_plot%spinor_mode`:
  `'total'` → sqrt(up²+dn²) — literally `sqrt(upspinor + dnspinor)`;
  `'up'`/`'down'` → sqrt of that component, optionally multiplied by
  `sign(1, Re w_component)` when `wannier_plot%spinor_phase` (default true).
  Result stored back into the real part of `wann_func`; **no phase fixing /
  reality check** for spinors (plot.F90:1973 comment).

There is **no translation of the WF data to the home cell**: the volumetric
grid is always the fixed supercell above. (`translate_home_cell` is an
unrelated keyword affecting `seedname_centres.xyz` output only —
wannier90_readwrite.F90:1372–1373, plot.F90:3174 `plot_write_xyz`.)
What `wannier_plot_mode` does instead is per-format (below).

Dispatch (plot.F90:2020–2028): `wannier_plot%format .eq. 'xcrysden'` →
`internal_xsf_format()`, `.eq. 'cube'` → `internal_cube_format(...)`, else
error. (Input validation accepts any string containing 'xcrys'/'cub',
wannier90_readwrite.F90:1260–1263, but only the exact strings dispatch.)

---

## C. Output formats

### C.1 XSF / XCrySDen (`internal_xsf_format`, plot.F90:2407–2483)

One file per WF: `format(a, '_', i5.5, '.xsf')` → `seedname_00001.xsf`
(plot.F90:2411, 2437). **All lengths in Angstrom** (Wannier90 internal units,
which matches the XSF standard). Layout:

1. Four comment lines starting `#` ("Generated by the Wannier90 code…", date).
2. If `index(wannier_plot%mode, 'mol') > 0`: the single line `ATOMS`
   (plot.F90:2445–2446); otherwise (crystal mode):

```
CRYSTAL
PRIMVEC
  <a1x a1y a1z>        ! (3f12.7), real_lattice(1,1:3), Angstrom
  <a2…>, <a3…>
CONVVEC
  <same three lines again>
PRIMCOORD
  <num_atoms>  1       ! format (i6,"  1")
```

   (plot.F90:2448–2458).
3. Atom lines, both modes, species-major order: `(a2,3x,3f12.7)` = 2-char
   species label + Cartesian position in Angstrom, exactly as in the input cell
   — **no replication, no translation** (plot.F90:2460–2464).
4. Blank line, then (plot.F90:2466–2472):

```
BEGIN_BLOCK_DATAGRID_3D
3D_field
BEGIN_DATAGRID_3D_UNKNOWN
  <ngs(1)*ngx  ngs(2)*ngy  ngs(3)*ngz>    ! (3i6)  — grid dims, NO +1
  <x_0ang y_0ang z_0ang>                  ! (3f12.6) origin, Angstrom
  <dirl(1,1:3)>                           ! (3f12.7) spanning vector 1
  <dirl(2,1:3)>, <dirl(3,1:3)>
```

   Rather than the periodic "+1 duplicated point" XSF convention, Wannier90
   writes N = ngs·ng points as a **general** datagrid whose spanning vectors
   are shortened by one grid spacing so the N points span them inclusively
   (plot.F90:2428–2433):

```fortran
fxcry(i) = real(ngs(i)*ng_i - 1, dp)/real(ng_i, dp)
dirl(:, j) = fxcry(:)*real_lattice(:, j)          ! row i = span vector i
```

   Origin = Cartesian position of grid point (nxx_lo, nyy_lo, nzz_lo), from
   r(n) = (n−1)/ng · a (plot.F90:2418–2426):

```fortran
x_0ang = -real(((ngs(1))/2)*ngx + 1, dp)/real(ngx, dp)*real_lattice(1, 1) - &
         real(((ngs(2))/2)*ngy + 1, dp)/real(ngy, dp)*real_lattice(2, 1) - &
         real(((ngs(3))/2)*ngz + 1, dp)/real(ngz, dp)*real_lattice(3, 1)
```

   Note the `+ 1`: origin index is nxx_lo = −(ngs/2)·ngx with coordinate
   (nxx_lo − 1)/ngx, hence the `(ngs/2)*ngx + 1` numerator — one grid spacing
   *below* −(ngs/2)·a.
5. Data (plot.F90:2473–2475): **real part only**, format `(6e13.5)` (6 values
   per line), **x fastest, then y, then z** (XSF column-major convention):

```fortran
write (file_unit, '(6e13.5)') &
  (((real(wann_func(nx, ny, nz, loop_b)), nx=nxx_lo, nxx_hi), ny=nyy_lo, nyy_hi), nz=nzz_lo, nzz_hi)
```

6. `END_DATAGRID_3D` / `END_BLOCK_DATAGRID_3D` (plot.F90:2476).

In molecule mode the only differences are `ATOMS` instead of the
CRYSTAL/PRIMVEC/PRIMCOORD block; the datagrid is identical.

### C.2 Gaussian cube (`internal_cube_format`, plot.F90:2075–2405)

One file per WF: `format(a, '_', i5.5, '.cube')` → `seedname_00001.cube`
(plot.F90:2151, 2181). **Everything divided by `bohr` on output → Bohr units**
(comment plot.F90:2325). `bohr` is `bohr_angstrom_internal` from
constants.F90:224; default build uses CODATA2006 → **0.52917720859 Å**
(constants.F90:92–100 default `#define CODATA2006`, value at :182; other
CODATA flags change it in the 10th decimal;
`USE_WANNIER90_V1_BOHR` → 0.5291772108, constants.F90:219).

Unlike xsf (whole supercell), cube writes a **cut-out box centred on the WF
centre** `wannier_data%centres(:, wann_index)` (the converged centre, Å):

- Box extent along each lattice direction i (plot.F90:2184–2194): the centre's
  fractional coordinate times |a_i|, ± the radius converted to a length along
  a_i:

```fortran
rstart(i) = (centre · recip_lattice(i,:))*moda(i)/twopi - twopi*wannier_plot%radius/(moda(i)*modb(i))
rend(i)   =  same + twopi*wannier_plot%radius/(moda(i)*modb(i))
```

  where `moda(i) = |a_i|` (Å), `modb(i) = |b_i|` (Å⁻¹) (plot.F90:2154–2161)
  and `dgrid(i) = moda(i)/ng_i` (plot.F90:2164).
- Integerization (plot.F90:2196–2201): `ilength = ceiling((rend-rstart)/dgrid)`,
  `istart = floor(rstart/dgrid) + 1`, `iend = istart + ilength - 1`.
- Cube origin, Cartesian Å (plot.F90:2204–2208):
  `orig(i) = Σ_j (istart(j)-1)*dgrid(j)*real_lattice(j,i)/moda(j)`.
  NOTE: origin corresponds to index istart−1, i.e. fractional (istart−1)/ng —
  this differs by one grid spacing from the supercell-array convention
  ((n−1)/ng); the data extraction below compensates via `qxx = nxx + istart - 1`.
- Data extraction (plot.F90:2235–2280): for each output point
  `q = n + istart - 1`, fold into the available supercell array by whole grid
  periods, `izz = int((abs(qzz) - 1)/ngz); if (qzz < nzz_lo) qzz = qzz + izz*ngz`;
  if q exceeds `n**_hi` → fatal error advising "(1) increase
  wannier_plot_supercell; (2) decrease wannier_plot_radius; (3) set
  wannier_plot_format=xcrysden" (plot.F90:2241–2248). Stored value =
  `real(wann_func(qxx, qyy, qzz, loop_w), dp)` — **real part only**
  (plot.F90:2277).
- `wannier_plot_mode` (plot.F90:2132–2135: 'mol' → lmol, 'crys' → lcrys):
  - **molecule mode** (plot.F90:2297–2301): the cube *origin* is translated by
    the lattice vector `irdiff = nint(comf - wcf)` (comf = fractional centre of
    mass of all atoms, plot.F90:2166–2175; wcf = fractional WF centre,
    plot.F90:2283) so the WF appears next to the atoms as given in the input;
    atoms written exactly as in the input (plot.F90:2353–2354).
  - **crystal mode** (default): atoms are replicated over image cells
    nxx,nyy,nzz ∈ [−ngs/2, (ngs+1)/2] and written iff their distance from the
    WF centre ≤ `wannier_plot%scale * wannier_plot%radius`
    (plot.F90:2302–2323 counts `icount` for the header; 2355–2373 writes them).

Cube file layout (plot.F90:2326–2385):

```
<comment line 1: "Generated by Wannier90 code http://www.wannier.org">
<comment line 2: "On <date> at <time>">
natoms  origx origy origz          ! (i4,3f13.5); natoms=num_atoms (mol) or icount (crys); orig/bohr
ilength(1)  a1x/(ngx*bohr) a1y/(ngx*bohr) a1z/(ngx*bohr)     ! (i4,3f13.5) voxel vectors
ilength(2)  a2/(ngy*bohr) ...
ilength(3)  a3/(ngz*bohr) ...
Z  1.00000  x/bohr y/bohr z/bohr   ! (i4,4f13.5) per atom; charge is dummy val_Q=1.0 (plot.F90:2137)
...
<data: (6E13.5), z fastest>
```

Positive grid counts ⇒ Bohr units per the cube convention. Atomic numbers come
from a lowercase 109-entry `periodic_table` lookup on `atom_data%symbol`
(plot.F90:2110–2121, 2140–2149; unmatched species get Z=0). Data loop
(plot.F90:2377–2385): x outer, then y, then **z fastest**, 6 values per line,
last line of each z-run may be short:

```fortran
do nxx = 1, ilength(1)
  do nyy = 1, ilength(2)
    do nzz = 1, ilength(3), 6
      nend = min(nzz + 5, ilength(3))
      write (file_unit, '(6E13.5)') wann_cube(nxx, nyy, nzz:nend)
```

---

## D. Keyword summary

| keyword | default | where parsed / defined | effect |
|---|---|---|---|
| `wannier_plot` | `.false.` | wannier90_readwrite.F90:414 | master switch; `plot_main` → `plot_wannier` (plot.F90:321–326) |
| `wannier_plot_list` | all WFs 1..num_wann | readwrite:1204–1250 (default fill at 1247–1249) | range vector (e.g. `1-4,6`); validated 1..num_wann (readwrite:1227) |
| `wannier_plot_supercell` | (2,2,2) | wannier90_types.F90:119; readwrite:1164–1188 (1 or 3 ints, all > 0) | supercell dims; grid bounds plot.F90:1674–1679 |
| `wannier_plot_format` | `'xcrysden'` | types:122; readwrite validation 1260–1263 | `'xcrysden'` → .xsf, `'cube'` → .cube (plot.F90:2020–2027) |
| `wannier_plot_mode` | `'crystal'` | types:123; readwrite:1194, validation 1264–1266 | xsf: ATOMS vs CRYSTAL header (plot.F90:2445); cube: origin shift vs atom replication (plot.F90:2132–2135) |
| `wannier_plot_radius` | 3.5 (Å) | types:120; readwrite:1252, must be ≥ 0 (:1275) | cube only: half-extent of cut-out box (plot.F90:2189/2193) |
| `wannier_plot_scale` | 1.0 | types:121; readwrite:1255, ≥ 0 | cube crystal mode only: atom-inclusion radius = scale*radius (plot.F90:2314, 2364) |
| `wannier_plot_spinor_mode` | `'total'` | types:124; readwrite:1197, validation 1267–1272 | spinor collapse: total / up / down (plot.F90:1952–1963) |
| `wannier_plot_spinor_phase` | `.true.` | types:125; readwrite:1200 | attach sign(Re) to up/down magnitudes (plot.F90:1946–1951) |
| `wvfn_formatted` | `.false.` | types:132; readwrite:1045 | UNK files formatted vs unformatted |
| `spin` | `'up'` (channel 1) | readwrite:1049–1060 | selects UNKppppp.1 vs .2 (non-spinor) |
| `translate_home_cell` | `.false.` | readwrite:1372 | NOT used by volumetric plotting; only `seedname_centres.xyz` (plot.F90:3174) |

## E. Implementation checklist / gotchas

1. Grid index n ↔ fractional coordinate (n−1)/ng; index 1 is the origin. All
   phases and origins carry this −1 (plot.F90:1880, 2418–2426).
2. UNK point ordering is x-fastest with `npoint = nx + (ny-1)*ngx + (nz-1)*ngx*ngy`.
3. Disentangled case: UNK contains all `num_bands` bands; keep only
   `lwindow(:,ik)` bands (in order), then apply `u_matrix_opt` (ndimwin ×
   num_wann), then `u_matrix` (num_wann × num_wann).
4. Divide by `num_kpts` exactly once, after summing all k.
5. Phase-fix each scalar WF: divide by unit-modulus value at the max-|w|²
   point → max point real positive; report max Im/Re over |Re| ≥ 0.01.
6. xsf: Å, whole supercell, N points (no +1) with spanning vectors scaled by
   (N−1)/ng, real part, x fastest, `(6e13.5)`.
7. cube: Bohr, radius-cut box around the WF centre, real part, z fastest,
   `(6E13.5)`, dummy atomic charge 1.0; fails if box leaves the computed
   supercell.
8. Spinor output is a magnitude (optionally signed), so no phase fix; both
   formats then plot it as the "real part".
