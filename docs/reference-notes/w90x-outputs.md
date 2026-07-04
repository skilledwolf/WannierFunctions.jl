# wannier90.x reference spec: small output features (cube / _r.dat / .bxsf / hr_diag / xyz)

Implementation-grade spec extracted from the reference Wannier90 source. All paths are relative
to `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`. Line numbers refer to those files.
All four features live in `src/plot.F90` (driver `plot_main`, plot.F90:48-330); keyword parsing in
`src/wannier90_readwrite.F90`; defaults in `src/wannier90_types.F90`.

Constants: reference build defaults to CODATA2006 (constants.F90:92-100 `#define CODATA2006`):
`bohr = bohr_angstrom_internal = 0.52917720859` Å (constants.F90:182, aliased at :224).

Cross-reference: the full WF-assembly pipeline feeding the cube writer (UNK reading, U/U_opt
application, Bloch sum, 1/N_k, phase fixing, spinor collapse) is specced in
`docs/reference-notes/wannier-plot.md` sections A-B; the cube layout there (C.2) agrees with
section 1 below.

---

## 1. Gaussian cube output (`wannier_plot = T`, `wannier_plot_format = cube`)

Writer: `internal_cube_format` (contained in `plot_wannier`), plot.F90:2075-2405. Dispatch at
plot.F90:2020-2028: `wannier_plot%format .eq. 'xcrysden'` → xsf, `.eq. 'cube'` → cube, else warn
`'wannier_plot_format not recognised in plot_wannier'`. Note the dispatch is an **exact** string
compare even though input validation only checks `index(format,'cub') > 0`
(wannier90_readwrite.F90:1260-1263) — only the literal value `cube` works.

### 1.1 Keywords and defaults

| keyword | default | file:line |
|---|---|---|
| `wannier_plot_format` | `'xcrysden'` | wannier90_types.F90:122; parsed wannier90_readwrite.F90:1191 |
| `wannier_plot_radius` | `3.5_dp` (Å) | wannier90_types.F90:120; parsed :1252; rejected if `< 0` (:1274) |
| `wannier_plot_scale` | `1.0_dp` | wannier90_types.F90:121; parsed :1255; rejected if `< 0` (:1278) |
| `wannier_plot_mode` | `'crystal'` | wannier90_types.F90:123; parsed :1194; must contain 'crys' or 'mol' (:1264) |
| `wannier_plot_supercell` | `(2,2,2)` | wannier90_types.F90:119; parsed :1164-1188 (1 int → all three, or 3 ints; all > 0) |
| `wannier_plot_list` | all WFs 1..num_wann | readwrite:1204-1250 (default fill 1247-1249) |
| `wannier_plot_spinor_mode` | `'total'` | wannier90_types.F90:124; parsed :1197 ('total'/'up'/'down', :1268-1273) |
| `wannier_plot_spinor_phase` | `.true.` | wannier90_types.F90:125; parsed :1200 |

Mode flags (plot.F90:2132-2135): `lmol = index(mode,'mol')>0`, `lcrys = index(mode,'crys')>0`
(only lmol is actually branched on).

### 1.2 Filename

```fortran
202 format(a, '_', i5.5, '.cube')                    ! plot.F90:2151
write (wancube, 202) trim(seedname), wann_index      ! plot.F90:2181
```

→ `seedname_00001.cube`, one file per entry of `wannier_plot_list` (`wann_index` is the WF
index from the list, not the loop counter).

### 1.3 Box geometry (all in Å until the final /bohr on write)

Inputs: converged centres `wannier_data%centres(:, wann_index)` (Cartesian Å);
`moda(i) = |a_i|`, `modb(i) = |b_i|` (plot.F90:2153-2161); `dgrid(i) = moda(i)/ng_i`
(plot.F90:2164) with `(ngx,ngy,ngz)` the UNK grid; `recip_lattice` includes 2π.

Extent along each lattice direction i (plot.F90:2184-2194) — centre projected onto direction
a_i, ± radius converted to a length along a_i:

```
proj(i)   = (Σ_c centre(c)*recip_lattice(i,c)) * moda(i)/twopi        ! = frac coord_i * |a_i|
rstart(i) = proj(i) - twopi*radius/(moda(i)*modb(i))
rend(i)   = proj(i) + twopi*radius/(moda(i)*modb(i))
```

Integerisation (plot.F90:2196-2201):

```
ilength(i) = ceiling((rend(i)-rstart(i))/dgrid(i))
istart(i)  = floor(rstart(i)/dgrid(i)) + 1
iend(i)    = istart(i) + ilength(i) - 1               ! only used for iprint>3 debug output
```

Cube origin, Cartesian Å (plot.F90:2203-2208): grid index n ↔ fractional (n−1)/ng, so

```
orig(c) = Σ_i (istart(i)-1) * dgrid(i) * real_lattice(i,c)/moda(i)    ! = Σ_i (istart(i)-1)/ng_i * a_i
```

`real_lattice(i,c)`: row i = lattice vector a_i, column c = Cartesian component.

### 1.4 Volumetric data extraction (plot.F90:2233-2280)

`wann_cube(nxx,nyy,nzz) = real(wann_func(qxx,qyy,qzz,loop_w), dp)` — **real part only**
(plot.F90:2277) of the assembled supercell WF (already /num_kpts and phase-fixed for scalars;
for spinors it is the signed/unsigned modulus per `spinor_mode`/`spinor_phase`, plot.F90:1936-1970 —
the cube then contains that non-negative (or sign(Re)-signed) magnitude).

Index mapping with fold-up-from-below only:

```fortran
qzz = nzz + istart(3) - 1
izz = int((abs(qzz) - 1)/ngz)
if (qzz .lt. nzz_lo) qzz = qzz + izz*ngz     ! nzz_lo = -((ngs(3))/2)*ngz  (plot.F90:1674-1679)
if (qzz .gt. nzz_hi) → warn-error            ! nzz_hi = ((ngs(3)+1)/2)*ngz - 1
```

(same for y, x; plot.F90:2235-2276). The error text (verbatim, 4 stdout lines then
`set_error_warn(error, 'Error plotting WF cube.', comm)`):

```
 Error plotting WF cube. Try one of the following:
    (1) increase wannier_plot_supercell;
    (2) decrease wannier_plot_radius;
    (3) set wannier_plot_format=xcrysden
```

### 1.5 Atoms

Atomic numbers: lowercase 109-entry `periodic_table` matched against `atom_data%symbol(isp)`
(plot.F90:2110-2121, 2140-2149); unmatched species keep `atomic_Z = 0`. Charge column is the
dummy `val_Q = 1.0_dp` (plot.F90:2137).

- **molecule mode** (`lmol`, plot.F90:2282-2301, 2353-2354): compute fractional WF centre `wcf`,
  fractional atomic centre-of-mass `comf` (plot.F90:2166-2175), `irdiff = nint(comf - wcf)`;
  translate the cube origin `orig += Σ_i irdiff(i)*a_i` so the box lands beside the input atoms;
  write all `num_atoms` atoms exactly as in the input.
- **crystal mode** (default, plot.F90:2302-2323 count pass, 2355-2373 write pass): replicate each
  atom over image cells `nxx,nyy,nzz ∈ [−ngs(i)/2, (ngs(i)+1)/2]` (integer division; inclusive —
  for ngs=2 that is −1..1, i.e. one more cell on top than the data supercell) and keep it iff its
  Cartesian distance from the **WF centre** satisfies
  `dist ≤ wannier_plot%scale * wannier_plot%radius` (plot.F90:2314, 2364). The count pass yields
  `icount` for the header; the write pass re-runs the identical loop.

### 1.6 File layout (plot.F90:2325-2385) — "everything in Bohr"

```
      Generated by Wannier90 code http://www.wannier.org        ! write(*, *) '     Generated by ...' (list-directed → 1 extra leading blank)
      On <cdate> at <ctime>                                     ! '     On ', cdate(9), ' at ', ctime(9)
natoms origx origy origz                                        ! '(i4,3f13.5)'; natoms = num_atoms (mol) | icount (crys); orig(:)/bohr
ilength(1) a1x/(ngx*bohr) a1y/(ngx*bohr) a1z/(ngx*bohr)         ! '(i4,3f13.5)' voxel vector 1
ilength(2) a2/(ngy*bohr) ...                                    ! '(i4,3f13.5)'
ilength(3) a3/(ngz*bohr) ...                                    ! '(i4,3f13.5)'
Z 1.00000 x/bohr y/bohr z/bohr                                  ! '(i4,4f13.5)' one line per (kept) atom
...
<volumetric data '(6E13.5)'>
```

`cdate/ctime` from `io_date` (io.F90:322-345): `cdate` = `(i2,a3,i4)` e.g. `29Aug2018`
(blank-padded day < 10), `ctime` = `(i2.2,":",i2.2,":",i2.2)` in a len-9 string (1 trailing blank).

Data ordering (plot.F90:2377-2385): x outer, y middle, **z fastest**, 6 values per line, short
last line per z-run:

```fortran
do nxx = 1, ilength(1)
  do nyy = 1, ilength(2)
    do nzz = 1, ilength(3), 6
      nend = min(nzz + 5, ilength(3))
      write (file_unit, '(6E13.5)') wann_cube(nxx, nyy, nzz:nend)
```

`E13.5` renders like `  0.26921E-01` / ` -0.14234E-01` (13 chars, leading `0.`, no E-field width).

### 1.7 Oracle: test-suite/tests/testw90_cube_format

`gaas.win` (verbatim keywords): `num_wann = 4`, `num_iter = 20`, `use_ws_distance = .false.`,
`search_shells=12`, `unit_cell_cart` in **bohr** (−5.367 0 5.367 / 0 5.367 5.367 / −5.367 5.367 0),
`atoms_frac` Ga (0,0,0), As (0.25,0.25,0.25), projections `As:sp3`, `mp_grid : 2 2 2` (8 explicit
kpoints), `wvfn_formatted=.true.`, and

```
wannier_plot = true
wannier_plot_list = 1
wannier_plot_supercell = 2
wannier_plot_format = cube
wannier_plot_radius = 2
wannier_plot_mode = crystal
```

UNK00001.1..UNK00008.1 are provided (formatted). jobconfig (tests/jobconfig:250-253):
`program = WANNIER90_CUBE`, `output = gaas_00001.cube`. userconfig (tests/userconfig:75-78):
parser `parse_cube.py` extracts **only line index 2** (the natoms+origin line) as 4 floats keyed
`'origin'`; tolerance `(1.0e-5, 1.0e-5, 'origin')` (abs, rel). Benchmark line 3:
`   2      3.22020     -1.07340     -1.07340` (2 atoms within 2 Å of the WF centre; 17³ box).
So the regression check is natoms + origin only — but produce the full file to match the
benchmark byte layout above.

---

## 2. `seedname_r.dat` (`write_rmn = T`) — position matrix elements

Writer: `plot_write_rmn`, plot.F90:2613-2711, called from `plot_main` (plot.F90:301-306).
Default `write_rmn = .false.` (wannier90_types.F90:68; parsed wannier90_readwrite.F90:988-989).
Setting `write_rmn` forces `hamiltonian_setup` (plot.F90:148-161) which supplies the
Wigner-Seitz R list `irvec(1:3, 1:nrpts)`/`nrpts` (hamiltonian.F90:147-148; same list as
`_hr.dat`) — but **ndegen is neither applied nor written** (see traps).

### 2.1 Quantity

For every R (all nrpts WS vectors) and every WF pair, with `M = m_matrix` the **final
Wannier-gauge overlaps** `M^(W)_{nm}(k,b) = <u^W_{nk}|u^W_{m,k+b}>` (first index = bra):

```
fac(k,R) = exp(-i*2π k·R) / N_k                                       ! plot.F90:2682-2683

n ≠ m (WYSV06 Eq. 44, linear form; plot.F90:2696-2698):
  <0n|r_α|Rm> = Σ_k fac(k,R) * i * Σ_b w_b b_α(k) * M^(W)_{nm}(k,b)

n = m (Im-ln form, all R; plot.F90:2686-2694):
  <0n|r_α|Rn> = − Σ_k fac(k,R) * Σ_b w_b b_α(k) * Im ln M^(W)_{nn}(k,b)
```

Comments in source: for R=0 the diagonal reduces to MV97 Eq.(32); otherwise Eq.(44) of
WYSV06 modified per MV97 Eqs.(27,29). `w_b = kmesh_info%wb(nn)` (Å²),
`b_α = kmesh_info%bk(ind,nn,nkp)` (Cartesian Å⁻¹) → positions in **Å**. Note the diagonal
result is real only at R=0 aggregate; in general `position` is complex and both parts are
written. The imaginary unit sits **outside** M (no conjugation anywhere); the diagonal uses
`aimag(log(M))` multiplied by the complex `fac`.

Index order trap: loops are `loop_rpt` (outer) → `m` → `n` (inner) (plot.F90:2673-2675), i.e.
**n varies fastest** in the file; the code accesses `m_matrix(n, m, nn, nkp)` and writes columns
`n, m` in that order (plot.F90:2705).

### 2.2 File layout

```
 written on <cdate> at <ctime>       ! character(len=33) header, list-directed write → leading blank (plot.F90:2666-2668)
           4                         ! write(unit,*) num_wann  → gfortran renders default integer in width 12 (plot.F90:2669)
          93                         ! write(unit,*) nrpts                                     (plot.F90:2670)
   -3    1    1    1    1    0.000949   -0.000000   -0.000949    0.000000   -0.002102   -0.000000
...
```

Data record (plot.F90:2705):

```fortran
write (file_unit, '(5I5,6F12.6)') irvec(:, loop_rpt), n, m, position(:)
```

`position(3)` is complex → 6 reals in order Re(x) Im(x) Re(y) Im(y) Re(z) Im(z). Total lines =
`3 + nrpts*num_wann**2` (diamond benchmark: 3 + 93*16 = 1491 ✓).

### 2.3 Oracle: test-suite/tests/testw90_rmn

`diamond.win`: `num_wann = 4`, `num_iter = 20`, `write_rmn = .true.`, `iprint = 0`, four s
projections at f=(0,0,0),(0,0,1/2),(0,1/2,0),(1/2,0,0), C atoms at ±(1/8,1/8,1/8) frac,
fcc cell ±1.61399 Å entries, `mp_grid = 4 4 4` (64 explicit kpoints). jobconfig
(tests/jobconfig:521-524): `program = WANNIER90_RMN_OK`, `output = diamond_r.dat`. userconfig
(tests/userconfig:58-73), parser `parse_rmn.py` reads `num_wann` (line 2), `nrpts` (line 3),
then all 11 columns of every data line. Tolerances: `irvec_a/b/c`, `index_i` (=col 4 = n),
`index_j` (=col 5 = m) at `(1.0e-10, 1.0e-10)` (exact integers); `real_x/imag_x/real_y/imag_y/
real_z/imag_z` at `(2.0e-6, 1.0)` — absolute 2e-6 (6-dp printout), relative effectively ignored.

---

## 3. `seedname.bxsf` (`fermi_surface_plot = T`) — XCrySDen Fermi surface

Writer: `plot_fermi_surface`, plot.F90:1360-1544 (module `w90_plot_mod` in src/plot.F90 — not a
separate module), called from `plot_main` at plot.F90:210-216 on root only. Requires
`hamiltonian_setup` + `hamiltonian_get_hr` (plot.F90:148-188; Γ-in-mesh warning at :175-178).

Keywords: `fermi_surface_plot` default `.false.` (wannier90_types.F90:55; parsed
wannier90_readwrite.F90:422-423); `fermi_surface_num_points` default `50`
(wannier90_types.F90:226; parsed :1301-1303, rejected only if `< 0`, :1314);
`fermi_surface_plot_format` default `'xcrysden'` (wannier90_types.F90:227; parsed :1305-1307;
must contain 'xcrys' — it is the only format, nothing else reads this string).

Fermi energy: `fermi_energy_list` is **always allocated** with n=1, value 0.0_dp when no
`fermi_energy` keyword is given (readwrite.F90:609-693, default n=1 at :634, fill at :690-692).
If a `fermi_energy_min/max` scan produced size > 1 → fatal
`"Error in plot: nfermi>1. Set the fermi level using the input parameter 'fermi_level'"`
(plot.F90:1414-1420; note the message names a keyword `fermi_level` that does not exist —
the real keyword is `fermi_energy`).

### 3.1 Grid and eigenvalues

`npts_plot = (num_points+1)**3` (plot.F90:1458). Loop order (plot.F90:1467-1501): `loop_x`
outer, `loop_y`, `loop_z` **innermost**; flattened index `ikp` increments in that order. The
k-point is fractional `k' = ((loop_x-1), (loop_y-1), (loop_z-1)) / num_points` — both endpoints
0 and 1 included (bxsf "general grid" with periodic duplicate).

H(k') interpolation (plot.F90:1473-1480) — note **+i** exponent and division by `ndegen` here:

```
rdotk = 2π * ((loop_x-1)*irvec(1,R) + (loop_y-1)*irvec(2,R) + (loop_z-1)*irvec(3,R)) / num_points
H(k') = Σ_R  (cos(rdotk) + i sin(rdotk)) / ndegen(R) * ham_r(:,:,R)
```

where `ham_r(:,:,R) = (1/N_k) Σ_k e^{-i2πk·R} H^(W)(k)` (hamiltonian.F90:422-428) and
`H^(W)(k) = U(k)† E(k) U(k)` (from `hamiltonian_get_hr`). No use_ws_distance anywhere in this
path. Eigenvalues via `ZHPEVX('N','A','U', ...)` on the packed upper triangle
(plot.F90:1482-1488), stored `eig_int(1:num_wann, ikp)` in eV; all `num_wann` bands written.

### 3.2 File layout (plot.F90:1503-1534) — all `write(bxsf_unit, *)` except the energies

```
  BEGIN_INFO
       #
       # this is a Band-XCRYSDEN-Structure-File
       # for Fermi Surface Visualisation
       #
       # Generated by the Wannier90 code http://www.wannier.org
       # On <cdate>  at <ctime>          ! list-directed: '      # On ', cdate, ' at ', ctime
       #
       Fermi Energy:  <fermi_energy_list(1)>   ! list-directed real (gfortran: ~17 sig digits)
  END_INFO
<blank line>
  BEGIN_BLOCK_BANDGRID_3D
 from_wannier_code
  BEGIN_BANDGRID_3D_fermi
 <num_wann>                              ! list-directed integer
 <np+1> <np+1> <np+1>                    ! np = fermi_surface_num_points
 0.0 0.0 0.0                             ! origin (literal string '0.0 0.0 0.0')
 <b1x b1y b1z>                           ! (recip_lattice(1,i), i=1,3)  — rows = b-vectors, Å⁻¹, includes 2π
 <b2x b2y b2z>
 <b3x b3y b3z>
 BAND:  <i>                              ! list-directed: 'BAND: ', i   — for i = 1..num_wann
 <eig>                                   ! write(unit,'(2E16.8)') eig_int(i, loop_kpt) — ONE value per line
 ...                                     ! loop_kpt = 1..npts_plot in ikp order (z fastest)
 END_BANDGRID_3D
  END_BLOCK_BANDGRID_3D
```

Every list-directed line carries one extra leading blank; the strings above show the in-code
literals (e.g. `' BEGIN_INFO'` → two leading spaces on disk). The energies use format
`(2E16.8)` but only one item is passed → exactly one 16-char field per line.
Stdout messages: `'Calculating Fermi surface'` (plot.F90:1411) and
`'Time to calculate interpolated Fermi surface ', f11.3, ' (sec)'` (plot.F90:1536-1537).
No test-suite oracle exists for .bxsf.

---

## 4. `write_hr_diag` — on-site H elements (stdout only)

plot.F90:250-261, inline in `plot_main`; default `.false.` (wannier90_types.F90:63; parsed
wannier90_readwrite.F90:973-974). Setting it forces `hamiltonian_setup` + `hamiltonian_get_hr`
(plot.F90:148-188). **No file is written** — output goes to `seedname.wout`, and only when
`print_output%iprint > 0` (default iprint = 1):

```fortran
write (stdout, *)
write (stdout, '(1x,a)') 'On-site Hamiltonian matrix elements'
write (stdout, '(3x,a)') '  n        <0n|H|0n> (eV)'
write (stdout, '(3x,a)') '-------------------------'
do i = 1, num_wann
  write (stdout, '(3x,i3,5x,f12.6)') i, real(ham_r(i, i, rpt_origin), kind=dp)
end do
write (stdout, *)
```

`rpt_origin` = index of R = (0,0,0) in the WS list (hamiltonian.F90:813); `ndegen(rpt_origin)`
is always 1, and no degeneracy division is applied (consistent: `_hr.dat` also stores raw
`ham_r`). Values are eV; real part only.

---

## 5. `write_xyz` / `translate_home_cell` — `seedname_centres.xyz`

Writer: `plot_write_xyz`, plot.F90:3174-3246, called from `plot_main` on root
(plot.F90:276-280). Defaults: `write_xyz = .false.` (wannier90_types.F90:70; parsed
wannier90_readwrite.F90:961-962); `translate_home_cell = .false.`
(wannier90_types.F90:95, member of `real_space_ham_type`; parsed wannier90_readwrite.F90:1372-1373).

`translate_home_cell = T` maps each centre to the home cell before writing via
`utility_translate_home` (utility.F90:609-649): Cartesian → fractional, then per component:
if `f < 0`: `f += ceiling(|f|)`; if `f > 1`: `f -= int(f)` (i.e. rationalised to [0,1] — a
component exactly in [0,1] is untouched, and f=1.0 stays 1.0), then back to Cartesian.
It affects **only this xyz file** (and, independently, transport's `tran_write_xyz`); it does
NOT affect volumetric plotting, `_hr.dat`, or interpolation (the `use_translation` machinery in
hamiltonian.F90:104-109 is driven by `bands_plot_mode`/`transport_mode`, not by this keyword).

File layout (plot.F90:3224-3240), standard xyz, Cartesian Å:

```
<num_wann + num_atoms>                       ! '(i6)'
 Wannier centres, written by Wannier90 on<cdate> at <ctime>    ! list-directed; NOTE: no space between "on" and cdate
X      <x> <y> <z>                           ! '("X",6x,3(f14.8,3x))' one line per WF centre
...
<label> <x> <y> <z>                          ! '(a2,5x,3(f14.8,3x))' per atom; label truncated to 2 chars
...
```

Atom loop: species outer, atoms inner (plot.F90:3235-3239), positions `atom_data%pos_cart`
as given in input (never translated). Comment-line date quirk: the string is
`'Wannier centres, written by Wannier90 on'//cdate//' at '//ctime` — for days ≥ 10 it reads
literally `...on18Jul2024 at ...`; for days < 10 the cdate's leading blank supplies the space.
Also writes to stdout: `'(/a)') ' Wannier centres written to file '//trim(seedname)//'_centres.xyz'`
(plot.F90:3242) and, if `iprint > 2`, a "Final centres" block scaled by `lenconfac`
(plot.F90:3214-3222). No test-suite oracle.

---

## Traps (all features)

1. **cube**: `wannier_plot_format` must be exactly `cube` (dispatch string-equality,
   plot.F90:2022) even though input validation passes anything containing `cub`.
2. **cube**: data is the real part of the phase-fixed WF for scalars but a (possibly
   sign(Re)-signed) spinor magnitude for `spinors=T`; division by num_kpts happens in
   different places (plot.F90:1964 spinor vs :1983 scalar) but exactly once either way.
3. **cube**: box indices fold **up only** — a box extending above the supercell errors out
   rather than wrapping (plot.F90:2240-2276). Atom images in crystal mode scan one cell
   further (+1 on the upper side) than the data supercell (plot.F90:2307-2309).
4. **cube**: `natoms` field is `i4` (overflows ≥ 10000); charge column is always `1.00000`;
   unknown species → `Z = 0` silently.
5. **_r.dat**: no `ndegen` anywhere — consumers must take the degeneracy list from
   `seedname_hr.dat` for R-sums. R list = same WS list as `_hr.dat` (so `use_ws_distance`
   never enters).
6. **_r.dat**: column order is `R1 R2 R3 n m` with **n fastest** (inner loop) — parser calls
   col-4 `index_i`, col-5 `index_j`. The diagonal (n=m) uses the Im-ln formula for **every** R,
   not just R=0, and its `fac` is complex → nonzero imaginary parts on diagonal lines are
   expected.
7. **_r.dat**: `m_matrix` here is the post-wannierisation (Wannier-gauge) overlap from the
   minimiser, not the raw `.mmn` — matches gauge-invariance validation policy: compare only
   physically meaningful aggregates when U differs.
8. **bxsf**: interpolation exponent is **+i2πk'·R** (forward, cos+isin) with 1/ndegen, on top of
   ham_r built with e^{−i2πk·R}/N_k. Grid has num_points+1 points per axis including the
   duplicated endpoint; z is the fastest flattened index. One energy per line despite the
   `(2E16.8)` descriptor.
9. **bxsf**: missing `fermi_energy` silently writes `Fermi Energy: 0.0` (default list [0.0]);
   a fermi scan (list size > 1) is a fatal error whose message suggests the nonexistent
   keyword `fermi_level`. `fermi_surface_num_points = 0` passes validation (`< 0` check) but
   divides by zero in `rdotk`.
10. **write_hr_diag**: stdout only, gated on `iprint > 0`; nothing to diff except `.wout`.
11. **xyz**: `translate_home_cell` rationalises to [0,1] (not [0,1)); components exactly 1.0
    are kept. The xyz atom label uses `label` (not `symbol`) truncated to 2 characters by the
    `a2` descriptor.
12. List-directed writes (`write(unit,*)`) add one leading blank before the first item; integer
    field widths are compiler-dependent (gfortran: default integer in width 12 — benchmarks
    were produced with gfortran, e.g. `           4` in `_r.dat`). To match file-precision,
    emulate gfortran list-directed formatting.
