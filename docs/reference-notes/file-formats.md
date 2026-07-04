# Wannier90 v3.1.0 — File Format Reference (implementation-grade)

Scope: exact on-disk layout of every Wannier90 file, from the reference source at
`reference/wannier90`. Every format string / record layout is cited `file:line`.
Confirmed against shipped example `tutorials/tutorial03/silicon.*`.

Primary source files:
- `src/readwrite.F90` — `.win` scalar/block parsing, `.eig` reader, units.
- `src/wannier90_readwrite.F90` — `.win` projections / kpath, `.wout` system/kmesh echo.
- `src/overlap.F90` — `.amn` / `.mmn` readers (`overlap_read`), `*_dump` writers.
- `src/kmesh.F90` — `.nnkp` writer (`kmesh_write`), b-shell report.
- `src/hamiltonian.F90` — `seedname_hr.dat`, `seedname_tb.dat` writers.
- `src/plot.F90` — `seedname_r.dat`, `seedname_band.{dat,kpt,gnu}`, `_band.labelinfo.dat`.
- `src/w90chk2chk.F90` — `.chk` (unformatted) and `.chk.fmt` (formatted) read+write, side by side. **The canonical reference for the `.chk` record layout.**
- `src/wannierise.F90` — `.wout` spread trace + centres/spread lines.
- `docs/docs/user_guide/wannier90/files.md` — official prose description.

---

## 0. CRITICAL cross-cutting gotchas (read first)

1. **Binary `.chk` uses Fortran sequential-unformatted records.** Each `read(unit) X`
   statement is ONE record. Fortran wraps every record with a record-length marker
   (typically a 4-byte integer before AND after each record; size/endianness is
   compiler/runtime dependent — gfortran uses 4-byte subrecord markers, may switch to
   8-byte for records >2 GiB). A multi-item read such as
   `read(chk_unit) ((real_lattice(i,j),i=1,3),j=1,3)` is a **single record of 9 reals**,
   not 9 records. A from-scratch Julia reader must consume these markers. **Strongly
   prefer `.chk.fmt`** (produce it with `w90chk2chk`) as the robust interchange path.
2. **`hr.dat` is lossy: matrix rows use `'(5I5,2F12.6)'` — only 6 decimals.** You cannot
   reproduce the reference to full precision *from* `hr.dat`. Precision by file:
   `hr.dat` = F12.6 (6 dp); `tb.dat`/`r.dat` = `E15.8` (8 sig figs); `.chk`/`.chk.fmt`
   = `G25.17` (full double). **The `.chk`/`.chk.fmt` is the only full-precision source.**
3. **`hr.dat`/`tb.dat`/`r.dat` header/count lines are list-directed** (`write(unit,*)`):
   leading whitespace, no fixed columns (date line, `num_wann`, `nrpts`, and in `tb.dat`
   the 3 lattice-vector lines). Only the degeneracy block (`15I5`) and matrix rows are
   fixed-format. **Parse by tokenizing, never by column slicing.**
4. **`ndegen` weighting:** `hr.dat`/`tb.dat` store `H(R)` **undivided**. The k-space
   interpolation consumer must divide each `H(R)` by `ndegen(R)`. WS sum rule:
   `Σ_R 1/ndegen(R) == mp_grid(1)*mp_grid(2)*mp_grid(3)` (checked to `eps8`,
   `hamiltonian.F90:850`).
5. **Im-ln branch cut.** Position matrix elements (`r.dat`, `tb.dat`) use
   `aimag(log(M))` = `atan2(Im,Re) ∈ (−π, π]`; `fac = exp(−i·rdotk)/num_kpts`,
   `rdotk = twopi·dot(kpt_latt(:,k), irvec(:,R))`. Reproduce the atan2 branch exactly.
6. **Lattice transpose.** `unit_cell_cart` is read as `real_lattice_tmp(row,col)` then
   `real_lattice = transpose(real_lattice_tmp)` (`readwrite.F90:1125`). Thus after read,
   **`real_lattice(i,:)` is lattice vector a_i** (row-major storage of the a-vectors).
   `.chk` stores `real_lattice` in this post-transpose orientation.
7. **Three different reader semantics for `.eig`/`.amn`/`.mmn`** — see §4.

Bohr constant (default `constants.F90`, v3 block): `bohr_angstrom_internal = 0.529177210544`
(older `bohr = 0.5291772108` also present; version-selected). `lenconfac = 1` for Ang,
`= 1/bohr` for Bohr (`readwrite.F90:181`). Energy always eV. Spreads always Å² internally,
scaled by `lenconfac**2` on output.

---

## 1. `seedname.win` (INPUT, formatted, free-form)

Master input. Case-insensitive keywords; `!` and `#` start comments; `:` and `=`
are token separators for scalars. Blocks are `begin <name>` … `end <name>`
(case-insensitive; `Begin`/`End` also accepted). Order of blocks/keywords is free.

### Units convention (`readwrite.F90:154-185`, `docs .../files.md`)
- `length_unit` = `ang` (default) or `bohr` — controls **`.wout` output** length unit only.
- Blocks `unit_cell_cart`, `atoms_cart`, `projections` may have an **optional first line**
  inside the block naming the unit: `ang` (default) or `bohr`. If present it is consumed
  before the numeric rows (`lunits` flag, `get_block_length(..., lunits)`).
- Energy windows: always eV. k-points: always fractional (crystallographic).
  `zona`: always Å⁻¹. Convergence thresholds: always Å².

### Scalar keywords (key ones)
- `num_bands` (int; defaults to `num_wann` if absent), `num_wann` (int, mandatory).
- `mp_grid = n1 n2 n3` (mandatory; `num_kpts = n1*n2*n3`, `readwrite.F90:418`).
- `num_iter`, `num_print_cycles`, `dis_num_iter`, `dis_mix_ratio`.
- `dis_win_min`, `dis_win_max`, `dis_froz_min`, `dis_froz_max` (eV).
- `bands_plot` (logical), `bands_num_points`, `write_hr`, `write_tb`, `write_rmn`,
  `write_xyz`, `write_bvec`, `use_ws_distance`, `guiding_centres`, `iprint`, `timing_level`.
Fortran logicals accept `.true./.false./T/F/true/false`; doubles accept `1.d0`, `17.0d0`.

### `begin unit_cell_cart … end unit_cell_cart`
Optional unit line (`ang`/`bohr`), then exactly **3 rows of 3 reals**: rows are the
Cartesian components of a_1, a_2, a_3. Stored transposed (see gotcha 6). Reciprocal
lattice built as `B_i = 2π/V (A_j × A_k)`, `V = A_1·(A_2×A_3)`.

### `begin atoms_frac … end atoms_frac` / `begin atoms_cart … end atoms_cart`
Each row: `<symbol> f1 f2 f3` (frac) or `<symbol> x y z` (cart, with optional unit line).
`atoms_frac` internally converted to Cartesian via `utility_cart_to_frac` with
`transpose(real_lattice)` (`readwrite.F90:1213`). Exactly one of the two may appear.
tutorial03 example:
```
Begin Atoms_Frac
Si  -0.25   0.75  -0.25
Si   0.00   0.00   0.00
End Atoms_Frac
```

### `begin projections … end projections`
Mini-language (parsed in the projections module; `wannier90_readwrite.F90:1534`
`..._read_projections`). Syntax per line (site : states):
```
<site> : <angular> [ : z = zx,zy,zz ] [ : x = xx,xy,xz ] [ : r=<n> ] [ : zona=<val> ]
```
- `<site>` = element symbol (e.g. `Si`), `f=fx,fy,fz` (fractional), `c=cx,cy,cz`
  (Cartesian), or `random`.
- `<angular>` = hydrogenic labels: `s`, `p`, `px`, `pz`, `sp3`, `sp2`, `d`, `dxy`,
  `f`, … or explicit `l=<L>,mr=<M>` (quantum numbers).
- `z=`/`x=` set the local axis frame (defaults z=(0,0,1), x=(1,0,0)); `r=` radial
  function index (default 1); `zona` = Z/a in Å⁻¹ (default 1.0).
The parsed result is emitted verbatim into the `projections` block of `.nnkp` (§6), where
each projector becomes: `site(3) l mr radial` + `z(3) x(3) zona`. tutorial03: `Si : sp3`
(expands to 4 sp3 projectors → 8 total for 2 Si atoms).

### `begin kpoints … end kpoints`
First (implicit) count = `num_kpts` (must equal `∏mp_grid`). Each row: `k1 k2 k3`
(fractional). Read via `get_keyword_block(..., 'kpoints', num_kpts, 3, ...)`
(`readwrite.F90:1028`). Order is the k-index order used everywhere downstream.

### `begin kpoint_path … end kpoint_path` (`readwrite.F90:479`)
Each line = **two** special points: `L1 f1 f2 f3  L2 g1 g2 g3`. `bands_num_spec_points
= 2 * (#lines)`. Labels+coords stored in `points(3, 2*npath)`, `labels(2*npath)`.
tutorial03:
```
begin kpoint_path
L 0.50000 0.50000 0.5000 G 0.00000 0.00000 0.0000
G 0.00000 0.00000 0.0000 X 0.50000 0.00000 0.5000
...
end kpoint_path
```
`bands_num_points` (default segment resolution) set separately as a scalar.

### `begin mp_grid` — no; `mp_grid` is a scalar keyword, not a block.

---

## 2. `seedname.eig` (INPUT, formatted)

Written by the DFT interface (pw2wannier90); W90 **reads** it (`readwrite.F90:759`).
One eigenvalue per line: `band_index  kpoint_index  eigenvalue_eV`. Loop nesting is
**k outer, band inner** → band index is fastest. Total lines = `num_bands * num_kpts`.
```fortran
do k = 1, num_kpts
  do n = 1, num_bands
    read(eig_unit,*) i, j, eigval(n,k)   ! aborts if i/=n or j/=k
```
**Positional AND enforced:** the reader hard-aborts on any `i≠n` or `j≠k` mismatch
(`readwrite.F90:768`). eigenvalue stored `eigval(band, kpt)`.
tutorial03 head: `1 1 -5.82184795595698` … total 768 = 12·64. Values in eV.

---

## 3. `seedname.amn` (INPUT, formatted) — projections A_{mn}(k)

Written by DFT interface; W90 reads in `overlap_read` (`overlap.F90:336-378`).
```
line 1: comment/date string (free text)
line 2: num_bands  num_kpts  num_proj      (read list-directed; example '(3i5)')
then num_bands*num_proj*num_kpts data lines:
        m  n  nkp   Re(A)  Im(A)
```
`m` = band index (1..num_bands), `n` = projection index (1..num_proj), `nkp` = k-index.
Placement: `au_matrix(m, proj2wann_map(n), nkp) = cmplx(Re,Im)` (`overlap.F90:376`).
**Index-placed / order-independent:** the explicit `(m,n,nkp)` on each line determines
placement; file order is conventional only. Conventional loop order (as emitted / in
tutorial03) is band `m` fastest, then projection `n`, then k. `proj2wann_map(n)<0` lines
are skipped (selective projections). Example header `12 64 8`; data lines 6144 = 12·8·64.
Sample:
```
Created on 25Feb2006 at 11:03:51
   12   64    8
    1    1    1    -315.325490567895 -503.972348172530
```
NOTE: W90's own `overlap_write` produces `seedname.amn_dump` (header literal `"header"`,
counts `'(3i5)'`, rows `'(3i5,2f18.12)' m ip ik au`) — a debug dump, **not** the standard
`.amn`. The de-facto standard `.amn` column layout comes from pw2wannier90 / the example
files, since W90 itself only reads it (list-directed, so exact spacing is irrelevant).

---

## 4. `seedname.mmn` (INPUT, formatted) — overlaps M_{mn}^{(k,b)}

Written by DFT interface; read in `overlap_read` (`overlap.F90:246-320`).
```
line 1: comment/date string
line 2: num_bands  num_kpts  nntot
then for each of (num_kpts*nntot) blocks:
   header line:  nkp  nkp2  g1 g2 g3        (nkp2 = neighbour k index; g = nncell shift)
   then num_bands*num_bands lines:  Re(M)  Im(M)
```
Inner data loop nesting: **n outer, m inner** → `mmn_tmp(m,n)` with `m` fastest
(`overlap.F90:282-287`). Stored `m_matrix_local(:,:,nn,kpt) = mmn_tmp`.
**Block-matched (not positional):** each block's `(nkp2,g1,g2,g3)` is searched against
`kmesh_info%nnlist(nkp,:)` / `nncell(:,nkp,:)` to find the internal neighbour index `nn`
(`overlap.F90:291-305`); an error is raised if 0 or >1 matches. So the `.mmn` block order
need not match W90's internal b-ordering — but the (nkp2,g) tuple must exactly match a
computed neighbour. Example header `12 64 3` (nntot=... per silicon example), first block
`1 2 0 0 0`. Column values are `f18.12`-ish reals from the DFT code (read list-directed).

---

## 5. `seedname.chk` / `seedname.chk.fmt` (INPUT/OUTPUT, binary / formatted)

Record layout (identical logical sequence for both; `w90chk2chk.F90:234-358` unformatted,
`642-732` formatted writer). **Each bullet = one Fortran record in the binary file.**
Formatted-writer format specifiers shown to fix precision.

| # | Content | binary read | `.chk.fmt` format |
|---|---------|-------------|-------------------|
| 1 | `header` (33 chars, date/time) | `read()` | `'(A33)'` |
| 2 | `num_bands` (int) | `read()` | `'(I0)'` |
| 3 | `num_exclude_bands` (int, ≥0) | `read()` | `'(I0)'` |
| 4 | `exclude_bands(1..num_exclude_bands)` — ONE binary record; in `.fmt` one int/line | `read()` | `'(I0)'` per line |
| 5 | `real_lattice(3,3)` = `((rl(i,j),i=1,3),j=1,3)`, 9 reals, one record | `read()` | `'(9G25.17)'` |
| 6 | `recip_lattice(3,3)` same ordering, 9 reals | `read()` | `'(9G25.17)'` |
| 7 | `num_kpts` (int) | `read()` | `'(I0)'` |
| 8 | `mp_grid(3)` | `read()` | `'(I0," ",I0," ",I0)'` |
| 9 | `kpt_latt(3,num_kpts)` — ONE binary record; in `.fmt` one k/line 3 reals | `read()` | `'(3G25.17)'` per k |
| 10 | `nntot` (int) | `read()` | `'(I0)'` |
| 11 | `num_wann` (int) | `read()` | `'(I0)'` |
| 12 | `checkpoint` (20-char tag, e.g. `postwann`/`postdis`) | `read()` | `'(A20)'` |
| 13 | `have_disentangled` (logical; in `.fmt` written as `'(I1)'` 0/1) | `read()` | `'(I1)'` |

If `have_disentangled` (binary records 14-17; formatted equivalents):
| 14 | `omega_invariant` (real) | `'(G25.17)'` |
| 15 | `lwindow(num_bands,num_kpts)` — ONE binary record of logicals; `.fmt` one 0/1 per line, nested `do nkp; do i=1,num_bands` | `'(I1)'` |
| 16 | `ndimwin(num_kpts)` — ONE binary record; `.fmt` one int/line | `'(I0)'` |
| 17 | `u_matrix_opt(num_bands,num_wann,num_kpts)` — ONE binary record; `.fmt` one complex `(Re,Im)` per line, nesting `nkp{ j{ i(num_bands) }}` | `'(2G25.17)'` |

Always (whether or not disentangled):
| 18 | `u_matrix(num_wann,num_wann,num_kpts)` — ONE binary record; `.fmt` nest `k{ j{ i(num_wann) }}` complex/line | `'(2G25.17)'` |
| 19 | `m_matrix(num_wann,num_wann,nntot,num_kpts)` — ONE binary record; `.fmt` nest `l(kpt){ k(nntot){ j{ i(num_wann) }}}` | `'(2G25.17)'` |
| 20 | `wannier_centres(3,num_wann)` — ONE binary record; `.fmt` 3 reals/line per WF | `'(3G25.17)'` |
| 21 | `wannier_spreads(num_wann)` — ONE binary record; `.fmt` one real/line | `'(G25.17)'` |

Notes:
- The stored `m_matrix` is the **rotated, `num_wann`-sized** overlap (allocated
  `(num_wann,num_wann,nntot,num_kpts)`, `w90chk2chk.F90:338`), NOT the raw `num_bands`
  `.mmn`. Do not expect `num_bands`.
- Fortran array element order inside a multi-item record is column-major with the written
  loop nesting; e.g. record 5 iterates `i` (rows) fastest.
- Both `real_lattice` and `recip_lattice` are in Ångström / Å⁻¹, in the post-transpose
  orientation of gotcha 6.

---

## 6. `seedname.nnkp` (OUTPUT, formatted) — written by `kmesh_write` (`kmesh.F90:958-1128`)

Written when `postproc_setup=.true.` (or `-pp`). Structure:
```
# File written on <date> at <time>            ('(4(a),/)')

calc_only_A  :  F                             ('(a,l2,/)')

begin real_lattice
  <a1x a1y a1z>                               each '(3(f12.7))'
  <a2 ...>
  <a3 ...>
end real_lattice

begin recip_lattice
  <b1 ...>   '(3f12.7)'
  <b2 ...>
  <b3 ...>
end recip_lattice

begin kpoints
  <num_kpts>                                  '(i6)'
  <k1 k2 k3>  (num_kpts lines)                '(3f14.8)' fractional
end kpoints

begin projections      (or begin spinor_projections if spinors)
  <num_proj>                                  '(i6)'
  for each proj:
     site(1) site(2) site(3)  l  mr  radial   '(3(f10.5,1x),2x,3i3)'
     z(1) z(2) z(3)  x(1) x(2) x(3)  zona      '(2x,3f11.7,1x,3f11.7,1x,f7.2)'
  (spinor variant adds a 3rd line: spin(1i3) spin_quant_axis(3f11.7))
end projections

[ begin auto_projections     (only if auto_projections requested)
    <num_proj>  '(i6)'
    0           '(i6)'      (reserved)
  end auto_projections ]

begin nnkpts
  <nntot>                                     '(i4)'
  for nkp=1..num_kpts, for nn=1..nntot:
     nkp  nnlist(nkp,nn)  nncell(1) nncell(2) nncell(3)   '(2i8,3x,3i4)'
end nnkpts

begin exclude_bands
  <num_exclude_bands>                         '(i4)'
  <band>  (one per line if >0)                '(i4)'
end exclude_bands
```
**nnkpts / kpb block is the load-bearing part:** for each k-point `nkp` and each of its
`nntot` neighbours, it gives the neighbour k-index `nnlist(nkp,nn)` and the reciprocal
lattice-vector shift `nncell(1:3)` (the "G-vector" placing the neighbour in the periodic
image). The DFT interface uses these tuples to produce matching `.mmn` blocks (§4). No
`ndegen`/weights here — b-vector weights `wb` are computed internally by W90, not stored
in `.nnkp` (they appear in `.wout` and optionally `.bvec`).

---

## 7. `seedname_hr.dat` (OUTPUT, formatted) — `hamiltonian_write_hr` (`hamiltonian.F90:631-692`)

```
line 1: <date/time header>          write(*,*)  (list-directed, leading space)
line 2: num_wann                    write(*,*)
line 3: nrpts                       write(*,*)
degeneracy block: ndegen(1..nrpts), 15 per line, '(15I5)'
then nrpts*num_wann*num_wann rows, nested do irpt{ do i(num_wann){ do j(num_wann) }}:
   irvec(1) irvec(2) irvec(3)  j  i  Re(H)  Im(H)     '(5I5,2F12.6)'
```
- `j` (fast) = left/row WF index, `i` = right/column WF index; value `ham_r(j,i,irpt)`.
- `irvec` = R in lattice-vector units (integers). Units eV. **6-decimal, lossy** (gotcha 2).
- H(R) stored undivided by ndegen (gotcha 4). Doc example header shows `num_wann`,
  `nrpts`, then a `15I5` degeneracy block, then the matrix rows.

## 8. `seedname_r.dat` (OUTPUT, formatted) — `plot_write_rmn` (`plot.F90:2613-2711`)

```
line 1: <date/time header>   write(*,*)
line 2: num_wann             write(*,*)
line 3: nrpts                write(*,*)
then nrpts*num_wann*num_wann rows, nested do rpt{ do m(num_wann){ do n(num_wann) }}:
   irvec(1) irvec(2) irvec(3)  n  m  <6 reals>        '(5I5,6F12.6)'
```
- 6 reals = `Re(r_x) Im(r_x) Re(r_y) Im(r_y) Re(r_z) Im(r_z)` for ⟨m0|r|nR⟩.
- Computed with Im-ln branch for diagonal (m==n) term, off-diagonal linear term
  (gotcha 5, `plot.F90:2686-2699`). `n` fast (row), `m` (col); value from `m_matrix(n,m,...)`.
- Units: Å (position). 6-decimal.

## 9. `seedname_tb.dat` (OUTPUT, formatted) — `hamiltonian_write_tb` (`hamiltonian.F90:862-994`)

Combined H + r + lattice, full-ish precision `E15.8`.
```
line 1: <date/time header>                 write(*,*)
line 2: a_1 (3 reals, Å)                    write(*,*) real_lattice(1,:)
line 3: a_2                                  write(*,*)
line 4: a_3                                  write(*,*)
line 5: num_wann                             write(*,*)
line 6: nrpts                                write(*,*)
degeneracy block: ndegen 15/line            '(15I5)'
H part: for irpt=1..nrpts:
   (blank line) irvec(1) irvec(2) irvec(3)  '(/,3I5)'
   for i(num_wann){ for j(num_wann) }:
      j  i  Re(H) Im(H)                      '(2I5,3x,2(E15.8,1x))'   ham_r(j,i,irpt)
r part: for irpt=1..nrpts:
   (blank line) irvec(1..3)                  '(/,3I5)'
   for i{ for j }:
      j  i  <6 reals rx,ry,rz Re/Im>         '(2I5,3x,6(E15.8,1x))'   pos_r(:)
```
Note the leading `/` in `'(/,3I5)'` writes a **blank line before each R block**. Lattice
vectors are `real_lattice(k,:)` = a_k (rows are a-vectors, gotcha 6), Å.

## 10. `seedname_wsvec.dat` (OUTPUT) — `ws_write_vec` (`ws_distance.F90:306`)

Written if `write_hr|write_rmn|write_tb`. Line 1: date/time + `use_ws_distance` value.
Per (R, m, n): a header `irvec(1..3)  m  n`, then the count of translation vectors T and
the T list that fold the WF back into the WS cell. Only relevant when
`use_ws_distance=.true.`; downstream H(R) matrix elements are shifted by these T.

## 11. `seedname.bvec` (OUTPUT) — if `write_bvec`
Line 1 date/time; line 2 `num_kpts  nntot`; then for each k and each neighbour the
b-vector Cartesian `(bx,by,bz)` and weight `wb`. (Weights otherwise only in `.wout`.)

---

## 12. Band-structure output (all if `bands_plot=.true.`)

### `seedname_band.dat` (`plot.F90:1157-1172`)
For each band `i`, for each path k-point `nkp`:
```
   xval(nkp)  eig_int(i,nkp) [ bands_proj(i,nkp) ]
```
Format `'(2E16.8)'` (or `'(3E16.8)'` with projection column if `write_proj`). A **blank
line** (`write(*,*) ' '`) separates bands. `xval` = cumulative path length (Å⁻¹ or the
recip-metric arc-length); `eig_int` in eV.

### `seedname_band.kpt` (`plot.F90:712-717`)
```
line 1: total_pts                    write(*,*)
then total_pts lines: k1 k2 k3  1.0   '(3f12.6,3x,a)'   (fractional coords + weight "1.0")
```
Directly reusable as a DFT band-path input.

### `seedname_band.gnu` (`plot.F90:1157-1217`)
Gnuplot script. Sets `xrange [0:xval(total_pts)]`, `yrange [emin:emax]`, xtics at the
special-point `xval`s with labels (`glabel`, `|`-joined at discontinuities), vertical
lines at each special point (format label 705), final line
`plot "seedname_band.dat"`. (Exact `70x` format numbers are cosmetic; structural
description suffices for reproduction.)

### `seedname_band.labelinfo.dat` (`plot.F90:721-741`) — companion to the above
One line per special point: `label  index  xval  k1 k2 k3` with
`'(a,3x,I10,3x,4f18.10)'`. `index` = position in the concatenated path, `xval` = its
x-coordinate, then the fractional k-coords. Load-bearing for reproducing gnu/agr tics.

---

## 13. `seedname.wout` (OUTPUT, formatted) — exact load-bearing lines

Verbosity via `iprint`. Length unit token `Ang`/`Bohr`; all lengths scaled by
`lenconfac`, spreads by `lenconfac**2`.

### System echo (`wannier90_readwrite.F90:1900-1994`)
```
 |   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |
 | k-point      Fractional Coordinate        Cartesian Coordinate (Ang^-1)    |
```
(Bohr variants swap the unit string.)

### b-vector shells report (`kmesh.F90:201-217`)
```
 |                    Distance to Nearest-Neighbour Shells                    |
 |                    ------------------------------------                    |
 |          Shell             Distance (Ang^-1)          Multiplicity         |
 |          -----             -----------------          ------------         |
 |             1                   0.398833                      6            |
```
Each shell row: `'(1x,a,11x,i3,17x,f10.6,19x,i4,12x,a)'` → `| <shell> <dist> <multi> |`,
`dist = dnn(shell)/lenconfac`. Followed by (`kmesh.F90:1388`)
`| The b-vectors are chosen automatically ...` and the used-shells / nearest-neighbour
count table, ending with the completeness-relation confirmation
`Completeness relation is fully satisfied [Eq. (B1), PRB 56, 12847 (1997)]`.

### Initial / Final State spread block (`wannierise.F90:492-500, 850-857`; formats 1000/1001)
```
 ------------------------------------------------------------------------------
 Initial State
  WF centre and spread    1  (  0.000000,  1.969243,  1.969243 )     1.52435832
  ...
  Sum of centres and spreads ( 11.815458, 11.815458, 11.815458 )    12.62663472
```
- Per-WF line format 1000: `2x,'WF centre and spread',i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8`
  → columns: WF index, (x,y,z centre), spread (⟨r²⟩−⟨r⟩²). Centre × lenconfac, spread × lenconfac².
- Sum line format 1001: `2x,'Sum of centres and spreads',1x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8`.

### Per-iteration CONV / SPRD trace (`wannierise.F90:522-529` and 700-765)
CONV line: `'(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)'` →
`<iter>  <ΔΩ>  <RMS grad>  <Ω_total>  <time>  <-- CONV`.
SPRD line (non-selective): `'(8x,a,F15.7,a,F15.7,a,F15.7,a)'` →
`        O_D=  <Ω_D>  O_OD=  <Ω_OD>  O_TOT=  <Ω_total>  <-- SPRD`.
DLTA line (from cycle ≥1): `Delta: O_D= ... O_OD= ... O_TOT= ... <-- DLTA`
(differences vs previous cycle). All values × lenconfac². `ΔΩ = Ω_tot(i) − Ω_tot(i−1)`;
`RMS grad = sqrt(|gcnorm1|)*lenconfac`.

### Final decomposition (`wannierise.F90:882-891`, non-selective branch)
```
         Spreads (Ang^2)       Omega I      =    12.480596753
        ================       Omega D      =     0.000000000
                               Omega OD     =     0.145856689
    Final Spread (Ang^2)       Omega Total  =    12.626453441
```
Formats: `'(3x,a21,a,f15.9)'` (Omega I and Total lines), `'(3x,a,f15.9)'` (D and OD).
Values: `om_i, om_d, om_od, om_tot` each × lenconfac². Identity
`Ω_tot = Ω_I + Ω_D + Ω_OD`. (Selective-localization mode replaces `Omega I` with
`Omega IOD`/`Omega Rest`; see `wannierise.F90:859-890`.)

### Timings (`docs .../files.md:575`)
`| <tag> : <Ncalls> <Time(s)>|`, then `All done: wannier90 exiting`.

---

## 14. `seedname.dmn` / UNK files (brief)
- `.dmn` (INPUT if `site_symmetry`): symmetry-adapted D-matrices, written by DFT interface.
- `UNKp.s` (INPUT if `wannier_plot`): unformatted-or-formatted (`wvfn_formatted`) periodic
  Bloch functions on a real-space grid. Record 1: `ngx ngy ngz ik nbnd`; then per band a
  record of `ngx*ngy*ngz` values (spinor case: 2 records per band). Filename
  `UNK<ik:5.5>.<spin:1>` or `UNK<ik:5.5>.NC` for spinors.

---

## 15. Reproduction checklist (Julia port)
- Read `.chk.fmt` (not binary) for U/U_dis/M/centres/spreads at full G25.17 precision.
- For H(R)/r(R): compute internally (Fourier of U-rotated M) rather than round-tripping
  the 6-dp `hr.dat`; use `hr.dat`/`tb.dat` only for cross-checking to their stated precision.
- Match `.eig` band-fast/k-slow ordering and the abort-on-mismatch invariant.
- Place `.amn` by explicit (m,n,k); match `.mmn` blocks by (nkp2, g) against nnkpts.
- Divide H(R) by `ndegen(R)`; verify Σ 1/ndegen == ∏mp_grid.
- Use atan2 branch for Im-ln; `fac = exp(−i·2π k·R)/num_kpts`.
- Keep `real_lattice(i,:) = a_i` (post-transpose); lengths Å, energies eV internally.
