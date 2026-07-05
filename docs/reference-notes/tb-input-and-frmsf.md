# TB-model input (`_hr.dat` / `_tb.dat` read-back) and FermiSurfer `.frmsf` output

Implementation-grade spec for two additions to WannierFunctions.jl:
- **A)** Read a `_hr.dat` or `_tb.dat` file back into an interpolable model (like
  WannierBerri's `System_tb`), so a user can interpolate directly from a tight-binding
  file with **no** `.chk`/`.mmn`/`.eig`.
- **B)** Write a FermiSurfer `.frmsf` 3D-tabulation file.

Sources:
- `reference/wannier90/src/hamiltonian.F90` — `hamiltonian_write_hr` (631-692),
  `hamiltonian_write_tb` (862-994). The reader must be the **exact inverse** of these.
- `reference/wannier90/src/plot.F90` — `plot_fermi_surface` (1360-1544): the k→bands
  interpolation engine (the physics for `.frmsf` is identical; only grid convention +
  output layout differ). NB wannier90's native FS file is **`.bxsf`** (XCrySDen), *not*
  `.frmsf` — `.frmsf` is not in the wannier90 tree.
- FermiSurfer official docs, [Input file page](https://mitsuaki1987.github.io/fermisurfer/en/_build/html/input.html)
  and [paper (arXiv:1811.06177)](https://arxiv.org/pdf/1811.06177), for the `.frmsf` layout.
- WannierBerri [`result/__tabresult.py`](https://docs.wannier-berri.org/en/latest/_modules/wannierberri/result/__tabresult.html)
  — cross-check of the `.frmsf` writer (shift flag, `nk` not `nk+1`, band-outer ordering,
  `E − E_F`).
- Existing project code this dovetails with: `src/output.jl` (`read_hr`/`write_hr`/
  `write_tb`), `src/operator.jl` (`TBOperator`), `src/berry.jl` (`BerryModel`),
  `src/interpolate.jl`. Companion note: `docs/reference-notes/file-formats.md` §7/§9/§12.

Constants: reference build defaults to CODATA2006
(`bohr_angstrom_internal = 0.52917720859`, `constants.F90:182`). Not load-bearing here —
H(R) is eV, r(R) is Å, reciprocal vectors are Å⁻¹, all pass through unchanged.

--------------------------------------------------------------------------------
## 0. Index / array conventions (critical for the reader)

Both writers loop **`i` (right/column WF) outer, `j` (left/row WF) inner** and print the
pair as `j i value`, writing `ham_r(j,i,irpt)` (`hamiltonian.F90:684, 945`). So on read:

```
row token = j   (left / "0j" WF, the bra)
col token = i   (right / "Ri" WF, the ket)
value     → O[j, i, irpt]          ! j is the FAST loop, matches src/output.jl read_hr
```

The physical element is ⟨**0**,j | O | **R**,i⟩ = "j in home cell, i at R" (the writer's
own comment, `hamiltonian.F90:448-449`). Storage is `O[j,i,ir]` = `(nw, nw, nrpts)`.

- **Fortran 1-based, column-major.** All indices below are 1-based.
- H(R) and r(R) are stored **UNDIVIDED by `ndegen`** (the wannier90 convention). The
  *consumer* (k-evaluation) divides: `O(k) = Σ_R e^{+i2π k·R} O(R)/ndegen(R)`
  (`operator.jl:9`, `TBOperator` call operator). Do **not** pre-divide on read.
- `irvec(:,irpt)` = R in **lattice-vector (integer) units**; `nrpts` = number of
  Wigner–Seitz vectors; `ndegen(irpt)` = WS degeneracy weight.
- WS sum rule (invariant a reader can assert): `Σ_R 1/ndegen(R) == ∏ mp_grid`
  (checked to `eps8`, `hamiltonian.F90:850`). But **`mp_grid` is NOT stored in either
  file** — see Traps.

--------------------------------------------------------------------------------
## 1. `_hr.dat` format — exact read-back (`hamiltonian_write_hr`, 631-692)

Writer (verbatim, `hamiltonian.F90:677-688`):
```fortran
write (file_unit, *) header                 ! list-directed: leading space + 33-char date string
write (file_unit, *) num_wann               ! list-directed int
write (file_unit, *) nrpts                  ! list-directed int
write (file_unit, '(15I5)') (ndegen(i), i=1, nrpts)
do irpt = 1, nrpts
  do i = 1, num_wann                         ! i = column (ket), OUTER
    do j = 1, num_wann                       ! j = row (bra), INNER (fast)
      write (file_unit, '(5I5,2F12.6)') irvec(:, irpt), j, i, ham_r(j, i, irpt)
```

On-disk layout:
```
line 1:  <date/time header>                 (free text, skip)
line 2:  num_wann                            (single int, tokenize)
line 3:  nrpts                               (single int, tokenize)
ndegen block: nrpts integers, 15 per line, format (15I5) → each width-5
then nrpts*num_wann*num_wann rows:
   R1 R2 R3  j  i  Re(H)  Im(H)              format (5I5,2F12.6)
```

Read algorithm (inverse; matches `src/output.jl:99-126` `read_hr`, already implemented):
1. Skip line 1.
2. `num_wann = parse(Int, line2)`, `nrpts = parse(Int, line3)`.
3. Accumulate `nrpts` integers for `ndegen` across wrapped `(15I5)` lines
   (**tokenize, never column-slice** — leading spaces; last line is partial).
4. For each of `nrpts*num_wann*num_wann` rows: split into 7 tokens
   `R1 R2 R3 j i re im`; set `Hr[j,i,irpt] = re + im*im_unit`. On the FIRST row of each
   irpt block (`i==1 && j==1`) record `irvec[irpt] = (R1,R2,R3)`.
5. Row-block ordering: `irpt` outer, then `i` (col), then `j` (row) fastest. The reader
   can rely on this ordering OR key off the explicit `(R1,R2,R3,j,i)` on each line
   (order-independent placement is safer — the file *is* self-describing per row).

**Precision: `F12.6` — LOSSY (6 decimals).** A model read from `_hr.dat` cannot reproduce
the reference to full precision. Use `_tb.dat` (`E15.8`) or `.chk.fmt` (`G25.17`) when
precision matters. (file-formats.md gotcha 2.)

Units: `Re/Im(H)` in **eV**.

--------------------------------------------------------------------------------
## 2. `_tb.dat` format — exact read-back (`hamiltonian_write_tb`, 862-994)

Combined lattice + H(R) + r(R), precision `E15.8`. Writer (`hamiltonian.F90:927-988`):
```fortran
write (file_unit, *) header                          ! date string
write (file_unit, *) real_lattice(1, :)              ! a_1 (Å), list-directed
write (file_unit, *) real_lattice(2, :)              ! a_2
write (file_unit, *) real_lattice(3, :)              ! a_3
write (file_unit, *) num_wann
write (file_unit, *) nrpts
write (file_unit, '(15I5)') (ndegen(i), i=1, nrpts)
! <0j|H|Ri> block
do irpt = 1, nrpts
  write (file_unit, '(/,3I5)') irvec(:, irpt)         ! LEADING BLANK LINE then R
  do i = 1, num_wann                                   ! col, outer
    do j = 1, num_wann                                 ! row, fast
      write (file_unit, '(2I5,3x,2(E15.8,1x))') j, i, ham_r(j, i, irpt)
! <0j|r|Ri> block  (three Cartesian components, complex)
do irpt = 1, nrpts
  write (file_unit, '(/,3I5)') irvec(:, irpt)         ! LEADING BLANK LINE then R
  do i = 1, num_wann
    do j = 1, num_wann
      write (file_unit, '(2I5,3x,6(E15.8,1x))') j, i, pos_r(:)   ! rx,ry,rz Re/Im interleaved
```

On-disk layout:
```
line 1:  <date/time header>                 (skip)
line 2:  a_1x a_1y a_1z                      (3 reals, Å) — real_lattice(1,:)
line 3:  a_2x a_2y a_2z
line 4:  a_3x a_3y a_3z
line 5:  num_wann
line 6:  nrpts
ndegen block: nrpts ints, (15I5)
H(R) block:  for each irpt:
   <blank line>
   R1 R2 R3                                  (3I5)
   then num_wann*num_wann rows (i outer, j fast):
      j  i  Re(H) Im(H)                       (2I5,3x,2(E15.8,1x))
r(R) block:  for each irpt:
   <blank line>
   R1 R2 R3                                  (3I5)
   then num_wann*num_wann rows (i outer, j fast):
      j  i  Re(rx) Im(rx) Re(ry) Im(ry) Re(rz) Im(rz)   (2I5,3x,6(E15.8,1x))
```

Read algorithm (inverse of `write_tb`, `src/output.jl:148-206`):
1. Skip line 1.
2. Read 3 lattice-vector lines. **`real_lattice(k,:) = a_k`** — rows are the a-vectors
   (post-transpose orientation, file-formats.md gotcha 6). In the project `Lattice`,
   `A` has **columns** = a-vectors, so build `A[:,k] = row_k` i.e.
   `A = permutedims(hcat(a1,a2,a3))'`… concretely: parse rows `a1,a2,a3`; set
   `A = hcat(a1, a2, a3)` (each parsed row becomes a column). Units Å.
3. `num_wann`, `nrpts` (tokenize).
4. `ndegen` block `(15I5)`, `nrpts` ints across wrapped lines.
5. H block: for each irpt, skip the blank line, read `R1 R2 R3` (record `irvec[irpt]`),
   then `num_wann*num_wann` rows `j i re im` → `Hr[j,i,irpt] = re + im*im_unit`.
6. r block: for each irpt, skip blank line, read `R1 R2 R3` (assert it matches
   `irvec[irpt]` — same order), then rows `j i re_x im_x re_y im_y re_z im_z` →
   `Ar[j,i,irpt,c] = re_c + im_c*im_unit` for c=1,2,3.

Precision: `E15.8` (8 sig figs). Units: H(R) **eV**; r(R) **Å**; lattice **Å**.

**List-directed lattice lines (byte detail, only matters for a byte-exact *writer*, not
the reader):** `write(*,*) real_lattice(k,:)` right-justifies each value in fields of
width 21/26/26 with a trailing 5 spaces (`src/output.jl:160-167`). The **reader must
tokenize**, so this is irrelevant on read.

--------------------------------------------------------------------------------
## 3. Building a `BerryModel`-equivalent from a TB file (no chk/mmn/eig)

The project already has the target struct: `BerryModel` (`src/berry.jl:34-44`) carries
`lattice, irvec, ndegen, Rcart, Hr(nw,nw,nr), Ar(nw,nw,nr,3), wsdist`. **A TB file supplies
exactly these** — no wannierisation, no overlaps needed. `TBOperator` (`src/operator.jl`)
is the same object for a single quantity (H → 1 component, r → 3 components).

### 3a. Constructor mapping (what to fill)

| BerryModel field | from `_tb.dat` | from `_hr.dat` |
|---|---|---|
| `irvec`, `ndegen` | read directly | read directly |
| `lattice` | from the 3 a-vector lines | **must be supplied separately** (not in file) |
| `Rcart` | `lattice.A * SVector(R...)` per irvec (matches `berry.jl:88`) | same |
| `Hr[j,i,ir]` | H block, undivided | matrix rows, undivided |
| `Ar[j,i,ir,c]` | **r block, read verbatim** (see 3b) | **unavailable → `Ar` empty (0 components)** |
| `wsdist` | `nothing` (a TB file already encodes the WS-folded H; `use_ws_distance` folding is baked in by whoever wrote it) | `nothing` |

### 3b. The r(R) block IS A(R) — read it verbatim, NO finite difference, NO Hermitisation

For a pure TB model there is **no `.mmn`**, so the postw90 finite-difference connection
`A_α(q) = Σ_b i w_b b_α M̃(q,b)` (`berry.jl:100-121`) is simply *unavailable*. The `_tb.dat`
r block is wannier90's own position operator ⟨0j|r|Ri⟩ (Wang–Yates–Souza–Vanderbilt
PRB 74 195118 Eq. 44, `hamiltonian.F90:951-988`). Reading it straight into `Ar` is exactly
what WannierBerri's `System_tb` does. So:
```
Ar[j,i,ir,c] = pos_r_c(j,i,ir)     ! verbatim, undivided by ndegen; the evaluator divides
```
Do **not** re-Hermitise and do **not** finite-difference. (The `TBOperator(:position,...)`
built here evaluates to A(k′) = Σ_R e^{+i2πk′·R} r(R)/ndegen(R), the position-operator
matrix in reciprocal space, i.e. the Berry connection.)

### 3c. Diagonal-convention caveat (a trap, NOT something to reconcile)

The `_tb.dat` r(R) diagonal uses the **Im-ln (log) convention**
(`pos_r = −Σ_b w_b b · aimag(log(M_nn)) · fac`, `hamiltonian.F90:973`), whereas a
`.chk`-sourced `BerryModel` uses **linear-for-all-elements + Hermitise**
(`berry.jl:60-62`, `transl_inv=F`). These are two valid discretisations of the *same*
physical A. A `_tb.dat`-sourced model therefore carries the Im-ln diagonal and its raw
`Ar[i,i,·,·]` differs on the diagonal from a chk-sourced one. **Gauge-covariant outputs
(Berry curvature, AHC, orbital magnetisation) still agree to interpolation precision** —
so validate those, never raw `Ar` elements. (Same gauge-invariance rule as the rest of the
project: compare bands/Ω/centres/interpolated observables, never raw operator elements.)

### 3d. Berry-quantity subset: H(R) only vs H(R)+r(R)

Directly from the project's own gating (`berry.jl:189-191`; AHC/Kubo/morb error on empty
`Ar`):

**Available with H(R) ONLY** (`_hr.dat`, or `_tb.dat` ignoring the r block):
- Interpolated band energies & eigenvectors (`interpolate.jl`).
- Band velocities / group velocities (dH/dk).
- Density of states (DOS), `geninterp` generic k-tabulation.
- BoltzWann transport (TDF/σ/Seebeck — needs only ε and v).
- **Fermi surface / `.frmsf`** (Part B — needs only H(R)).

**Requires H(R) + r(R)** (`_tb.dat` with r block, or a `.chk`+`.mmn` build):
- Anomalous Hall conductivity (AHC) and Berry curvature.
- Kubo optical conductivity, orbital magnetisation (morb), gyrotropic tensors.
- Shift current, spin Hall conductivity (SHC — also needs the spin operator S(R), a
  further ingredient beyond r(R)).

(Rule of thumb: anything whose integrand contains the Berry connection A or its curl needs
r(R); anything that is a functional of ε(k) and ∂ε/∂k alone needs only H(R).)

--------------------------------------------------------------------------------
## 4. FermiSurfer `.frmsf` output

**wannier90 itself writes `.bxsf` (XCrySDen), not `.frmsf`.** The physics engine, however,
is identical to `plot_fermi_surface` (`plot.F90:1360-1544`): interpolate H(R)→H(k) on a
regular grid, diagonalise, tabulate eigenvalues. Only the **grid convention** and the
**file layout** differ. The `.frmsf` layout below is from the FermiSurfer docs, cross-checked
against WannierBerri's writer.

### 4a. Exact `.frmsf` file layout

```
line 1:  nk1 nk2 nk3          ! grid points per reciprocal direction (integers)
line 2:  1                    ! shift flag (see 4c)
line 3:  nbnd                 ! number of bands
line 4:  b1x b1y b1z          ! reciprocal lattice vector 1 (Å⁻¹)
line 5:  b2x b2y b2z          ! reciprocal lattice vector 2
line 6:  b3x b3y b3z          ! reciprocal lattice vector 3
then nbnd*nk1*nk2*nk3 energy values, ONE per line
[optional] then nbnd*nk1*nk2*nk3 colour/quantity values, ONE per line, SAME ordering
```

### 4b. Value ordering (the #1 correctness detail)

The energy array is `eig(nk3, nk2, nk1, nbnd)` — Fortran-flattened, i.e.:
```
do ibnd = 1, nbnd           ! band OUTERMOST
  do ik1 = 1, nk1
    do ik2 = 1, nk2
      do ik3 = 1, nk3       ! ik3 INNERMOST / FASTEST
        write value
```
The optional colour block follows in **exactly the same order**. WannierBerri writes both
via `flatten(order='F')` over `Enk.data[:, iband]` (band-outer, k-inner) — confirming this.

In Julia (column-major), if you store `E[ik3, ik2, ik1, ibnd]` then the natural
`vec(E)` (or `E[:]`) already yields this Fortran order — write `vec(E[:,:,:,b])` per band,
or `vec(E)` for the whole array. If you store `E[ibnd, ik1, ik2, ik3]` you must permute:
`vec(permutedims(E, (4,3,2,1)))`.

### 4c. The four load-bearing details (each a silent-bug axis)

1. **Grid = `nk` points, NO boundary doubling.** frmsf uses exactly `nk1 nk2 nk3` points
   and does *not* repeat the BZ-boundary point. This is the #1 divergence from `.bxsf`,
   which writes `num_points+1` points and loops `(num_points+1)^3`
   (`plot.F90:1458,1468-1470,1520`). **Do not copy the bxsf `+1` loop.** Sample the grid
   `k = (i1-1)/nk1, (i2-1)/nk2, (i3-1)/nk3` for `i∈1..nk` (fractional, k-type 0). WannierBerri
   writes `grid[0] grid[1] grid[2]` directly (no `+1`).

2. **Shift flag line = `1`.** WannierBerri hardcodes `"1 \n"`. The flag selects the k-grid
   type / fractional-coordinate convention (docs allow 0/1/2). Use `1` to match WannierBerri
   output for a Γ-inclusive Monkhorst grid; `0` is the unshifted convention. For a
   Wannier-interpolation grid starting at Γ, `1` is the safe default.

3. **Energies are absolute; FermiSurfer assumes E_F = 0 by default.** FermiSurfer draws the
   isosurface at 0. WannierBerri therefore writes **`E − E_F`** so the surface lands at the
   Fermi level with no menu shift. **Recommended: write `ε_int(k) − E_F`** (accept a
   `fermi_energy` argument; default 0.0 → absolute energies, then the user shifts in the
   FermiSurfer "Shift Fermi Energy" menu). Contrast `.bxsf`, which writes absolute energies
   *and* a `Fermi Energy:` header field (`plot.F90:1513`); `.frmsf` has **no Fermi-energy
   field**, so the shift must be baked into the values.

4. **Reciprocal vectors, one per line, `14.8f` per component (Å⁻¹).** WannierBerri:
   `"  ".join("{:14.8f}".format(x) for x in v)` per row. Rows are b1,b2,b3.
   Project convention: `lattice.B` has **columns** = b_i (`chk.jl:15`), so write
   `B[:,1], B[:,2], B[:,3]` as the three rows (each b_i's 3 Cartesian components).

### 4d. The interpolation engine (reuse `plot_fermi_surface` / `interpolate.jl`)

Per grid k-point `k' = ((i1-1)/nk1, (i2-1)/nk2, (i3-1)/nk3)`:
```
H(k') = Σ_R e^{+i 2π k'·R} H(R) / ndegen(R)          ! op(k') in operator.jl / eval_hk
ε(k') = eigvals(Hermitian(H(k')))                     ! ascending; nbnd = num_wann
```
This is byte-for-byte the bxsf loop (`plot.F90:1473-1488`: `fac = exp(+i·rdotk)/ndegen`,
`rdotk = 2π[(i1-1)R1+(i2-1)R2+(i3-1)R3]/num_points`, pack upper triangle, `ZHPEVX('N',...)`),
except the grid has `nk` (not `nk+1`) points and the k-spacing divides by `nk` not
`num_points`. `nbnd = num_wann`. Sort eigenvalues ascending (ZHPEVX does; `eigvals` does).

For a **colour** overlay (e.g. Berry curvature magnitude, band velocity, orbital character),
compute the per-band scalar at each k in the same loop and emit the second block in the same
(ibnd, ik1, ik2, ik3) order. Berry-curvature colour needs H(R)+r(R) (a `BerryModel` with
non-empty `Ar`); a plain band-velocity or spin colour needs only H(R) (+ S(R) for spin).

--------------------------------------------------------------------------------
## 5. Traps / silent-mismatch checklist

**TB-file reader:**
1. **`j i` order, `j` fast.** Row token = row/bra WF `j`; store `O[j,i,ir]`. Both writers
   loop `i` (col) outer, `j` (row) inner. Swapping them transposes every H(R)/r(R) block —
   silently wrong bands unless H is Hermitian-symmetric in a way that hides it.
2. **Never divide by `ndegen` on read.** Store undivided; divide only in k-evaluation
   (`operator.jl:9,69`). Double-dividing under-weights degenerate R.
3. **Tokenize, don't column-slice.** Header, `num_wann`, `nrpts`, and (in `_tb.dat`) the 3
   lattice lines are list-directed (leading whitespace, no fixed columns). Only `(15I5)` and
   the matrix rows are fixed-format. The `(15I5)` degeneracy block wraps at 15 and the last
   line is partial — accumulate until you have `nrpts`.
4. **`_tb.dat` has a LEADING blank line before every R block** (`'(/,3I5)'`). Skip it in
   both the H and r sections. Missing this shifts the whole parse by one line.
5. **Lattice orientation:** `_tb.dat` rows are `real_lattice(k,:) = a_k` (post-transpose,
   file-formats.md gotcha 6). Build the project `Lattice.A` with **columns** = a-vectors:
   `A = hcat(a1, a2, a3)` where `a_k` is parsed row k. Get this wrong and every Rcart /
   reciprocal vector is transposed.
6. **`mp_grid` is NOT in either file.** You cannot recover it from `_hr.dat`/`_tb.dat`
   alone. If a downstream routine needs `mp_grid` (e.g. to rebuild a matching k-grid), it
   must be supplied by the caller or inferred from `Σ_R 1/ndegen(R) = ∏mp_grid` plus the
   R-vector extents — not uniquely determined in general, so require it as an argument.
7. **`_hr.dat` is 6-dp LOSSY** (`F12.6`). A model read from `_hr.dat` will not reproduce
   full-precision bands. Prefer `_tb.dat` (`E15.8`) or `.chk.fmt` (`G25.17`).
8. **`_hr.dat` has no r(R) and no lattice.** An `_hr.dat`-sourced model is **H(R)-only**
   (bands/DOS/velocities/BoltzWann/Fermi-surface) and needs the lattice supplied separately;
   Berry-connection outputs (AHC/Kubo/morb) must error, matching `berry.jl:190-191`.
9. **`_tb.dat` r(R) diagonal is Im-ln, not linear+Hermitise.** Read verbatim into `Ar`; do
   NOT finite-difference or re-Hermitise. Validate gauge-covariant outputs, never raw `Ar`
   (§3c).
10. **`wsdist = nothing` for a TB-file model.** A `_tb.dat`/`_hr.dat` already bakes in
    whatever `use_ws_distance` folding the producer used; do not re-apply WS-distance
    minimal-image folding on top.

**`.frmsf` writer:**
11. **`nk`, not `nk+1`.** No BZ-boundary doubling (unlike `.bxsf`). Grid line and value
    count both use `nk1*nk2*nk3` per band. Copying the bxsf `+1` loop is the #1 silent bug.
12. **Value order: band OUTER, ik3 INNER** (`eig(nk3,nk2,nk1,nbnd)`, Fortran flatten). The
    colour block follows in the identical order. Wrong nesting scrambles the surface.
13. **Write `ε − E_F`** (FermiSurfer isosurface is at 0; no Fermi field in the file). Default
    `fermi_energy=0` → absolute energies (user shifts in-app). `.bxsf` differs (absolute +
    header field).
14. **Shift flag `1`** on line 2 (matches WannierBerri). k-spacing divides by `nk`
    (fractional `k=(i-1)/nk`), not by `num_points`.
15. **Reciprocal vectors as rows** b1/b2/b3, `14.8f` components, Å⁻¹; project `lattice.B`
    columns are the b_i, so emit `B[:,i]` as row i.
