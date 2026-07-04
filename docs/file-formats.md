# File formats

A user-facing guide to the Wannier90-compatible files this package reads and writes. These are
the standard Wannier90 v3.1.0 formats; your DFT interface (e.g. `pw2wannier90`) produces the
inputs, and they are consumed here unchanged. For the exhaustive, implementation-grade spec with
reference source citations, see `docs/reference-notes/file-formats.md`.

Throughout: the **seedname** is the common prefix (e.g. `diamond`), so the files are
`diamond.win`, `diamond.amn`, and so on. Energies are always in **eV**; k-points are always
**fractional** (crystallographic); lengths default to **Ångström**.

Quick map:

| File | Direction | Contents |
|------|-----------|----------|
| `.win` | input | master input: cell, k-mesh, projections, parameters |
| `.eig` | input | DFT eigenvalues `ε_{mk}` (eV) |
| `.amn` | input | trial projections `A_{mn}(k)` |
| `.mmn` | input | neighbour overlaps `M_{mn}^{(k,b)}` |
| `_hr.dat` | output | real-space Hamiltonian `H(R)` (tight-binding model) |

All of these except the binary `.chk` are plain text and **list-directed / tokenised** — exact
column spacing is irrelevant; parse by tokens, never by fixed columns.

---

## `seedname.win` — master input (input)

Free-form, case-insensitive text. `!` and `#` start comments; `:` and `=` separate a key from
its value. Blocks are delimited by `begin <name>` … `end <name>`. Order is free.

Key scalar keywords:

- `num_wann` — number of Wannier functions (**mandatory**).
- `num_bands` — number of Bloch bands (defaults to `num_wann` ⇒ isolated case, no
  disentanglement).
- `mp_grid = n1 n2 n3` — Monkhorst–Pack subdivisions; `num_kpts = n1·n2·n3`.
- `num_iter` — number of MLWF minimisation iterations.
- Disentanglement (only when `num_bands > num_wann`): `dis_win_min/max`,
  `dis_froz_min/max` (eV), `dis_num_iter`.
- Interpolation / output: `bands_plot`, `bands_num_points`, `write_hr`, `use_ws_distance`.

Key blocks:

```
begin unit_cell_cart          # optional first line: 'ang' (default) or 'bohr'
  a1x a1y a1z                 # three rows: Cartesian components of a1, a2, a3
  a2x a2y a2z
  a3x a3y a3z
end unit_cell_cart

begin atoms_frac              # (or atoms_cart) — one per line
  C  -0.125 -0.125 -0.125
end atoms_frac

begin projections            # trial orbitals: site : angular-type
  f=0.0,0.0,0.0 : s          # site as f=frac / c=cart / element symbol
end projections

begin kpoints                # num_kpts rows of fractional k
  0.0 0.0 0.0
  ...
end kpoints

begin kpoint_path            # for band interpolation: two labelled points per line
  L 0.5 0.5 0.5  G 0.0 0.0 0.0
  G 0.0 0.0 0.0  X 0.5 0.0 0.5
end kpoint_path
```

Note: the cell may be given in `bohr`, but Wannier centres/spreads are always reported in
Å / Å² regardless.

---

## `seedname.eig` — eigenvalues (input)

One eigenvalue per line, with **band index fastest** (k outer, band inner):

```
  band_index   kpoint_index   eigenvalue_eV
  1   1   -5.82184795595698
  2   1    ...
```

Total lines = `num_bands · num_kpts`. The indices are positional and validated — a mismatch is
a hard error. Stored as `ε_{mk}`. Required for interpolation (`build_hr` / `interpolate_bands`);
localisation itself does not use eigenvalues.

---

## `seedname.amn` — trial projections (input)

The overlaps `A_{mn}(k) = ⟨ψ_{mk} | g_n⟩` of Bloch states with the trial orbitals `g_n`, used to
seed the initial gauge.

```
line 1: comment / date string
line 2: num_bands  num_kpts  num_proj
then num_bands·num_proj·num_kpts data lines:
    m   n   k    Re(A)   Im(A)
```

`m` = band, `n` = projection, `k` = k-index. Placement is by the explicit `(m,n,k)` on each
line, so file order is conventional only (band `m` fastest, then `n`, then `k`).

---

## `seedname.mmn` — overlaps (input)

The neighbour overlaps `M_{mn}^{(k,b)} = ⟨u_{mk} | u_{n,k+b}⟩`, the core quantity of the whole
method.

```
line 1: comment / date string
line 2: num_bands  num_kpts  nntot
then, for each of (num_kpts · nntot) blocks:
    k   k'   g1 g2 g3          # k' = neighbour k-index; g = reciprocal-lattice shift
    Re(M) Im(M)                # num_bands·num_bands lines, band m fastest (n outer, m inner)
```

Each block's `(k', g)` header identifies which neighbour `b` it is; the package matches these
against the b-vectors it derives from the mesh, so block order need not match the internal
ordering, but every `(k', g)` tuple must correspond to a computed neighbour.

---

## `seedname_hr.dat` — real-space Hamiltonian (output)

The interpolated tight-binding model `H(R)` — the primary product of Wannier interpolation.

```
line 1: date/time header
line 2: num_wann
line 3: nrpts
degeneracy block: ndegen(1..nrpts), 15 integers per line, format (15I5)
then nrpts·num_wann·num_wann rows:
    Rx  Ry  Rz   j   i   Re(H)   Im(H)          format (5I5,2F12.6)
```

Conventions:

- `R` (`Rx Ry Rz`) is a lattice vector in integer lattice-vector units; the row value is
  `H_{j,i}(R)` in **eV**.
- The inner index `j` (the row / left WF) varies fastest.
- `H(R)` is stored **undivided by `ndegen`**; a consumer interpolating `H(k)` must apply the
  weight `1/ndegen(R)` per R-vector. The set obeys `Σ_R 1/ndegen(R) = ∏ mp_grid`.
- The `F12.6` format is **lossy** (6 decimals): `_hr.dat` is a portable tight-binding export,
  not a full-precision checkpoint. This package computes `H(R)` internally from the localised
  gauge rather than round-tripping through the 6-decimal file.

---

## A note on `.chk` (checkpoint)

Wannier90's binary `.chk` (and its formatted `.chk.fmt`) is the only **full-precision** carrier
of the final gauge and overlaps. Round-tripping the `.chk`/`.chk.fmt` with `wannier90.x` is on
this package's roadmap and is not yet supported; the current pipeline reconstructs everything
from `.amn/.mmn/.eig` + the localisation it runs itself.
