# Wannier90 v3.1.0 — `.win` Parameters & Feature Inventory

Scope: `src/readwrite.F90`, `src/wannier90_readwrite.F90`, `src/types.F90`,
`src/wannier90_types.F90`. Reference tree:
`/Users/wolft/Dev/wannier90_greenfield/reference/wannier90`.

This document is the input-parameter spec for a from-scratch Julia reimplementation of
`wannier90.x` (wannierise + Wannier interpolation). All defaults and validation rules are
cited to `file:line`. Math conventions that matter for numerical reproduction are flagged
with **GOTCHA**.

---

## 0. Global conventions established here (read first)

- **Precision**: `dp = kind(1.0d0)` — IEEE double. `constants.F90:52`.
- **Internal length unit is Angstrom.** `length_unit` default `'ang'`
  (`types.F90:58`). Bohr input is converted *to* Angstrom on read.
- **`bohr` conversion constant = `0.52917721092`** (CODATA 2010, the default build;
  `constants.F90:210,224`). A legacy `USE_WANNIER90_V1_BOHR` build uses
  `0.5291772108` (`constants.F90:219`). **GOTCHA**: to match reference numerics you
  must use the *same* Bohr constant the reference binary was compiled with. Default =
  CODATA 2010. Julia should hard-code `0.52917721092`.
- **`length_unit` semantics** (`readwrite.F90:154-189`): valid values `'ang'` /
  `'bohr'`. If `'bohr'`, `lenconfac = 1/bohr` (this factor is used *only for output
  formatting*, converting internal Angstrom back to Bohr on print). It does **not**
  change internal storage. Input coordinate blocks carry their *own* units line (see
  §2). After reading, `length_unit(1:1)` is upper-cased for printout ('Ang'/'Bohr').
- **`energy_unit`** is read (`readwrite.F90:170`) but is effectively unused
  (comment `wannier90_readwrite.F90:254` "is this not used???"). Energies are eV
  throughout. **Do not apply any energy unit conversion.**
- **Fortran arrays are 1-based, column-major.** All `(3, N)` arrays store the 3-vector
  contiguously (fast index = Cartesian/component). Julia should keep the same layout for
  1:1 index correspondence.
- Input parsing is **case-insensitive** and lowercases keywords; block matching uses
  `begin <kw>` / `end <kw>` (`readwrite.F90:3099-3100`).
- Reading is *destructive*: once a keyword/block is consumed it is blanked from
  `settings%in_data` (e.g. `readwrite.F90:3225`). First occurrence wins for library
  `entries`; for `.win` text, a duplicate `begin`/`end` is an error
  (`readwrite.F90:3130,3148`).

---

## 1. Mandatory & core scalar parameters

Defaults come from the type declaration unless a read routine overrides. Column
"Default" = value if the keyword is absent.

### System size

| Keyword | Type | Default | Meaning / validation | Source |
|---|---|---|---|---|
| `num_wann` | int | **none (mandatory)** | Number of Wannier functions. Error if absent or `<= 0`. | `readwrite.F90:191-212` |
| `num_bands` | int | `num_wann` | Number of Bloch bands entering the calc. If absent → `num_wann`. Error if `num_bands < num_wann`. | `readwrite.F90:326-359` |
| `mp_grid` | int(3) | **none (mandatory)** | Monkhorst–Pack grid dims. Error if absent or any `< 1`. Sets `num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)`. | `readwrite.F90:388-420` |
| `exclude_bands` | int range list | (none) | Bands to drop before wannierisation (1-based). Range syntax `n1-n2,n3,...`. Must be positive. | `readwrite.F90:285-324` |
| `total_bands` | int | 0 | Convenience: total DFT bands, used with `exclude_bands` to derive `num_bands`. | `readwrite.F90:214-229` |
| `gamma_only` | logical | `.false.` | Gamma-point-only branch. Error if `.true.` and `num_kpts /= 1`. | `readwrite.F90:361-386` |

**Disentanglement trigger**: `disentanglement = (num_bands > num_wann)`
(`wannier90_readwrite.F90:137,1851,2506`). There is **no** separate `disentanglement`
keyword — it is purely derived. When `num_bands == num_wann` the disentangle step is
skipped entirely.

### System physics (`read_system`, `readwrite.F90:422-477`)

| Keyword | Type | Default | Meaning | Notes |
|---|---|---|---|---|
| `spinors` | logical | `.false.` | WFs are spinors. | Sets `num_elec_per_state = 1` when true. |
| `num_elec_per_state` | int | 2 (1 if spinors) | Spin degeneracy per state. Only 1 or 2 allowed; must be 1 if `spinors`. | `readwrite.F90:453-468` |
| `num_valence_bands` | int | (unset) | # valence bands (used by some spread/DOS paths). Must be `> 0` if given. | `readwrite.F90:470-476` |

### Verbosity / control (`read_verbosity`, `read_algorithm_control`)

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `iprint` | int | 1 | Output verbosity. `>= 2` also enables `svd_omega` printing. On non-root MPI ranks forced to 0. | `readwrite.F90:129-135`, `types.F90:53` |
| `timing_level` | int | 1 | Timing verbosity. | `readwrite.F90:116`, `types.F90:55` |
| `optimisation` | int | (module default) | Algorithm/optimisation level flag. | `readwrite.F90:149` |
| `length_unit` | char | `'ang'` | `'ang'` or `'bohr'` (output formatting only). | `readwrite.F90:173`, `types.F90:58` |
| `energy_unit` | char | (unused) | Read but ignored. | `readwrite.F90:170` |

---

## 2. Geometry blocks (units subtlety)

### `unit_cell_cart` block — `read_lattice` (`readwrite.F90:1110-1130`)

- 3×3 block of lattice vectors in **rows** of the input.
- **GOTCHA (transpose)**: internal storage is
  `real_lattice = transpose(real_lattice_tmp)` (`readwrite.F90:1125`). So the *i*-th
  input **row** (a lattice vector) becomes the *i*-th **column** of `real_lattice`. i.e.
  `real_lattice(:, i)` is lattice vector *i*. When you form reciprocal vectors / do
  cart↔frac you must respect this: `real_lattice` columns = a1,a2,a3.
- **Units line**: the block may have an optional first line `ang` or `bohr`. Handled in
  `get_keyword_block` (`readwrite.F90:3180-3223`): if `blen == rows+1` and the extra
  line says `bohr`, all real values are multiplied by `bohr` (→ Angstrom). `ang` =
  no-op. Only `unit_cell_cart` is allowed to have the extra units line among the 3-row
  blocks (`readwrite.F90:3175-3184`). Missing/unrecognised unit → error.
- Mandatory: error "Did not find the cell information" if absent (`readwrite.F90:1127`).

### `atoms_cart` / `atoms_frac` blocks — `read_atoms` (`readwrite.F90:1132-…`)

- Two mutually-relevant blocks:
  - `atoms_cart`: label + Cartesian position (Angstrom, or Bohr via units line).
  - `atoms_frac`: label + fractional (lattice) coordinates.
- Positions stored in `atom_data%pos_cart(:, species_atom, species)` (Angstrom),
  grouped by species. `species_num`, `label`, `symbol`, `num_atoms`, `num_species`
  populated (`types.F90:238-250`). Cart↔frac conversion uses `real_lattice`
  (`utility_cart_to_frac`).
- **GOTCHA**: atoms are optional for a pure wannierise/interpolation run — they only
  feed `.xyz` output, projection site defaulting, and plotting. Do not make them
  mandatory.

### `kpoints` block — `read_kpoints` (`readwrite.F90:991-1056`)

- Block of `num_kpts` rows × 3 columns of **fractional** k-coordinates
  (`kpt_latt(3, num_kpts)`).
- **GOTCHA (default k-mesh)**: if the `kpoints` block is absent, w90 *generates* a
  regular MP grid itself with triple loop ordering (`readwrite.F90:1035-1047`):
  ```
  ik=1
  do ia=1,mp_grid(1); do ib=1,mp_grid(2); do ic=1,mp_grid(3)
    kpt(:,ik) = ((ia-1)/mp1, (ib-1)/mp2, (ic-1)/mp3); ik++
  ```
  i.e. `ic` (3rd dim) is the fastest-varying index. If you generate k-points yourself,
  match this ordering exactly or the `.mmn`/`.amn`/`.eig` row order will mismatch.
- The variable is named `kpt_cart` internally but holds fractional coords (misnomer);
  no units conversion is applied to k-points.

---

## 3. Wannierise minimisation controls (`read_wannierise`, `wannier90_readwrite.F90:625-788`)

| Keyword | Type | Default | Validation | Source |
|---|---|---|---|---|
| `num_iter` | int | 100 | `>= 0` | `wannier90_readwrite.F90:649-656`; `wannier90_types.F90:192` |
| `num_print_cycles` | int | 1 | `>= 0` | `:640-647`; `types:190` |
| `num_dump_cycles` | int | 100 | `>= 0` (checkpoint write interval) | `:631-638`; `types:188` |
| `num_cg_steps` | int | 5 | `>= 0` | `:658-665`; `types:194` |
| `conv_tol` | real | `1.0e-10` | `>= 0` (spread convergence, Å²) | `:667-674`; `types:196` |
| `conv_window` | int | **-1** (disabled) | see below | `:681-685`; `types:197` |
| `conv_noise_amp` | real | -1.0 (disabled) | random-noise kick amplitude | `:676-678`; `types:204` |
| `conv_noise_num` | int | 3 | `>= 0` | `:687-694`; `types:205` |
| `trial_step` | real | 2.0 | line-search trial step; error if both `trial_step` and `fixed_step` set | `:732-739`; `types:201` |
| `fixed_step` | real | -999.0 (off) | if `> 0` sets `lfixstep=.true.`; must be `> 0` if given | `:722-730`; `types:200` |
| `precond` | logical | `.false.` | preconditioned CG | `:741-743`; `types:202` |
| `guiding_centres` | logical | `.false.` | enable guiding centres | `:696-698`; `wannier90_types:178` |
| `num_guide_cycles` | int | 1 | `>= 0` | `:704-711`; `types:179` |
| `num_no_guide_iter` | int | 0 | `>= 0` | `:713-720`; `types:180` |
| `use_ss_functional` | logical | `.false.` | Stengel–Spalding spread functional variant | `:700-702`; `types:199` |

**GOTCHA — `conv_window` default is coupled to `conv_noise_amp`**
(`wannier90_readwrite.F90:681-682`):
```
wann_control%conv_window = -1
if (wann_control%conv_noise_amp > 0.0_dp) wann_control%conv_window = 5
```
The `= 3` initializer in the type (`wannier90_types:197` has *no* initializer for
`wann_control%conv_window`; the value 3 in the doc comment is for `dis_control`) is
overridden. Effective default: **-1 (no convergence-window check)**, becoming **5** if
noise is enabled, unless the user sets it explicitly. Reproduce this exactly.

**Selective localisation / constrained centres** (`slwf_*`):

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `slwf_num` | int | `num_wann` | # objective WFs; `1..num_wann`; `< num_wann` ⇒ `selective_loc=.true.` | `:745-758` |
| `slwf_constrain` | logical | `.false.` | constrain centres; ignored unless selective_loc | `:760-776` |
| `slwf_lambda` | real | 1.0 | Lagrange multiplier, `>= 0` | `:778-787`; `types:172` |
| `slwf_centres` block | real(4×N) | — | per-WF centre constraints | `readwrite.F90:1739` |

---

## 4. Disentanglement controls (`read_disentangle`, `wannier90_readwrite.F90:791-880`; `read_dis_manifold`, `readwrite.F90:794-877`)

### Iteration controls

| Keyword | Type | Default | Validation | Source |
|---|---|---|---|---|
| `dis_num_iter` | int | 200 | `>= 0` | `wannier90_readwrite.F90:808-814`; `wannier90_types:143` |
| `dis_mix_ratio` | real | 0.5 | `0 < r <= 1` (else error) | `:816-822`; `types:145` |
| `dis_conv_tol` | real | `1.0e-10` | `>= 0` | `:824-830`; `types:147` |
| `dis_conv_window` | int | 3 | `>= 0` | `:832-838`; `types:149` |

### Energy windows (`read_dis_manifold`, `readwrite.F90:794-877`)

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `dis_win_min` | real | `-huge(dp)` | outer window lower bound (eV) | `readwrite.F90:809`; `types:213` |
| `dis_win_max` | real | `+huge(dp)` | outer window upper bound. Error if `win_max < win_min`. | `:813-821`; `types:215` |
| `dis_froz_min` | real | `-huge(dp)` | frozen (inner) window lower | `:829`; `types:217` |
| `dis_froz_max` | real | `+huge(dp)` | frozen (inner) window upper | `:823`; `types:219` |

**GOTCHA (frozen window logic)** (`readwrite.F90:823-840`):
- `dis_froz_max` present ⇒ `frozen_states = .true.` (energy frozen window active).
- Error if `dis_froz_min` present but `dis_froz_max` absent.
- Error if `froz_max < froz_min`.
- Default outer-window bounds are `±huge(0.0_dp)` **not** finite; effectively "all
  bands in window". Match this sentinel.

### Projectability-based freezing (`readwrite.F90:844-876`) — advanced, off by default

| Keyword | Type | Default | Meaning |
|---|---|---|---|
| `dis_froz_proj` | logical | `.false.` | use projectability frozen window instead of energy |
| `dis_proj_min` | real | 0.01 | lower projectability threshold, `∈ [0,1]` |
| `dis_proj_max` | real | 0.95 | upper projectability threshold, `∈ [0,1]`, `>= proj_min` |

### Disentanglement spheres (`dis_spheres_*`) — advanced

| Keyword | Type | Default | Notes |
|---|---|---|---|
| `dis_spheres_first_wann` | int | 1 | `>= 1`, `<= num_bands-num_wann+1` |
| `dis_spheres_num` | int | 0 | `>= 0`; if `>0` read `dis_spheres` block (4×N: kx,ky,kz,radius; radius `> 1e-15`) |

---

## 5. k-mesh / finite-difference b-vectors (`read_kmesh_data`, `readwrite.F90:879-989`)

These control the b-vector shell selection used to build the finite-difference operator
(the heart of the Marzari–Vanderbilt gradient). Critical for numerical match.

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `search_shells` | int | 36 | # shells to scan for satisfying B1. `>= 0`. | `readwrite.F90:897-903`; `types:149` |
| `search_supcell_size` | int | 5 | recip-cell supercell size to search for shells. `>= 0`. | `:904-910`; `types:150` |
| `kmesh_tol` | real | `1.0e-6` | tolerance for shell degeneracy / B1 tests. `>= 0`. | `:936-942`; `types:151` |
| `shell_list` | int list | (auto) | explicit shells to use (bypasses auto). `1..max_shells_h`. | `:944-971` |
| `num_shells` | int | (obsolete) | must equal len(shell_list) if given, else error. | `:973-979` |
| `skip_b1_tests` | logical | `.false.` | skip B1 condition check (Marzari–Vanderbilt PRB 56,12847). For Z2PACK etc. | `:986-988`; `types:146` |
| `kmesh_shell_from_file` | logical | `.false.` | read b-vector shells from file (complex cases) | `:932-934`; `types:144` |
| `higher_order_n` | int | 1 | higher-order finite-difference order. `>= 0`. Sets `max_shells_h = n(4n²+15n+17)/6`. | `:911-922`; `types:141` |
| `higher_order_nearest_shells` | logical | `.false.` (experimental) | if false, `max_shells_aux = 6` | `:924-929`; `types:142` |

**GOTCHA (b-vector ordering)**: b-vectors, weights `wb`, and neighbour lists
(`nnlist`, `nncell`) are derived in `kmesh.F90` (out of this scope file) but the
`.mmn`/`.nnkp` file b-vector ordering is fixed by that algorithm. `use_ss_functional`
adds reordering arrays (`nnord/nninv/nnrev`, `types:197-199`). For a faithful port,
reproduce the shell-selection and b-vector enumeration order from `kmesh.F90` — the
overlap-matrix column order depends on it. `max_shells = 6`, `num_nnmax = 12`
(`types.F90:128-129`).

**`nnkpts` block** (`wannier90_readwrite.F90:1434-1457`): explicit neighbour list,
allowed **only** in post-processing setup mode (`kmesh_info%explicit_nnkpts`).

---

## 6. Wigner–Seitz / real-space distance (`read_ws_data`, `readwrite.F90:694-737`)

Controls the interpolation R-vector treatment (crucial for band interpolation match).

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `use_ws_distance` | logical | `.true.` | use minimal-image WS distance when placing WFs for H(R) interpolation | `readwrite.F90:706`; `types:88` |
| `ws_distance_tol` | real | `1.0e-5` | absolute tol for "equivalent" (degenerate) WS distances | `:710`; `types:89` |
| `ws_search_size` | int(3) | `(2,2,2)` | supercell extent (each dir) searched for WS-cell points. 1 scalar → replicated to 3; or 3 ints. All `> 0`. | `:714-736`; `types:91` |

**GOTCHA (WS weighting/degeneracy)**: `ws_distance_type` (`types:97-113`) stores
`irdist/crdist/ndeg` — when a WF shift is on the WS-cell boundary, several equivalent
images are kept and the H(R) contribution is split with weight `1/ndeg`. Getting
`ndeg` and the tie-breaking tolerance right is essential to reproduce interpolated
bands. `use_ws_distance=.true.` is the **default** in v3.x; the older
`translate_home_cell`-only behaviour differs. See separate ws_distance notes.

---

## 7. Interpolation / plotting / output flags

### k-path for band interpolation (`read_kpath`, `readwrite.F90:479-542`)

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `kpoint_path` block | — | — | pairs of special points (label + frac coords). Stored as `2×nseg` labels/points. | `readwrite.F90:497-531`; `types:253-264` |
| `bands_num_points` | int | 100 | # points in first segment (`num_points_first_segment`); `>= 0` if bands_plot. | `:533-541`; `types:259` |
| `bands_plot` | logical | `.false.` | enable band-structure interpolation output | `wannier90_readwrite.F90:418`; `wannier90_types:53` |
| `bands_plot_format` | char | `'gnuplot'` | output format | `wannier90_readwrite.F90:1085`; `wannier90_types:110` |
| `bands_plot_mode` | char | `'s-k'` | slater-koster ('s-k') vs cut-off | `:1089`; `types:109` |
| `bands_plot_project` | int range | — | project onto WF subset | `:1093-1113` |
| `bands_plot_dim` | int | 3 (`system_dim`) | dimensionality for plotting | `:604`; `types:90` |
| `explicit_kpath` / `explicit_kpath_labels` blocks | — | — | user-supplied explicit k list for bands (cannot combine with `kpoint_path`) | `readwrite.F90:544-608`, `1058-1108` |

### H(R) output — `write_hr` / `hr_plot` (`wannier90_readwrite.F90:978-988`)

- `write_hr` (logical, default `.false.`) writes `seedname_hr.dat`
  (`types.F90` output_file, `wannier90_types:60`).
- **GOTCHA (alias)**: legacy keyword `hr_plot` (`:978`) is read into a temp and mapped
  to `write_hr` (`:984`). Support both spellings; they set the same flag.

### Wannier function real-space plotting (`read_wann_plot`, block reads at `wannier90_readwrite.F90:1164-1255`)

| Keyword | Type | Default | Meaning | Source (type: `wannier90_types.F90`) |
|---|---|---|---|---|
| `wannier_plot` | logical | `.false.` | enable WF plotting | `:54` |
| `wannier_plot_supercell` | int(1 or 3) | `(2,2,2)` | supercell for plotting | `:119` |
| `wannier_plot_format` | char | `'xcrysden'` | | `:122` |
| `wannier_plot_mode` | char | `'crystal'` | crystal vs molecule | `:123` |
| `wannier_plot_list` | int range | — | which WFs to plot | `readwrite ...:1205` |
| `wannier_plot_radius` | real | 3.5 | | `:121` |
| `wannier_plot_scale` | real | 1.0 | | `:121` |
| `wannier_plot_spinor_mode` | char | `'total'` | | `:124` |
| `wannier_plot_spinor_phase` | logical | `.true.` | | `:125` |
| `wvfn_formatted` | logical | `.false.` | read formatted UNKp files | `:132` |
| `spin` | char | (→ `spin_channel=1`) | up/down channel for plotting | `readwrite ...:1049` |

### Home-cell translation (`wannier90_readwrite.F90:1372-1400`)

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `translate_home_cell` | logical | `.false.` | translate WF centres into home cell (affects `.xyz`; future: H(R)) | `wannier90_types:95` |
| `translation_centre_frac` | real(3) | `(0,0,0)` | centre for translation (frac) | `wannier90_types:100` |
| `automatic_translation` | logical | `.true.` | (internal, derived) | `wannier90_types:102` |

### Initial guess / projections

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `projections` block | — | — | trial orbital definitions (site;l,m;radial;zona;z/x axes). | `readwrite.F90:3926-…` |
| `auto_projections` | logical | `.false.` | code supplies count only; external `.amn` (e.g. SCDM) | `wannier90_readwrite.F90:1574`; `wannier90_types:270` |
| `select_projections` | int range | — | subset of projections → WFs | `:1623` |
| `use_bloch_phases` | logical | `.false.` | use Bloch phases as initial guess (no `.amn`) | `wannier90_readwrite.F90:1400` |
| `spinors` | logical | `.false.` | (see §1) affects projection spin handling | |

**Projection per-orbital defaults** (`readwrite.F90:3926+`, `types.F90:174-181`):
`radial = 1`, `zona = 1.0`, quantisation axis `z = (0,0,1)`, `x = (1,0,0)`,
`s_qaxis = (0,0,1)`. `site`, `l`, `m` have no default (must be given).
`proj_zona_def = 1.0`, `proj_radial_def = 1` (`readwrite.F90:3975-3980`).

### Execution / restart

| Keyword | Type | Default | Meaning | Source |
|---|---|---|---|---|
| `postproc_setup` | logical | `.false.` | write `.nnkp` and stop (the `-pp` mode) | `wannier90_readwrite.F90:898`; `wannier90_types:51` |
| `restart` | char | `' '` | `''`/`default`/`wannierise`/`plot`/`transport`; requires `.chk`. Forced `''` if postproc_setup. | `wannier90_readwrite.F90:923-940` |
| `write_xyz` | logical | `.false.` | write WF centres `.xyz` | `wannier90_readwrite.F90:961` |
| `write_u_matrices` | logical | `.false.` | write `_u.mat`/`_u_dis.mat` | `:1019` |
| `write_bvec` | logical | `.false.` | write b-vectors/weights | `:1023` |
| `write_tb` | logical | `.false.` | write tight-binding `_tb.dat` | `:992` |
| `write_rmn` / `write_r2mn` / `write_proj` / `write_hr_diag` | logical | `.false.` | various matrix dumps | `:965-988` |
| `guiding_centres` | logical | `.false.` | (see §3) | |

---

## 8. Out-of-scope for greenfield core (present in these files but for postw90/transport/FS)

Read here but **not needed** for wannierise+interpolation:
- **Transport** (`tran_*`, `transport`, `transport_mode`): `wannier90_readwrite.F90:410-514`,
  `wannier90_types.F90:232-254`. Landauer–Büttiker — out of scope.
- **Fermi surface** (`fermi_surface_plot`, `fermi_surface_num_points` default 50,
  `fermi_surface_plot_format` `'xcrysden'`): `wannier90_types:222-228`. Low priority.
- **Site symmetry** (`site_symmetry`, `symmetrize_eps` default `1.d-3`):
  `wannier90_readwrite.F90:386-391`, `wannier90_types:275-283`. Advanced; skip initially.
- `dist_cutoff*`, `hr_cutoff`, `one_dim_axis`, `system_dim` — mostly plot/transport.
- `fermi_energy`, `write_vdw_data`, `dump_inputs`, `distk` (library-only MPI dist),
  `cp_pp`, `calc_only_A`.

---

## 9. FEATURE INVENTORY for greenfield `wannier90.x` core (drives architecture)

### MUST-HAVE (minimum to reproduce a standard wannierise + interpolation run)

1. **`.win` parser**: keyword scalars (int/real/logical/char), vectors, ranges
   (`n1-n2,n3`), and blocks with optional `ang`/`bohr` units line. Case-insensitive,
   destructive-consume, duplicate detection. (`readwrite.F90:2751-3646`).
2. **Core parameters**: everything in §1, §3, §4, §5, §6. Especially the derived
   `disentanglement = num_bands > num_wann` and the coupled `conv_window` default.
3. **Geometry**: `unit_cell_cart` (with transpose!), `kpoints` (with auto-MP-grid
   fallback + exact loop ordering), `mp_grid` → `num_kpts`.
4. **Reading external matrices**: `.eig` (`read_eigvals`, `readwrite.F90:739-792`, with
   band/kpoint index-match validation), `.amn`, `.mmn` (in overlap.F90 — separate scope).
5. **k-mesh / b-vectors**: shell search, B1 condition, weights `wb`, neighbour lists —
   the finite-difference backbone (kmesh.F90 scope).
6. **Disentanglement** (energy-window path only to start): outer/frozen windows,
   `dis_num_iter`, `dis_mix_ratio`, `dis_conv_tol/window`.
7. **Wannierise (MLWF)**: MV gradient, CG (`num_cg_steps`), line search
   (`trial_step`/`fixed_step`), `conv_tol`/`conv_window`, spread decomposition.
8. **WS-distance module** (§6) + **band interpolation** (`kpoint_path`,
   `bands_num_points`, `bands_plot`) — the "interpolation" half of the task.
9. **`write_hr`/`hr_plot`** (H(R) output) and **`.chk`** checkpoint read/write
   (`readwrite.F90:2232-2668`) for restart and interop.
10. Projections initial guess (`projections` block) OR external `.amn` ingestion; plus
    `exclude_bands`.

### NICE-TO-HAVE (add after core reproduces reference)

- `wannier_plot` (real-space WF cubes/xsf) and UNKp reading.
- `guiding_centres` + `num_guide_cycles`/`num_no_guide_iter`.
- `precond` preconditioned CG, `conv_noise_*` noise-kick escape.
- `use_bloch_phases`, `auto_projections`/SCDM interop, `select_projections`.
- Selective localisation / constrained centres (`slwf_*`).
- `spinors` support (spinor WFs).
- `gamma_only` optimised branch.
- `write_tb`, `write_u_matrices`, `write_bvec`, `write_xyz`, `write_rmn`.
- Projectability disentanglement (`dis_froz_proj`, `dis_proj_*`), `dis_spheres_*`.
- Higher-order finite differences (`higher_order_n`).

### OUT-OF-SCOPE (postw90 / niche — do not build into wannier90.x core)

- All `tran_*` transport, `transport_mode`, quantum conductance (transport.F90).
- `fermi_surface_plot`.
- Site symmetry (`site_symmetry`, `symmetrize_eps`, sitesym.F90).
- `write_vdw_data`, vdW C6.
- Library/MPI-distribution plumbing (`distk`), `cp_pp`, `calc_only_A`, `dump_inputs`.
- Everything under `src/postw90/`.

---

## 10. Silent-mismatch checklist (port carefully)

1. **Bohr constant**: `0.52917721092` (CODATA 2010 default). A different value silently
   shifts all Bohr-input geometries.
2. **`real_lattice = transpose(input)`** — lattice vectors are **columns** internally.
3. **Auto k-mesh loop order**: dim-3 fastest (`readwrite.F90:1035-1047`).
4. **`conv_window` default = -1** (not 3), → 5 if `conv_noise_amp>0`.
5. **`disentanglement` is derived**, not an input flag.
6. **Frozen window sentinels** `±huge(dp)`; `dis_froz_max` presence toggles
   `frozen_states`.
7. **Units line on blocks** multiplies by `bohr`; only `unit_cell_cart` may carry it
   among 3-row blocks; k-points never converted.
8. **`length_unit`/`lenconfac` are output-only**; internal storage stays Angstrom.
9. **`energy_unit` ignored** — energies are eV; do not convert.
10. **b-vector / neighbour ordering** from kmesh must match, or `.mmn`/overlap column
    order diverges (verify against kmesh.F90 notes).
11. **WS `ndeg` weighting** (`1/ndeg` on boundary images) governs interpolated H(R).
12. **`num_elec_per_state`** flips to 1 under `spinors` — affects any electron-count /
    occupation logic.
