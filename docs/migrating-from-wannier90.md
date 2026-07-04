# Migrating from Wannier90

If you already run `wannier90.x`, this package slots into the same place in your workflow. It
reads the **same input files** your DFT interface produces — `.win`, `.amn`, `.mmn`, `.eig` — so
no changes to the DFT → Wannier90 setup are required. What changes is that you drive the
localisation and interpolation from Julia and get the results back as native data structures
instead of parsing a `.wout`.

## The standard Wannier90 pipeline vs. this package

A typical Wannier90 run:

```
1. SCF + NSCF (DFT)                         → wavefunctions
2. wannier90.x -pp seedname                 → seedname.nnkp (b-vectors requested)
3. pw2wannier90 (or your interface)         → seedname.amn/.mmn/.eig
4. wannier90.x seedname                     → localise + interpolate → seedname.wout, _hr.dat, _band.dat
```

Steps 1–3 are unchanged: keep producing `.amn/.mmn/.eig` exactly as before. This package
replaces **step 4**:

```julia
using WannierFunctions

model = read_model("seedname")               # reads .win/.amn/.mmn/.eig
res   = run_wannier(model)                   # ← replaces the wannier90.x minimiser
                                             #   (add win_max=…, froz_max=… for entangled bands)
H = hamiltonian_operator(model, res)         # ← the H(R) that goes into _hr.dat
E = bands(H, kpath)                          # ← replaces bands_plot
```

(And step 2 is covered too: `bin/wannier90.jl -pp seedname` writes the `.nnkp`.)

## Where your `.win` keywords map

| `.win` keyword | Effect here |
|----------------|-------------|
| `num_wann`, `num_bands` | read from `.win`; `num_bands > num_wann` triggers disentanglement automatically |
| `mp_grid` | drives the k-mesh, b-vectors, and the Wigner–Seitz set |
| `num_iter` | passed to `wannierise(model; num_iter = …)` |
| `unit_cell_cart`, `atoms_frac/cart` | build the lattice; cell may be in `bohr` |
| `projections` | consumed to build the initial gauge from the `.amn` |
| `kpoints` | the mesh; `kpoint_path` gives the interpolation path |
| `dis_win_*`, `dis_froz_*`, `dis_num_iter`, `dis_mix_ratio` | disentanglement — **supported**; `run_wannier` auto-selects it when `num_bands > num_wann` |
| `bands_plot`, `bands_num_points`, `write_hr`, `write_tb` | honoured by the `bin/wannier90.jl` CLI to write the band / `_hr.dat` / `_tb.dat` files |

## What you get back

Instead of scraping the `.wout`, you read fields directly:

- `res.spread.Ω, .ΩI, .ΩOD, .ΩD` — the spread decomposition (the `Omega I/D/OD/Total` lines of a
  `.wout`), in Å².
- `res.spread.centres` (3 × num_wann, Cartesian Å) and `res.spread.spreads` (Å²) — the per-WF
  centre-and-spread lines.
- `res.U` — the final gauge, per k-point.
- `res.omega_trace`, `res.niter`, `res.converged` — the convergence trace.

These reproduce the reference `.wout` numbers to the test-suite tolerances (~1e-6 on the Omega
components, ~1e-5 Å on centres), for both the isolated-bands and the disentanglement cases.

## Supported vs. not-yet-supported

**Supported and validated:**

- Reading `.win/.amn/.mmn/.eig`.
- b-vector shells and B1 weights from the mesh.
- Initial (Löwdin-projected) gauge, centres, and spread.
- MLWF localisation for `num_bands == num_wann` (isolated bands).
- **Disentanglement** (`num_bands > num_wann`) with outer + frozen energy windows
  (Souza–Marzari–Vanderbilt Ω_I minimisation) — validated on silicon and copper.
- Wannier interpolation: `H(k) → H(R)` on the Wigner–Seitz set, and band interpolation on a
  k-path, including the `use_ws_distance` minimal-image refinement (the reference default;
  validated against the reference binary to ~2e-5 eV).
- **Output writers**: `.wout`, `_hr.dat`, `_tb.dat`, `_band.dat/.kpt/.labelinfo.dat`, driven by
  the `bin/wannier90.jl` command-line front end (a drop-in for `wannier90.x`).

**Also supported:**

- **`-pp` mode**: `bin/wannier90.jl -pp seed` (or `postproc_setup = .true.`) generates the k-mesh
  from the `.win` alone and writes `seed.nnkp`, byte-identical to `wannier90.x -pp` — so a new
  DFT → Wannier workflow can start here, not just finish here.
- **Position operator** `⟨0m|r|Rn⟩`: `position_operator(model, res)`, and `_tb.dat` is written
  with real r-blocks (validated against the reference binary to the E15.8 file precision).
- **Strict input validation**: unknown `.win` keywords error with a did-you-mean suggestion
  (checked against the reference parser's own keyword catalogue); recognised-but-unsupported
  keywords warn once. Pass `read_win(path; strict=false)` to downgrade to warnings.
- **Two optimisers**: the CLI/.win path uses `:w90` (reference-exact); the Julia-native API
  defaults to `:rcg` (Riemannian CG with real convergence). Same minima, verified by tests.

**Also supported (Phase 1):** the binary `.chk` checkpoint in both directions (write ours,
`wannier90.x restart=plot` / `postw90.x` consume it), `gamma_only` inputs (half b-set expanded
exactly), `spinors` in `-pp`, `exclude_bands`, and `wannier_plot` (formatted UNK → `.xsf`).

**Not yet supported:**

- **`.chk.fmt`** (the formatted transport variant; binary `.chk` works).
- **`guiding_centres`** branch selection for the `Im ln` sheet (default off; the CLI warns if you
  set it and falls back to the principal branch).
- `_r.dat` and Berry-phase observables (the position operator itself is implemented; `_tb.dat`
  carries real r-blocks).
- **Γ-only real-orthogonal parity**: Γ inputs run on the (exact) complex path and may converge
  marginally below the reference's real-gauge value; a parity mode is on the roadmap.
- Projectability (`dis_froz_proj`) / symmetry-adapted (SAWF) variants, and `postw90.x`
  post-processing (out of scope for the core).

## Behavioural notes to expect

- **Convergence semantics differ by optimiser.** The `.win`/CLI path uses `:w90`, where —
  matching Wannier90 — the convergence check is off by default (`conv_window = -1`) and the loop
  runs the full `num_iter`; `converged = false` there just means the optional early stop wasn't
  enabled. The Julia-native default `:rcg` has a real convergence criterion and `converged`
  means what it says.
- **Units.** Centres and spreads come back in Å / Å² even if your `.win` cell is in `bohr`,
  matching the `.wout` convention. Energies are eV.
- **Branch cut.** The `Im ln` principal branch is used (as in the default, no-guiding-centres
  Wannier90 path); results match provided your starting projections are reasonable.
