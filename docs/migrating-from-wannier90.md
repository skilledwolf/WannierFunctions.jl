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
using Wannier90, StaticArrays

model = read_model("seedname")          # reads .win/.amn/.mmn/.eig
res   = wannierise(model; num_iter = 20)   # ← replaces the wannier90.x minimiser

irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
Hr, _ = build_hr(res.U, model.eig, model.kgrid, irvec) # ← the H(R) that goes into _hr.dat
bands = interpolate_bands(Hr, irvec, ndegen, kpath)    # ← replaces bands_plot
```

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

**Not yet supported:**

- **`.chk` / `.chk.fmt`** read/write for full-precision interchange with `wannier90.x`.
- **`guiding_centres`** branch selection for the `Im ln` sheet (default off; the CLI warns if you
  set it and falls back to the principal branch).
- **Position operator** `r(R)` / `_r.dat` and Berry-phase observables (`_tb.dat` is written with a
  zero r-block placeholder).
- **Γ-only** real-gauge minimiser (a distinct algorithm in the reference); use the general path.
- Projectability (`dis_froz_proj`) / symmetry-adapted (SAWF) variants, and `postw90.x`
  post-processing (out of scope for the core).

## Behavioural notes to expect

- **`converged` is often `false` at the end.** Wannier90's convergence-window check is off by
  default (`conv_window = -1`), so the loop runs the full `num_iter`. This package matches that;
  `converged = false` simply means the optional early-stop criterion was not enabled, not that
  the run failed. Inspect `res.omega_trace` to confirm the spread has plateaued.
- **Units.** Centres and spreads come back in Å / Å² even if your `.win` cell is in `bohr`,
  matching the `.wout` convention. Energies are eV.
- **Branch cut.** The `Im ln` principal branch is used (as in the default, no-guiding-centres
  Wannier90 path); results match provided your starting projections are reasonable.
