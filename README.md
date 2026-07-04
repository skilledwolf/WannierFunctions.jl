# Wannier90.jl

A modern Julia reimplementation of the Wannier90 core (`wannier90.x`): construction of
**maximally-localised Wannier functions** (MLWF) and **Wannier interpolation** of band
structures.

This is an **independent, from-scratch reimplementation**. It is not affiliated with, nor
derived from the source of, the official [Wannier90](https://wannier.org) project. It aims for
**drop-in compatibility** with the standard Wannier90 input files — `.win`, `.amn`, `.mmn`,
`.eig` — so that an existing DFT → Wannier90 workflow can be pointed at this package with no
changes to the files produced by your DFT interface (e.g. `pw2wannier90`).

## What it does

Given the projections `A_{mn}(k)`, overlaps `M_{mn}^{(k,b)}`, and eigenvalues `ε_{mk}` that a
DFT code writes for Wannier90, this package:

- builds the finite-difference **b-vector** shells and B1 weights from the k-mesh;
- computes the **Marzari–Vanderbilt spread functional** `Ω = Ω_I + Ω_OD + Ω_D` and its
  Wannier centres;
- **localises** the gauge by conjugate-gradient minimisation of `Ω` (unitary update via the
  matrix exponential of an anti-Hermitian generator);
- Fourier-transforms the Wannier-gauge Hamiltonian `H(k) → H(R)` on the Wigner–Seitz
  supercell and **interpolates** band energies on an arbitrary k-path.

## Validation status

The full pipeline is validated against the **reference Wannier90 v3.1.0 test suite** (the golden
`.wout` benchmarks shipped with the reference) and, for band interpolation, against a locally
built reference `wannier90.x` binary:

| Milestone | Scope | Status |
|-----------|-------|--------|
| **M0** | I/O, k-mesh + b-vectors, initial gauge, spread | ✅ validated |
| **M1** | Maximal localisation (MLWF gauge optimisation) | ✅ validated |
| **M2** | Wannier interpolation `H(R) → H(k)`, band structure | ✅ validated |
| **M3** | Disentanglement (Souza–Marzari–Vanderbilt Ω_I minimisation) | ✅ validated |

Spread components (`Ω_I / Ω_D / Ω_OD / Ω_Total`) reproduce the reference to **~1e-8** (well inside
the test-suite tolerances of 1e-6 absolute on the Ω components and 1e-5 Å on centres, which
require both the absolute and relative bounds to hold). Validated cases:

| Case | System | Bands → WF | Disentangle | Final Ω (mine / ref) |
|------|--------|-----------|-------------|----------------------|
| `testw90_example01` | GaAs | 4 → 4 | no | 4.466881 / 4.466880976 |
| `testw90_example05` | diamond | 4 → 4 | no | 2.320904915 / 2.320904915 |
| `testw90_example03` | silicon | 12 → 8 | **yes** (frozen) | 14.499574503 / 14.499574503 |
| `testw90_example04` | copper (metal) | 12 → 7 | **yes** (frozen) | 4.028040094 / 4.028040058 |

For silicon the disentanglement Ω_I convergence trace matches the reference **iteration by
iteration** (12.70775084 → 11.99932157 → … → 11.849193709). Interpolated band structures match
the reference `wannier90.x` to **~1e-5 eV** across the k-path (diamond and silicon).

The tiny residuals on GaAs/copper (~4e-8) are the CODATA bohr-radius constant at its last digit
(their `.win` cells are in bohr); diamond/silicon (cells in Å) match to machine precision. No
benchmark numbers are invented — see `docs/reference-notes/` for provenance.

The four reference cases are all FCC single-shell meshes, so the test suite additionally covers the
multi-shell B1 weight solve on a synthetic tetragonal mesh (two shells, completeness to 1e-16), the
`use_ws_distance` grid invariant, and the k-path discontinuity handling — 209 tests in total.

Not yet implemented: the Γ-only real-gauge minimiser, `.chk` interchange, guiding-centre branch
selection, and `postw90.x` post-processing — see the roadmap. (The `use_ws_distance` minimal-image
interpolation refinement — the reference default — *is* implemented and validated against the
reference binary to ~2e-5 eV.)

## Installation

Requires **Julia ≥ 1.10** and **StaticArrays** (the only runtime dependency beyond the
standard library `LinearAlgebra` / `Printf`).

```julia
pkg> add StaticArrays
pkg> dev /path/to/wannier90_greenfield     # or: add https://…  once published
```

Or run scripts directly against the project:

```bash
julia --project=/path/to/wannier90_greenfield yourscript.jl
```

## Quickstart

Using the diamond example (`num_wann = 4`, 4×4×4 mesh; `.win/.amn/.mmn/.eig` present):

```julia
using Wannier90
using StaticArrays

seed  = "scratch/diamond_bands/diamond"           # seedname, no extension
model = read_model(seed)                           # reads .win/.amn/.mmn/.eig

res = wannierise(model; num_iter = 20)             # MLWF gauge optimisation
sp  = res.spread
println("Ω = ", sp.Ω, "  (ΩI=", sp.ΩI, " ΩOD=", sp.ΩOD, " ΩD=", sp.ΩD, ")")
# → Ω = 2.320904914…  (ΩI=1.954619860  ΩOD=0.366285055  ΩD≈0)

irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
Hr, _ = build_hr(res.U, model.eig, model.kgrid, irvec)   # H(k) → H(R)

kpath = [SVector(0.0, 0.0, 0.0), SVector(0.5, 0.0, 0.5)] # Γ, X (fractional)
bands = interpolate_bands(Hr, irvec, ndegen, kpath)      # (num_wann × npts), ascending
println("bands at Γ = ", bands[:, 1])
```

`SpreadResult` carries `centres::Matrix{Float64}` (3 × num_wann, Cartesian Å),
`spreads::Vector{Float64}` (Å²), and `Ω, ΩI, ΩOD, ΩD`. `wannierise` returns the per-k gauge
`U` (num_wann × num_wann × nkpt), the rotated overlaps `Mrot`, the `omega_trace`, `niter`, and
`converged`.

> Note: Wannier90's default convergence check is off (`conv_window = -1`), so the loop runs the
> full `num_iter`; `converged` reflects that and is not an error.

### One call for either case: `run_wannier`

`run_wannier` auto-selects the isolated-bands path or disentanglement (`num_bands > num_wann`,
using the `dis_*` energy windows in the `.win`) and returns a uniform result you can interpolate:

```julia
model, win, res = run_wannier("scratch/silicon/silicon")   # silicon: 12 bands → 8 WF, frozen window
println(res.disentangled, "  Ω = ", res.spread.Ω)          # true   Ω = 14.499574503

kpts, xvals, labels, idx = generate_kpath(win, model.lattice; bands_num_points = 100)
bands = interpolate(model, res, kpts)                       # (num_wann × npts) interpolated energies
```

### Command line (drop-in for `wannier90.x`)

```bash
julia --project=/path/to/wannier90_greenfield bin/wannier90.jl silicon
```

Reads `silicon.win` (+ `.amn/.mmn/.eig`), runs the full pipeline, and writes `silicon.wout` plus,
when the `.win` requests them, `silicon_hr.dat` (`write_hr`/`hr_plot`), `silicon_tb.dat`
(`write_tb`), and the band files `silicon_band.dat/.kpt/.labelinfo.dat` (`bands_plot`) — the same
outputs, in the same formats, as `wannier90.x silicon`.

## Roadmap

- **`.chk` / `.chk.fmt` interchange** for full-precision round-tripping with `wannier90.x`.
- **Position operator** `r(R)`, `_r.dat`, and Wannier-gauge observables (Berry-phase quantities).
- Γ-only (real-gauge) minimiser path (a distinct algorithm from the general complex minimiser).
- Projectability-based (`dis_froz_proj`) and symmetry-adapted (SAWF) variants; `postw90.x`
  post-processing (BoltzWann, AHC, …) is out of scope for the core.

## Documentation

- [`docs/theory.md`](docs/theory.md) — the math: spread functional, b-vectors, gauge gradient,
  disentanglement, and interpolation.
- [`docs/file-formats.md`](docs/file-formats.md) — the `.win/.amn/.mmn/.eig/_hr.dat` formats.
- [`docs/migrating-from-wannier90.md`](docs/migrating-from-wannier90.md) — mapping an existing
  Wannier90 workflow onto this package.

## References

- N. Marzari and D. Vanderbilt, *Maximally localized generalized Wannier functions for
  composite energy bands*, Phys. Rev. B **56**, 12847 (1997).
- I. Souza, N. Marzari, and D. Vanderbilt, *Maximally localized Wannier functions for entangled
  energy bands*, Phys. Rev. B **65**, 035109 (2001).
- N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, and D. Vanderbilt, *Maximally localized
  Wannier functions: Theory and applications*, Rev. Mod. Phys. **84**, 1419 (2012).
