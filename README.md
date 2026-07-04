# WannierFunctions.jl

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

Spread components (`Ω_I / Ω_D / Ω_OD / Ω_Total`) reproduce the reference to **print precision
(≤5e-10)** — far inside the test-suite tolerances of 1e-6 absolute on the Ω components and 1e-5 Å
on centres. Validated cases:

| Case | System | Bands → WF | Disentangle | Final Ω (mine / ref) |
|------|--------|-----------|-------------|----------------------|
| `testw90_example01` | GaAs | 4 → 4 | no | 4.466880976 / 4.466880976 |
| `testw90_example05` | diamond | 4 → 4 | no | 2.320904915 / 2.320904915 |
| `testw90_example03` | silicon | 12 → 8 | **yes** (frozen) | 14.499574503 / 14.499574503 |
| `testw90_example04` | copper (metal) | 12 → 7 | **yes** (frozen) | 4.028040058 / 4.028040058 |

For silicon the disentanglement Ω_I convergence trace matches the reference **iteration by
iteration** (12.70775084 → 11.99932157 → … → 11.849193709). Against a locally built
`wannier90.x`: interpolated bands match to ~1e-5 eV (the `_band.dat` file precision), `_hr.dat`
is byte-identical, the `_tb.dat` position-operator block matches to 1.3e-8 Å (the E15.8 file
floor), and the `.nnkp` written by `-pp` mode is **byte-identical** (modulo the date header) on
all four systems. The bohr constant is the reference's default **CODATA2006**, so bohr-specified
cells reproduce reference numbers exactly.

The four reference cases are all FCC single-shell meshes, so the test suite additionally covers
the multi-shell B1 weight solve on a synthetic tetragonal mesh (two shells, completeness to
1e-16), the `use_ws_distance` grid invariant, k-path discontinuity handling, optimizer parity,
and the operator-API invariants — **250 tests**, passing identically at 1 and 8 threads.

**Checkpoint interchange works in both directions**: we read `wannier90.x`-written `.chk` files
exactly, and `wannier90.x restart=plot` consumes ours — for the disentangled silicon case it
reproduces its own band structure to 1e-6 eV from our checkpoint. The CLI writes `seedname.chk`
after every run, so `postw90.x` can post-process our results directly.

Not yet implemented: the Γ-only real-gauge minimiser, guiding-centre branch selection, spinors,
WF plotting, and `postw90.x` post-processing — see the roadmap.

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

One call runs the whole pipeline — it auto-selects disentanglement when `num_bands > num_wann`
and returns a result you can interpolate. Using the shipped diamond example:

```julia
using WannierFunctions

model = read_model("examples/data/diamond")   # .win/.amn/.mmn/.eig; rich show in the REPL
res   = run_wannier(model)                    # :rcg optimiser, true convergence
res.spread                                    # Ω = 2.320904915 Å², centres, Ω_I/Ω_OD/Ω_D

H = hamiltonian_operator(model, res)          # interpolable operator H(R)
E = bands(H, [[0.0, 0.0, 0.0], [0.5, 0.0, 0.5]])   # energies at Γ and X (fractional k)

r = position_operator(model, res)             # ⟨0m|r|Rn⟩, 3 components — same TBOperator type
```

Everything the `.wout` would tell you is a field: `res.spread.Ω/.ΩI/.ΩOD/.ΩD`,
`res.spread.centres` (3 × num_wann, Å), `res.converged`, `res.niter`, and for disentangled runs
`res.dis.omega_I_trace`. For entangled bands pass the windows as plain keywords:

```julia
model = read_model("scratch/silicon/silicon")           # 12 bands → 8 WF
res = run_wannier(model; win_max = 17.0, froz_max = 6.4, dis_num_iter = 120, dis_mix_ratio = 1.0)
res.disentangled, res.spread.Ω                          # true, 14.499574503
```

Two optimisers are built in: `algorithm = :rcg` (default — Riemannian Polak–Ribière+ conjugate
gradient on the product-of-unitaries manifold, with a real convergence criterion; typically
converges in fewer iterations than a fixed sweep count) and `algorithm = :w90` (the reference
optimiser, reproduced exactly for drop-in parity — it is what the CLI uses). Both land on the
same minima; the test suite asserts it.

Input handling is strict: an unknown `.win` keyword is an **error with a did-you-mean
suggestion** (`num_itre` → "did you mean `num_iter`?"), checked against a catalogue of 278
keywords generated from the reference parser's source. Recognised-but-unsupported keywords warn
once and are ignored.

k-loops are threaded (gated by problem size, so small systems don't pay scheduling overhead):
start Julia with `julia -t auto` for dense interpolation workloads. `benchmark/run.jl` has the
numbers.

### Command line (drop-in for `wannier90.x`)

```bash
julia --project=/path/to/wannier90_greenfield bin/wannier90.jl -pp silicon   # setup: writes silicon.nnkp
julia --project=/path/to/wannier90_greenfield bin/wannier90.jl silicon       # full run
```

`-pp` generates the k-mesh from `silicon.win` alone and writes `silicon.nnkp` for the DFT
interface (pw2wannier90 etc.) — byte-identical to `wannier90.x -pp` output. The full run reads
`silicon.win` (+ `.amn/.mmn/.eig`), uses the reference-faithful `:w90` optimiser, and writes
`silicon.wout` plus, when the `.win` requests them, `silicon_hr.dat` (`write_hr`/`hr_plot`),
`silicon_tb.dat` with real position-operator blocks (`write_tb`), and the band files
`silicon_band.dat/.kpt/.labelinfo.dat` (`bands_plot`) — the same outputs, in the same formats, as
`wannier90.x silicon`.

**Γ-only calculations are supported** (silane, benzene validated): the half-stored b-vector set
is expanded to the closed full set on load using M(−b) = M(b)†, so the general machinery is exact.
One documented difference: we minimise over the full unitary group, while the reference Γ code
restricts to real orthogonal gauges — our converged spread can come out marginally *lower*
(benzene: 12.9583353 vs the reference's stationary 12.9583380; Ω_I identical to 3e-12).

**Spinor projections** (`spinors = .true.`): `-pp` writes the `spinor_projections` block —
byte-identical to `wannier90.x -pp` on the Pt (SOC) test input.

**Wannier-function plotting** (`wannier_plot = .true.`): reads formatted UNK files and writes
XCrySDen `.xsf` volumetric grids — all four GaAs Wannier functions match the reference output to
the file's e13.5 precision, including the global phase convention.

## Post-processing: Berry curvature and the anomalous Hall conductivity

The first slice of `postw90.x` functionality is implemented on the operator layer
(`berry_task = ahc`): `BerryModel(seedname)` consumes any completed run (checkpoint + `.eig` +
`.mmn`), builds H(R) and the Berry connection A(R) with the exact postw90 conventions, and
`anomalous_hall(bm; fermi_energy, kmesh)` integrates the occupied-manifold Berry curvature
(WYSV J0/J1/J2 decomposition). Validated on the bcc-Fe reference case — spinor, 28 bands
disentangled to 18 WFs — reproducing `postw90.x` **to every printed digit**:
σ = (0.0334, 0.0572, 1222.1510) S/cm on the benchmark 10³ mesh, in under a second on 8 threads.
A time-reversal-symmetric insulator (diamond) integrates to zero as it must. The CLI honours
`berry = true` / `berry_task = ahc` / `berry_kmesh` / `fermi_energy` from the `.win`.

**SCDM automatic projections**: `scdm_projections(model; dir)` computes initial projections
directly from the UNK wavefunctions (QRCP column selection, isolated/erfc/gaussian smearing) —
no `projections` block needed. On GaAs, an SCDM start converges to the identical
gauge-invariant minimum (Ω = 4.466880976) as the hand-chosen sp³ projections.

Both `.chk` (binary) and `.chk.fmt` (formatted) checkpoints read and write, validated bit-exact
against `w90chk2chk.x` conversions.

## Roadmap

- More postw90 surface on the same layer: orbital magnetisation, spin Hall, optical (Kubo),
  Boltzmann transport; adaptive k-mesh refinement.
- `guiding_centres` branch selection; Γ-only real-orthogonal parity mode; `.cube` plot output.
- Berry-phase observables on top of `TBOperator` (the position operator is already in place).
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
