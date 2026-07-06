# WannierFunctions.jl

[![CI](https://github.com/skilledwolf/WannierFunctions.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/skilledwolf/WannierFunctions.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://tobiaswolf.net/WannierFunctions.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A modern, from-scratch Julia implementation of the complete **Wannier90 tool chain** —
maximally-localised Wannier functions, disentanglement, Wannier interpolation, and the full
`postw90` physics surface — **drop-in compatible** with the standard `.win/.amn/.mmn/.eig`
files, and also usable entirely **in memory** from a Julia DFT code (DFTK.jl), with no files
and no external binaries.

It is an independent reimplementation, not affiliated with nor derived from the source of the
official [Wannier90](https://wannier.org) project. Every feature is validated against
reference oracles — see [Validation](docs/src/validation.md) for the methodology and numbers
(454 tests; spread components typically reproduce the reference to print precision).

## What's inside

| | |
|---|---|
| **Core** | Marzari–Vanderbilt localisation (`:w90` reference-exact and `:rcg` Riemannian CG optimisers), Souza–Marzari–Vanderbilt disentanglement, Wigner–Seitz Wannier interpolation, Γ-only real-orthogonal algorithm, spinors |
| **Drop-in binaries** | [`bin/wannier90.jl`](bin/wannier90.jl) (incl. `-pp`), [`bin/postw90.jl`](bin/postw90.jl), [`bin/w90chk2chk.jl`](bin/w90chk2chk.jl) — same inputs, same output files, same formats |
| **postw90 physics** | Berry curvature + AHC (adaptive mesh, Fermi scans), orbital magnetisation, Kubo optical conductivity + JDOS, spin Hall (Qiao + Ryoo, tetrahedron method), shift current, k·p expansion, gyrotropic tensors, DOS (spin-decomposed, projected), BoltzWann, kpath/kslice, geninterp, spin moments |
| **Advanced wannierisation** | SCDM automatic projections, PDWF projectability disentanglement, `dis_spheres`, guiding centres, preconditioned CG, SLWF+C selective localisation, symmetry-adapted WFs (incl. constrained disentanglement), Stengel–Spaldin functional, higher-order finite differences |
| **Beyond parity** | Ballistic (Landauer) transport, tight-binding model input (`_hr.dat`/`_tb.dat` → the whole physics stack), irreducible-BZ symmetrised integration, circular injection current, FermiSurfer export, in-memory **DFTK.jl bridge** with projection-free SCDM |

The one reference feature deliberately out of scope is lead–conductor–lead transport
(`tran_lcr`); see [compatibility notes](docs/src/wannier90-compat.md) for that and every other
known behavioural difference.

## Installation

Requires **Julia ≥ 1.10**; the only runtime dependency is StaticArrays.

```julia
pkg> add https://github.com/skilledwolf/WannierFunctions.jl   # (not yet registered)
```

To also get the drop-in command-line binaries, either clone the repository and run them from
`bin/`, or install launchers from the package:

```julia
julia> using WannierFunctions
julia> install_cli()    # writes wannier90.jl / postw90.jl / w90chk2chk.jl to ~/.julia/bin
```

## Quick start

**Drop-in, from files** (any DFT interface that writes Wannier90 inputs):

```bash
wannier90.jl -pp silicon   # writes silicon.nnkp for pw2wannier90
wannier90.jl silicon       # localise → silicon.wout, _hr.dat, bands, .chk
postw90.jl silicon         # AHC/DOS/kpath/… per the .win keywords
```

(shown with the `install_cli()` launchers on `PATH`; from a clone, use
`julia --project=. bin/wannier90.jl …` instead)

**As a library** (results are data, not a `.wout` to parse):

```julia
using WannierFunctions

model = read_model("silicon")                     # .win/.amn/.mmn/.eig
res   = run_wannier(model; win_max=17.0, froz_max=6.4)   # disentangle + localise
res.spread.Ω, res.spread.centres                  # 14.499574503 Å², 3×num_wann

H = hamiltonian_operator(model, res)              # interpolable H(R)
E = bands(H, [[0,0,0], [0.5,0,0.5]])              # energies at Γ and X
```

**All-Julia, no files** — a DFTK SCF handed straight to the wannieriser; for an isolated
manifold, `num_wann` is the *only* Wannier-specific input (SCDM automatic projections):

```julia
using DFTK, WannierFunctions
scfres = self_consistent_field(basis; tol=1e-10)  # symmetries=false, full k-grid
wmodel = wannier_model(scfres; num_wann=4)        # ← the entire Wannier specification
res    = wannierise(wmodel; algorithm=:w90)
```

More in [Getting started](docs/src/getting-started.md), and nine runnable
[examples](examples/README.md) from GaAs localisation up to **magic-angle twisted bilayer
graphene from first principles** (figures in [examples/output](examples/output)).

## Documentation

| Page | For |
|------|-----|
| [Getting started](docs/src/getting-started.md) | New users: both workflows end to end, reading results, choosing windows/projections, troubleshooting |
| [How-to guides](docs/src/howto.md) | Experienced users: recipes for every module — disentanglement schemes, postw90 tasks, transport, symmetry, TB input, DFTK |
| [Wannier90 compatibility](docs/src/wannier90-compat.md) | Binary/keyword mapping, conventions, and every known behavioural difference |
| [Validation](docs/src/validation.md) | The oracle methodology and the complete parity results |
| [Theory](docs/src/theory.md) | The mathematics: spread functional, b-vectors, gradients, disentanglement, interpolation |
| [File formats](docs/src/file-formats.md) | `.win/.amn/.mmn/.eig/_hr.dat/…` specifications |
| [Migrating from Wannier90](docs/src/migrating-from-wannier90.md) | Mapping an existing workflow |
| [Examples](examples/README.md) | Nine runnable scripts with reference numbers and plots |
| [docs/reference-notes/](docs/reference-notes/) | Implementation-grade notes per algorithm (conventions, oracle anchors, upstream quirks) |

Build the HTML docs locally with `julia --project=docs docs/make.jl` (works offline).

## Design notes for the impatient

- **Strict input validation**: an unknown `.win` keyword is an error with a did-you-mean
  suggestion, checked against a catalogue generated from the reference parser.
- **Two optimisers, same minima**: `:rcg` (native default, true convergence criterion) and
  `:w90` (reference-exact trajectory, used by the CLI). The test suite asserts they agree.
- **Everything is an operator**: `TBOperator` carries H(R), r(R), or any Wannier-gauge
  operator; the Berry-physics stack consumes them uniformly, whether they came from a
  checkpoint, a `_tb.dat`, or a live DFTK run.
- **Threaded** k/R loops (size-gated); start Julia with `-t auto` for dense interpolation.
- **Honesty first**: validation compares only gauge-invariant quantities; caveats
  (e.g. cross-code tolerances, upstream bugs found) are documented, not hidden — see
  [Validation](docs/src/validation.md).

## Citing

See [`CITATION.cff`](CITATION.cff). The physics references:

- N. Marzari and D. Vanderbilt, Phys. Rev. B **56**, 12847 (1997).
- I. Souza, N. Marzari, and D. Vanderbilt, Phys. Rev. B **65**, 035109 (2001).
- N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, and D. Vanderbilt,
  Rev. Mod. Phys. **84**, 1419 (2012).
- G. Pizzi *et al.*, *Wannier90 as a community code*, J. Phys.: Condens. Matter **32**,
  165902 (2020) — the reference implementation this package is validated against.

MIT license.
