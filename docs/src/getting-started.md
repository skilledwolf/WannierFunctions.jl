# Getting started

This page walks through both ways of using the package — the drop-in file workflow and the
all-Julia in-memory workflow — and explains how to read the results and choose the few
parameters that actually matter.

## Installation

```julia
pkg> add https://github.com/skilledwolf/WannierFunctions.jl   # (not yet registered)
```

Julia ≥ 1.10. For the in-memory DFT workflow you also need `pkg> add DFTK`; the plotting
examples additionally use `Plots`.

## Workflow A: from Wannier90 input files

If a DFT interface (e.g. Quantum ESPRESSO's `pw2wannier90`) already produces
`seedname.amn/.mmn/.eig` for you, nothing upstream changes:

```bash
# 1. write the b-vector request for the DFT interface (replaces `wannier90.x -pp`)
julia --project=. bin/wannier90.jl -pp seedname

# 2. run pw2wannier90 as usual → seedname.amn/.mmn/.eig

# 3. localise + interpolate (replaces `wannier90.x`)
julia --project=. bin/wannier90.jl seedname

# 4. post-process (replaces `postw90.x`): AHC, DOS, kpath, BoltzWann, … per the .win
julia --project=. bin/postw90.jl seedname
```

The same run as a library call, with results as data instead of a `.wout` to parse:

```julia
using WannierFunctions

model = read_model("seedname")        # rich REPL display: bands, WFs, k-mesh, shells
res   = run_wannier(model)            # isolated manifold: just localise
res.spread.Ω                          # total spread (Å²)
res.spread.centres                    # 3 × num_wann Wannier centres (Cartesian Å)
res.converged, res.niter
```

For entangled bands (`num_bands > num_wann`) pass the windows as keywords —
`run_wannier` selects disentanglement automatically:

```julia
res = run_wannier(model; win_max = 17.0, froz_max = 6.4)
res.disentangled                      # true
res.omega_I                           # the gauge-invariant part Ω_I
res.dis.omega_I_trace                 # per-iteration convergence trace
```

Band interpolation from the result:

```julia
H = hamiltonian_operator(model, res)              # H(R) on the Wigner–Seitz set
E = bands(H, [[0,0,0], [0.5,0.5,0.5], [0.5,0,0.5]])   # fractional k-points
```

## Workflow B: all-Julia, no files (DFTK.jl)

With DFTK installed, the package extension hands an SCF result straight to the wannieriser.
Two requirements on the DFT side: run on the **full (symmetry-unreduced) Monkhorst–Pack
grid** (`symmetries = false`), and converge a few empty bands if you need them.

The minimal-input path uses **SCDM automatic projections** — for an isolated group of bands
(e.g. the four valence bands of silicon), `num_wann` is the only Wannier-specific input:

```julia
using DFTK, WannierFunctions

model  = model_DFT(lattice, atoms, positions; functionals=LDA(), symmetries=false)
basis  = PlaneWaveBasis(model; Ecut=14, kgrid=(4,4,4))
scfres = self_consistent_field(basis; tol=1e-10)

wmodel = wannier_model(scfres; num_wann=4)        # SCDM: no projections, no chemistry input
res    = wannierise(wmodel; algorithm=:w90, num_iter=500, conv_tol=1e-10, conv_window=5)
```

If you prefer explicit trial orbitals (they document the physics and are needed for
projectability-based disentanglement):

```julia
projs  = [DFTK.HydrogenicWannierProjection(center, 2, 1, 0, 4.0) for center in centers]  # 2p_z
wmodel = wannier_model(scfres, projs; num_wann=4, num_bands=20)
```

Runnable, validated versions of both:
[`examples/06`](../../examples/06_dftk_end_to_end.jl),
[`examples/07`](../../examples/07_dftk_scdm_minimal.jl), and the bilayer-graphene /
twisted-bilayer showcases [`examples/08`](../../examples/08_dftk_bilayer_graphene.jl) and
[`examples/09`](../../examples/09_tbg_local_stacking.jl).

## Reading the results

Everything a `.wout` reports is a field:

| Field | Meaning |
|-------|---------|
| `res.spread.Ω`, `.ΩI`, `.ΩOD`, `.ΩD` | the spread decomposition Ω = Ω_I + Ω_OD + Ω_D (Å²) |
| `res.spread.centres`, `.spreads` | per-WF centres (3 × n, Cartesian Å) and spreads |
| `res.converged`, `res.niter` | convergence status |
| `res.U` | the final gauge U(k) |
| `res.disentangled`, `res.omega_I`, `res.dis` | disentanglement results (when active) |

A good run looks like: `Ω_I` converged and *frozen* thereafter (it is gauge-invariant),
`Ω` decreasing monotonically to convergence, WF spreads of order 1–3 Å² for covalent
orbitals, and centres where chemistry says they should be (bond midpoints, atoms).

## The three decisions that matter

**1. Initial projections.** Ordered from most to least automatic:

- *SCDM* (`wannier_model(scfres; num_wann)` or `scdm_projections(model; dir)`): selected from
  the wavefunctions themselves by column-pivoted QR — zero chemistry input. For entangled
  manifolds add an energy window: `scdm_mode=:erfc, scdm_mu, scdm_sigma`.
- *Atomic orbitals* (a `projections` block, or `HydrogenicWannierProjection` with DFTK):
  explicit and self-documenting; enables projectability-based disentanglement. Note that
  atomic orbitals on different sites overlap — Löwdin-orthonormalise the projection columns
  before using projectabilities as weights (see `examples/08`).

**2. Disentanglement (only when `num_bands > num_wann`).** Ordered by how much you need to
know about the band structure:

- *Energy windows* (`win_min/win_max` outer, `froz_min/froz_max` frozen): the classic scheme;
  requires looking at the bands once.
- *PDWF projectability freezing* (`froz_proj=true, proj_min=0.02, proj_max=0.95`): freeze by
  orbital character instead of energy — robust when unwanted states (surface/vacuum states in
  slabs, semicore) intrude into any energy window.
- *SCDM-erfc*: the smearing function does the selecting; combine with a generous outer window.

**3. The optimiser.** `algorithm = :rcg` (default) has a true convergence criterion and
usually needs fewer iterations; `algorithm = :w90` reproduces the reference trajectory exactly
(what the CLI uses, and what you want when comparing gauge-dependent outputs against
`wannier90.x`). Both reach the same minima.

## Troubleshooting

| Symptom | Likely cause → fix |
|---------|--------------------|
| Huge spreads (≫ 5 Å²/WF), no convergence | Wrong states in the frozen manifold — in slabs an energy window catches vacuum states at Γ: use PDWF projectability freezing instead |
| `projectability ∉ [0,1]` error | Non-orthogonal atomic projectors (e.g. overlapping pz): Löwdin-orthonormalise `model.A` per k first |
| `gamma_only` errors about k-points | Γ-only mode needs exactly one k-point and the half b-vector `.mmn` convention |
| Metals: SCF or windows behave oddly | Use smearing in the DFT step; for the wannierisation give the frozen window a margin below E_F |
| Unknown `.win` keyword error | Deliberate — the parser is strict, and suggests the correct spelling; see [compatibility](wannier90-compat.md) |
| Interpolated bands wiggle between k-points | k-mesh too coarse for the WF spread; densify `mp_grid` (and keep `use_ws_distance` on — it is the default) |

Next: the [How-to guides](howto.md) for every module, or
[Migrating from Wannier90](migrating-from-wannier90.md) if you have an existing workflow.
