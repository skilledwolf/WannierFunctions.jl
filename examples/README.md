# Examples

Runnable examples for `WannierFunctions.jl`. Run each from the repository root:

```bash
julia --project=. examples/01_gaas_localization.jl
julia --project=. examples/02_diamond_interpolation.jl
julia --project=. examples/03_silicon_disentanglement.jl
```

| Script | System | What it shows |
|--------|--------|---------------|
| `01_gaas_localization.jl` | GaAs, 4 → 4 WF | Maximally-localised Wannier functions (no disentanglement); centres and spread. |
| `02_diamond_interpolation.jl` | diamond, 4 → 4 WF | Localisation **and** band-structure interpolation along L–Γ–X–K–Γ; writes `_band.dat`. |
| `03_silicon_disentanglement.jl` | silicon, 12 → 8 WF | **Disentanglement** with an outer + frozen energy window; the Ω_I convergence trace. |

The GaAs and diamond inputs (`.win/.amn/.mmn/.eig`) are shipped under `data/`. Silicon's overlap
file is ~2.7 MB and is not shipped; example 3 stages it from the reference Wannier90 tree if you
have it under `reference/wannier90` (otherwise it prints instructions).

Each script prints its result next to the corresponding reference Wannier90 benchmark number, so
you can see the agreement directly.

| `05_berry_ahc.jl` | bcc Fe, 18 WF (SOC) | **Berry curvature + anomalous Hall conductivity** from a finished run's checkpoint; reproduces `postw90.x` digit-for-digit. |
| `06_dftk_end_to_end.jl` | silicon, 4 → 4 WF | **All-Julia DFT → Wannier** pipeline: a DFTK SCF handed to the wannieriser in memory (no files, no binaries); needs `] add DFTK`. |
| `07_dftk_scdm_minimal.jl` | silicon, 4 → 4 WF | The **minimal-input workflow**: DFTK + **SCDM** automatic projections — `num_wann` is the only Wannier-specific input; needs `] add DFTK`. |
| `08_dftk_bilayer_graphene.jl` | AB bilayer graphene, 20 → 4 WF | DFTK slab → four **pz Wannier functions** via ortho-atomic projections + **PDWF projectability freezing** (vacuum states never intrude); band plot with the Bernal quadratic touching and γ₁ splitting at K; needs `] add DFTK Plots`. |
| `08_dftk_bilayer_graphene_minimal.jl` | AB bilayer graphene, 20 → 4 WF | The **lean** companion: same pz + PDWF model in far less code, one figure. Uses **`scdm_auto`** to show *why* the graphene π model cannot be reduced to an energy-only recipe — the projectability-vs-energy fit residual is large because π/σ overlap in energy (an instructive boundary of SCDM-erfc). |
| `09_tbg_local_stacking.jl` | twisted bilayer graphene | **Magic-angle physics from first principles** via the local-stacking approximation: 9 DFTK stackings → interlayer Dirac coupling T(d) → Bistritzer–MacDonald w₀/w₁ (C₃-exact trio) → moiré flat bands and the magic-angle dip; needs `] add DFTK Plots`, ~15 min (cached). |

## Command-line equivalent

Anything the scripts do, the drop-in CLI does from the shell — e.g. the diamond example is
equivalent to:

```bash
julia --project=. bin/wannier90.jl examples/data/diamond
```

which writes `diamond.wout` and (because `diamond.win` sets `bands_plot`) the band files, exactly
as `wannier90.x diamond` would.
