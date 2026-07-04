# Examples

Runnable examples for `Wannier90.jl`. Run each from the repository root:

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

## Command-line equivalent

Anything the scripts do, the drop-in CLI does from the shell — e.g. the diamond example is
equivalent to:

```bash
julia --project=. bin/wannier90.jl examples/data/diamond
```

which writes `diamond.wout` and (because `diamond.win` sets `bands_plot`) the band files, exactly
as `wannier90.x diamond` would.
