# Example 1 — Maximally-localised Wannier functions for GaAs valence bands.
#
# The simplest case: 4 isolated valence bands → 4 Wannier functions, no disentanglement.
# Mirrors Wannier90 tutorial/example 1. Run from the repository root:
#
#   julia --project=. examples/01_gaas_localization.jl
#
using Wannier90
using Printf

seed = joinpath(@__DIR__, "data", "gaas")

# read_model parses seedname.win/.amn/.mmn (GaAs ships no .eig — pure localisation).
model = read_model(seed)
@printf("GaAs: %d bands → %d Wannier functions, %d k-points (%d×%d×%d)\n",
        model.num_bands, model.num_wann, length(model.kgrid.frac), model.kgrid.mp_grid...)

# Minimise the Marzari–Vanderbilt spread (20 iterations, as in the reference input).
res = wannierise(model; num_iter = 20)
s = res.spread

println("\nWannier centres (Å) and spreads (Å²):")
for n in 1:model.num_wann
    @printf("  WF %d  (% .6f, % .6f, % .6f)   %.6f\n",
            n, s.centres[1, n], s.centres[2, n], s.centres[3, n], s.spreads[n])
end
@printf("\nΩ = %.9f Å²   (Ω_I = %.9f,  Ω_OD = %.9f,  Ω_D = %.9f)\n", s.Ω, s.ΩI, s.ΩOD, s.ΩD)

# The reference Wannier90 benchmark for this case:
@printf("reference Ω = 4.466880976 Å²  →  difference = %.2e\n", abs(s.Ω - 4.466880976))
