# Example 3 — Disentanglement: 12 bands → 8 Wannier functions for silicon (valence + conduction).
#
# Silicon's sp³ + antibonding manifold is entangled with higher conduction bands, so an outer
# energy window (dis_win_max) and a frozen inner window (dis_froz_max) are used to extract an
# 8-dimensional subspace before localisation. Mirrors Wannier90 example 3.
#
# The silicon overlap file (.mmn, ~2.7 MB) is not shipped with this repo; this script stages the
# inputs from the reference test suite if it is present. Run from the repository root:
#
#   julia --project=. examples/03_silicon_disentanglement.jl
#
using Wannier90
using Printf

# Locate/stage the silicon inputs.
datadir = joinpath(@__DIR__, "data")
seed = joinpath(datadir, "silicon")
ref = joinpath(@__DIR__, "..", "reference", "wannier90")
if !isfile(seed * ".mmn")
    tests = joinpath(ref, "test-suite", "tests", "testw90_example03")
    mmnbz = joinpath(ref, "test-suite", "checkpoints", "si_geninterp", "silicon.mmn.bz2")
    if isfile(joinpath(tests, "silicon.win")) && isfile(mmnbz)
        @info "staging silicon inputs from the reference test suite"
        for f in ("silicon.win", "silicon.amn", "silicon.eig")
            cp(joinpath(tests, f), joinpath(datadir, f); force = true)
        end
        run(pipeline(`bunzip2 -kc $mmnbz`, stdout = seed * ".mmn"))
    else
        error("""
              silicon inputs not found. Provide silicon.{win,amn,mmn,eig} in $(datadir),
              or clone the reference Wannier90 tree under reference/wannier90 so this script
              can stage them (silicon.mmn is ~2.7 MB and is not shipped here).""")
    end
end

model, win, res = run_wannier(seed)
@printf("silicon: %d bands → %d WF   (dis_win_max=%.1f eV, dis_froz_max=%.1f eV)\n",
        model.num_bands, model.num_wann, win.dis_win_max, win.dis_froz_max)
@printf("disentangled = %s\n", res.disentangled)

# The disentanglement Ω_I convergence trace (matches the reference iteration by iteration).
println("\nΩ_I convergence (first 6 iterations):")
for (it, wi1, wi, d) in res.dis.omega_I_trace[1:min(6, end)]
    @printf("  iter %2d:  Ω_I = %.8f   (Δ = %.2e)\n", it, wi, d)
end

s = res.spread
@printf("\nFinal spread: Ω = %.9f Å²  (Ω_I = %.9f, Ω_OD = %.9f, Ω_D = %.9f)\n",
        s.Ω, s.ΩI, s.ΩOD, s.ΩD)
@printf("reference Ω = 14.499574503 Å²  →  difference = %.2e\n", abs(s.Ω - 14.499574503))
