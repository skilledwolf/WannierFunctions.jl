# Example 2 — Wannier interpolation of the diamond band structure.
#
# 4 isolated valence bands → 4 Wannier functions, then interpolate the bands along L–Γ–X–K–Γ.
# Mirrors Wannier90 example 5. Run from the repository root:
#
#   julia --project=. examples/02_diamond_interpolation.jl
#
using WannierFunctions
using Printf

seed = joinpath(@__DIR__, "data", "diamond")

# run_wannier reads the seed, runs the full pipeline, and returns a result ready to interpolate.
model, win, res = run_wannier(seed)
s = res.spread
@printf("diamond: Ω = %.9f Å²  (Ω_I = %.9f, Ω_OD = %.9f, Ω_D = %.9f)\n",
        s.Ω, s.ΩI, s.ΩOD, s.ΩD)
@printf("reference Ω = 2.320904915 Å²  →  difference = %.2e\n\n", abs(s.Ω - 2.320904915))

# Build the k-path from the .win `kpoint_path` block and interpolate.
kpts, xvals, labels, idx = generate_kpath(win, model.lattice; bands_num_points = 50)
bands = interpolate(model, res, kpts)                 # (num_wann × npoints), ascending per k
@printf("interpolated %d bands at %d k-points along: %s\n",
        size(bands, 1), length(kpts), join(labels, " – "))

# Write the standard band-structure files (band.dat/.kpt/.labelinfo), as wannier90.x would.
write_band_dat(seed * "_band.dat", xvals, bands)
write_band_kpt(seed * "_band.kpt", kpts)
write_labelinfo(seed * "_band.labelinfo.dat", labels, idx, xvals[idx], kpts[idx])
println("wrote diamond_band.dat / .kpt / .labelinfo.dat next to the input files")

# Show the valence-band top at Γ (the last k-point of the path).
@printf("\nband energies at Γ (eV): %s\n",
        join((@sprintf("%.4f", e) for e in sort(bands[:, end])), "  "))
