# All-Julia DFT → Wannier pipeline: a DFTK.jl silicon LDA calculation handed to the
# wannieriser entirely in memory — no .amn/.mmn/.eig files, no external binaries.
#
# Requires DFTK in the active environment:  ] add DFTK
#
# The package extension (ext/WannierFunctionsDFTKExt.jl) provides
# `wannier_model(scfres, projections; num_wann)`: the b-vector neighbour list comes from
# WannierFunctions' own kmesh search and the ⟨u|u⟩ / ⟨ψ|g⟩ matrix elements from DFTK's
# plane-wave routines. Everything downstream (spread, interpolation, Berry-phase physics)
# is then available on the resulting Model.

using DFTK
using WannierFunctions
using LinearAlgebra, Printf

# --- 1. silicon LDA ground state on the full (unreduced) 4×4×4 grid -------------------
a = 10.26                                     # lattice constant (Bohr)
lattice = a / 2 * [[0 1 1.0]; [1 0 1.0]; [1 1 0.0]]
Si = ElementPsp(:Si; psp=load_psp("hgh/lda/si-q4"))
atoms = [Si, Si]
positions = [ones(3) / 8, -ones(3) / 8]

model = model_DFT(lattice, atoms, positions; functionals=LDA(), symmetries=false)
basis = PlaneWaveBasis(model; Ecut=14, kgrid=(4, 4, 4))
scfres = self_consistent_field(basis; tol=1e-10)

# --- 2. hand off to the wannieriser in memory -----------------------------------------
# Trial projections: Gaussians near the four bonds (the minimiser refines the centres).
centers = [[1, 1, 1] / 8, [-3, 1, 1] / 8, [1, -3, 1] / 8, [1, 1, -3] / 8]
projs = [DFTK.GaussianWannierProjection(c) for c in centers]

wmodel = wannier_model(scfres, projs; num_wann=4)
res = wannierise(wmodel; num_iter=500, algorithm=:w90, conv_tol=1e-10, conv_window=5)
@printf("Ω = %.6f Å²   converged = %s\n", res.spread.Ω, res.converged)

# The four maximally-localised WFs sit on the Si–Si bond centres:
Binv = inv(Matrix(wmodel.lattice.A))
for n in 1:4
    println("  WF $n centre (frac): ", round.(Binv * res.spread.centres[:, n]; digits=4))
end

# --- 3. Wannier interpolation reproduces the ab-initio bands --------------------------
irvec, ndegen = wigner_seitz(wmodel.lattice, wmodel.kgrid.mp_grid)
Hr, _ = build_hr(res.U, wmodel.eig, wmodel.kgrid, irvec)
maxdev = 0.0
for (ik, kf) in enumerate(wmodel.kgrid.frac)
    E = interpolate_bands(Hr, irvec, ndegen, [kf])[:, 1]
    global maxdev = max(maxdev, maximum(abs.(E .- wmodel.eig[:, ik])))
end
@printf("max |E_interp − E_scf| on the SCF mesh = %.2e eV\n", maxdev)
