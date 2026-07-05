# The minimal-input Wannier workflow: DFTK + SCDM automatic projections.
#
# No `projections` block, no orbital guesses, no trial centres, no files, no external
# binaries. For an isolated band manifold (here: the four valence bands of silicon) the ONLY
# Wannier-specific input is `num_wann`. The initial gauge is selected from the wavefunctions
# themselves by SCDM (column-pivoted QR on the real-space orbitals — Damle, Lin & Ying,
# JCTC 11, 1463 (2015)), and the MV minimiser polishes it to the maximally-localised gauge.
#
# Requires DFTK in the active environment:  ] add DFTK

using DFTK
using WannierFunctions
using LinearAlgebra, Printf

# --- 1. any converged DFTK calculation on a full (symmetry-unreduced) MP grid -----------
a = 10.26                                     # Si lattice constant (Bohr)
lattice = a / 2 * [[0 1 1.0]; [1 0 1.0]; [1 1 0.0]]
Si = ElementPsp(:Si; psp=load_psp("hgh/lda/si-q4"))
model = model_DFT(lattice, [Si, Si], [ones(3) / 8, -ones(3) / 8];
                  functionals=LDA(), symmetries=false)
basis = PlaneWaveBasis(model; Ecut=14, kgrid=(4, 4, 4))
scfres = self_consistent_field(basis; tol=1e-10)

# --- 2. wannierise. This line is the entire Wannier-specific specification. -------------
wmodel = wannier_model(scfres; num_wann=4)

res = wannierise(wmodel; num_iter=500, algorithm=:w90, conv_tol=1e-10, conv_window=5)
@printf("Ω = %.6f Å²   converged = %s\n", res.spread.Ω, res.converged)

# SCDM found the Si–Si bond centres without being told any chemistry:
Binv = inv(Matrix(wmodel.lattice.A))
for n in 1:4
    println("  WF $n centre (frac): ", round.(Binv * res.spread.centres[:, n]; digits=4))
end

# --- 3. and the interpolation reproduces the ab-initio bands ----------------------------
irvec, ndegen = wigner_seitz(wmodel.lattice, wmodel.kgrid.mp_grid)
Hr, _ = build_hr(res.U, wmodel.eig, wmodel.kgrid, irvec)
maxdev = 0.0
for (ik, kf) in enumerate(wmodel.kgrid.frac)
    E = interpolate_bands(Hr, irvec, ndegen, [kf])[:, 1]
    global maxdev = max(maxdev, maximum(abs.(E .- wmodel.eig[:, ik])))
end
@printf("max |E_interp − E_scf| on the SCF mesh = %.2e eV\n", maxdev)

# Entangled manifolds (metals, valence+conduction) need two more numbers, not more
# chemistry: an energy window for the SCDM smearing —
#   wannier_model(scfres; num_wann, num_bands, scdm_mode=:erfc, scdm_mu=εF_eV, scdm_sigma=2.0)
# followed by the usual disentanglement windows in run_wannier/disentangle.
