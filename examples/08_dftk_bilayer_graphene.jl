# AB (Bernal) bilayer graphene, all in Julia: DFTK ground state → four pz Wannier functions
# (Hydrogenic pz projections + disentanglement + maximal localisation, no files) → Wannier
# bands vs the ab-initio bands along Γ–K–M–Γ, with the parabolic band touching at K that is
# the hallmark of Bernal stacking.
#
# Requires:  ] add DFTK Plots
#
# Runtime: a few minutes (Ecut 18 Ha, 6×6×1 k-grid, ~20 Å vacuum).

using DFTK
using WannierFunctions
using LinearAlgebra, Printf, StaticArrays
using Plots
ENV["GKSwstype"] = "100"                      # headless GR

# --- 1. geometry: AB bilayer, interlayer 3.35 Å, 20 Å cell height ----------------------
a = 2.46 / WannierFunctions.BOHR              # in-plane lattice constant (Bohr)
c = 20.0 / WannierFunctions.BOHR
lattice = [a  -a/2       0;
           0  a*√3/2     0;
           0  0          c]
dz = (3.35 / 20.0)                            # interlayer separation (fractional)
z1, z2 = 0.5 - dz / 2, 0.5 + dz / 2
# layer 1: A₁ (0,0), B₁ (1/3,2/3); layer 2 shifted by the bond vector: A₂ above B₁ — Bernal
positions = [[0.0, 0.0, z1], [1/3, 2/3, z1],
             [1/3, 2/3, z2], [2/3, 4/3 - 1, z2]]
C = ElementPsp(:C; psp=load_psp("hgh/lda/c-q4"))
atoms = [C, C, C, C]

model = model_DFT(lattice, atoms, positions; functionals=LDA(),
                  temperature=1e-3, smearing=DFTK.Smearing.Gaussian(),
                  symmetries=false)
basis = PlaneWaveBasis(model; Ecut=18, kgrid=(6, 6, 1))
scfres = self_consistent_field(basis; tol=1e-8,
                               nbandsalg=DFTK.AdaptiveBands(model; n_bands_converge=20))
εF = DFTK.auconvert(DFTK.Unitful.eV, scfres.εF).val
@printf("SCF done. εF = %.4f eV\n", εF)

# --- 2. wannierise: four pz orbitals (one per carbon) ----------------------------------
projs = [DFTK.HydrogenicWannierProjection(p, 2, 1, 0, 4.0) for p in positions]   # 2p_z
wmodel = wannier_model(scfres, projs; num_wann=4, num_bands=20)

# Löwdin-orthonormalise the atomic projectors in the band space ("ortho-atomic", as in
# pw2wannier90): neighbouring pz orbitals overlap (the π bond!), so raw projectabilities can
# exceed 1; orthonormal columns make them proper weights ∈ [0, 1] for the PDWF freeze below.
for k in 1:size(wmodel.A, 3)
    Ak = wmodel.A[:, :, k]
    wmodel.A[:, :, k] = Ak * inv(sqrt(Hermitian(Ak' * Ak)))
end

# π/π* bands entangle with σ states below and nearly-free vacuum states of the slab above —
# an energy-window freeze would catch the vacuum states at Γ. Select and freeze by pz
# PROJECTABILITY instead (the PDWF scheme, Qiao–Pizzi–Marzari): bands with > 95% pz character
# are frozen wherever they sit in energy, bands with 2–95% enter the disentanglement pool,
# and the rest are discarded.
dis = disentangle(wmodel; win_min=εF - 20.0, win_max=εF + 10.0,
                  froz_proj=true, proj_min=0.02, proj_max=0.95, num_iter=500)
res = localize(dis.U0, dis.Mrot0, wmodel.bvectors;
               num_iter=1000, algorithm=:w90, conv_tol=1e-10, conv_window=5)
@printf("Ω_I = %.4f Å²,  Ω = %.4f Å²,  converged = %s\n",
        dis.omega_I, res.spread.Ω, res.converged)
sr = res.spread
for n in 1:4
    @printf("  WF %d: centre z = %7.3f Å, spread = %.3f Å²\n", n, sr.centres[3, n], sr.spreads[n])
end

# --- 3. Wannier-interpolated bands vs the ab-initio bands ------------------------------
irvec, ndegen = wigner_seitz(wmodel.lattice, wmodel.kgrid.mp_grid)
Hr, _ = build_hr(res.U, dis.eigval_opt, wmodel.kgrid, irvec)

Γ, K, M = [0.0, 0.0, 0.0], [1/3, 1/3, 0.0], [0.5, 0.0, 0.0]
function path(pts, n)
    ks = SVector{3,Float64}[]
    for i in 1:length(pts)-1
        for t in range(0, 1; length=n + 1)[1:end-1]
            push!(ks, SVector{3,Float64}((pts[i] .+ t .* (pts[i+1] .- pts[i]))...))
        end
    end
    push!(ks, SVector{3,Float64}(pts[end]...))
    ks
end
kpts = path([Γ, K, M, Γ], 60)
B = Matrix(wmodel.lattice.B)
xs = pushfirst!(cumsum([norm(B * (kpts[i+1] - kpts[i])) for i in 1:length(kpts)-1]), 0.0)
Ew = interpolate_bands_ws(Hr, irvec, ndegen, res.spread.centres, wmodel.lattice,
                          wmodel.kgrid.mp_grid, kpts) .- εF

# ab-initio reference on a sparser version of the same path
kref = path([Γ, K, M, Γ], 12)
bands_ref = compute_bands(scfres, DFTK.ExplicitKpoints([DFTK.Vec3(k...) for k in kref]);
                          n_bands=16)
Eref = hcat([DFTK.auconvert.(DFTK.Unitful.eV, εk) .|> x -> x.val
             for εk in bands_ref.eigenvalues]...) .- εF
xref = pushfirst!(cumsum([norm(B * (kref[i+1] - kref[i])) for i in 1:length(kref)-1]), 0.0)

lbl = ["Γ", "K", "M", "Γ"]
lx = [xs[1 + 60*(i-1)] for i in 1:3]
push!(lx, xs[end])
p1 = plot(; xticks=(lx, lbl), ylabel="E − E_F (eV)", legend=:topright,
          title="AB bilayer graphene: pz Wannier model vs DFT", ylims=(-11, 7))
vline!(p1, lx; color=:gray, alpha=0.4, label="")
hline!(p1, [0.0]; color=:gray, style=:dash, alpha=0.6, label="")
for b in 1:size(Eref, 1)
    scatter!(p1, xref, Eref[b, :]; color=:black, ms=2.2, msw=0,
             label=(b == 1 ? "DFT (DFTK)" : ""))
end
for b in 1:4
    plot!(p1, xs, Ew[b, :]; color=:crimson, lw=2, label=(b == 1 ? "Wannier (4 × pz)" : ""))
end

# zoom at K: the quadratic band touching of Bernal stacking (massive chiral fermions)
kz = path([[1/3 - 0.05, 1/3 - 0.05, 0.0], K, [1/3 + 0.035, 1/3 + 0.07, 0.0]], 80)
xz = pushfirst!(cumsum([norm(B * (kz[i+1] - kz[i])) for i in 1:length(kz)-1]), 0.0)
Ez = interpolate_bands_ws(Hr, irvec, ndegen, res.spread.centres, wmodel.lattice,
                          wmodel.kgrid.mp_grid, kz) .- εF
p2 = plot(; xlabel="k around K (Å⁻¹)", ylabel="E − E_F (eV)", legend=false,
          title="Quadratic touching at K", ylims=(-1.2, 1.2))
for b in 1:4
    plot!(p2, xz .- xz[findmin(abs.(Ez[2, :] .- Ez[3, :]))[2]], Ez[b, :]; color=:crimson, lw=2)
end
hline!(p2, [0.0]; color=:gray, style=:dash, alpha=0.6)

out = joinpath(@__DIR__, "output")
mkpath(out)
savefig(plot(p1, p2; layout=(1, 2), size=(1100, 420), margin=5Plots.mm),
        joinpath(out, "08_bilayer_graphene_bands.png"))
println("wrote ", joinpath(out, "08_bilayer_graphene_bands.png"))
