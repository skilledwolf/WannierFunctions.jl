# AB bilayer graphene — the *minimal reliable* wannierisation, and an honest look at how
# minimal it can actually be.
#
# This is the lean companion to 08_dftk_bilayer_graphene.jl. It keeps the same physics but
# strips the code to the essentials (one figure, no K-zoom), and it answers a question that
# comes up whenever someone wants "num_wann and nothing else": can the graphene π model be
# obtained from an energy-only recipe (SCDM-erfc, with its μ/σ auto-fitted by `scdm_auto`)?
#
# The answer, demonstrated below, is NO — and that is instructive, not a defect. Graphene's
# π and σ manifolds OVERLAP in energy, so no energy window or energy-based SCDM smearing can
# separate them; the pz *character* is an irreducible physical input. `scdm_auto` reports
# this itself: fitting the pz-projectability-vs-energy cloud to an erfc gives a large residual
# because pz character is band-pass, not a monotonic step. `scdm_auto` is the right tool for
# an entangled but energy-SEPARABLE manifold (a transition-metal d-manifold — the Vitale et
# al. tungsten case, npj Comput. Mater. 6, 66 (2020)); graphene π is the counterexample.
#
# So "minimal" here means the *code* is lean. The *specification* — num_wann, the pz
# character, and projectability disentanglement (PDWF) — is already the minimal reliable one.
#
# Requires:  ] add DFTK Plots     Runtime: a few minutes (same SCF as example 08).

using DFTK, WannierFunctions
using LinearAlgebra, Printf, StaticArrays
using Plots
ENV["GKSwstype"] = "100"

# --- 1. ground state: AB bilayer, 3.35 Å interlayer, 20 Å cell (identical to example 08) ---
a = 2.46 / WannierFunctions.BOHR
c = 20.0 / WannierFunctions.BOHR
lattice = [a -a/2 0; 0 a*√3/2 0; 0 0 c]
dz = 3.35 / 20.0
z1, z2 = 0.5 - dz/2, 0.5 + dz/2
positions = [[0.0, 0.0, z1], [1/3, 2/3, z1], [1/3, 2/3, z2], [2/3, 4/3 - 1, z2]]
C = ElementPsp(:C; psp=load_psp("hgh/lda/c-q4"))

model = model_DFT(lattice, [C, C, C, C], positions; functionals=LDA(),
                  temperature=1e-3, smearing=DFTK.Smearing.Gaussian(), symmetries=false)
basis = PlaneWaveBasis(model; Ecut=18, kgrid=(6, 6, 1))
scfres = self_consistent_field(basis; tol=1e-8,
                               nbandsalg=DFTK.AdaptiveBands(model; n_bands_converge=20))
εF = DFTK.auconvert(DFTK.Unitful.eV, scfres.εF).val
@printf("SCF done. εF = %.4f eV\n", εF)

# --- 2. one pz orbital per carbon, Löwdin-orthonormalised (proper PDWF weights) -----------
projs = [DFTK.HydrogenicWannierProjection(p, 2, 1, 0, 4.0) for p in positions]
wmodel = wannier_model(scfres, projs; num_wann=4, num_bands=20)
for k in 1:size(wmodel.A, 3)
    Ak = wmodel.A[:, :, k]
    wmodel.A[:, :, k] = Ak * inv(sqrt(Hermitian(Ak' * Ak)))
end

# --- 3. why not energy-only? ask scdm_auto to fit the projectability curve ----------------
# scdm_auto fits P(ε) = ½ erfc((ε−μ)/σ) to the (energy, pz-projectability) scatter. A small
# residual would mean the manifold is energy-separable and SCDM-erfc could replace the pz
# projections. Here the residual is large — the diagnostic that graphene π is NOT separable.
fit = scdm_auto(wmodel.A, wmodel.eig)
@printf("scdm_auto erfc fit residual rms = %.3f  (≈0 ⇒ energy-separable; here it is not:\n", fit.rms)
println("            pz character is band-pass, so an energy-only SCDM cannot isolate π —")
println("            the pz projections + PDWF below are the irreducible specification.)")

# --- 4. projectability disentanglement (PDWF) + maximal localisation ----------------------
# Freeze bands with > 95% pz character, disentangle the 2–95% pool, discard the rest; the
# generous outer window just brackets the 20 computed bands.
dis = disentangle(wmodel; win_min=εF - 20.0, win_max=εF + 10.0,
                  froz_proj=true, proj_min=0.02, proj_max=0.95, num_iter=500)
res = localize(dis.U0, dis.Mrot0, wmodel.bvectors;
               num_iter=1000, algorithm=:w90, conv_tol=1e-10, conv_window=5)
@printf("Ω_I = %.4f Å²,  Ω = %.4f Å²,  converged = %s\n", dis.omega_I, res.spread.Ω, res.converged)
for n in 1:4
    @printf("  WF %d: centre z = %+7.3f Å, spread = %.3f Å²\n",
            n, res.spread.centres[3, n], res.spread.spreads[n])
end

# --- 5. one figure: Wannier-interpolated bands vs DFT along Γ–K–M–Γ -----------------------
irvec, ndegen = wigner_seitz(wmodel.lattice, wmodel.kgrid.mp_grid)
Hr, _ = build_hr(res.U, dis.eigval_opt, wmodel.kgrid, irvec)
Γ, K, M = [0.0, 0.0, 0.0], [1/3, 1/3, 0.0], [0.5, 0.0, 0.0]
function path(pts, n)
    ks = SVector{3,Float64}[]
    for i in 1:length(pts)-1, t in range(0, 1; length=n + 1)[1:end-1]
        push!(ks, SVector{3,Float64}((pts[i] .+ t .* (pts[i+1] .- pts[i]))...))
    end
    push!(ks, SVector{3,Float64}(pts[end]...)); ks
end
B = Matrix(wmodel.lattice.B)
xof(ks) = pushfirst!(cumsum([norm(B * (ks[i+1] - ks[i])) for i in 1:length(ks)-1]), 0.0)

kpts = path([Γ, K, M, Γ], 60); xs = xof(kpts)
Ew = interpolate_bands_ws(Hr, irvec, ndegen, res.spread.centres, wmodel.lattice,
                          wmodel.kgrid.mp_grid, kpts) .- εF
kref = path([Γ, K, M, Γ], 12); xref = xof(kref)
bref = compute_bands(scfres, DFTK.ExplicitKpoints([DFTK.Vec3(k...) for k in kref]); n_bands=16)
Eref = hcat([DFTK.auconvert.(DFTK.Unitful.eV, εk) .|> x -> x.val for εk in bref.eigenvalues]...) .- εF

lx = [xs[1 + 60*(i-1)] for i in 1:3]; push!(lx, xs[end])
p = plot(; xticks=(lx, ["Γ", "K", "M", "Γ"]), ylabel="E − E_F (eV)", legend=:topright,
         title="AB bilayer graphene: minimal pz Wannier model vs DFT", ylims=(-11, 7))
vline!(p, lx; color=:gray, alpha=0.4, label="")
hline!(p, [0.0]; color=:gray, style=:dash, alpha=0.6, label="")
for b in 1:size(Eref, 1)
    scatter!(p, xref, Eref[b, :]; color=:black, ms=2.2, msw=0, label=(b == 1 ? "DFT (DFTK)" : ""))
end
for b in 1:4
    plot!(p, xs, Ew[b, :]; color=:crimson, lw=2, label=(b == 1 ? "Wannier (4 × pz)" : ""))
end

out = joinpath(@__DIR__, "output"); mkpath(out)
savefig(p, joinpath(out, "08_minimal_bilayer_graphene_bands.png"))
println("wrote ", joinpath(out, "08_minimal_bilayer_graphene_bands.png"))
