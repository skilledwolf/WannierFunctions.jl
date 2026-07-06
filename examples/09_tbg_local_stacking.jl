# Twisted bilayer graphene from first principles, via the LOCAL-STACKING approximation
# (Jung, Raoux, Qiao & MacDonald, PRB 89, 205414; Bistritzer & MacDonald, PNAS 108, 12233):
#
#   At small twist angles the moiré pattern looks locally like an UNTWISTED bilayer with a
#   slowly varying stacking shift d(r). So:
#     A. run DFTK on the untwisted bilayer at a 3×3 grid of stacking shifts d (with the
#        first-harmonic corrugation of the interlayer distance, AA 3.60 Å ↔ AB 3.35 Å),
#     B. build the four-orbital projected-pz Wannier Hamiltonian H(k; d) for each shift
#        (Löwdin projection — deterministic, gauge-smooth in d),
#     C. read off the interlayer Dirac-point coupling T(d) = H(K)[layer1, layer2] and
#        Fourier-decompose it over the stacking cell: the C₃-degenerate dominant trio of
#        harmonics gives the Bistritzer–MacDonald tunnelings w₀ (AA) and w₁ (AB),
#     D. solve the BM continuum model with these first-principles parameters: moiré flat
#        bands at θ ≈ 1.05° and the flat-band width vs angle (the magic angle).
#
# Requires:  ] add DFTK Plots      Runtime: ~10–15 minutes (9 small SCFs).

using DFTK
using WannierFunctions
using LinearAlgebra, Printf, StaticArrays
using Plots
ENV["GKSwstype"] = "100"

const A0 = 2.46                                # graphene lattice constant (Å)
const HBAR_VF_GUESS = 5.25                     # only used to seed axis ranges

# --- A. DFT sweep over stacking shifts --------------------------------------------------
a = A0 / WannierFunctions.BOHR
c = 20.0 / WannierFunctions.BOHR
lattice = [a  -a/2     0;
           0  a*√3/2   0;
           0  0        c]
C = ElementPsp(:C; psp=load_psp("hgh/lda/c-q4"))

# corrugation: d_z(d) = c0 + c1 Σ_n cos(G_n·d), calibrated to AA = 3.60 Å, AB = 3.35 Å
dz_of(d1, d2) = begin
    # first-star stacking harmonics of the 120° lattice convention: b₁, b₂, b₁−b₂
    s = cos(2π * d1) + cos(2π * d2) + cos(2π * (d1 - d2))
    c1 = (3.60 - 3.35) / 4.5
    c0 = 3.60 - 3 * c1
    c0 + c1 * s
end

"Projected-pz Bloch Hamiltonian H_W(k) on the SCF grid for stacking shift (d1, d2)."
function pz_model(d1, d2)
    dz = dz_of(d1, d2) / 20.0                  # fractional
    z1, z2 = 0.5 - dz / 2, 0.5 + dz / 2
    positions = [[0.0, 0.0, z1], [1/3, 2/3, z1],
                 [0.0 + d1, 0.0 + d2, z2], [1/3 + d1, 2/3 + d2, z2]]
    model = model_DFT(lattice, [C, C, C, C], positions; functionals=LDA(),
                      temperature=1e-3, smearing=DFTK.Smearing.Gaussian(),
                      symmetries=false)
    basis = PlaneWaveBasis(model; Ecut=16, kgrid=(6, 6, 1))
    scfres = self_consistent_field(basis; tol=1e-7,
                                   nbandsalg=DFTK.AdaptiveBands(model; n_bands_converge=12))
    projs = [DFTK.HydrogenicWannierProjection(p, 2, 1, 0, 4.0) for p in positions]
    wm = wannier_model(scfres, projs; num_wann=4, num_bands=12)
    εF = DFTK.auconvert(DFTK.Unitful.eV, scfres.εF).val
    nk = length(wm.kgrid.frac)
    Hk = Array{ComplexF64,3}(undef, 4, 4, nk)
    for k in 1:nk
        A = wm.A[:, :, k]
        V = A * inv(sqrt(Hermitian(A' * A)))   # Löwdin isometry: deterministic pz subspace
        Hk[:, :, k] = V' * Diagonal(wm.eig[:, k] .- εF) * V
    end
    return wm, Hk
end

shifts = [(i / 3, j / 3) for i in 0:2 for j in 0:2]
cache = joinpath(@__DIR__, "output", "09_tbg_dft_cache.jls")
using Serialization
if isfile(cache)
    println("loading cached DFT sweep from $cache (delete it to recompute)")
    T, kfrac0, Hk0 = deserialize(cache)
else
    T = Dict{Tuple{Float64,Float64},Matrix{ComplexF64}}()
    kfrac0 = nothing
    Hk0 = nothing
    for (n, (d1, d2)) in enumerate(shifts)
        @printf("[%d/9] stacking d = (%.3f, %.3f), d_z = %.3f Å … ", n, d1, d2, dz_of(d1, d2))
        wm, Hk = pz_model(d1, d2)
        (d1, d2) == (0.0, 0.0) && (global kfrac0 = wm.kgrid.frac; global Hk0 = Hk)
        iK = findfirst(k -> maximum(abs.(mod.(k .- [1/3, 1/3, 0.0] .+ 0.5, 1.0) .- 0.5)) < 1e-8,
                       wm.kgrid.frac)
        T[(d1, d2)] = Hk[1:2, 3:4, iK]         # interlayer block at the Dirac momentum
        @printf("|T| = %.0f meV\n", 1e3 * sqrt(sum(abs2, T[(d1, d2)]) / 4))
    end
    mkpath(dirname(cache))
    serialize(cache, (T, kfrac0, Hk0))
end

# --- B/C. BM parameters from the stacking Fourier transform ----------------------------
# In the wannier90 Bloch convention the interlayer block carries a global e^{∓iK·d} phase
# (the layer-2 orbitals ride with the shifted atoms while the Bloch sums use lattice vectors
# only). Unwind it, then the two-centre expansion T(d) e^{±iK·d} = Σ_G c_G e^{−iG·d} is
# dominated by the C₃-related trio of smallest |K+G|: m ∈ {(0,0), (-1,0), (0,-1)}. The sign
# convention is fixed empirically by whichever choice makes the trio C₃-degenerate.
trio = ((0, 0), (-1, 0), (0, -1))
function harmonics(sgn)
    h(m1, m2, α, β) = sum(T[(d1, d2)][α, β] * cis(2π * (sgn * (d1 + d2) / 3 + m1 * d1 + m2 * d2))
                          for (d1, d2) in shifts) / 9
    cAA = [h(m..., 1, 1) for m in trio]
    cAB = [h(m..., 1, 2) for m in trio]
    w0s = abs.(cAA)
    w1s = abs.(cAB)
    spread = (maximum(w0s) - minimum(w0s)) + (maximum(w1s) - minimum(w1s))
    leak = maximum(abs(h(m1, m2, 1, 2)) for m1 in -1:1 for m2 in -1:1
                   if !((m1, m2) in trio))
    (; cAA, cAB, w0s, w1s, spread, leak, sgn)
end
hp, hm = harmonics(+1), harmonics(-1)
hbest = hp.spread <= hm.spread ? hp : hm
@printf("\nphase convention: e^{%+diK·d}  (trio spreads %+d: %.1f meV, %+d: %.1f meV)\n",
        hp.spread <= hm.spread ? 1 : -1, 1, 1e3 * hp.spread, -1, 1e3 * hm.spread)
w0 = sum(hbest.w0s) / 3
w1 = sum(hbest.w1s) / 3
@printf("BM tunnelings from DFT:  w0 = %.1f meV  (trio spread %.1f)   w1 = %.1f meV  (trio spread %.1f)\n",
        1e3 * w0, 1e3 * (maximum(hbest.w0s) - minimum(hbest.w0s)),
        1e3 * w1, 1e3 * (maximum(hbest.w1s) - minimum(hbest.w1s)))
@printf("largest non-trio harmonic: %.1f meV  (should be small)\n", 1e3 * hbest.leak)

# Fermi velocity from the layer-1 intralayer Dirac cone, evaluated directly on the SCF grid
# (K and its nearest grid neighbour — no interpolation, so the σ/vacuum entanglement far
# from K cannot pollute the estimate).
Bm = 2π * inv(Matrix(lattice * WannierFunctions.BOHR))'   # reciprocal lattice (Å⁻¹)
onmesh(target) = findfirst(k -> maximum(abs.(mod.(k .- target .+ 0.5, 1.0) .- 0.5)) < 1e-8,
                           kfrac0)                        # DFTK stores k in (−1/2, 1/2]
iK = onmesh([1/3, 1/3, 0.0])
iK1 = onmesh([1/3 + 1/6, 1/3, 0.0])
EK = eigvals(Hermitian(Hk0[1:2, 1:2, iK]))
EK1 = eigvals(Hermitian(Hk0[1:2, 1:2, iK1]))
dk = norm(Bm * [1/6, 0.0, 0.0])
ħvF = ((EK1[2] - EK1[1]) - (EK[2] - EK[1])) / (2 * dk)    # cone splitting slope (eV·Å)
@printf("ħ v_F = %.2f eV·Å   (α = w1/(ħv_F kθ) sets the magic angle)\n", ħvF)

# --- D. Bistritzer–MacDonald continuum model --------------------------------------------
"BM moiré bands at twist θ (rad) along a path in the moiré BZ; returns (x, E)."
function bm_bands(θ, w0, w1, ħvF; nshell=4, npath=40)
    kD = 4π / (3 * A0)
    kθ = 2 * kD * sin(θ / 2)
    q = [kθ .* (sin(2π * (j - 1) / 3), -cos(2π * (j - 1) / 3)) for j in 1:3]
    Tj = [w0 * [1 0; 0 1] .+ w1 * [0 cis(-2π * (j - 1) / 3); cis(2π * (j - 1) / 3) 0]
          for j in 1:3]
    # bipartite momentum lattice (layer 1 / layer 2), BFS from the two Dirac points
    sites = Dict{NTuple{3,Float64},Int}()      # (Qx, Qy, layer)
    key(Q, l) = (round(Q[1] / kθ; digits=6), round(Q[2] / kθ; digits=6), l)
    frontier = [((0.0, 0.0), 1.0)]
    sites[key((0.0, 0.0), 1.0)] = 1
    Qs = [((0.0, 0.0), 1)]
    while !isempty(frontier)
        nf = []
        for (Q, l) in frontier, j in 1:3
            Q2 = l == 1 ? (Q[1] + q[j][1], Q[2] + q[j][2]) : (Q[1] - q[j][1], Q[2] - q[j][2])
            l2 = 3.0 - l
            hypot(Q2...) > nshell * kθ && continue
            haskey(sites, key(Q2, l2)) && continue
            sites[key(Q2, l2)] = length(Qs) + 1
            push!(Qs, (Q2, Int(l2)))
            push!(nf, (Q2, l2))
        end
        frontier = nf
    end
    nQ = length(Qs)
    # moiré path Km → Γm → Mm → K'm: the two BZ corners are the layer Dirac momenta
    # (Km = 0 for layer 1, K'm = q1 for layer 2); Γm completes the equilateral hexagon
    # centre and Mm is the midpoint of the Km–K'm edge.
    Km = (0.0, 0.0)
    Kp = q[1]
    Γm = (√3 * kθ / 2, -kθ / 2)
    Mm = (0.0, -kθ / 2)
    nodes = [Km, Γm, Mm, Kp]
    ks = Tuple{Float64,Float64}[]
    for i in 1:3, t in range(0, 1; length=npath + 1)[1:end-1]
        push!(ks, nodes[i] .+ t .* (nodes[i+1] .- nodes[i]))
    end
    push!(ks, nodes[end])
    xs = pushfirst!(cumsum([hypot((ks[i+1] .- ks[i])...) for i in 1:length(ks)-1]), 0.0)
    E = zeros(2nQ, length(ks))
    H = zeros(ComplexF64, 2nQ, 2nQ)
    for (ik, k) in enumerate(ks)
        H .= 0
        for (i, (Q, l)) in enumerate(Qs)
            κ = (k[1] - Q[1], k[2] - Q[2])
            H[2i-1:2i, 2i-1:2i] = ħvF * (κ[1] * [0 1; 1 0] + κ[2] * [0 -im; im 0])
            l == 1 || continue
            for j in 1:3
                Q2 = (Q[1] + q[j][1], Q[2] + q[j][2])
                i2 = get(sites, key(Q2, 2.0), 0)
                i2 == 0 && continue
                H[2i-1:2i, 2i2-1:2i2] = Tj[j]'
                H[2i2-1:2i2, 2i-1:2i] = Tj[j]
            end
        end
        E[:, ik] = eigvals(Hermitian(H))
    end
    return xs, E, nQ
end

θdeg = 1.05
xs, E, nQ = bm_bands(θdeg * π / 180, w0, w1, ħvF)
mid = nQ                                        # the two flat bands are E[nQ] and E[nQ+1]
p1 = plot(; ylabel="E (meV)", legend=false, ylims=(-160, 160),
          title=@sprintf("TBG moiré bands, θ = %.2f° (w0 = %.0f, w1 = %.0f meV)",
                         θdeg, 1e3 * w0, 1e3 * w1),
          xticks=(xs[[1, 41, 81, 121]], ["Kₘ", "Γₘ", "Mₘ", "K'ₘ"]))
for b in max(1, mid - 4):min(2nQ, mid + 5)
    plot!(p1, xs, 1e3 .* E[b, :]; color=:gray, lw=1)
end
plot!(p1, xs, 1e3 .* E[mid, :]; color=:crimson, lw=2.5)
plot!(p1, xs, 1e3 .* E[mid+1, :]; color=:crimson, lw=2.5)

# flat-band width vs angle: the magic angle
θs = 0.80:0.05:1.60
widths = Float64[]
for θ in θs
    _, Eθ, nQθ = bm_bands(θ * π / 180, w0, w1, ħvF; nshell=5, npath=16)
    push!(widths, 1e3 * (maximum(Eθ[nQθ+1, :]) - minimum(Eθ[nQθ, :])))
end
θmagic = θs[argmin(widths)]
θpred = 2 * asind(w1 / (0.586 * ħvF * 2 * (4π / (3 * A0))))     # α = w1/(ħv_F kθ) = 0.586
@printf("magic angle from the scan: %.2f°  (BM condition α = 0.586 predicts %.2f°)\n",
        θmagic, θpred)
# NB: the LDA Fermi velocity is ~10%% below the GW/experimental one, which pushes the magic
# angle above the observed 1.05–1.1°; the BM condition α ≈ 0.586 is what the model tests.
p2 = plot(θs, widths; marker=:circle, lw=2, color=:navy, legend=false,
          xlabel="twist angle θ (deg)", ylabel="flat-band width (meV)",
          title=@sprintf("magic angle ≈ %.2f° (LDA v_F)", θmagic))
vline!(p2, [θmagic]; color=:crimson, style=:dash)

# interlayer coupling over the stacking cell (first-harmonic model through the 9 DFT points)
fine = range(0, 1; length=60)
TAAmap = [1e3 * abs(sum(hbest.cAA[i] * cis(-2π * (trio[i][1] * d1 + trio[i][2] * d2))
                        for i in 1:3))
          for d2 in fine, d1 in fine]
p3 = heatmap(fine, fine, TAAmap; xlabel="d₁", ylabel="d₂", color=:viridis,
             title="|T_AA(d)| (meV) over the stacking cell", aspect_ratio=1)
scatter!(p3, [d1 for (d1, _) in shifts], [d2 for (_, d2) in shifts];
         color=:white, ms=3, label="DFT samples")

out = joinpath(@__DIR__, "output")
mkpath(out)
savefig(plot(p3, p1, p2; layout=(1, 3), size=(1500, 420), margin=6Plots.mm),
        joinpath(out, "09_tbg_local_stacking.png"))
println("wrote ", joinpath(out, "09_tbg_local_stacking.png"))
