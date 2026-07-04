# BoltzWann — Boltzmann transport in the relaxation-time approximation on Wannier-interpolated
# bands (Pizzi, Volja, Kozinsky, Fornari & Marzari, CPC 185, 422 (2014); boltzwann.F90).
#
#   TDF_ij(ε) = τ/(V N_k) Σ_{n,k} g_s v_i v_j δ(ε − ε_nk)   [1/ħ² · eV·fs/Å; g_s = 2 non-spinor]
#   σ_ij(μ,T)      =  e³/ħ²·10⁻⁵ · Σ_ε TDF_ij (−∂f/∂ε) Δε              [1/Ω/m]
#   (σS)_ij(μ,T)   =  e³/ħ²·10⁻⁵ · Σ_ε TDF_ij (−∂f/∂ε)(ε−μ) Δε / T    [A/m/K]
#   S              =  −σ⁻¹·(σS)                                        [V/K]
#   K_ij(μ,T)      =  e³/ħ²·10⁻⁵ · Σ_ε TDF_ij (−∂f/∂ε)(ε−μ)² Δε / T   [W/m/K]
# δ uses the same fixed-width Gaussian / histogram-fallback machinery as the DOS.

using LinearAlgebra
using StaticArrays

const KB_SI = 1.3806504e-23        # CODATA2006

"Symmetric-tensor packing order used by BoltzWann files: xx, xy, yy, xz, yz, zz."
const TDF_PAIRS = ((1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3))

"""
    BoltzWannResult

TDF energy grid + tensor (6 × nE), and per (μ, T): electrical conductivity σ (1/Ω/m),
σ·S (A/m/K), Seebeck S (V/K), and K (W/m/K), each 3×3 stored as `[iμ, iT]` matrices.
"""
struct BoltzWannResult
    energies::Vector{Float64}
    tdf::Matrix{Float64}
    mus::Vector{Float64}
    temps::Vector{Float64}
    elcond::Array{Float64,4}       # 3×3×nμ×nT
    sigmas::Array{Float64,4}
    seebeck::Array{Float64,4}
    kappa::Array{Float64,4}
end

"−∂f/∂ε of the Fermi function (1/eV), with the reference's overflow cutoff."
function minus_fermi_derivative(E::Float64, mu::Float64, KT::Float64)
    x = (E - mu) / KT
    abs(x) > 36.0 && return 0.0
    return exp(x) / (KT * (exp(x) + 1.0)^2)
end

"""
    boltzwann(bm; kmesh, relax_time, mus, temps, tdf_energy_step=0.001,
              tdf_smr_width=0.0, win=nothing, elec_per_state=2) -> BoltzWannResult

Transport distribution function and RTA transport tensors. `relax_time` in fs; `win` is the
(emin, emax) band-energy window (defaults to the eigenvalue range of H(R)'s source bands as
postw90 uses the disentanglement window — pass it explicitly for parity).
"""
function boltzwann(bm::BerryModel; kmesh::NTuple{3,Int}=(25, 25, 25), relax_time::Float64=10.0,
                   mus::Vector{Float64}, temps::Vector{Float64},
                   tdf_energy_step::Float64=0.001, tdf_smr_width::Float64=0.0,
                   win::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                   elec_per_state::Int=2)
    # TDF energy grid (boltzwann.F90:270-290)
    exceed = max(10.0 * tdf_smr_width, 0.2)
    if win === nothing
        error("boltzwann: pass win=(win_min, win_max) — the reference uses the disentanglement window")
    end
    emin = win[1] - exceed
    nE = Int(floor((win[2] - win[1] + 2 * exceed) / tdf_energy_step)) + 1
    nE == 1 && (nE = 2)
    es = [emin + (i - 1) * tdf_energy_step for i in 1:nE]

    nktot = prod(kmesh)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    acc = [zeros(6, nE) for _ in 1:nktot]
    Threads.@threads for idx in 1:nktot
        _tdf_kpoint!(acc[idx], bm, kl[idx], es, tdf_smr_width, elec_per_state)
    end
    tdf = sum(acc) .* (relax_time / (nktot * cell_volume(bm.lattice)))

    # transport integrals per (μ, T)
    nμ, nT = length(mus), length(temps)
    elcond = zeros(3, 3, nμ, nT); sigmas = zeros(3, 3, nμ, nT)
    seebeck = zeros(3, 3, nμ, nT); kappa = zeros(3, 3, nμ, nT)
    unit = ELEM_CHARGE_SI^3 / HBAR_SI^2 * 1e-5
    for iT in 1:nT, iμ in 1:nμ
        KT = temps[iT] * KB_SI / ELEM_CHARGE_SI
        w = [minus_fermi_derivative(es[ie], mus[iμ], KT) for ie in 1:nE]
        c0 = zeros(6); c1 = zeros(6); c2 = zeros(6)
        for ie in 1:nE
            dμ = es[ie] - mus[iμ]
            for c in 1:6
                t = tdf[c, ie] * w[ie]
                c0[c] += t
                c1[c] += t * dμ
                c2[c] += t * dμ^2
            end
        end
        σ = zeros(3, 3); σS = zeros(3, 3); K = zeros(3, 3)
        for (c, (i, j)) in enumerate(TDF_PAIRS)
            σ[i, j] = σ[j, i] = c0[c] * tdf_energy_step * unit
            σS[i, j] = σS[j, i] = c1[c] * tdf_energy_step / temps[iT] * unit
            K[i, j] = K[j, i] = c2[c] * tdf_energy_step / temps[iT] * unit
        end
        elcond[:, :, iμ, iT] = σ
        sigmas[:, :, iμ, iT] = σS
        kappa[:, :, iμ, iT] = K
        seebeck[:, :, iμ, iT] = abs(det(σ)) > 0 ? -σ \ σS : zeros(3, 3)
    end
    return BoltzWannResult(es, tdf, mus, temps, elcond, sigmas, seebeck, kappa)
end

function _tdf_kpoint!(out::Matrix{Float64}, bm::BerryModel, kf::SVector{3,Float64},
                      es::Vector{Float64}, smr_width::Float64, elec_per_state::Int)
    E, dE = eig_deleig(bm, kf; deriv=true)
    nE = length(es)
    binwidth = es[2] - es[1]
    for n in 1:length(E)
        vv = SVector(dE[1, n], dE[2, n], dE[3, n])
        if smr_width / binwidth < 2.0                       # histogram fallback
            ie = round(Int, (E[n] - es[1]) / binwidth) + 1
            if 1 <= ie <= nE
                w = elec_per_state / binwidth
                for (c, (i, j)) in enumerate(TDF_PAIRS)
                    out[c, ie] += w * vv[i] * vv[j]
                end
            end
        else
            for ie in 1:nE
                abs(es[ie] - E[n]) > 10.0 * smr_width && continue
                δ = elec_per_state * _w0gauss((es[ie] - E[n]) / smr_width) / smr_width
                for (c, (i, j)) in enumerate(TDF_PAIRS)
                    out[c, ie] += δ * vv[i] * vv[j]
                end
            end
        end
    end
    return out
end
