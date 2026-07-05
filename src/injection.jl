# Circular photogalvanic injection current (Lihm & Park, PRB 105, 045201 (2022), Eq. 10;
# the WannierBerri `InjectionCurrent` calculator). The injection-current rate tensor is
#
#   η_abc(ω) = C/(N_k V) Σ_k Σ_{m≠n} (f_n − f_m) (v^a_m − v^a_n) A^b_mn A^c_nm δ(E_m − E_n − ω)
#
# with A the covariant (Hermitian) Berry connection A^b_mn = [U†A^W_b U]_mn + i (∂_b H̄)_mn/(E_n−E_m)
# for m≠n, v^a the band velocity, f the T=0 occupation, and δ a Gaussian of fixed width. The
# constant C = −0.0011617929… reproduces WannierBerri's units. Validated against WannierBerri on
# a shared tight-binding model.

using LinearAlgebra
using StaticArrays

const INJ_FACTOR = -0.001161792936599752    # WannierBerri factor_injection_current

"""
    injection_current(bm; freqs, fermi_energy, kmesh, smr_width=0.1) -> Array{Float64,4}

Circular injection-current tensor η_abc(ω) (3×3×3 × nfreq), on a Γ-centred `kmesh`, matching
the WannierBerri `InjectionCurrent` calculator (Gaussian smearing `smr_width`). Needs a
`BerryModel` with the Berry connection A(R) (from `.mmn` or a `_tb.dat`).
"""
function injection_current(bm::BerryModel; freqs::Vector{Float64}, fermi_energy::Float64,
                           kmesh::NTuple{3,Int}=(12, 12, 12), smr_width::Float64=0.1)
    size(bm.Ar, 4) == 3 || error("injection_current needs A(R) (.mmn or _tb.dat r(R))")
    nfreq = length(freqs)
    nktot = prod(kmesh)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    per_k = Vector{Array{Float64,4}}(undef, nktot)
    Threads.@threads for idx in 1:nktot
        per_k[idx] = _inj_kpoint(bm, kl[idx], freqs, fermi_energy, smr_width)
    end
    η = reduce(+, per_k)
    η .*= INJ_FACTOR / (nktot * cell_volume(bm.lattice))
    return η
end

function _inj_kpoint(bm::BerryModel, kf::SVector{3,Float64}, freqs::Vector{Float64},
                     ef::Float64, w::Float64)
    kd = _berry_kdata(bm, kf)
    E, U = kd.E, kd.U
    nw = length(E)
    nfreq = length(freqs)
    # covariant Berry connection A^b_mn = U†A_b U + i (U†∂_bH U)_mn/(E_n − E_m), and velocities
    AH = [U' * kd.A[b] * U for b in 1:3]
    for b in 1:3, m in 1:nw, n in 1:nw
        n == m && continue
        AH[b][m, n] += im * kd.dHh[b][m, n] / (E[n] - E[m])
    end
    v = [real(kd.dHh[a][m, m]) for a in 1:3, m in 1:nw]
    occ = [E[m] < ef ? 1.0 : 0.0 for m in 1:nw]
    gnorm = 1.0 / (w * sqrt(pi))
    out = zeros(3, 3, 3, nfreq)
    for m in 1:nw, n in 1:nw
        n == m && continue
        occf = occ[n] - occ[m]
        occf == 0.0 && continue
        Emn = E[m] - E[n]
        for a in 1:3, b in 1:3, c in 1:3
            imn = occf * (v[a, m] - v[a, n]) * real(AH[b][m, n] * AH[c][n, m])
            imn == 0.0 && continue
            for (iw, ω) in enumerate(freqs)
                out[a, b, c, iw] += imn * gnorm * exp(-((Emn - ω) / w)^2)
            end
        end
    end
    return out
end
