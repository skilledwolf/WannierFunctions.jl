# Stengel–Spaldin functional (`use_ss_functional`, wannier90 v4-dev; Stengel & Spaldin,
# PRB 73, 075121 (2006)). Instead of the MV sum of per-k Berry phases, the objective uses the
# k-AVERAGED diagonal overlap in a uniform b-vector order,
#
#   M̄_nn(b) = (1/N_k) Σ_k M_nn(k, b),      Ω_SS = Σ_n Σ_b w_b (1 − |M̄_nn(b)|²),
#
# whose minimisation gives Wannier functions with real-space-invariant ("single-point")
# centres  r̄_n = −Σ_b w_b b·Im ln M̄_nn(b).  Port of the use_ss_functional branches of
# wann_omega / wann_domega (wannierise.F90); the final REPORTED spread is the standard MV
# decomposition at the SS-optimal gauge, matching the reference's Final State block.

using LinearAlgebra

"""
    SSData

Uniform b-vector order maps for the Stengel–Spaldin functional: `nnord[nn, k]` is the index at
k of the b-vector equal to `bvec[:, nn, 1]`, `nnrev[nn, k]` the index of its negative
(the reference's `kmesh_info%nnord` / `%nnrev`).
"""
struct SSData
    nnord::Matrix{Int}
    nnrev::Matrix{Int}
end

"Build the uniform-order index maps from the b-vector table."
function ss_data(bv::BVectors)
    nk = size(bv.bvec, 3)
    nntot = bv.nntot
    nnord = Matrix{Int}(undef, nntot, nk)
    nnrev = Matrix{Int}(undef, nntot, nk)
    for k in 1:nk, nn in 1:nntot
        v = SVector{3,Float64}(bv.bvec[1, nn, 1], bv.bvec[2, nn, 1], bv.bvec[3, nn, 1])
        io = ir = 0
        for c in 1:nntot
            u = SVector{3,Float64}(bv.bvec[1, c, k], bv.bvec[2, c, k], bv.bvec[3, c, k])
            norm(u - v) < 1e-8 && (io = c)
            norm(u + v) < 1e-8 && (ir = c)
        end
        (io == 0 || ir == 0) && error("ss_data: b-vector sets differ across k-points")
        nnord[nn, k] = io
        nnrev[nn, k] = ir
    end
    return SSData(nnord, nnrev)
end

"k-averaged diagonal overlaps M̄_nn(b) in the uniform b-order (num_wann × nntot)."
function _ss_summnn(Mrot::Array{ComplexF64,4}, ss::SSData)
    nw, _, nntot, nk = size(Mrot)
    s = zeros(ComplexF64, nw, nntot)
    for k in 1:nk, nn in 1:nntot
        cnn = ss.nnord[nn, k]
        for n in 1:nw
            s[n, nn] += Mrot[n, n, cnn, k]
        end
    end
    return s ./ nk
end

"""
    ss_spread(Mrot, bv, ss) -> SpreadResult

The SS objective and single-point centres. Only `Ω` (the objective) and `centres` are
meaningful during minimisation; the ΩI/OD/D fields are zeroed — the standard MV decomposition
is recomputed at the converged gauge for reporting.
"""
function ss_spread(Mrot::Array{ComplexF64,4}, bv::BVectors, ss::SSData)
    nw = size(Mrot, 1)
    smn = _ss_summnn(Mrot, ss)
    centres = zeros(3, nw)
    spreads = zeros(nw)
    for nn in 1:bv.nntot
        w = bv.wb[nn, 1]
        for n in 1:nw
            lnv = imag(log(smn[n, nn]))
            for i in 1:3
                centres[i, n] -= w * bv.bvec[i, nn, 1] * lnv
            end
            spreads[n] += w * (1.0 - abs2(smn[n, nn]))
        end
    end
    return SpreadResult(centres, spreads, sum(spreads), 0.0, 0.0, 0.0)
end

"""
    ss_gradient(Mrot, bv, ss) -> G

Gradient of Ω_SS in the reference's `cdodq` convention (drop-in for `omega_gradient`):
four M̄-weighted terms over the ±b pairs, divided by N_k.
"""
function ss_gradient(Mrot::Array{ComplexF64,4}, bv::BVectors, ss::SSData)
    nw, _, nntot, nk = size(Mrot)
    smn = _ss_summnn(Mrot, ss)
    G = zeros(ComplexF64, nw, nw, nk)
    for k in 1:nk, nn in 1:nntot
        cnn = ss.nnord[nn, k]
        cnn2 = ss.nnrev[nn, k]
        w = bv.wb[nn, 1]
        for n in 1:nw, m in 1:nw
            G[m, n, k] += w * (Mrot[m, n, cnn, k] * conj(smn[n, nn]) -
                               conj(Mrot[n, m, cnn2, k]) * conj(smn[m, nn]) -
                               conj(Mrot[n, m, cnn, k]) * smn[m, nn] +
                               Mrot[m, n, cnn2, k] * smn[n, nn])
        end
    end
    G ./= nk
    return G
end
