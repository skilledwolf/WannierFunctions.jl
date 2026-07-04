# Wannier spread functional and its Marzari–Vanderbilt decomposition.
#
# Formulas (PRB 56, 12847 (1997), single-k-point-per-neighbour finite differences), evaluated on
# the gauge-rotated overlaps M̃ = M̃_{k,b}[m,n]:
#
#   centre    r̄_n   = −(1/N_k) Σ_{k,b} w_b b Im ln M̃_{nn}
#   ⟨r²⟩_n         = (1/N_k) Σ_{k,b} w_b [ (1 − |M̃_{nn}|²) + (Im ln M̃_{nn})² ]
#   spread    Ω_n   = ⟨r²⟩_n − |r̄_n|²                 ;   Ω = Σ_n Ω_n
#   Ω_I            = (1/N_k) Σ_{k,b} w_b ( N_w − Σ_{m,n} |M̃_{mn}|² )   (gauge-invariant)
#   Ω_OD           = (1/N_k) Σ_{k,b} w_b Σ_{m≠n} |M̃_{mn}|²
#   Ω_D            = (1/N_k) Σ_{k,b} w_b Σ_n ( −Im ln M̃_{nn} − b·r̄_n )²
# with Ω = Ω_I + Ω_OD + Ω_D.
#
# The k-sums are threaded (per-k partials, reduced after the loop).

using LinearAlgebra
using StaticArrays

"Container for the results of a spread evaluation."
struct SpreadResult
    centres::Matrix{Float64}   # (3 × num_wann), Cartesian Å
    spreads::Vector{Float64}   # (num_wann), Ų
    Ω::Float64
    ΩI::Float64
    ΩOD::Float64
    ΩD::Float64
end

"""
    compute_spread(Mrot, bv) -> SpreadResult

Evaluate Wannier centres, per-function spreads, and the Ω_I/Ω_OD/Ω_D decomposition from the
gauge-rotated overlaps `Mrot[n,n',b,k]` and the neighbour geometry `bv::BVectors`.
"""
function compute_spread(Mrot::Array{ComplexF64,4}, bv::BVectors)
    nw = size(Mrot, 1)
    nntot = size(Mrot, 3)
    nk = size(Mrot, 4)
    invNk = 1.0 / nk

    # Pass 1 — centres (needed by Ω_D): per-k partials, then reduce.
    rk = zeros(3, nw, nk)
    @maybe_threads (nk >= THREAD_MIN) for k in 1:nk
        @inbounds for b in 1:nntot
            w = bv.wb[b, k]
            bx, by, bz = bv.bvec[1, b, k], bv.bvec[2, b, k], bv.bvec[3, b, k]
            for n in 1:nw
                imln = imag(log(Mrot[n, n, b, k]))
                f = invNk * w * imln
                rk[1, n, k] -= f * bx
                rk[2, n, k] -= f * by
                rk[3, n, k] -= f * bz
            end
        end
    end
    r = dropdims(sum(rk; dims=3); dims=3)

    # Pass 2 — spreads and the decomposition: per-k partials, then reduce.
    r2k = zeros(nw, nk)
    ΩIk = zeros(nk); ΩODk = zeros(nk); ΩDk = zeros(nk)
    @maybe_threads (nk >= THREAD_MIN) for k in 1:nk
        @inbounds for b in 1:nntot
            w = bv.wb[b, k]
            bx, by, bz = bv.bvec[1, b, k], bv.bvec[2, b, k], bv.bvec[3, b, k]
            sall = 0.0
            for n in 1:nw, m in 1:nw
                sall += abs2(Mrot[m, n, b, k])
            end
            sdiag = 0.0
            for n in 1:nw
                mnn = Mrot[n, n, b, k]
                a2 = abs2(mnn)
                sdiag += a2
                imln = imag(log(mnn))
                r2k[n, k] += invNk * w * ((1.0 - a2) + imln^2)
                q = -imln - (bx * r[1, n] + by * r[2, n] + bz * r[3, n])
                ΩDk[k] += invNk * w * q^2
            end
            ΩIk[k] += invNk * w * (nw - sall)
            ΩODk[k] += invNk * w * (sall - sdiag)
        end
    end
    r2 = dropdims(sum(r2k; dims=2); dims=2)
    ΩI, ΩOD, ΩD = sum(ΩIk), sum(ΩODk), sum(ΩDk)

    spreads = [r2[n] - (r[1, n]^2 + r[2, n]^2 + r[3, n]^2) for n in 1:nw]
    Ω = sum(spreads)
    return SpreadResult(r, spreads, Ω, ΩI, ΩOD, ΩD)
end
