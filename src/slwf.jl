# Selective localisation of Wannier functions with constrained centres (SLWF+C;
# Wang–Jia–Yuan–Bernevig PRB 90, 165125 (2014)). The localiser minimises the spread of only
# the first `slwf_num` Wannier functions (the reported objective is "Omega Total_C", NOT the
# total Ω), optionally with a penalty λ·Σ|r̄_n − c_n|² pinning their centres to target sites
# `c_n`. The unitary rotation is still full U(N) over all WFs — only the objective and its
# gradient are restricted to the selected subset; the non-selected WFs absorb the
# delocalisation. Exact conventions in docs/reference-notes/slwf-constrained.md.
#
# With λ = 1 the objective reduces to Ω_C = Σ_{n≤slwf_num} [spread_n + |r̄_n − c_n|²].

using LinearAlgebra
using StaticArrays

"""
    SLWF

Selective-localisation controls: minimise the spread of Wannier functions `1..num`, with an
optional centre constraint (`constrain`) of strength `lambda` toward the Cartesian sites
`centres` (3×num).
"""
struct SLWF
    num::Int
    lambda::Float64
    constrain::Bool
    centres::Matrix{Float64}     # 3 × num (Cartesian Å)
end

"""
    slwf_omega(Mrot, bv, slwf) -> (; ΩC, centres)

The SLWF+C objective Ω_C and the (standard) Wannier centres, from the gauge-rotated overlaps.
"""
function slwf_omega(Mrot::Array{ComplexF64,4}, bv::BVectors, slwf::SLWF)
    nw = size(Mrot, 1); nntot = size(Mrot, 3); nk = size(Mrot, 4)
    invNk = 1.0 / nk
    S = slwf.num
    λ = slwf.constrain ? slwf.lambda : 0.0
    # centres (all WFs, standard MV definition)
    r = zeros(3, nw)
    for k in 1:nk, b in 1:nntot
        w = bv.wb[b, k]
        bb = SVector(bv.bvec[1, b, k], bv.bvec[2, b, k], bv.bvec[3, b, k])
        for n in 1:nw
            r[:, n] .-= (invNk * w * imag(log(Mrot[n, n, b, k]))) .* bb
        end
    end
    om_iod = 0.0; om_d = 0.0; om_nu = 0.0
    for k in 1:nk, b in 1:nntot
        w = bv.wb[b, k]
        bb = SVector(bv.bvec[1, b, k], bv.bvec[2, b, k], bv.bvec[3, b, k])
        summ = 0.0
        for n in 1:S
            mnn = Mrot[n, n, b, k]
            iml = imag(log(mnn))
            summ += abs2(mnn) - λ * iml^2
            brn = bb[1] * r[1, n] + bb[2] * r[2, n] + bb[3] * r[3, n]
            om_d += invNk * (1.0 - λ) * w * (iml + brn)^2
            if slwf.constrain
                bcn = bb[1] * slwf.centres[1, n] + bb[2] * slwf.centres[2, n] +
                      bb[3] * slwf.centres[3, n]
                om_nu += invNk * 2.0 * λ * w * iml * bcn
            end
        end
        om_iod += invNk * w * (S - summ)
    end
    if slwf.constrain
        for n in 1:S
            om_nu += λ * (slwf.centres[1, n]^2 + slwf.centres[2, n]^2 + slwf.centres[3, n]^2)
        end
    end
    return (; ΩC=om_iod + om_d + om_nu, centres=r)
end

"""
    slwf_gradient(Mrot, bv, slwf) -> G

The SLWF+C gradient field (same convention as `omega_gradient`, ×4/N_k), a verbatim port of
wann_domega's `selective_loc` branch: the standard MV gradient restricted to the selected WF
block plus the λ centre-constraint terms. Only elements with m ≤ slwf_num or n ≤ slwf_num are
nonzero (the non-selected block is free).
"""
function slwf_gradient(Mrot::Array{ComplexF64,4}, bv::BVectors, slwf::SLWF)
    nw = size(Mrot, 1); nntot = size(Mrot, 3); nk = size(Mrot, 4)
    S = slwf.num
    λ = slwf.constrain ? slwf.lambda : 0.0
    im2 = complex(0.0, -0.5)
    # centres r̄_n (standard MV) from Mrot
    r = slwf_omega(Mrot, bv, slwf).centres
    G = zeros(ComplexF64, nw, nw, nk)
    for k in 1:nk
        for b in 1:nntot
            w = bv.wb[b, k]
            bb = SVector(bv.bvec[1, b, k], bv.bvec[2, b, k], bv.bvec[3, b, k])
            Mk = @view Mrot[:, :, b, k]
            lnt = [w * imag(log(Mk[n, n])) for n in 1:nw]
            rnkb = [bb[1] * r[1, n] + bb[2] * r[2, n] + bb[3] * r[3, n] for n in 1:nw]
            r0kb = zeros(nw)
            if slwf.constrain
                for n in 1:S
                    r0kb[n] = bb[1] * slwf.centres[1, n] + bb[2] * slwf.centres[2, n] +
                              bb[3] * slwf.centres[3, n]
                end
            end
            crt = [Mk[m, n] / Mk[n, n] for m in 1:nw, n in 1:nw]
            cr = [Mk[m, n] * conj(Mk[n, n]) for m in 1:nw, n in 1:nw]
            for n in 1:nw, m in 1:nw
                if m <= S && n <= S
                    g = w * 0.5 * (cr[m, n] - conj(cr[n, m]))
                    g -= (crt[m, n] * lnt[n] + conj(crt[n, m] * lnt[m])) * im2
                    g -= (crt[m, n] * rnkb[n] + conj(crt[n, m] * rnkb[m])) * im2
                    if slwf.constrain
                        g += λ * (crt[m, n] * lnt[n] + conj(crt[n, m] * lnt[m])) * im2
                        g += w * λ * (crt[m, n] * rnkb[n] + conj(crt[n, m] * rnkb[m])) * im2
                        g -= λ * (crt[m, n] * lnt[n] + conj(crt[n, m]) * lnt[m]) * im2
                        g -= w * λ * (r0kb[n] * crt[m, n] + r0kb[m] * conj(crt[n, m])) * im2
                    end
                    G[m, n, k] += g
                elseif m <= S            # n > S
                    g = -w * 0.5 * conj(cr[n, m])
                    g -= conj(crt[n, m] * (lnt[m] + w * rnkb[m])) * im2
                    if slwf.constrain
                        g += λ * conj(crt[n, m] * (lnt[m] + w * rnkb[m])) * im2 -
                             λ * conj(crt[n, m]) * lnt[m] * im2
                        g -= w * λ * r0kb[m] * conj(crt[n, m]) * im2
                    end
                    G[m, n, k] += g
                elseif n <= S            # m > S
                    g = w * cr[m, n] * 0.5 - crt[m, n] * (lnt[n] + w * rnkb[n]) * im2
                    if slwf.constrain
                        g += λ * crt[m, n] * (lnt[n] + w * rnkb[n]) * im2 -
                             λ * crt[m, n] * lnt[n] * im2
                        g -= w * λ * r0kb[n] * crt[m, n] * im2
                    end
                    G[m, n, k] += g
                end
            end
        end
    end
    G .*= (4.0 / nk)
    return G
end
