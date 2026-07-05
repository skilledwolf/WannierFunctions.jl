# Γ-only real-orthogonal localisation (wann_main_gamma, wannierise.F90): at Γ with real
# wavefunctions the overlaps are real-symmetric up to the ±b pairing, and the MV spread is
# minimised over REAL orthogonal gauges by Jacobi-like 2×2 rotation sweeps that jointly
# "diagonalise" the weighted real and imaginary parts of the half-set M(b) (Gygi/Silvestrelli).
# The result is a real gauge — real Wannier functions — matching wannier90's gamma_only output.
#
# Conventions: the model M holds the CLOSED expanded set (read_model doubles the file's half
# set with adjoint partners); this driver works on the file half with doubled weights
# (w90's gamma kmesh solves B1 on the half set, so w_half = 2·w_full).

using LinearAlgebra

"One Jacobi sweep over all (nw1 < nw2) pairs (internal_new_u_and_m_gamma, ported verbatim)."
function _gamma_sweep!(m_w::Array{Float64,3}, ur::Matrix{Float64})
    nw = size(ur, 1)
    tn = size(m_w, 3)
    for nw1 in 1:nw, nw2 in nw1+1:nw
        a11 = 0.0; a12 = 0.0; a22 = 0.0
        for nn in 1:tn
            d = m_w[nw1, nw1, nn] - m_w[nw2, nw2, nn]
            a11 += d^2
            a12 += m_w[nw1, nw2, nn] * d
            a22 += m_w[nw1, nw2, nn]^2
        end
        a12 *= 2.0
        a22 *= 4.0
        a21 = a22 - a11
        local θ
        if abs(a12) > 1e-10
            twoθ = 0.5 * (a21 + sqrt(a21^2 + 4.0 * a12^2)) / a12
            θ = 0.5 * atan(twoθ)
        elseif a21 < 1e-10
            θ = 0.0
        else
            θ = 0.25 * π
        end
        cc, ss = cos(θ), sin(θ)
        for nn in 1:tn
            for i in 1:nw                       # M ← M·R
                r1 = m_w[i, nw1, nn] * cc + m_w[i, nw2, nn] * ss
                r2 = -m_w[i, nw1, nn] * ss + m_w[i, nw2, nn] * cc
                m_w[i, nw1, nn] = r1
                m_w[i, nw2, nn] = r2
            end
            for j in 1:nw                       # M ← Rᵀ·M
                r1 = cc * m_w[nw1, j, nn] + ss * m_w[nw2, j, nn]
                r2 = -ss * m_w[nw1, j, nn] + cc * m_w[nw2, j, nn]
                m_w[nw1, j, nn] = r1
                m_w[nw2, j, nn] = r2
            end
        end
        for i in 1:nw                           # U ← U·R
            r1 = ur[i, nw1] * cc + ur[i, nw2] * ss
            r2 = -ur[i, nw1] * ss + ur[i, nw2] * cc
            ur[i, nw1] = r1
            ur[i, nw2] = r2
        end
    end
    return m_w
end

"Centres and spread from the weighted half-set matrices (wann_omega_gamma). Returns
(rave, spreads_n, om_d, om_od, om_tot); `om_i` is the frozen first-pass invariant."
function _gamma_omega(m_w::Array{Float64,3}, wb::Vector{Float64}, bk::Matrix{Float64},
                      om_i::Float64)
    nw = size(m_w, 1)
    nn2 = length(wb)
    wbtot = sum(wb)
    lntmp = Matrix{Float64}(undef, nw, nn2)
    for nn in 1:nn2, n in 1:nw
        # orthorhombic (3 half-b's): pure atan2 branch, no sheet
        lntmp[n, nn] = atan(m_w[n, n, 2nn], m_w[n, n, 2nn-1])
    end
    rave = zeros(3, nw)
    for n in 1:nw, nn in 1:nn2, i in 1:3
        rave[i, n] -= wb[nn] * bk[i, nn] * lntmp[n, nn]
    end
    mnn2 = zeros(nw)
    spreads = fill(wbtot, nw)
    for n in 1:nw
        for nn in 1:nn2
            mnn2[n] += m_w[n, n, 2nn-1]^2 + m_w[n, n, 2nn]^2
            spreads[n] += wb[nn] * lntmp[n, nn]^2
        end
        spreads[n] -= mnn2[n]
    end
    om_od = wbtot * nw - sum(mnn2) - om_i
    om_d = 0.0
    if nn2 != 3
        for nn in 1:nn2, n in 1:nw
            brn = bk[1, nn] * rave[1, n] + bk[2, nn] * rave[2, n] + bk[3, nn] * rave[3, n]
            om_d += wb[nn] * (lntmp[n, nn] + brn)^2
        end
    end
    return rave, spreads, om_d, om_od, om_i + om_d + om_od
end

"""
    _localize_gamma(U0, Mrot0, bv; num_iter, conv_tol, conv_window) -> WannieriseResult

wannier90's `gamma_only` localiser: real-orthogonal Jacobi sweeps on the half b-set. `U0` must
be real-valued (the Löwdin projection of real Γ data is); `Mrot0` is the expanded closed set,
of which the first half is the file half. The returned gauge is real.
"""
function _localize_gamma(U0::Array{ComplexF64,3}, Mrot0::Array{ComplexF64,4}, bv::BVectors;
                         num_iter::Int=100, conv_tol::Float64=CONV_TOL_DEFAULT,
                         conv_window::Int=-1, verbose::Bool=false,
                         guides::Union{Nothing,Matrix{Float64}}=nothing, kwargs...)
    size(U0, 3) == 1 || error("gamma_only localisation requires a single k-point")
    maximum(abs.(imag.(U0))) < 1e-8 ||
        error("gamma_only: the initial gauge is not real — realify the subspace first")
    nw = size(U0, 1)
    nn2 = bv.nntot ÷ 2
    wb = 2.0 .* bv.wb[1:nn2, 1]                       # half-set convention: w_half = 2 w_full
    bk = Matrix{Float64}(bv.bvec[:, 1:nn2, 1])

    m_w = Array{Float64,3}(undef, nw, nw, 2nn2)
    for nn in 1:nn2
        sq = sqrt(wb[nn])
        m_w[:, :, 2nn-1] = sq .* real.(@view Mrot0[:, :, nn, 1])
        m_w[:, :, 2nn] = sq .* imag.(@view Mrot0[:, :, nn, 1])
    end
    ur = Matrix{Float64}(I, nw, nw)

    # first-pass gauge-invariant part: Ω_I = wbtot·nw − Σ |m_w|²
    om_i = sum(wb) * nw - sum(abs2, m_w)
    _, _, _, _, om = _gamma_omega(m_w, wb, bk, om_i)
    omega_trace = Float64[om]
    history = fill(Inf, max(conv_window, 0))
    converged = false
    niter = 0
    for iter in 1:num_iter
        niter = iter
        _gamma_sweep!(m_w, ur)
        _, _, _, _, om_new = _gamma_omega(m_w, wb, bk, om_i)
        push!(omega_trace, om_new)
        verbose && @info "gamma iter $iter" Ω=om_new
        Δ = om_new - om
        om = om_new
        if conv_window > 1
            history = circshift(history, -1)
            history[end] = Δ
            if iter >= conv_window && all(h -> abs(h) <= conv_tol, history)
                converged = true
                break
            end
        end
    end

    # Final gauge and expanded overlaps; report the standard full-set MV decomposition.
    U = similar(U0)
    U[:, :, 1] = U0[:, :, 1] * ur
    Mrot = similar(Mrot0)
    for nn in 1:nn2
        sq = 1.0 / sqrt(wb[nn])
        Mhalf = sq .* complex.(m_w[:, :, 2nn-1], m_w[:, :, 2nn])
        Mrot[:, :, nn, 1] = Mhalf
        Mrot[:, :, nn2+nn, 1] = Mhalf'
    end
    # Reporting: for a non-orthorhombic half-set the reference keeps guiding centres on, whose
    # sheets unwrap the centre log branches (essential for elongated cells, e.g. chains); for
    # the 3-b case it disables them and om_d ≡ 0 by the B1 identity.
    sr = compute_spread(Mrot, bv; guides=(nn2 == 3 ? nothing : guides))
    return WannieriseResult(U, Mrot, sr, omega_trace, niter, converged)
end

"Real orthonormal basis of the column span of a (numerically conjugation-closed) complex
embedding; errors if the span is not real. Used to realify the Γ disentangled subspace."
function realify_subspace(Uopt::Matrix{ComplexF64})
    nb, nw = size(Uopt)
    F = svd(hcat(real.(Uopt), imag.(Uopt)))
    B = F.U[:, 1:nw]                                  # nb × nw real orthonormal
    resid = maximum(abs.(Uopt .- B * (B' * Uopt)))
    resid < 1e-8 || error("gamma_only: disentangled subspace is not conjugation-closed " *
                          "(realification residual $resid)")
    return B
end
