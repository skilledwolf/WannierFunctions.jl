# Maximal-localisation iteration (Marzari–Vanderbilt gauge optimisation) for the
# num_bands == num_wann case. Mirrors `wann_main` in the reference `src/wannierise.F90`:
# analytic Ω gradient, Fletcher–Reeves conjugate-gradient search direction, parabolic line
# search, and a unitary update U ← U·exp(ΔW) via the matrix exponential of an anti-Hermitian
# generator. See docs/reference-notes/localization.md for the exact conventions.

using LinearAlgebra

"exp of an anti-Hermitian matrix `X`, via the Hermitian eigendecomposition of i·X."
function expm_antiherm(X::AbstractMatrix{ComplexF64})
    F = eigen(Hermitian(im .* X))        # i·X = V diag(ev) V†,  ev real
    return F.vectors * Diagonal(cis.(-F.values)) * F.vectors'
end

"Per-k anti-Hermitian generators → per-k unitary rotations."
function expm_all(gen::Array{ComplexF64,3})
    nw, _, nk = size(gen)
    R = Array{ComplexF64,3}(undef, nw, nw, nk)
    for k in 1:nk
        R[:, :, k] = expm_antiherm(@view gen[:, :, k])
    end
    return R
end

"Apply gauge rotation `R[:,:,k]` in place: U ← U·R and M̃_{k,b} ← R_k† M̃_{k,b} R_{k+b}."
function apply_rotation!(U::Array{ComplexF64,3}, Mrot::Array{ComplexF64,4},
                        kpb::Matrix{Int}, R::Array{ComplexF64,3})
    nw, _, nk = size(U)
    nntot = size(Mrot, 3)
    for k in 1:nk
        U[:, :, k] = (@view U[:, :, k]) * (@view R[:, :, k])
    end
    for k in 1:nk, b in 1:nntot
        kb = kpb[b, k]
        Mrot[:, :, b, k] = (@view R[:, :, k])' * (@view Mrot[:, :, b, k]) * (@view R[:, :, kb])
    end
    return nothing
end

"""
    omega_gradient(Mrot, bv, centres) -> G

Analytic gradient dΩ/dW as the anti-Hermitian matrix `G[:,:,k]`, following `wann_domega`:
`G = (4/N_k) Σ_{k,b} w_b ( A[R] − S[T] )` with `R_{mn}=M_{mn}·conj(M_nn)`,
`R̃_{mn}=M_{mn}/M_nn`, `q_n = Im ln M_nn + b·r_n`.
"""
function omega_gradient(Mrot::Array{ComplexF64,4}, bv::BVectors, centres::Matrix{Float64})
    nw = size(Mrot, 1); nntot = size(Mrot, 3); nk = size(Mrot, 4)
    G = zeros(ComplexF64, nw, nw, nk)
    lnt = Vector{Float64}(undef, nw)      # w_b · Im ln M_nn
    rnkb = Vector{Float64}(undef, nw)     # b · r_n
    mnn = Vector{ComplexF64}(undef, nw)
    for k in 1:nk, b in 1:nntot
        w = bv.wb[b, k]
        bx, by, bz = bv.bvec[1, b, k], bv.bvec[2, b, k], bv.bvec[3, b, k]
        @inbounds for n in 1:nw
            mnn[n] = Mrot[n, n, b, k]
            lnt[n] = w * imag(log(mnn[n]))
            rnkb[n] = bx * centres[1, n] + by * centres[2, n] + bz * centres[3, n]
        end
        @inbounds for n in 1:nw, m in 1:nw
            mmn = Mrot[m, n, b, k]
            mnm = Mrot[n, m, b, k]
            crt_mn = mmn / mnn[n]
            crt_nm = mnm / mnn[m]
            cr_mn = mmn * conj(mnn[n])
            cr_nm = mnm * conj(mnn[m])
            # A[R] = (R − R†)/2
            G[m, n, k] += w * 0.5 * (cr_mn - conj(cr_nm))
            # −S[T], with T split into the w_b-carrying ln part and the explicit-w_b rnkb part
            G[m, n, k] -= (crt_mn * lnt[n] + conj(crt_nm * lnt[m])) * complex(0.0, -0.5)
            G[m, n, k] -= w * (crt_mn * rnkb[n] + conj(crt_nm * rnkb[m])) * complex(0.0, -0.5)
        end
    end
    G .*= (4.0 / nk)
    return G
end

# Parabolic line search (internal_optimal_step): returns (alphamin, falphamin, lquad).
function optimal_step(wann_om::Float64, trial_om::Float64, doda0::Float64, trial_step::Float64)
    fac = trial_om - wann_om
    local shift
    if abs(fac) > floatmin(Float64)
        fac = 1.0 / fac
        shift = 1.0
    else
        fac = 1.0e6
        shift = fac * trial_om - fac * wann_om
    end
    eqb = fac * doda0
    eqa = shift - eqb * trial_step
    local alphamin, falphamin, lquad
    if abs(eqa / (fac * wann_om)) > eps(Float64)
        lquad = true
        alphamin = -0.5 * eqb / eqa * trial_step^2
        falphamin = wann_om - 0.25 * eqb^2 / (fac * eqa) * trial_step^2
    else
        lquad = false
        alphamin = trial_step
        falphamin = trial_om
    end
    if doda0 * alphamin > 0.0
        lquad = false
        alphamin = trial_step
        falphamin = trial_om
    end
    return alphamin, falphamin, lquad
end

"Result of a wannierisation run."
struct WannieriseResult
    U::Array{ComplexF64,3}         # final gauge (num_wann × num_wann × nkpt, applied on projections)
    Mrot::Array{ComplexF64,4}      # final gauge-rotated overlaps
    spread::SpreadResult           # final centres / spreads / Ω decomposition
    omega_trace::Vector{Float64}   # total Ω at iteration 0,1,…,num_iter
    niter::Int
    converged::Bool
end

"""
    wannierise(model; num_iter, trial_step=2.0, num_cg_steps=5,
               conv_tol=1e-10, conv_window=-1, verbose=false) -> WannieriseResult

Minimise the Wannier spread starting from the Löwdin-projected gauge. With the defaults
(`conv_window = -1`) the loop runs the full `num_iter` iterations, matching Wannier90.
"""
function wannierise(model::Model; num_iter::Int=model_num_iter(model),
                    trial_step::Float64=2.0, num_cg_steps::Int=5,
                    conv_tol::Float64=CONV_TOL_DEFAULT, conv_window::Int=-1,
                    verbose::Bool=false)
    bv = model.bvectors
    kpb = bv.kpb
    nw = model.num_wann
    nk = nkpt(model.kgrid)
    wbtot = sum(@view bv.wb[:, 1])

    U = initial_gauge(model.A)
    Mrot = rotate_overlaps(model.M, U, kpb)
    sr = compute_spread(Mrot, bv)
    omega_trace = Float64[sr.Ω]
    verbose && @info "iter 0" Ω=sr.Ω ΩI=sr.ΩI ΩOD=sr.ΩOD ΩD=sr.ΩD

    cdqkeep = zeros(ComplexF64, nw, nw, nk)
    ncg = 0
    gcnorm0 = 0.0
    converged = false
    history = fill(Inf, max(conv_window, 0))
    iter = 0
    while iter < num_iter
        iter += 1
        G = omega_gradient(Mrot, bv, sr.centres)

        # Fletcher–Reeves conjugate-gradient coefficient.
        gcnorm1 = sum(abs2, G)
        local gcfac
        if iter == 1 || ncg >= num_cg_steps
            gcfac = 0.0; ncg = 0
        elseif gcnorm0 > eps(Float64)
            gcfac = gcnorm1 / gcnorm0
            if gcfac > 3.0
                gcfac = 0.0; ncg = 0
            else
                ncg += 1
            end
        else
            gcfac = 0.0; ncg = 0
        end
        gcnorm0 = gcnorm1

        cdq = G .+ gcfac .* cdqkeep
        doda0 = -real(dot(G, cdq)) / (4.0 * wbtot)     # dot(a,b)=Σ conj(a)·b
        if doda0 > 0.0
            if ncg > 0
                cdq = copy(G); ncg = 0; gcfac = 0.0
                doda0 = -real(dot(G, cdq)) / (4.0 * wbtot)
                if doda0 > 0.0
                    cdq = -cdq; doda0 = -doda0
                end
            else
                cdq = -cdq; doda0 = -doda0
            end
        end
        cdqkeep = copy(cdq)

        # Trial step, then parabolic optimal step.
        Rtrial = expm_all(cdq .* (trial_step / (4.0 * wbtot)))
        Ut = copy(U); Mt = copy(Mrot)
        apply_rotation!(Ut, Mt, kpb, Rtrial)
        srt = compute_spread(Mt, bv)

        alphamin, _, lquad = optimal_step(sr.Ω, srt.Ω, doda0, trial_step)

        old_om = sr.Ω
        if lquad
            Ropt = expm_all(cdq .* (alphamin / (4.0 * wbtot)))
            apply_rotation!(U, Mrot, kpb, Ropt)
            sr = compute_spread(Mrot, bv)
        else
            U, Mrot, sr = Ut, Mt, srt
        end
        push!(omega_trace, sr.Ω)
        verbose && @info "iter $iter" Ω=sr.Ω α=alphamin

        if conv_window > 1
            history = circshift(history, -1)
            history[end] = sr.Ω - old_om
            if all(h -> abs(h) <= conv_tol, history)
                converged = true
                break
            end
        end
    end

    return WannieriseResult(U, Mrot, sr, omega_trace, iter, converged)
end

model_num_iter(::Model) = 100    # fallback; callers usually pass num_iter from the .win
