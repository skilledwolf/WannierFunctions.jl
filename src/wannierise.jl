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
    @maybe_threads (nk >= THREAD_MIN) for k in 1:nk
        R[:, :, k] = expm_antiherm(@view gen[:, :, k])
    end
    return R
end

"Apply gauge rotation `R[:,:,k]` in place: U ← U·R and M̃_{k,b} ← R_k† M̃_{k,b} R_{k+b}."
function apply_rotation!(U::Array{ComplexF64,3}, Mrot::Array{ComplexF64,4},
                        kpb::Matrix{Int}, R::Array{ComplexF64,3})
    nw, _, nk = size(U)
    nntot = size(Mrot, 3)
    @maybe_threads (nk >= THREAD_MIN) for k in 1:nk
        U[:, :, k] = (@view U[:, :, k]) * (@view R[:, :, k])
    end
    # NB: Mrot[:, :, b, k] reads R at both k and its neighbour, but writes only slice (b, k),
    # and R is not mutated here — safe to thread over k.
    @maybe_threads (nk >= THREAD_MIN) for k in 1:nk
        for b in 1:nntot
            kb = kpb[b, k]
            Mrot[:, :, b, k] = (@view R[:, :, k])' * (@view Mrot[:, :, b, k]) * (@view R[:, :, kb])
        end
    end
    return nothing
end

"""
    omega_gradient(Mrot, bv, centres) -> G

Analytic gradient dΩ/dW as the anti-Hermitian matrix `G[:,:,k]`, following `wann_domega`:
`G = (4/N_k) Σ_{k,b} w_b ( A[R] − S[T] )` with `R_{mn}=M_{mn}·conj(M_nn)`,
`R̃_{mn}=M_{mn}/M_nn`, `q_n = Im ln M_nn + b·r_n`.
"""
function omega_gradient(Mrot::Array{ComplexF64,4}, bv::BVectors, centres::Matrix{Float64};
                        guides::Union{Nothing,Matrix{Float64}}=nothing)
    nw = size(Mrot, 1); nntot = size(Mrot, 3); nk = size(Mrot, 4)
    G = zeros(ComplexF64, nw, nw, nk)
    @maybe_threads (nk >= THREAD_MIN) for k in 1:nk
      # thread-local scratch
      lnt = Vector{Float64}(undef, nw)      # w_b · Im ln M_nn
      rnkb = Vector{Float64}(undef, nw)     # b · r_n
      mnn = Vector{ComplexF64}(undef, nw)
      for b in 1:nntot
        w = bv.wb[b, k]
        bx, by, bz = bv.bvec[1, b, k], bv.bvec[2, b, k], bv.bvec[3, b, k]
        @inbounds for n in 1:nw
            mnn[n] = Mrot[n, n, b, k]
            lnt[n] = w * _guided_imln(mnn[n], guides, n, bx, by, bz)
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
    end
    G .*= (4.0 / nk)
    return G
end

# Preconditioner: real-space Lorentzian low-pass filter of the gradient (precond-cg.md).
# g_R = (1/N_k) Σ_k e^{-i2πk·R} G_k ; g_R *= 1/(1+|R|²/α), α = 10·Ω/nw ; g̃_k = Σ_R e^{+i2πk·R} g_R/ndegen.
# With the filter ≡ 1 the forward+backward pair is the identity, so it is purely the R-weighting.
function _precond_filter(G::Array{ComplexF64,3}, pc::NamedTuple, om_tot::Float64, nw::Int)
    kfrac, irvec, ndegen, Rcart = pc.kfrac, pc.irvec, pc.ndegen, pc.Rcart
    nk = size(G, 3); nr = length(irvec)
    α = 10.0 * om_tot / nw
    gR = zeros(ComplexF64, nw, nw, nr)
    for ir in 1:nr
        R = SVector{3,Float64}(irvec[ir]...)
        acc = @view gR[:, :, ir]
        for k in 1:nk
            fac = cis(-TWOPI * dot(kfrac[k], R)) / nk
            @views acc .+= fac .* G[:, :, k]
        end
        acc .*= 1.0 / (1.0 + dot(Rcart[ir], Rcart[ir]) / α)
    end
    Gp = zeros(ComplexF64, nw, nw, nk)
    for k in 1:nk
        acc = @view Gp[:, :, k]
        for ir in 1:nr
            fac = cis(TWOPI * dot(kfrac[k], SVector{3,Float64}(irvec[ir]...))) / ndegen[ir]
            @views acc .+= fac .* gR[:, :, ir]
        end
    end
    return Gp
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
    wannierise(model; num_iter=100, algorithm=:rcg, ...) -> WannieriseResult

Minimise the Wannier spread starting from the Löwdin-projected gauge (isolated-bands case,
num_bands == num_wann). `algorithm = :rcg` (default) uses Riemannian conjugate gradient with a
true convergence criterion (`num_iter` is a maximum); `algorithm = :w90` reproduces the reference
Wannier90 optimiser exactly (fixed `num_iter` sweeps unless `conv_window > 1`).
"""
function wannierise(model::Model; num_iter::Int=model_num_iter(model), sitesym=nothing, kwargs...)
    U0 = initial_gauge(model.A)
    # site_symmetry: symmetrise the initial gauge so the whole trajectory stays symmetry-adapted.
    # The reconstruction U(Rk) = d_band·U(k)·d_wann† uses the band representation on the left
    # (for the isolated case d_band is the original .dmn band-rep; the rotation/gradient use
    # d_wann on both sides, keeping the relation invariant under U ← U·R).
    sitesym === nothing || symmetrize_u!(U0, sitesym, sitesym.d_band, sitesym.d_wann)
    Mrot0 = rotate_overlaps(model.M, U0, model.bvectors.kpb)
    return localize(U0, Mrot0, model.bvectors; num_iter=num_iter, sitesym=sitesym, kwargs...)
end

"""
    localize(U0, Mrot0, bv; num_iter=100, algorithm=:rcg, kwargs...) -> WannieriseResult

Run the spread minimisation from an initial square gauge `U0` (num_wann × num_wann × nkpt) and
its gauge-rotated overlaps `Mrot0`. Shared by the isolated-bands path (`wannierise`) and the
post-disentanglement handoff. Ω_I is invariant under this step.

`algorithm`:
- `:rcg` — Riemannian Polak–Ribière+ conjugate gradient on the product-of-unitaries manifold,
  parabolic-model line search with Armijo backtracking safeguard, convergence when |ΔΩ| stays
  below `conv_tol` (default 1e-10 Ų) for `conv_window` (default 3) iterations. Modern default.
- `:w90` — the reference Wannier90 minimiser, reproduced exactly (Fletcher–Reeves CG, parabolic
  line search, fixed sweep count with convergence checking off unless `conv_window > 1`). Use for
  bit-faithful parity with `wannier90.x`.
"""
function localize(U0::Array{ComplexF64,3}, Mrot0::Array{ComplexF64,4}, bv::BVectors;
                  algorithm::Symbol=:rcg, kwargs...)
    algorithm === :w90 && return _localize_w90(U0, Mrot0, bv; kwargs...)
    algorithm === :rcg && return _localize_rcg(U0, Mrot0, bv; kwargs...)
    error("unknown localisation algorithm $algorithm (expected :rcg or :w90)")
end

function _localize_w90(U0::Array{ComplexF64,3}, Mrot0::Array{ComplexF64,4}, bv::BVectors;
                  num_iter::Int=100, trial_step::Float64=2.0, num_cg_steps::Int=5,
                  conv_tol::Float64=CONV_TOL_DEFAULT, conv_window::Int=-1,
                  verbose::Bool=false,
                  guides::Union{Nothing,Matrix{Float64}}=nothing,
                  precond::Union{Nothing,NamedTuple}=nothing, slwf=nothing, sitesym=nothing)
    slwf === nothing || error("SLWF+C requires algorithm = :rcg (the Ω_C objective path)")
    sitesym === nothing || error("site_symmetry requires algorithm = :rcg")
    kpb = bv.kpb
    nw = size(U0, 1)
    nk = size(U0, 3)
    wbtot = sum(@view bv.wb[:, 1])

    U = copy(U0)
    Mrot = copy(Mrot0)
    sr = compute_spread(Mrot, bv; guides=guides)
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
        G = omega_gradient(Mrot, bv, sr.centres; guides=guides)
        # Preconditioned CG: replace the steepest-descent gradient by the Lorentzian
        # real-space-filtered gradient M⁻¹G (the FR direction and mixed inner product use it;
        # the line-search slope keeps the TRUE gradient — see precond-cg.md).
        Gpre = precond === nothing ? G : _precond_filter(G, precond, sr.Ω, nw)

        # Fletcher–Reeves conjugate-gradient coefficient.
        gcnorm1 = precond === nothing ? sum(abs2, G) : real(dot(Gpre, G))
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

        # Steepest-descent component uses the (possibly filtered) gradient; the CG-memory term
        # β·d_prev is untouched. The line-search slope doda0 always uses the TRUE gradient G.
        cdq = Gpre .+ gcfac .* cdqkeep
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
        srt = compute_spread(Mt, bv; guides=guides)

        alphamin, _, lquad = optimal_step(sr.Ω, srt.Ω, doda0, trial_step)

        old_om = sr.Ω
        if lquad
            Ropt = expm_all(cdq .* (alphamin / (4.0 * wbtot)))
            apply_rotation!(U, Mrot, kpb, Ropt)
            sr = compute_spread(Mrot, bv; guides=guides)
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

# ---------------------------------------------------------------------------
# :rcg — Riemannian Polak–Ribière+ conjugate gradient on ∏_k U(N).
#
# The localisation problem is optimisation over the product of unitary groups; anti-Hermitian
# matrices are the Lie algebra (tangent space at identity), and U ← U·exp(s·X) moves along the
# geodesic generated by X. `omega_gradient` returns the field G with dΩ/ds along exp(s·X/(4w_tot))
# equal to −Re⟨G,X⟩/(4w_tot), i.e. +G is the steepest-descent generator (Wannier90's convention,
# reused here so step scales are comparable across the two algorithms).
# ---------------------------------------------------------------------------

function _localize_rcg(U0::Array{ComplexF64,3}, Mrot0::Array{ComplexF64,4}, bv::BVectors;
                       num_iter::Int=1000, trial_step::Float64=2.0,
                       conv_tol::Float64=CONV_TOL_DEFAULT, conv_window::Int=3,
                       verbose::Bool=false, num_cg_steps::Int=0,
                       guides::Union{Nothing,Matrix{Float64}}=nothing,
                       precond::Union{Nothing,NamedTuple}=nothing,
                       slwf=nothing, sitesym=nothing)   # num_cg_steps/guides/precond: :w90-only
    (guides === nothing && precond === nothing) ||
        error("guiding_centres / precond require algorithm = :w90 (the reference path)")
    kpb = bv.kpb
    wbtot = sum(@view bv.wb[:, 1])
    fourw = 4.0 * wbtot

    U = copy(U0)
    Mrot = copy(Mrot0)
    sr = compute_spread(Mrot, bv; slwf=slwf)
    omega_trace = Float64[sr.Ω]

    inner(A, B) = real(dot(A, B))            # Re Σ_k tr(A_k† B_k)

    G_prev = zeros(ComplexF64, size(U))
    d = zeros(ComplexF64, size(U))
    step = trial_step
    history = fill(Inf, max(conv_window, 1))
    converged = false
    iter = 0
    while iter < num_iter
        iter += 1
        G = slwf === nothing ? omega_gradient(Mrot, bv, sr.centres) : slwf_gradient(Mrot, bv, slwf)
        sitesym === nothing || symmetrize_gradient!(G, sitesym)   # project onto symmetric tangent
        gnorm2 = inner(G, G)
        if gnorm2 < 1e-24                     # stationary point
            converged = true
            iter -= 1
            break
        end

        # Polak–Ribière+ with automatic reset.
        if iter == 1
            d = copy(G)
        else
            β = max(0.0, inner(G, G .- G_prev) / inner(G_prev, G_prev))
            d = G .+ β .* d
            # reset to steepest descent if not a descent direction
            inner(G, d) <= 0 && (d = copy(G))
        end
        G_prev = G
        # site_symmetry: propagate the (representative-k) rotation to the star so U stays
        # symmetry-adapted after the step.
        sitesym === nothing || symmetrize_rotation!(d, sitesym)

        slope = -inner(G, d) / fourw          # dΩ/ds at s=0 (negative)

        # Parabolic-model step (one trial evaluation), Armijo backtracking as safeguard.
        s = step
        local Ut, Mt, srt
        accepted = false
        for bt in 1:12
            R = expm_all(d .* (s / fourw))
            Ut = copy(U); Mt = copy(Mrot)
            apply_rotation!(Ut, Mt, kpb, R)
            srt = compute_spread(Mt, bv; slwf=slwf)
            if srt.Ω <= sr.Ω + 1.0e-4 * s * slope    # Armijo sufficient decrease
                # one parabolic refinement through (0, sr.Ω), slope, (s, srt.Ω)
                denom = srt.Ω - sr.Ω - slope * s
                if denom > 0
                    s_opt = clamp(-slope * s^2 / (2 * denom), 0.1 * s, 3.0 * s)
                    if abs(s_opt - s) > 1e-3 * s
                        Ro = expm_all(d .* (s_opt / fourw))
                        Uo = copy(U); Mo = copy(Mrot)
                        apply_rotation!(Uo, Mo, kpb, Ro)
                        sro = compute_spread(Mo, bv; slwf=slwf)
                        if sro.Ω < srt.Ω
                            Ut, Mt, srt, s = Uo, Mo, sro, s_opt
                        end
                    end
                end
                accepted = true
                break
            end
            s *= 0.5
        end
        if !accepted
            # Line search cannot decrease Ω. With a vanishing slope this *is* convergence
            # (we are at the minimum to within line-search resolution), not a failure.
            converged = abs(slope) < 1e-9
            break
        end

        old_om = sr.Ω
        U, Mrot, sr = Ut, Mt, srt
        step = min(2.0 * s, 20.0)             # carry an adaptive initial step
        push!(omega_trace, sr.Ω)
        verbose && @info "rcg iter $iter" Ω=sr.Ω step=s

        history = circshift(history, -1)
        history[end] = sr.Ω - old_om
        if iter >= conv_window && all(h -> abs(h) <= conv_tol, history)
            converged = true
            break
        end
    end

    return WannieriseResult(U, Mrot, sr, omega_trace, iter, converged)
end

model_num_iter(::Model) = 100    # fallback; callers usually pass num_iter from the .win
