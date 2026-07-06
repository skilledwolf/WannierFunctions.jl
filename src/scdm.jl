# SCDM — automatic initial projections from the wavefunctions themselves
# (Damle, Lin & Ying, J. Chem. Theory Comput. 11, 1463 (2015); SCDM-k variant as used by
# pw2wannier90's scdm_proj). Removes the need to guess a `projections` block:
#
#   1. At the anchor k-point (Γ), form W[m, r] = f(ε_m) ψ*_m(r) on the real-space grid and run
#      column-pivoted QR; the first num_wann pivot columns select grid points r_j where the
#      occupied manifold is well represented ("selected columns of the density matrix").
#   2. At every k, the projection is the wavefunction value at those points,
#      A_mn(k) = f(ε_mk) · conj(ψ_mk(r_n)) = f(ε_mk) · e^{-i k·r_n} · conj(u_mk(r_n)).
#
# The smearing f selects the manifold: isolated bands use f ≡ 1; entangled cases use
# erfc ((μ, σ) as in pw2wannier90) or a Gaussian window.

using LinearAlgebra
using Printf

"""
    scdm_projections(model; dir, mode=:isolated, mu=0.0, sigma=1.0) -> A

Compute SCDM initial projections from the UNK files in `dir`. Returns `A`
(num_bands × num_wann × nkpt), ready to replace `model.A` (or write with [`write_amn`](@ref)).

`mode`:
- `:isolated`  — f ≡ 1 (num_bands == num_wann or a clean isolated group).
- `:erfc`      — f(ε) = erfc((ε − μ)/σ)/2, for entangled valence-like manifolds.
- `:gaussian`  — f(ε) = exp(−(ε − μ)²/σ²), for selecting bands around μ.

The pivot points are chosen at the k-point closest to Γ (the k-list must contain one).
`:erfc`/`:gaussian` require band energies (`model.eig`).
"""
function scdm_projections(model::Model; dir::AbstractString=".",
                          mode::Symbol=:isolated, mu::Float64=0.0, sigma::Float64=1.0)
    nb, nw, nk = model.num_bands, model.num_wann, nkpt(model.kgrid)
    ng0, _, _ = read_unk(joinpath(dir, @sprintf("UNK%05d.%1d", 1, 1)))
    getu = function (k)
        ngk, ikk, u = read_unk(joinpath(dir, @sprintf("UNK%05d.%1d", k, 1)))
        (ngk == ng0 && ikk == k) || error("UNK file $k: header mismatch")
        size(u, 2) >= nb || error("UNK has $(size(u, 2)) bands, model needs $nb")
        u
    end
    return scdm_amn(getu, ng0, model.kgrid.frac, nb, nw;
                    eig=model.eig, mode=mode, mu=mu, sigma=sigma)
end

"""
    scdm_amn(getu, ng, kfrac, num_bands, num_wann;
             eig=nothing, mode=:isolated, mu=0.0, sigma=1.0) -> A

Array-level SCDM core: `getu(k)` returns the periodic parts `u_mk(r)` on the `ng` real-space
grid as an `(npts × ≥num_bands)` matrix with x fastest (the UNK layout; `vec` of an
`ng₁×ng₂×ng₃` array). Used by [`scdm_projections`](@ref) (UNK files) and by the DFTK
extension (in-memory wavefunctions).
"""
function scdm_amn(getu::Function, ng::NTuple{3,Int}, kfrac::AbstractVector,
                  num_bands::Int, num_wann::Int;
                  eig::Union{Nothing,Matrix{Float64}}=nothing,
                  mode::Symbol=:isolated, mu::Float64=0.0, sigma::Float64=1.0)
    nb, nw, nk = num_bands, num_wann, length(kfrac)

    f = if mode === :isolated
        (m, k) -> 1.0
    elseif mode === :erfc
        eig !== nothing || error("scdm :erfc needs band energies")
        (m, k) -> 0.5 * erfc_((eig[m, k] - mu) / sigma)
    elseif mode === :gaussian
        eig !== nothing || error("scdm :gaussian needs band energies")
        (m, k) -> exp(-((eig[m, k] - mu) / sigma)^2)
    else
        error("scdm mode $mode (expected :isolated, :erfc, or :gaussian)")
    end

    # Anchor k-point: closest to Γ.
    kΓ = argmin([sum(abs2, k) for k in kfrac])
    sum(abs2, kfrac[kΓ]) < 1e-8 ||
        @warn "no Γ point in the k-list; using the closest k for SCDM pivots" k = kfrac[kΓ]

    uΓ = getu(kΓ)
    npts = prod(ng)

    # Column-pivoted QR on W[m, r] = f(ε_mΓ) ψ*_m(r). At Γ, ψ = u.
    W = Matrix{ComplexF64}(undef, nb, npts)
    for m in 1:nb
        fm = f(m, kΓ)
        @views W[m, :] .= fm .* conj.(uΓ[:, m])
    end
    piv = qr(W, ColumnNorm()).p[1:nw]

    # Pivot fractional coordinates (grid index n ↔ (n−1)/ng).
    rfrac = Vector{NTuple{3,Float64}}(undef, nw)
    for (j, p) in enumerate(piv)
        p0 = p - 1
        nx = p0 % ng[1]
        ny = (p0 ÷ ng[1]) % ng[2]
        nz = p0 ÷ (ng[1] * ng[2])
        rfrac[j] = (nx / ng[1], ny / ng[2], nz / ng[3])
    end

    A = Array{ComplexF64,3}(undef, nb, nw, nk)
    for k in 1:nk
        u = getu(k)
        kf = kfrac[k]
        for (n, p) in enumerate(piv)
            phase = cis(-TWOPI * (kf[1] * rfrac[n][1] + kf[2] * rfrac[n][2] + kf[3] * rfrac[n][3]))
            for m in 1:nb
                A[m, n, k] = f(m, k) * phase * conj(u[p, m])
            end
        end
    end
    return A
end

# erfc via the standard Abramowitz–Stegun 7.1.26 rational approximation (|err| < 1.5e-7),
# adequate for a smearing window and avoids a SpecialFunctions dependency.
function erfc_(x::Float64)
    z = abs(x)
    t = 1.0 / (1.0 + 0.5 * z)
    e = t * exp(-z^2 - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
        t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
        t * (-0.82215223 + t * 0.17087277)))))))))
    return x >= 0 ? e : 2.0 - e
end

# --------------------------------------------------------------------------------------
# scdm_auto — fit the SCDM erfc smearing (μ, σ) from a projectability-vs-energy curve,
# the Vitale et al. high-throughput protocol (npj Comput. Mater. 6, 66 (2020)).
#
# The projectability of Kohn–Sham state |ψ_mk⟩ onto the trial manifold is
# P_mk = ⟨ψ_mk|P̂|ψ_mk⟩ = [A_k (A_k†A_k)⁻¹ A_k†]_mm ∈ [0,1] (diagonal of the orthogonal
# projector built from the .amn columns — the same quantity the PDWF window uses, made
# gauge/normalisation-proof by the (A†A)⁻¹). The cloud {(ε_mk, P_mk)} follows a
# complementary error function; fitting it gives (μ_fit, σ_fit), and the SCDM smearing
# f(ε) = ½ erfc((ε−μ)/σ) is then set to μ = μ_fit − k·σ_fit, σ = σ_fit (k = 3 by default,
# as in aiida-wannier90-workflows), shifting the window down so the manifold is kept with
# weight ≈ 1.

"""
    scdm_auto(proj, eig; sigma_factor=3.0) -> (; mu, sigma, mu_fit, sigma_fit, rms)
    scdm_auto(A,    eig; sigma_factor=3.0) -> (; mu, sigma, mu_fit, sigma_fit, rms)

Fit the SCDM erfc smearing parameters from a projectability-vs-energy curve (Vitale et al.,
npj Comput. Mater. **6**, 66 (2020)), removing the last hand-set numbers from an SCDM-erfc
wannierisation.

The first form takes a projectability matrix `proj` and band energies `eig`, both
`num_bands × nkpt` (energies in eV). The second computes the projectability itself from an
`.amn` array `A` (`num_bands × num_wann × nkpt`) as the diagonal of the orthogonal projector
onto each k-point's trial columns, so it works for non-orthonormal trial orbitals.

Returns a named tuple: `mu`/`sigma` are the SCDM parameters ready to pass to
[`scdm_projections`](@ref) / `wannier_model` (`mu = mu_fit − sigma_factor·sigma_fit`),
`mu_fit`/`sigma_fit` are the raw erfc fit, and `rms` is the fit residual (a large value
warns that the projectability is not erfc-like — e.g. a manifold not separable in energy).

The fit is only as clean as the projectability. It is designed for **pseudo-atomic-orbital**
projectabilities — the ones Quantum ESPRESSO's `atomic_proj`/`projwfc` produce, and what the
Vitale protocol assumes — where `P(ε)` drops monotonically from ≈1 to ≈0. Crude trial
orbitals (e.g. a hydrogenic guess that misses the true radial shape, or a valence shell an
`n ≤ 3` hydrogenic cannot represent) can give a non-erfc cloud and a degenerate fit with a
large `rms`; treat that `rms` as the honest signal that this projectability is not a good
basis for the fit, not as a value to feed onward.
"""
function scdm_auto(proj::AbstractMatrix{<:Real}, eig::AbstractMatrix{<:Real};
                   sigma_factor::Real=3.0)
    size(proj) == size(eig) || error("scdm_auto: proj $(size(proj)) and eig $(size(eig)) differ")
    μf, σf, rms = _fit_erfc(vec(Float64.(eig)), vec(Float64.(proj)))
    return (; mu = μf - sigma_factor * σf, sigma = σf, mu_fit = μf, sigma_fit = σf, rms = rms)
end

function scdm_auto(A::AbstractArray{<:Complex,3}, eig::AbstractMatrix{<:Real}; kwargs...)
    nb, nw, nk = size(A)
    (nb, nk) == size(eig) || error("scdm_auto: A is $(nb)×$(nw)×$(nk); eig must be $(nb)×$(nk)")
    P = Matrix{Float64}(undef, nb, nk)
    for k in 1:nk
        Ak = A[:, :, k]
        # Projectability = diagonal of the orthogonal projector onto the trial column space.
        # Build it from the left singular vectors (Uᵣ Uᵣ†): ∈ [0,1] by construction and robust
        # to rank-deficient / near-dependent columns, unlike a Cholesky solve of A†A.
        F = svd(Ak)
        tol = maximum(F.S; init = 0.0) * max(size(Ak)...) * eps(Float64)
        r = count(>(tol), F.S)
        Q = @view F.U[:, 1:r]
        P[:, k] = clamp.(vec(sum(abs2, Q; dims = 2)), 0.0, 1.0)
    end
    return scdm_auto(P, eig; kwargs...)
end

# 2-parameter Levenberg–Marquardt fit of P(ε) = ½ erfc((ε−μ)/σ) to a scatter {(ε_i, P_i)}.
function _fit_erfc(ε::Vector{Float64}, p::Vector{Float64}; maxiter::Int=300)
    length(ε) >= 3 || error("scdm_auto: need at least 3 (energy, projectability) points")
    invsqrtπ = 1 / sqrt(π)

    # Initial guess: weight each point by p(1−p) (peaks in the transition region) to locate μ,
    # its spread to size σ. Falls back to the energy median/spread if there is no transition.
    w = p .* (1 .- p)
    sw = sum(w)
    μ, σ = if sw > 1e-6
        m = sum(w .* ε) / sw
        v = sum(w .* (ε .- m) .^ 2) / sw
        m, max(sqrt(v), 1e-3)
    else
        me = sum(ε) / length(ε)
        me, max(sqrt(sum((ε .- me) .^ 2) / length(ε)) / 4, 1e-3)
    end

    resid(μ, σ) = 0.5 .* erfc_.((ε .- μ) ./ σ) .- p
    cost(μ, σ) = sum(abs2, resid(μ, σ))

    λ = 1e-3
    c = cost(μ, σ)
    for _ in 1:maxiter
        t = (ε .- μ) ./ σ
        e = exp.(-t .^ 2)
        gμ = e .* (invsqrtπ / σ)              # ∂/∂μ [½ erfc(t)]
        gσ = (t .* e) .* (invsqrtπ / σ)       # ∂/∂σ [½ erfc(t)]
        r = 0.5 .* erfc_.(t) .- p
        # normal equations JᵀJ δ = −Jᵀr  (2×2)
        a11 = sum(gμ .^ 2); a12 = sum(gμ .* gσ); a22 = sum(gσ .^ 2)
        b1 = sum(gμ .* r); b2 = sum(gσ .* r)
        improved = false
        for _ in 1:30
            d = (a11 + λ) * (a22 + λ) - a12^2
            abs(d) < 1e-30 && (λ *= 10; continue)
            δμ = -((a22 + λ) * b1 - a12 * b2) / d
            δσ = -(-a12 * b1 + (a11 + λ) * b2) / d
            μn, σn = μ + δμ, max(σ + δσ, 1e-4)
            cn = cost(μn, σn)
            if cn < c
                μ, σ, c = μn, σn, cn
                λ = max(λ / 3, 1e-12)
                improved = true
                break
            else
                λ *= 5
            end
        end
        improved || break
    end
    rms = sqrt(c / length(ε))
    return μ, σ, rms
end

"""
    write_amn(path, A; header)

Write projections in the `.amn` format (`m n k Re Im`, band index fastest).
"""
function write_amn(path::AbstractString, A::Array{ComplexF64,3};
                   header::AbstractString="written by WannierFunctions.jl (SCDM)")
    nb, nw, nk = size(A)
    open(path, "w") do io
        println(io, header)
        @printf(io, "%12d%12d%12d\n", nb, nk, nw)
        for k in 1:nk, n in 1:nw, m in 1:nb
            @printf(io, "%5d%5d%5d%18.12f%18.12f\n", m, n, k, real(A[m, n, k]), imag(A[m, n, k]))
        end
    end
    return path
end
