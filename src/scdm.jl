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
