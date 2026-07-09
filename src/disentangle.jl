# Disentanglement — Souza–Marzari–Vanderbilt subspace selection (PRB 65, 035109 (2001)).
#
# For num_bands > num_wann, pick a num_wann-dimensional subspace at every k that minimises the
# gauge-invariant spread Ω_I. Produces, per k, a rectangular embedding U_opt (ndimwin × num_wann)
# and hands off num_wann×num_wann overlaps + an initial square gauge to the MV localisation.
# Conventions follow src/disentangle.F90; see docs/reference-notes/disentanglement.md.
#
# Ragged storage: the outer-window dimension ndimwin varies per k, so per-k quantities are held as
# Vectors of Matrices rather than dense 4-D arrays.

using LinearAlgebra
using StaticArrays

"Per-k energy-window bookkeeping (all indices window-local unless noted)."
struct WindowData
    nfirstwin::Vector{Int}          # global band index of window bottom (imin), per k
    ndimwin::Vector{Int}            # window dimension, per k
    ndimfroz::Vector{Int}           # number of frozen bands, per k
    indxfroz::Vector{Vector{Int}}   # window-local indices of frozen bands, per k
    indxnfroz::Vector{Vector{Int}}  # window-local indices of non-frozen bands, per k
    lfrozen::Vector{Vector{Bool}}   # frozen mask (length ndimwin), per k
    frozen::Bool                    # any frozen states anywhere
    winbands::Vector{Vector{Int}}   # absolute band indices in the window, per k (length ndimwin)
end
# Contiguous-window compat constructor: winbands = nfirstwin : nfirstwin+ndimwin-1
WindowData(nfirstwin, ndimwin, ndimfroz, indxfroz, indxnfroz, lfrozen, frozen) =
    WindowData(nfirstwin, ndimwin, ndimfroz, indxfroz, indxnfroz, lfrozen, frozen,
               [collect(nfirstwin[k]:nfirstwin[k]+ndimwin[k]-1) for k in 1:length(nfirstwin)])

"""
    dis_windows(eig, num_wann; win_min=-Inf, win_max=Inf, froz_min=-Inf, froz_max=nothing)
        -> WindowData

Select the outer (`win_min`/`win_max`) and frozen (`froz_min`/`froz_max`) energy windows per
k-point (all in eV). A frozen window exists iff `froz_max !== nothing`. Assumes `eig[:,k]`
ascending.
"""
function dis_windows(eig::Matrix{Float64}, num_wann::Int;
                     win_min::Float64=-Inf, win_max::Float64=Inf,
                     froz_min::Float64=-Inf, froz_max::Union{Nothing,Float64}=nothing,
                     spheres::Union{Nothing,Vector{<:Tuple}}=nothing,
                     kfrac::Union{Nothing,Vector}=nothing,
                     recip::Union{Nothing,AbstractMatrix}=nothing,
                     sphere_first_wann::Int=1)
    nb, nk = size(eig)
    frozen = froz_max !== nothing
    nfirstwin = zeros(Int, nk); ndimwin = zeros(Int, nk); ndimfroz = zeros(Int, nk)
    indxfroz = [Int[] for _ in 1:nk]
    indxnfroz = [Int[] for _ in 1:nk]
    lfrozen = [Bool[] for _ in 1:nk]
    winbands = [Int[] for _ in 1:nk]
    for k in 1:nk
        e = @view eig[:, k]
        imin = findfirst(i -> e[i] >= win_min && e[i] <= win_max, 1:nb)
        imax = findlast(i -> e[i] <= win_max && e[i] >= win_min, 1:nb)
        (imin === nothing || imax === nothing) && error("empty outer window at k=$k")
        # dis_spheres: outside every sphere, replace the window with the fixed absolute band
        # range [first_wann, first_wann+num_wann-1] (no disentanglement freedom).
        if spheres !== nothing && !_in_any_sphere(kfrac[k], spheres, recip)
            imin = sphere_first_wann
            imax = sphere_first_wann + num_wann - 1
        end
        nd = imax - imin + 1
        nd >= num_wann || error("ndimwin ($nd) < num_wann ($num_wann) at k=$k")
        nfirstwin[k] = imin; ndimwin[k] = nd
        winbands[k] = collect(imin:imax)
        lf = falses(nd)
        if frozen
            for i in imin:imax
                (e[i] >= froz_min && e[i] <= froz_max) && (lf[i-imin+1] = true)
            end
        end
        nf = count(lf)
        nf <= num_wann || error("ndimfroz ($nf) > num_wann at k=$k")
        ndimfroz[k] = nf
        indxfroz[k] = findall(lf)
        indxnfroz[k] = findall(.!lf)
        lfrozen[k] = lf
    end
    return WindowData(nfirstwin, ndimwin, ndimfroz, indxfroz, indxnfroz, lfrozen, frozen, winbands)
end

"""
    dis_windows_proj(eig, A, num_wann; win_min, win_max, froz_min=-Inf, froz_max=nothing,
                     proj_min, proj_max) -> WindowData

PDWF (projectability-disentangled) window selection (Qiao–Pizzi–Marzari 2023): partition the
bands inside the outer energy window by projectability `p_nk = Σ_w |A_nw(k)|²`. States with
`p ≥ proj_max` (or in the energy frozen window, if given) are frozen; `proj_min ≤ p < proj_max`
are the disentanglement pool; `p < proj_min` are discarded from the window entirely (so the
window is generally non-contiguous). Expects (quasi-)orthonormal `.amn` columns (0 ≤ p ≤ 1).
"""
function dis_windows_proj(eig::Matrix{Float64}, A::Array{ComplexF64,3}, num_wann::Int;
                          win_min::Float64=-Inf, win_max::Float64=Inf,
                          froz_min::Float64=-Inf, froz_max::Union{Nothing,Float64}=nothing,
                          proj_min::Float64=0.0, proj_max::Float64=1.0)
    nb, nk = size(eig)
    frozen = froz_max !== nothing
    nfirstwin = zeros(Int, nk); ndimwin = zeros(Int, nk); ndimfroz = zeros(Int, nk)
    indxfroz = [Int[] for _ in 1:nk]
    indxnfroz = [Int[] for _ in 1:nk]
    lfrozen = [Bool[] for _ in 1:nk]
    winbands = [Int[] for _ in 1:nk]
    for k in 1:nk
        e = @view eig[:, k]
        bands = Int[]; lf = Bool[]
        for i in 1:nb
            (e[i] < win_min || e[i] > win_max) && continue          # discard: outside energy win
            p = sum(abs2(A[i, w, k]) for w in 1:num_wann)
            # Non-orthogonal atomic projectors (overlapping orbitals on different sites) can
            # push the projectability slightly above 1; tolerate a few percent and clamp.
            (p < -1e-8 || p > 1.05) &&
                error("dis_windows_proj: projectability $p ∉ [0,1] at band $i, k=$k — " *
                      ".amn columns far from orthonormal")
            p = clamp(p, 0.0, 1.0)
            isfroz = p >= proj_max || (frozen && froz_min <= e[i] <= froz_max)
            if isfroz
                push!(bands, i); push!(lf, true)
            elseif p >= proj_min                                    # disentangle pool
                push!(bands, i); push!(lf, false)
            end                                                     # else discard
        end
        nd = length(bands)
        nd >= num_wann || error("PDWF ndimwin ($nd) < num_wann ($num_wann) at k=$k")
        nf = count(lf)
        nf <= num_wann || error("PDWF ndimfroz ($nf) > num_wann at k=$k")
        nfirstwin[k] = isempty(bands) ? 1 : bands[1]
        ndimwin[k] = nd
        ndimfroz[k] = nf
        winbands[k] = bands
        lfrozen[k] = lf
        indxfroz[k] = findall(lf)
        indxnfroz[k] = findall(.!lf)
    end
    # PDWF freezes by projectability even with no energy frozen window; the frozen-locking
    # path must run whenever ANY state is frozen anywhere.
    any_frozen = frozen || any(>(0), ndimfroz)
    return WindowData(nfirstwin, ndimwin, ndimfroz, indxfroz, indxnfroz, lfrozen, any_frozen, winbands)
end

"True iff fractional k is inside any (centre_frac, radius) sphere (2π reciprocal metric)."
function _in_any_sphere(kf, spheres::Vector{<:Tuple}, recip::AbstractMatrix)
    kfv = SVector{3,Float64}(kf...)
    for (c, r) in spheres
        df = kfv .- SVector{3,Float64}(c...)
        dff = round.(df) .- df                     # nearest periodic image (anint(df) − df)
        dk = recip * dff                           # Cartesian Å⁻¹ (recip columns = b_i, 2π)
        dot(dk, dk) < r^2 && return true
    end
    return false
end

"Window-slim the eigenvalues, projections A, and overlaps M to window-local indices."
function slim_data(eig, A, M, kpb, wd::WindowData)
    nb, nk = size(eig)
    nntot = size(M, 3)
    num_wann = size(A, 2)
    # winbands[k] holds the absolute band indices in the window (contiguous for the energy
    # path, possibly non-contiguous for PDWF where mid-window low-projectability states are
    # discarded); all window-local extraction goes through it.
    wb = wd.winbands
    eigwin = [Float64[eig[wb[k][i], k] for i in 1:wd.ndimwin[k]] for k in 1:nk]
    Awin = Vector{Matrix{ComplexF64}}(undef, nk)
    for k in 1:nk
        Awin[k] = ComplexF64[A[wb[k][i], j, k] for i in 1:wd.ndimwin[k], j in 1:num_wann]
    end
    # Mwin[k][nn] is ndimwin[k] × ndimwin[k2]
    Mwin = Vector{Vector{Matrix{ComplexF64}}}(undef, nk)
    for k in 1:nk
        Mwin[k] = Vector{Matrix{ComplexF64}}(undef, nntot)
        for nn in 1:nntot
            k2 = kpb[nn, k]
            Mwin[k][nn] = ComplexF64[M[wb[k][i], wb[k2][j], nn, k]
                                     for i in 1:wd.ndimwin[k], j in 1:wd.ndimwin[k2]]
        end
    end
    return eigwin, Awin, Mwin
end

"SVD/polar orthonormalisation of a rectangular projection block: U = Z·V† (drop singular values)."
function svd_orthonormalize(A::AbstractMatrix{ComplexF64})
    F = svd(A)
    return F.U * F.Vt
end

"""
    dis_project(Awin, wd) -> U_opt

Initial optimal-subspace embedding per k from the trial projections (SMV Sec. III.D), then lock
the frozen states (SMV Eq. 27) when a frozen window is present.
"""
function dis_project(Awin::Vector{Matrix{ComplexF64}}, wd::WindowData, num_wann::Int)
    nk = length(Awin)
    Uopt = Vector{Matrix{ComplexF64}}(undef, nk)
    for k in 1:nk
        Uopt[k] = svd_orthonormalize(Awin[k])           # ndimwin × num_wann, orthonormal columns
    end
    wd.frozen && dis_proj_froz!(Uopt, wd, num_wann)
    return Uopt
end

"Lock frozen bands into U_opt (SMV Eq. 27): project the trial subspace onto the non-frozen space."
function dis_proj_froz!(Uopt, wd::WindowData, num_wann::Int)
    for k in 1:length(Uopt)
        nf = wd.ndimfroz[k]
        nf == 0 && continue
        nd = wd.ndimwin[k]
        if num_wann > nf
            U = Uopt[k]
            Ps = U * U'                                  # projector onto trial subspace
            Q = Diagonal([wd.lfrozen[k][n] ? 0.0 : 1.0 for n in 1:nd])
            CQPQ = Hermitian(Q * Ps * Q)
            F = eigen(CQPQ)                              # ascending eigenvalues
            nsel = num_wann - nf
            # Take the nsel largest eigenvectors of CQPQ = Q Pₛ Q (support on the non-frozen
            # subspace; frozen directions are eigenvalue-0).
            all(v -> -1e-8 <= v <= 1 + 1e-8, F.values) ||
                error("dis_proj_froz: QPQ eigenvalues outside [0,1] at k=$k")
            # Ortho-fix (reference default, disentangle.F90:2123-2227): when a *required*
            # eigenvalue is ≈0 it is degenerate with the frozen null space of QPQ, and the
            # returned eigenvector may point into the frozen span. Re-select those vectors by
            # explicit orthogonality to the frozen states (a frozen state m is the unit vector
            # e_{indxfroz[m]}, so orthogonality is |v[ifz]| ≤ eps8 for every frozen index).
            nzero = count(j -> F.values[j] < 1e-8, nd-nsel+1:nd)
            Unew = zeros(ComplexF64, nd, num_wann)
            if nzero == 0
                for (col, l) in enumerate(nf+1:num_wann)          # ascending, as the reference
                    Unew[:, l] = F.vectors[:, nd-nsel+col]
                end
            else
                goods = nsel - nzero
                vmap = Int[nd - c + 1 for c in 1:goods]           # top eigenvectors, descending
                for _ in 1:nzero
                    found = 0
                    for v in nd:-1:1
                        v in vmap && continue
                        if all(ifz -> abs(F.vectors[ifz, v]) <= 1e-8, wd.indxfroz[k])
                            found = v
                            break
                        end
                    end
                    found == 0 && error("dis_proj_froz: ortho-fix failed to find enough " *
                                        "frozen-orthogonal eigenvectors at k=$k")
                    push!(vmap, found)
                end
                for (col, l) in enumerate(nf+1:num_wann)
                    Unew[:, l] = F.vectors[:, vmap[col]]
                end
            end
        else
            Unew = zeros(ComplexF64, nd, num_wann)
        end
        # frozen columns 1..nf: unit vectors on the frozen bands
        for (l, ifz) in enumerate(wd.indxfroz[k])
            Unew[:, l] .= 0
            Unew[ifz, l] = 1
        end
        Uopt[k] = Unew
    end
    return Uopt
end

"Z-matrix (SMV Eq. 21) restricted to the non-frozen rows, for a single k."
function zmatrix(Mwin_k::Vector{Matrix{ComplexF64}}, Uopt, kpb, wd::WindowData, bv::BVectors,
                 k::Int, num_wann::Int)
    nfz = wd.indxnfroz[k]
    ndimk = length(nfz)
    Z = zeros(ComplexF64, ndimk, ndimk)
    cbw = Matrix{ComplexF64}(undef, size(Mwin_k[1], 1), num_wann)   # ndimwin[k] × num_wann
    for nn in 1:bv.nntot
        w = bv.wb[nn, k]
        k2 = kpb[nn, k]
        mul!(cbw, Mwin_k[nn], Uopt[k2])
        for n in 1:ndimk
            q = nfz[n]
            for m in 1:n
                p = nfz[m]
                csum = zero(ComplexF64)
                @inbounds for l in 1:num_wann
                    csum += cbw[p, l] * conj(cbw[q, l])
                end
                Z[m, n] += w * csum
            end
        end
    end
    for n in 1:ndimk, m in 1:n-1
        Z[n, m] = conj(Z[m, n])
    end
    return Hermitian(Z)
end

"Ω_I of the current subspaces: (1/N_k) Σ_k Σ_b w_b (num_wann − Σ_{mn}|⟨w_m,k|w_n,k+b⟩|²)."
function omega_invariant(Mwin, Uopt, kpb, bv::BVectors, num_wann::Int)
    nk = length(Uopt)
    wbtot = sum(@view bv.wb[:, 1])
    ndmax = maximum(size(u, 1) for u in Uopt)
    t1 = Matrix{ComplexF64}(undef, num_wann, ndmax)
    cww = Matrix{ComplexF64}(undef, num_wann, num_wann)
    tot = 0.0
    for k in 1:nk
        s = num_wann * wbtot
        for nn in 1:bv.nntot
            k2 = kpb[nn, k]
            t1v = @view t1[:, 1:size(Uopt[k2], 1)]
            mul!(t1v, Uopt[k]', Mwin[k][nn])
            mul!(cww, t1v, Uopt[k2])                     # num_wann × num_wann
            s -= bv.wb[nn, k] * sum(abs2, cww)
        end
        tot += s
    end
    return tot / nk
end

"Result of disentanglement, ready for MV localisation and interpolation."
struct DisentangleResult
    U0::Array{ComplexF64,3}        # initial square gauge (num_wann × num_wann × nkpt)
    Mrot0::Array{ComplexF64,4}     # num_wann×num_wann gauge overlaps for the localiser
    eigval_opt::Matrix{Float64}    # subspace eigenvalues (num_wann × nkpt), for interpolation
    Uopt::Vector{Matrix{ComplexF64}}  # optimal-subspace embeddings (ndimwin × num_wann) per k
    omega_I::Float64               # converged Ω_I (invariant through localisation)
    omega_I_trace::Vector{Tuple{Int,Float64,Float64,Float64}}  # (iter, ΩI(i-1), ΩI(i), Δ)
    niter::Int
end

"""
    disentangle(model; win_min=-Inf, win_max=Inf, froz_min=-Inf, froz_max=nothing,
                num_iter=200, mix_ratio=0.5, conv_tol=1e-10, conv_window=3,
                verbose=false) -> DisentangleResult

Full SMV disentanglement: window selection, projection + frozen locking, Z-matrix subspace
iteration, subspace-Hamiltonian diagonalisation, and handoff to the localiser. Windows are in eV;
a frozen window exists iff `froz_max` is given. The defaults (no windows) disentangle over all
`num_bands` states.

The `disentangle(model, win::WinInput)` method pulls all of these from a parsed `.win` instead.
"""
function disentangle(model::Model;
                     win_min::Float64=-Inf, win_max::Float64=Inf,
                     froz_min::Float64=-Inf, froz_max::Union{Nothing,Float64}=nothing,
                     num_iter::Int=200, mix_ratio::Float64=0.5,
                     conv_tol::Float64=1e-10, conv_window::Int=3, verbose::Bool=false,
                     spheres::Union{Nothing,Vector{<:Tuple}}=nothing, sphere_first_wann::Int=1,
                     froz_proj::Bool=false, proj_min::Float64=0.0, proj_max::Float64=1.0,
                     sitesym=nothing, gamma::Bool=false)
    model.eig !== nothing || error("disentanglement requires band energies (.eig)")
    nw = model.num_wann
    bv = model.bvectors
    kpb = bv.kpb
    nk = nkpt(model.kgrid)
    dis_num_iter, dis_mix_ratio = num_iter, mix_ratio
    dis_conv_tol, dis_conv_window = conv_tol, conv_window
    # The reference assumes the .eig band order is ascending per k (window selection relies on it).
    # Do NOT sort here: sorting the energies without also permuting the A/M rows would desync them.
    eig = model.eig
    for k in 1:nk
        issorted(@view eig[:, k]) ||
            @warn "eig at k=$k is not ascending; disentanglement windowing assumes ascending band " *
                  "order (matching the reference). Tiny inversions among degenerate bands are benign." maxlog=1
    end

    wd = if froz_proj
        dis_windows_proj(eig, model.A, nw; win_min=win_min, win_max=win_max,
                         froz_min=froz_min, froz_max=froz_max,
                         proj_min=proj_min, proj_max=proj_max)
    else
        dis_windows(eig, nw; win_min=win_min, win_max=win_max,
                    froz_min=froz_min, froz_max=froz_max,
                    spheres=spheres, kfrac=(spheres === nothing ? nothing : model.kgrid.frac),
                    recip=(spheres === nothing ? nothing : model.lattice.B),
                    sphere_first_wann=sphere_first_wann)
    end
    eigwin, Awin, Mwin = slim_data(eig, model.A, model.M, kpb, wd)
    Uopt = dis_project(Awin, wd, nw)
    # Symmetry-adapted mode (site_symmetry with num_bands > num_wann): the reference's
    # constrained Ω_I minimiser. Requires no frozen states (matches the reference's
    # 'not implemented in symmetry-adapted mode' guard) and symmetric initial embeddings.
    if sitesym !== nothing
        all(==(0), wd.ndimfroz) ||
            error("site_symmetry: a frozen window is not implemented in symmetry-adapted " *
                  "mode (matches the reference)")
        _symmetrize_uopt!(Uopt, sitesym, wd)
    end

    wbtot = sum(@view bv.wb[:, 1])
    trace = Tuple{Int,Float64,Float64,Float64}[]
    history = fill(Inf, dis_conv_window)
    Zin = Vector{Matrix{ComplexF64}}(undef, nk)
    # gemm scratch for the per-iteration frozen-overlap Ω_I sums
    froz_t1 = Matrix{ComplexF64}(undef, nw, maximum(wd.ndimwin))
    froz_cww = Matrix{ComplexF64}(undef, nw, nw)
    niter = 0
    for iter in 1:dis_num_iter
        niter = iter
        # Build/mix the Z matrices (only for k with non-frozen states).
        Zout = Vector{Matrix{ComplexF64}}(undef, nk)
        @maybe_threads (nk >= THREAD_MIN) for k in 1:nk
            if wd.ndimfroz[k] == nw
                Zout[k] = zeros(ComplexF64, 0, 0)
            else
                Zout[k] = Matrix(zmatrix(Mwin[k], Uopt, kpb, wd, bv, k, nw))
            end
        end
        sitesym === nothing || symmetrize_zmatrix!(Zout, sitesym)
        if iter == 1
            Zin = Zout
        else
            for k in 1:nk
                isempty(Zin[k]) && continue
                # In symmetry-adapted mode only the irreducible representatives are mixed
                # (only they feed the constrained update); the rest are never read.
                sitesym === nothing || sitesym.ir2ik[sitesym.ik2ir[k]] == k || continue
                Zin[k] .= dis_mix_ratio .* Zout[k] .+ (1 - dis_mix_ratio) .* Zin[k]
                Zout[k] .= Zin[k]'          # Zout[k] is rebuilt next iteration — free scratch
                Zin[k] .= (Zin[k] .+ Zout[k]) ./ 2
            end
        end

        local womegai1
        if sitesym !== nothing
            # Constrained Ω_I update at each irreducible representative; each carries its
            # star's weight nsym/|G_k|, and −tr Re λ replaces the Z-eigenvalue sum.
            wk = zeros(nk)
            for ir in 1:sitesym.nkptirr
                ik = sitesym.ir2ik[ir]
                ngk = count(==(ik), @view sitesym.kptsym[:, ir])
                wk[ik] = nw * wbtot * sitesym.nsym / ngk
                λ = dis_extract_symmetry!(Uopt[ik], Zin[ik], sitesym, ik, wd.ndimwin[ik])
                wk[ik] -= sum(real(λ[j, j]) for j in 1:nw)
            end
            _symmetrize_uopt!(Uopt, sitesym, wd)       # re-project reps + propagate the stars
            womegai1 = sum(wk) / nk
        else
            # Ω_I(i-1): num_wann·wbtot per k, minus frozen overlap, minus the largest Z
            # eigenvalues. The frozen contribution is evaluated for ALL k against the
            # previous-iteration neighbour subspaces BEFORE any k is updated (matches the
            # reference's Ω_I(i-1) definition).
            wk = fill(nw * wbtot, nk)
            for k in 1:nk
                nf = wd.ndimfroz[k]
                nf == 0 && continue
                for nn in 1:bv.nntot
                    k2 = kpb[nn, k]
                    t1v = @view froz_t1[:, 1:wd.ndimwin[k2]]
                    mul!(t1v, Uopt[k]', Mwin[k][nn])
                    mul!(froz_cww, t1v, Uopt[k2])
                    for m in 1:nf, n in 1:nw
                        wk[k] -= bv.wb[nn, k] * abs2(froz_cww[m, n])
                    end
                end
            end
            # Diagonalise Z (fixed from iteration start), subtract eigenvalues, update U_opt.
            for k in 1:nk
                nf = wd.ndimfroz[k]
                nf < nw || continue
                F = eigen(Hermitian(Zin[k]))                # ascending
                nsel = nw - nf
                wk[k] -= sum(@view F.values[end-nsel+1:end])
                vecs = @view F.vectors[:, end-nsel+1:end]
                for (col, l) in enumerate(nf+1:nw)
                    Uopt[k][:, l] .= 0
                    for (i, ir) in enumerate(wd.indxnfroz[k])
                        Uopt[k][ir, l] = vecs[i, col]
                    end
                end
            end
            womegai1 = sum(wk) / nk
        end

        womegai = omega_invariant(Mwin, Uopt, kpb, bv, nw)
        delta = womegai1 / womegai - 1
        push!(trace, (iter, womegai1, womegai, delta))
        verbose && @info "DIS iter $iter" ΩI_prev=womegai1 ΩI=womegai Δ=delta

        history = circshift(history, -1); history[end] = delta
        if iter >= dis_conv_window && all(h -> abs(h) < dis_conv_tol, history)
            break
        end
    end

    omega_I = omega_invariant(Mwin, Uopt, kpb, bv, nw)

    # Γ-only: the optimal subspace of real Γ data is conjugation-closed — rotate the
    # embedding to a real orthonormal basis so the real-orthogonal localiser can take over.
    if gamma
        nk == 1 || error("gamma_only disentanglement requires a single k-point")
        Uopt[1] = ComplexF64.(realify_subspace(Uopt[1]))
    end

    # Post-iteration: diagonalise H in the optimal subspace, rotate U_opt to the eigenbasis.
    # In symmetry-adapted mode only the eigenvalues are taken — U_opt must keep its
    # star-relation U(Rk) = d·U(k)·D†, so the eigenbasis rotation is skipped (reference: the
    # ceamp replacement is guarded by `.not. lsitesymmetry`).
    eigval_opt = Matrix{Float64}(undef, nw, nk)
    for k in 1:nk
        Hsub = Uopt[k]' * Diagonal(eigwin[k]) * Uopt[k]        # num_wann × num_wann
        if gamma                                               # keep the embedding real
            Fr = eigen(Symmetric(real.(Hsub)))
            eigval_opt[:, k] = Fr.values
            Uopt[k] = Uopt[k] * Fr.vectors
        else
            F = eigen(Hermitian(Hsub))                         # ascending
            eigval_opt[:, k] = F.values
            sitesym === nothing && (Uopt[k] = Uopt[k] * F.vectors)
        end
    end

    # Handoff: num_wann×num_wann overlaps + initial square gauge from ⟨ψ̃|g⟩.
    Mrot0 = Array{ComplexF64,4}(undef, nw, nw, bv.nntot, nk)
    U0 = Array{ComplexF64,3}(undef, nw, nw, nk)
    if sitesym === nothing
        for k in 1:nk
            caa = Uopt[k]' * Awin[k]                           # num_wann × num_wann
            if gamma                                           # real Löwdin, real gauge
                Fr = svd(real.(caa))
                U0[:, :, k] = ComplexF64.(Fr.U * Fr.Vt)
            else
                U0[:, :, k] = svd_orthonormalize(caa)
            end
        end
    else
        # Square gauge at the representatives only, then symmetrise in the Wannier
        # representation (d_band := d_wann) and reconstruct the stars.
        stloc = replace_d_matrix_band(sitesym)
        for ir in 1:sitesym.nkptirr
            ik = sitesym.ir2ik[ir]
            U0[:, :, ik] = svd_orthonormalize(Uopt[ik]' * Awin[ik])
        end
        symmetrize_u!(U0, stloc, stloc.d_band, stloc.d_wann)
    end
    for k in 1:nk, nn in 1:bv.nntot
        k2 = kpb[nn, k]
        Mrot0[:, :, nn, k] = Uopt[k]' * Mwin[k][nn] * Uopt[k2]
    end
    # Fold the initial square gauge into the overlaps for the localiser.
    Mrot0 = rotate_overlaps(Mrot0, U0, kpb)

    return DisentangleResult(U0, Mrot0, eigval_opt, Uopt, omega_I, trace, niter)
end

"Symmetrise the disentanglement subspace embeddings (uniform window) using the band
representation d_band and the Wannier representation d_wann of a `Sitesym`."
function _symmetrize_uopt!(Uopt::Vector{Matrix{ComplexF64}}, sitesym, wd::WindowData)
    nd = wd.ndimwin[1]
    all(==(nd), wd.ndimwin) ||
        error("site_symmetry disentanglement currently requires a uniform window")
    nw = size(Uopt[1], 2)
    nk = length(Uopt)
    arr = Array{ComplexF64,3}(undef, nd, nw, nk)
    for k in 1:nk
        arr[:, :, k] = Uopt[k]
    end
    symmetrize_u!(arr, sitesym, sitesym.d_band, sitesym.d_wann; n=nd)
    for k in 1:nk
        Uopt[k] = arr[:, :, k]
    end
    return Uopt
end

"Compat method: pull the windows and iteration controls from a parsed `.win`."
function disentangle(model::Model, win::WinInput; kwargs...)
    # dis_spheres block: `kx ky kz radius` rows (fractional k, radius Å⁻¹ 2π convention)
    spheres = nothing
    sphere_first_wann = _getint(win.raw, "dis_spheres_first_wann", 1)
    if _getint(win.raw, "dis_spheres_num", 0) > 0 && haskey(win.blocks, "dis_spheres")
        spheres = Tuple{SVector{3,Float64},Float64}[]
        for ln in win.blocks["dis_spheres"]
            t = split(ln)
            length(t) >= 4 || continue
            push!(spheres, (SVector{3,Float64}(parse_f64.(t[1:3])...), parse_f64(t[4])))
        end
    end
    # PDWF: projectability frozen window
    froz_proj = _getbool(win.raw, "dis_froz_proj", false)
    proj_min = _getfloat(win.raw, "dis_proj_min", 0.0)
    proj_max = _getfloat(win.raw, "dis_proj_max", 1.0)
    return disentangle(model;
        win_min=win.dis_win_min, win_max=win.dis_win_max,
        froz_min=win.dis_froz_min,
        froz_max=(win.dis_froz_max == -Inf ? nothing : win.dis_froz_max),
        num_iter=win.dis_num_iter, mix_ratio=win.dis_mix_ratio,
        conv_tol=_getfloat(win.raw, "dis_conv_tol", 1e-10),
        conv_window=_getint(win.raw, "dis_conv_window", 3),
        spheres=spheres, sphere_first_wann=sphere_first_wann,
        froz_proj=froz_proj, proj_min=proj_min, proj_max=proj_max, kwargs...)
end
