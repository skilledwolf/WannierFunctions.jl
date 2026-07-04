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
end

"""
    dis_windows(eig, win) -> WindowData

Select the outer (`dis_win_min/max`) and frozen (`dis_froz_min/max`) energy windows per k-point.
Assumes `eig[:,k]` ascending.
"""
function dis_windows(eig::Matrix{Float64}, win::WinInput, num_wann::Int)
    nb, nk = size(eig)
    frozen = win.dis_froz_max != -Inf
    nfirstwin = zeros(Int, nk); ndimwin = zeros(Int, nk); ndimfroz = zeros(Int, nk)
    indxfroz = [Int[] for _ in 1:nk]
    indxnfroz = [Int[] for _ in 1:nk]
    lfrozen = [Bool[] for _ in 1:nk]
    for k in 1:nk
        e = @view eig[:, k]
        imin = findfirst(i -> e[i] >= win.dis_win_min && e[i] <= win.dis_win_max, 1:nb)
        imax = findlast(i -> e[i] <= win.dis_win_max && e[i] >= win.dis_win_min, 1:nb)
        (imin === nothing || imax === nothing) && error("empty outer window at k=$k")
        nd = imax - imin + 1
        nd >= num_wann || error("ndimwin ($nd) < num_wann ($num_wann) at k=$k")
        nfirstwin[k] = imin; ndimwin[k] = nd
        lf = falses(nd)
        if frozen
            for i in imin:imax
                (e[i] >= win.dis_froz_min && e[i] <= win.dis_froz_max) && (lf[i-imin+1] = true)
            end
        end
        nf = count(lf)
        nf <= num_wann || error("ndimfroz ($nf) > num_wann at k=$k")
        ndimfroz[k] = nf
        indxfroz[k] = findall(lf)
        indxnfroz[k] = findall(.!lf)
        lfrozen[k] = lf
    end
    return WindowData(nfirstwin, ndimwin, ndimfroz, indxfroz, indxnfroz, lfrozen, frozen)
end

"Window-slim the eigenvalues, projections A, and overlaps M to window-local indices."
function slim_data(eig, A, M, kpb, wd::WindowData)
    nb, nk = size(eig)
    nntot = size(M, 3)
    num_wann = size(A, 2)
    eigwin = [Float64[eig[wd.nfirstwin[k]+i-1, k] for i in 1:wd.ndimwin[k]] for k in 1:nk]
    Awin = Vector{Matrix{ComplexF64}}(undef, nk)
    for k in 1:nk
        f = wd.nfirstwin[k]
        Awin[k] = ComplexF64[A[f+i-1, j, k] for i in 1:wd.ndimwin[k], j in 1:num_wann]
    end
    # Mwin[k][nn] is ndimwin[k] × ndimwin[k2]
    Mwin = Vector{Vector{Matrix{ComplexF64}}}(undef, nk)
    for k in 1:nk
        fk = wd.nfirstwin[k]
        Mwin[k] = Vector{Matrix{ComplexF64}}(undef, nntot)
        for nn in 1:nntot
            k2 = kpb[nn, k]
            fk2 = wd.nfirstwin[k2]
            Mwin[k][nn] = ComplexF64[M[fk+i-1, fk2+j-1, nn, k]
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
            # take the nsel largest eigenvectors
            vecs = F.vectors[:, nd-nsel+1:nd]
            Unew = zeros(ComplexF64, nd, num_wann)
            for (col, l) in enumerate(nf+1:num_wann)
                Unew[:, l] = vecs[:, col]
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
    for nn in 1:bv.nntot
        w = bv.wb[nn, k]
        k2 = kpb[nn, k]
        cbw = Mwin_k[nn] * Uopt[k2]                      # ndimwin[k] × num_wann
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
    tot = 0.0
    for k in 1:nk
        s = num_wann * wbtot
        for nn in 1:bv.nntot
            k2 = kpb[nn, k]
            cww = Uopt[k]' * Mwin[k][nn] * Uopt[k2]      # num_wann × num_wann
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
    disentangle(model; dis_num_iter, dis_mix_ratio, dis_conv_tol=1e-10, dis_conv_window=3,
                verbose=false) -> DisentangleResult

Full SMV disentanglement: window selection, projection + frozen locking, Z-matrix subspace
iteration, subspace-Hamiltonian diagonalisation, and handoff to the localiser.
"""
function disentangle(model::Model, win::WinInput;
                     dis_num_iter::Int=win.dis_num_iter, dis_mix_ratio::Float64=win.dis_mix_ratio,
                     dis_conv_tol::Float64=1e-10, dis_conv_window::Int=3, verbose::Bool=false)
    model.eig !== nothing || error("disentanglement requires band energies (.eig)")
    nw = model.num_wann
    bv = model.bvectors
    kpb = bv.kpb
    nk = nkpt(model.kgrid)
    eig = sort(model.eig, dims=1)                      # ensure ascending per k

    wd = dis_windows(eig, win, nw)
    eigwin, Awin, Mwin = slim_data(eig, model.A, model.M, kpb, wd)
    Uopt = dis_project(Awin, wd, nw)

    wbtot = sum(@view bv.wb[:, 1])
    trace = Tuple{Int,Float64,Float64,Float64}[]
    history = fill(Inf, dis_conv_window)
    Zin = Vector{Matrix{ComplexF64}}(undef, nk)
    niter = 0
    for iter in 1:dis_num_iter
        niter = iter
        # Build/mix the Z matrices (only for k with non-frozen states).
        Zout = Vector{Matrix{ComplexF64}}(undef, nk)
        for k in 1:nk
            wd.ndimfroz[k] == nw && (Zout[k] = zeros(ComplexF64, 0, 0); continue)
            Zout[k] = Matrix(zmatrix(Mwin[k], Uopt, kpb, wd, bv, k, nw))
        end
        if iter == 1
            Zin = Zout
        else
            for k in 1:nk
                isempty(Zin[k]) && continue
                Zin[k] = dis_mix_ratio .* Zout[k] .+ (1 - dis_mix_ratio) .* Zin[k]
                Zin[k] = (Zin[k] + Zin[k]') ./ 2
            end
        end

        # Ω_I(i-1): num_wann·wbtot per k, minus frozen overlap, minus the largest Z eigenvalues.
        # The frozen contribution is evaluated for ALL k against the previous-iteration neighbour
        # subspaces BEFORE any k is updated (matches the reference's Ω_I(i-1) definition).
        wk = fill(nw * wbtot, nk)
        for k in 1:nk
            nf = wd.ndimfroz[k]
            nf == 0 && continue
            for nn in 1:bv.nntot
                k2 = kpb[nn, k]
                cww = Uopt[k]' * Mwin[k][nn] * Uopt[k2]
                for m in 1:nf, n in 1:nw
                    wk[k] -= bv.wb[nn, k] * abs2(cww[m, n])
                end
            end
        end
        # Now diagonalise Z (fixed from iteration start), subtract eigenvalues, update U_opt.
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

    # Post-iteration: diagonalise H in the optimal subspace, rotate U_opt to the eigenbasis.
    eigval_opt = Matrix{Float64}(undef, nw, nk)
    for k in 1:nk
        Hsub = Uopt[k]' * Diagonal(eigwin[k]) * Uopt[k]        # num_wann × num_wann
        F = eigen(Hermitian(Hsub))                             # ascending
        eigval_opt[:, k] = F.values
        Uopt[k] = Uopt[k] * F.vectors
    end

    # Handoff: num_wann×num_wann overlaps + initial square gauge from ⟨ψ̃|g⟩.
    Mrot0 = Array{ComplexF64,4}(undef, nw, nw, bv.nntot, nk)
    U0 = Array{ComplexF64,3}(undef, nw, nw, nk)
    for k in 1:nk
        caa = Uopt[k]' * Awin[k]                               # num_wann × num_wann
        U0[:, :, k] = svd_orthonormalize(caa)
    end
    for k in 1:nk, nn in 1:bv.nntot
        k2 = kpb[nn, k]
        Mrot0[:, :, nn, k] = Uopt[k]' * Mwin[k][nn] * Uopt[k2]
    end
    # Fold the initial square gauge into the overlaps for the localiser.
    Mrot0 = rotate_overlaps(Mrot0, U0, kpb)

    return DisentangleResult(U0, Mrot0, eigval_opt, Uopt, omega_I, trace, niter)
end
