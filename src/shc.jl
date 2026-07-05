# Spin Hall conductivity — Qiao–Zhou–Yuan–Zhao method (PRB 98, 214402 (2018) [QZYZ18]) as in
# postw90's berry_get_shc_klist / get_SHC_R; exact conventions in
# docs/reference-notes/shc-dos-boltz-kslice.md §1. Computes the k-resolved Berry-curvature-like
# term of σ^{spin γ}_{αβ} and its Fermi-scan integral in (ħ/e)·S/cm.

using LinearAlgebra
using StaticArrays

"""
    read_spn(path; num_bands, num_kpts) -> spn

Read a formatted `.spn`: `spn[n, m, ik, s] = ⟨ψ_n|σ_s|ψ_m⟩` (s = x,y,z), stored lower-triangular
(n ≤ m) with Hermitian completion.
"""
function read_spn(path::AbstractString; num_bands::Int, num_kpts::Int)
    open(path, "r") do io
        readline(io)
        nb, nk = parse.(Int, split(readline(io)))
        (nb == num_bands && nk == num_kpts) ||
            error(".spn dims ($nb,$nk) ≠ model ($num_bands,$num_kpts)")
        spn = Array{ComplexF64,4}(undef, nb, nb, nk, 3)
        for ik in 1:nk, m in 1:nb, n in 1:m
            for s in 1:3
                t = split(readline(io))
                spn[n, m, ik, s] = complex(parse(Float64, t[1]), parse(Float64, t[2]))
                spn[m, n, ik, s] = conj(spn[n, m, ik, s])
            end
        end
        return spn
    end
end

"""
    ShcModel

`BerryModel` plus the QZYZ operators: SS(R) (σ, 3 comps), SH(R) (σH), SR(R) (σ(r−R)_α) and
SHR(R) (σH(r−R)_α), each for the chosen spin component γ built on demand — here all three γ
are stored (γ, α indices: [.,.,R,γ] and [.,.,R,γ,α]).
"""
struct ShcModel
    bm::BerryModel
    SSr::Array{ComplexF64,4}       # (nw, nw, nr, 3)
    SHr::Array{ComplexF64,4}       # (nw, nw, nr, 3)
    SRr::Array{ComplexF64,5}       # (nw, nw, nr, 3, 3)   [γ, α]
    SHRr::Array{ComplexF64,5}
end

Base.show(io::IO, ::MIME"text/plain", m::ShcModel) =
    print(io, "ShcModel: ", num_wann(m.bm), " WF, ", length(m.bm.irvec), " R-vectors (H, A, S, SH, SR, SHR)")

"""
    ShcModel(seedname; scissors_shift=0.0, num_valence_bands=0) -> ShcModel

Assemble from `seedname.{win,mmn,eig,chk,spn}` (Qiao method: no `.uHu` needed). A nonzero
`scissors_shift` (eV) rigidly shifts the ab-initio bands above `num_valence_bands` before
H(R) and σH are built (postw90's scissors correction — Qiao method only).
"""
function ShcModel(seedname::AbstractString; scissors_shift::Float64=0.0,
                  num_valence_bands::Int=0)
    win = read_win(seedname * ".win")
    chk = isfile(seedname * ".chk") ? read_chk(seedname * ".chk") :
          read_chk_fmt(seedname * ".chk.fmt")
    eig = read_eig(seedname * ".eig")
    scissors_shift != 0.0 && (eig[num_valence_bands+1:end, :] .+= scissors_shift)
    M, kpb, gpb, nb, nk, nntot = read_mmn(seedname * ".mmn")
    lattice = Lattice(win.unit_cell)
    kgrid = KGrid(win.kpoints, win.mp_grid)
    bv = build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol=win.kmesh_tol)
    bm = BerryModel(chk, eig, bv, kgrid, lattice; use_ws_distance=win.use_ws_distance)
    spn = read_spn(seedname * ".spn"; num_bands=nb, num_kpts=nk)

    nw = num_wann(chk)
    vs, winidx = gauge_v_windows(chk, nb)

    SSq = zeros(ComplexF64, nw, nw, nk, 3)
    SHq = zeros(ComplexF64, nw, nw, nk, 3)
    SRq = zeros(ComplexF64, nw, nw, nk, 3, 3)
    SHRq = zeros(ComplexF64, nw, nw, nk, 3, 3)
    for q in 1:nk
        wq = winidx[q]
        for s in 1:3
            σ = @view spn[:, :, q, s]
            σH = σ * Diagonal(@view eig[:, q])                # SH_o = σ·diag(ε)  (QZYZ Eq. 48)
            SSq[:, :, q, s] = vs[q]' * σ[wq, wq] * vs[q]
            SHq[:, :, q, s] = vs[q]' * σH[wq, wq] * vs[q]
            for b in 1:nntot
                qb = kpb[b, q]
                So = @view M[:, :, b, q]
                SM = (σ * So)[wq, winidx[qb]]                 # Eq. 50
                SHM = (σH * So)[wq, winidx[qb]]               # Eq. 51
                SMq = vs[q]' * SM * vs[qb]
                SHMq = vs[q]' * SHM * vs[qb]
                w = bv.wb[b, q]
                for α in 1:3
                    f = w * bv.bvec[α, b, q]
                    @views SRq[:, :, q, s, α] .+= f .* (SMq .- SSq[:, :, q, s])
                    @views SHRq[:, :, q, s, α] .+= f .* (SHMq .- SHq[:, :, q, s])
                end
            end
        end
    end

    nr = length(bm.irvec)
    SSr = zeros(ComplexF64, nw, nw, nr, 3)
    SHr = zeros(ComplexF64, nw, nw, nr, 3)
    SRr = zeros(ComplexF64, nw, nw, nr, 3, 3)
    SHRr = zeros(ComplexF64, nw, nw, nr, 3, 3)
    Threads.@threads for ir in 1:nr
        R = SVector{3,Float64}(bm.irvec[ir]...)
        for q in 1:nk
            fac = cis(-TWOPI * dot(kgrid.frac[q], R)) / nk
            @views for s in 1:3
                SSr[:, :, ir, s] .+= fac .* SSq[:, :, q, s]
                SHr[:, :, ir, s] .+= fac .* SHq[:, :, q, s]
                for α in 1:3
                    SRr[:, :, ir, s, α] .+= (im * fac) .* SRq[:, :, q, s, α]
                    SHRr[:, :, ir, s, α] .+= (im * fac) .* SHRq[:, :, q, s, α]
                end
            end
        end
    end
    return ShcModel(bm, SSr, SHr, SRr, SHRr)
end

"""
    read_shu(path; num_bands, num_kpts, nntot) -> X

Read a Fortran-unformatted `.sHu` or `.sIu` (identical layout):
`X[m, n, s, b, ik] = ⟨u_mk|σ_s H_k|u_n,k+b⟩` (`.sHu`) or `⟨u_mk|σ_s|u_n,k+b⟩` (`.sIu`),
with b running in `.mmn` neighbour order. On disk each record is a ket-fastest nb×nb matrix
that must be transposed (pw2wannier90 convention, as postw90 does on read).
"""
function read_shu(path::AbstractString; num_bands::Int, num_kpts::Int, nntot::Int)
    open(path, "r") do io
        _frec(io)                                            # header (60 chars)
        dims = Int.(_frec(io, Int32))
        dims == [num_bands, num_kpts, nntot] ||
            error("$(basename(path)) dims $dims ≠ model [$num_bands, $num_kpts, $nntot]")
        X = Array{ComplexF64,5}(undef, num_bands, num_bands, 3, nntot, num_kpts)
        for ik in 1:num_kpts, b in 1:nntot, s in 1:3
            M = reshape(Vector(_frec(io, ComplexF64)), num_bands, num_bands)
            X[:, :, s, b, ik] = transpose(M)
        end
        return X
    end
end

"""
    ShcRyooModel

`BerryModel` plus the Ryoo–Park–Souza operators (RPS19, PRB 99, 235113): S(R) (from `.spn`),
SAA(R) = ⟨0|σ_γ(r−R)_α|R⟩ (from `.sIu`) and SBB(R) = ⟨0|σ_γ H(r−R)_α|R⟩ (from `.sHu`).
"""
struct ShcRyooModel
    bm::BerryModel
    SSr::Array{ComplexF64,4}       # (nw, nw, nr, 3)
    SAAr::Array{ComplexF64,5}      # (nw, nw, nr, 3, 3)   [γ, α]
    SBBr::Array{ComplexF64,5}
end

Base.show(io::IO, ::MIME"text/plain", m::ShcRyooModel) =
    print(io, "ShcRyooModel: ", num_wann(m.bm), " WF, ", length(m.bm.irvec),
          " R-vectors (H, A, S, SAA, SBB)")

"""
    ShcRyooModel(seedname; transl_inv_full=false) -> ShcRyooModel

Assemble from `seedname.{win,mmn,eig,chk,spn,sHu,sIu}`. With `transl_inv_full = true` the
one-shell translation-invariant scheme is used for A(R), SAA(R) and SBB(R) (phases
e^{ib·r₀} and e^{−ib·R/2} with r₀ = (r̄_i+r̄_j)/2, plus the (r₀−R/2)·S / ·σH diagonal
corrections).

NB: the reference implementation accumulates the windowed σ-rotation across the three spin
components (get_oper.F90:2680), so its SAA/SBB for γ = y hold σx+σy and for γ = z hold
σx+σy+σz. We compute the clean per-component operators; results agree for γ = x (which is
what the shipped oracle tests use) and differ — deliberately — for γ = y, z.
"""
function ShcRyooModel(seedname::AbstractString; transl_inv_full::Bool=false)
    win = read_win(seedname * ".win")
    chk = isfile(seedname * ".chk") ? read_chk(seedname * ".chk") :
          read_chk_fmt(seedname * ".chk.fmt")
    eig = read_eig(seedname * ".eig")
    M, kpb, gpb, nb, nk, nntot = read_mmn(seedname * ".mmn")
    lattice = Lattice(win.unit_cell)
    kgrid = KGrid(win.kpoints, win.mp_grid)
    bv = build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol=win.kmesh_tol)
    bm = BerryModel(chk, eig, bv, kgrid, lattice; use_ws_distance=win.use_ws_distance,
                    transl_inv_full=transl_inv_full)
    spn = read_spn(seedname * ".spn"; num_bands=nb, num_kpts=nk)
    shu = read_shu(seedname * ".sHu"; num_bands=nb, num_kpts=nk, nntot=nntot)
    siu = read_shu(seedname * ".sIu"; num_bands=nb, num_kpts=nk, nntot=nntot)

    nw = num_wann(chk)
    vs, winidx = gauge_v_windows(chk, nb)

    # Original WS set (bm.irvec is the EXPANDED set when transl_inv_full && use_ws_distance)
    irvec, ndegen = wigner_seitz(lattice, kgrid.mp_grid)
    ws0 = win.use_ws_distance ?
          ws_translate_dist(irvec, chk.centres, lattice, kgrid.mp_grid) : nothing

    # S(q) and σH(q) (σH needed for the transl_inv_full correction of SBB)
    SSq = zeros(ComplexF64, nw, nw, nk, 3)
    SHq = zeros(ComplexF64, nw, nw, nk, 3)
    for q in 1:nk, s in 1:3
        σ = @view spn[:, :, q, s]
        σH = σ * Diagonal(@view eig[:, q])
        SSq[:, :, q, s] = vs[q]' * σ[winidx[q], winidx[q]] * vs[q]
        SHq[:, :, q, s] = vs[q]' * σH[winidx[q], winidx[q]] * vs[q]
    end
    nr = length(bm.irvec)
    SSr = zeros(ComplexF64, nw, nw, nr, 3)
    SHr = zeros(ComplexF64, nw, nw, nr, 3)
    if !transl_inv_full
        for s in 1:3
            SSr[:, :, :, s] = fourier_q_to_R((@view SSq[:, :, :, s]), kgrid, bm.irvec)
            SHr[:, :, :, s] = fourier_q_to_R((@view SHq[:, :, :, s]), kgrid, bm.irvec)
        end
    else
        # bm lives on the expanded set — scatter S/σH onto it (pre-divided convention)
        irvec90, idx90 = _pw90_rset(irvec, ws0)
        @assert irvec90 == bm.irvec
        for s in 1:3
            SSr[:, :, :, s] = _scatter_ws(fourier_q_to_R((@view SSq[:, :, :, s]), kgrid, irvec),
                                          irvec, ndegen, ws0, idx90, nr)
            SHr[:, :, :, s] = _scatter_ws(fourier_q_to_R((@view SHq[:, :, :, s]), kgrid, irvec),
                                          irvec, ndegen, ws0, idx90, nr)
        end
    end

    SAAr = _shu_to_R(siu, SSr, vs, winidx, bv, kpb, kgrid, bm, chk,
                     irvec, ndegen, ws0, transl_inv_full)
    SBBr = _shu_to_R(shu, SHr, vs, winidx, bv, kpb, kgrid, bm, chk,
                     irvec, ndegen, ws0, transl_inv_full)
    return ShcRyooModel(bm, SSr, SAAr, SBBr)
end

# SXX(R) from raw .sHu/.sIu matrices: gauge-rotate on the window, accumulate the b-sum with
# +i·w_b·b_α. With transl_inv_full: e^{ib·r₀} phases at accumulation, canonical b-slots
# scattered onto the expanded minimal-image R-set, per-slot e^{−ib·R̃/2}, then the
# (r₀ − R̃/2)·partner diagonal correction — all on the expanded set, matching the reference
# order (get_oper.F90:2710-2769).
function _shu_to_R(X::Array{ComplexF64,5}, partner_R::Array{ComplexF64,4},
                   vs, winidx, bv::BVectors, kpb::Matrix{Int}, kgrid::KGrid,
                   bm::BerryModel, chk::Checkpoint, irvec::Vector{NTuple{3,Int}},
                   ndegen::Vector{Int}, ws0, transl_inv_full::Bool)
    nw = size(partner_R, 1)
    nk = nkpt(kgrid)
    nr = length(bm.irvec)
    nntot = size(X, 4)
    out = zeros(ComplexF64, nw, nw, nr, 3, 3)
    if !transl_inv_full
        SXq = zeros(ComplexF64, nw, nw, nk, 3, 3)
        for q in 1:nk, b in 1:nntot
            qb = kpb[b, q]
            w = bv.wb[b, q]
            for s in 1:3
                core = vs[q]' * X[winidx[q], winidx[qb], s, b, q] * vs[qb]
                for a in 1:3
                    @views SXq[:, :, q, s, a] .+= (im * w * bv.bvec[a, b, q]) .* core
                end
            end
        end
        for s in 1:3, a in 1:3
            out[:, :, :, s, a] = fourier_q_to_R((@view SXq[:, :, :, s, a]), kgrid, bm.irvec)
        end
        return out
    end
    r0 = [(chk.centres[a, i] + chk.centres[a, j]) / 2 for i in 1:nw, j in 1:nw, a in 1:3]
    slot = _bvec_canonical_slots(bv)
    SXq_b = zeros(ComplexF64, nw, nw, nk, nntot, 3, 3)
    for q in 1:nk, b in 1:nntot
        qb = kpb[b, q]
        w = bv.wb[b, q]
        bvq = SVector(bv.bvec[1, b, q], bv.bvec[2, b, q], bv.bvec[3, b, q])
        b0 = slot[b, q]
        for s in 1:3
            core = vs[q]' * X[winidx[q], winidx[qb], s, b, q] * vs[qb]
            for j in 1:nw, i in 1:nw
                ph = cis(dot(bvq, SVector(r0[i, j, 1], r0[i, j, 2], r0[i, j, 3])))
                core[i, j] *= ph
            end
            for a in 1:3
                @views SXq_b[:, :, q, b0, s, a] .+= (im * w * bvq[a]) .* core
            end
        end
    end
    irvec90, idx90 = _pw90_rset(irvec, ws0)
    @assert irvec90 == bm.irvec
    for s in 1:3, a in 1:3, b0 in 1:nntot
        Xb = fourier_q_to_R((@view SXq_b[:, :, :, b0, s, a]), kgrid, irvec)
        Xb90 = _scatter_ws(Xb, irvec, ndegen, ws0, idx90, nr)
        bk1 = SVector(bv.bvec[1, b0, 1], bv.bvec[2, b0, 1], bv.bvec[3, b0, 1])
        for ir in 1:nr
            ph2 = cis(-dot(bk1, bm.Rcart[ir]) / 2)
            @views out[:, :, ir, s, a] .+= ph2 .* Xb90[:, :, ir]
        end
    end
    for ir in 1:nr, a in 1:3, s in 1:3, j in 1:nw, i in 1:nw
        out[i, j, ir, s, a] += (r0[i, j, a] - bm.Rcart[ir][a] / 2) * partner_R[i, j, ir, s]
    end
    return out
end

"Canonical b-slot map: slot[b,q] = index b₀ with bvec(:,b₀,1) == bvec(:,b,q) (tol 1e-7)."
function _bvec_canonical_slots(bv::BVectors)
    nntot, nk = size(bv.wb, 1), size(bv.wb, 2)
    slot = Matrix{Int}(undef, nntot, nk)
    for q in 1:nk, b in 1:nntot
        found = 0
        for b0 in 1:nntot
            if abs(bv.bvec[1, b0, 1] - bv.bvec[1, b, q]) < 1e-7 &&
               abs(bv.bvec[2, b0, 1] - bv.bvec[2, b, q]) < 1e-7 &&
               abs(bv.bvec[3, b0, 1] - bv.bvec[3, b, q]) < 1e-7
                found = b0
                break
            end
        end
        found == 0 && error("no canonical b-vector slot for b=$b at q=$q")
        slot[b, q] = found
    end
    return slot
end

"ws-aware operator interpolation O(k) from O(R) (shared by all SHC operators)."
function _ft_op(bm::BerryModel, Or::AbstractArray{ComplexF64,3}, kf::SVector{3,Float64})
    nw = size(Or, 1)
    O = zeros(ComplexF64, nw, nw)
    if bm.wsdist === nothing
        for ir in 1:length(bm.irvec)
            fac = cis(TWOPI * dot(kf, SVector{3,Float64}(bm.irvec[ir]...))) / bm.ndegen[ir]
            @views O .+= fac .* Or[:, :, ir]
        end
    else
        for ir in 1:length(bm.irvec)
            nd0 = bm.ndegen[ir]
            for j in 1:nw, i in 1:nw
                dl = bm.wsdist[i, j, ir]
                w = 1.0 / (nd0 * length(dl))
                o = Or[i, j, ir]
                for Rt in dl
                    O[i, j] += w * cis(TWOPI * dot(kf, SVector{3,Float64}(Rt...))) * o
                end
            end
        end
    end
    return O
end

"""
    shc_fermiscan(sm; fermi_energies, kmesh, γ=3, α=1, β=2, adaptive=true,
                  adpt_fac=√2, adpt_max=1.0, smr_width=0.0, eigval_max=Inf)
        -> Vector (per Fermi energy, (ħ/e)·S/cm)

Fermi-scan spin Hall conductivity σ^{spin γ}_{αβ} (QZYZ18 Berry-curvature-like term).
"""
function shc_fermiscan(sm::ShcModel; fermi_energies::Vector{Float64},
                       kmesh::NTuple{3,Int}=(25, 25, 25), γ::Int=3, α::Int=1, β::Int=2,
                       adaptive::Bool=true, adpt_fac::Float64=sqrt(2.0), adpt_max::Float64=1.0,
                       smr_width::Float64=0.0, eigval_max::Float64=Inf)
    bm = sm.bm
    nw = num_wann(bm)
    nf = length(fermi_energies)
    nktot = prod(kmesh)
    Δk = kmesh_spacing(bm.lattice, kmesh)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    per_k = Vector{Vector{Float64}}(undef, nktot)
    Threads.@threads for idx in 1:nktot
        per_k[idx] = _shc_kpoint(sm, kl[idx], fermi_energies, Δk, γ, α, β,
                                 adaptive, adpt_fac, adpt_max, smr_width, eigval_max)
    end
    fac = 1.0e8 * ELEM_CHARGE_SI^2 / (HBAR_SI * cell_volume(bm.lattice)) / 2.0
    return (fac / nktot) .* sum(per_k)
end

function _shc_kpoint(sm::ShcModel, kf::SVector{3,Float64}, efs::Vector{Float64}, Δk::Float64,
                     γ::Int, α::Int, β::Int, adaptive::Bool, adpt_fac::Float64,
                     adpt_max::Float64, smr_width::Float64, eigval_max::Float64)
    E, ω = _shc_k_band(sm, kf, Δk, γ, α, β, adaptive, adpt_fac, adpt_max, smr_width, eigval_max)
    nw = length(E)
    return [sum(Float64(E[n] < ef) * ω[n] for n in 1:nw) for ef in efs]
end

"Common per-k SHC data: kd, band energies/velocities, D_h, AA_β, and the js matrix."
function _shc_k_setup(sm, kf::SVector{3,Float64}, γ::Int, α::Int, β::Int)
    bm = sm.bm
    nw = num_wann(bm)
    kd = _berry_kdata(bm, kf)
    E = kd.E
    dE = [real(kd.dHh[c][n, n]) for c in 1:3, n in 1:nw]
    Dh = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    for c in 1:3, m in 1:nw, n in 1:nw
        (n != m && abs(E[m] - E[n]) > 1e-7) && (Dh[c][n, m] = kd.dHh[c][n, m] / (E[m] - E[n]))
    end
    AAβ = kd.U' * kd.A[β] * kd.U .+ im .* Dh[β]              # Eq. (25) WYSV06
    js = _shc_js(sm, kd, kf, γ, α, Dh)
    return E, dE, AAβ, js
end

"js (Qiao/QZYZ18 Eq. 23): B = ∂ε_m·S + ε_m·K − L (column m), js = (B + B†)/2."
function _shc_js(sm::ShcModel, kd, kf::SVector{3,Float64}, γ::Int, α::Int, Dh)
    bm = sm.bm
    E, U = kd.E, kd.U
    nw = length(E)
    Sk = U' * _ft_op(bm, (@view sm.SSr[:, :, :, γ]), kf) * U
    SRk = -im .* (U' * _ft_op(bm, (@view sm.SRr[:, :, :, γ, α]), kf) * U)
    K = SRk .+ Sk * Dh[α]
    SHRk = -im .* (U' * _ft_op(bm, (@view sm.SHRr[:, :, :, γ, α]), kf) * U)
    SHk = U' * _ft_op(bm, (@view sm.SHr[:, :, :, γ]), kf) * U
    L = SHRk .+ SHk * Dh[α]
    B = zeros(ComplexF64, nw, nw)
    for m in 1:nw, n in 1:nw
        B[n, m] = real(kd.dHh[α][m, m]) * Sk[n, m] + E[m] * K[n, m] - L[n, m]
    end
    return (B .+ B') ./ 2
end

"js (Ryoo/RPS19 Eqs. 21, 26, 37–40): ½{σ,v} from S, SAA, SBB and the full bar velocity."
function _shc_js(sm::ShcRyooModel, kd, kf::SVector{3,Float64}, γ::Int, α::Int, Dh)
    bm = sm.bm
    E, U = kd.E, kd.U
    nw = length(E)
    Sk = U' * _ft_op(bm, (@view sm.SSr[:, :, :, γ]), kf) * U
    SAA = U' * _ft_op(bm, (@view sm.SAAr[:, :, :, γ, α]), kf) * U
    SBB = U' * _ft_op(bm, (@view sm.SBBr[:, :, :, γ, α]), kf) * U
    VV0 = kd.dHh[α]                                          # U†(∂_α H^W)U, full matrix
    sv = VV0 * Sk .+ Sk * VV0
    js = zeros(ComplexF64, nw, nw)
    for m in 1:nw, n in 1:nw
        js[n, m] = (sv[n, m] - im * (E[m] * SAA[n, m] - SBB[n, m])
                    + im * (E[n] * conj(SAA[m, n]) - conj(SBB[m, n]))) / 2
    end
    return js
end

"Band-resolved SHC k-term ω_n = Ω^{spin γ}_{n,αβ}(k) (Ų, no occupation factor) and energies."
function _shc_k_band(sm::Union{ShcModel,ShcRyooModel}, kf::SVector{3,Float64}, Δk::Float64,
                     γ::Int, α::Int, β::Int, adaptive::Bool, adpt_fac::Float64,
                     adpt_max::Float64, smr_width::Float64, eigval_max::Float64)
    E, dE, AAβ, js = _shc_k_setup(sm, kf, γ, α, β)
    nw = length(E)
    ω = zeros(nw)
    for n in 1:nw
        E[n] > eigval_max && continue
        for m in 1:nw
            (m == n || E[m] > eigval_max) && continue
            η = adaptive ?
                min(adpt_fac * norm(SVector(dE[1, m] - dE[1, n], dE[2, m] - dE[2, n],
                                            dE[3, m] - dE[3, n])) * Δk, adpt_max) : smr_width
            prod = js[n, m] * im * (E[m] - E[n]) * AAβ[m, n]
            ω[n] += -2.0 * imag(prod) / ((E[m] - E[n])^2 + η^2)
        end
    end
    return E, ω
end

"""
    shc_imjv(sm, kf; γ=3, α=1, β=2) -> (E, imjv)

Band energies and the raw SHC integrand matrix `imjv[n,m] = Im[j^{spin γ}_{α,nm}·v_{β,mn}]`
(no energy denominator, no smearing) at one k — the quantity the tetrahedron method integrates.
Equals the per-pair numerator of the Gaussian [`shc_fermiscan`](@ref).
"""
function shc_imjv(sm::Union{ShcModel,ShcRyooModel}, kf::SVector{3,Float64};
                  γ::Int=3, α::Int=1, β::Int=2)
    E, dE, AAβ, js = _shc_k_setup(sm, kf, γ, α, β)
    nw = length(E)
    imjv = zeros(nw, nw)
    for m in 1:nw, n in 1:nw
        n == m && continue
        imjv[n, m] = imag(js[n, m] * im * (E[m] - E[n]) * AAβ[m, n])
    end
    return E, imjv
end

"""
    shc_freqscan(sm; freqs, fermi_energy, kmesh, γ=3, α=1, β=2, adaptive=true,
                 adpt_fac=√2, adpt_max=1.0, smr_width=0.0, eigval_max=Inf)
        -> Vector{ComplexF64} (per frequency, (ħ/e)·S/cm)

ac spin Hall conductivity σ^{spin γ}_{αβ}(ω) on a frequency list (eV), with T = 0
occupations at the single `fermi_energy`. `sm` is an `ShcModel` (Qiao) or `ShcRyooModel`.
"""
function shc_freqscan(sm::Union{ShcModel,ShcRyooModel}; freqs::Vector{Float64},
                      fermi_energy::Float64, kmesh::NTuple{3,Int}=(25, 25, 25),
                      γ::Int=3, α::Int=1, β::Int=2,
                      adaptive::Bool=true, adpt_fac::Float64=sqrt(2.0), adpt_max::Float64=1.0,
                      smr_width::Float64=0.0, eigval_max::Float64=Inf)
    bm = sm.bm
    nktot = prod(kmesh)
    Δk = kmesh_spacing(bm.lattice, kmesh)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    per_k = Vector{Vector{ComplexF64}}(undef, nktot)
    Threads.@threads for idx in 1:nktot
        per_k[idx] = _shc_k_freq(sm, kl[idx], freqs, fermi_energy, Δk, γ, α, β,
                                 adaptive, adpt_fac, adpt_max, smr_width, eigval_max)
    end
    fac = 1.0e8 * ELEM_CHARGE_SI^2 / (HBAR_SI * cell_volume(bm.lattice)) / 2.0
    return (fac / nktot) .* sum(per_k)
end

function _shc_k_freq(sm, kf::SVector{3,Float64}, freqs::Vector{Float64}, ef::Float64,
                     Δk::Float64, γ::Int, α::Int, β::Int, adaptive::Bool, adpt_fac::Float64,
                     adpt_max::Float64, smr_width::Float64, eigval_max::Float64)
    E, dE, AAβ, js = _shc_k_setup(sm, kf, γ, α, β)
    nw = length(E)
    ωl = zeros(ComplexF64, length(freqs))
    for n in 1:nw
        (E[n] > eigval_max || E[n] >= ef) && continue        # occupation on n only
        for m in 1:nw
            (m == n || E[m] > eigval_max) && continue
            η = adaptive ?
                min(adpt_fac * norm(SVector(dE[1, m] - dE[1, n], dE[2, m] - dE[2, n],
                                            dE[3, m] - dE[3, n])) * Δk, adpt_max) : smr_width
            rfac = E[m] - E[n]
            p = -2.0 * imag(js[n, m] * im * rfac * AAβ[m, n])
            for (i, w) in enumerate(freqs)
                ωl[i] += p / (rfac^2 - (w + im * η)^2)
            end
        end
    end
    return ωl
end

"""
    write_shc(path, x, shc; freq_scan=false) -> path

Write a postw90 `-shc-fermiscan.dat` (real `shc`) or `-shc-freqscan.dat`
(`freq_scan = true`, complex `shc`) file with the exact reference formats.
"""
function write_shc(path::AbstractString, x::Vector{Float64}, shc::Vector; freq_scan::Bool=false)
    open(path, "w") do io
        if freq_scan
            println(io, "#No.   Frequency(eV)   Re(sigma)((hbar/e)*S/cm)   Im(sigma)((hbar/e)*S/cm)")
            for i in 1:length(x)
                @printf(io, "%4d %12.6f  %s %s \n", i, x[i],
                        strip(fortran_e(real(shc[i]), 17, 8)), strip(fortran_e(imag(shc[i]), 17, 8)))
            end
        else
            println(io, "#No.   Fermi energy(eV)   SHC((hbar/e)*S/cm)")
            for i in 1:length(x)
                @printf(io, "%4d %12.6f %s\n", i, x[i], strip(fortran_e(real(shc[i]), 17, 8)))
            end
        end
    end
    return path
end
