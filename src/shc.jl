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
    ShcModel(seedname) -> ShcModel

Assemble from `seedname.{win,mmn,eig,chk,spn}` (Qiao method: no `.uHu` needed).
"""
function ShcModel(seedname::AbstractString)
    win = read_win(seedname * ".win")
    chk = isfile(seedname * ".chk") ? read_chk(seedname * ".chk") :
          read_chk_fmt(seedname * ".chk.fmt")
    eig = read_eig(seedname * ".eig")
    M, kpb, gpb, nb, nk, nntot = read_mmn(seedname * ".mmn")
    lattice = Lattice(win.unit_cell)
    kgrid = KGrid(win.kpoints, win.mp_grid)
    bv = build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol=win.kmesh_tol)
    bm = BerryModel(chk, eig, bv, kgrid, lattice; use_ws_distance=win.use_ws_distance)
    spn = read_spn(seedname * ".spn"; num_bands=nb, num_kpts=nk)

    nw = num_wann(chk)
    vs = Vector{Matrix{ComplexF64}}(undef, nk)
    winidx = Vector{Vector{Int}}(undef, nk)
    for q in 1:nk
        if chk.have_disentangled
            nd = chk.ndimwin[q]
            vs[q] = chk.u_matrix_opt[1:nd, :, q] * chk.u_matrix[:, :, q]
            winidx[q] = findall(@view chk.lwindow[:, q])
        else
            vs[q] = chk.u_matrix[:, :, q]
            winidx[q] = collect(1:nb)
        end
    end

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
    bm = sm.bm
    nw = num_wann(bm)
    kd = _berry_kdata(bm, kf)
    E, U = kd.E, kd.U
    dE = [real(kd.dHh[c][n, n]) for c in 1:3, n in 1:nw]
    Dh = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    for c in 1:3, m in 1:nw, n in 1:nw
        (n != m && abs(E[m] - E[n]) > 1e-7) && (Dh[c][n, m] = kd.dHh[c][n, m] / (E[m] - E[n]))
    end
    AAβ = U' * kd.A[β] * U .+ im .* Dh[β]                    # Eq. (25) WYSV06

    # js (QZYZ18 Eq. 23): B = ∂ε_m·S + ε_m·K − L (column m), js = (B + B†)/2
    Sk = U' * _ft_op(bm, (@view sm.SSr[:, :, :, γ]), kf) * U
    SRk = -im .* (U' * _ft_op(bm, (@view sm.SRr[:, :, :, γ, α]), kf) * U)
    K = SRk .+ Sk * Dh[α]
    SHRk = -im .* (U' * _ft_op(bm, (@view sm.SHRr[:, :, :, γ, α]), kf) * U)
    SHk = U' * _ft_op(bm, (@view sm.SHr[:, :, :, γ]), kf) * U
    L = SHRk .+ SHk * Dh[α]
    B = zeros(ComplexF64, nw, nw)
    for m in 1:nw, n in 1:nw
        B[n, m] = dE[α, m] * Sk[n, m] + E[m] * K[n, m] - L[n, m]
    end
    js = (B .+ B') ./ 2

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
    return [sum(Float64(E[n] < ef) * ω[n] for n in 1:nw) for ef in efs]
end
