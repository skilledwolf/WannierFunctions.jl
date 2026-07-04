# Orbital magnetisation (postw90's `berry_task = morb`), following Lopez–Vanderbilt–
# Thonhauser–Souza PRB 85, 014435 (2012) [LVTS12] as implemented in berry.F90/get_oper.F90;
# exact conventions in docs/reference-notes/kubo-morb-geninterp.md §3.
#
# Needs two operators beyond the AHC set, built from the ab-initio matrix elements:
#   B_a(R)  = <0n|H(r−R)_a|Rm> :  B_a(q) = i Σ_b w_b b_a  V(q)† diag(ε_q) M(q,b) V(q+b)
#   C_ab(R) = <0n|r_a H (r−R)_b|Rm> :
#             C_ab(q) = Σ_{b1,b2} w_b1 b1_a w_b2 b2_b  V(q+b1)† uHu(q;b1,b2) V(q+b2)
# with uHu(m,n) = <u_{m,q+b1}|H_q|u_{n,q+b2}> from the .uHu file. The g/h traces (LVTS12
# Eqs. 66/56) then give M_orb = −(eV_au/bohr²)·(−2Im[g] + −2Im[h] − 2E_F·(−2Im[f])) in μ_B/cell.

using LinearAlgebra
using StaticArrays

const EV_AU = 3.674932540e-2      # eV → Hartree (CODATA2006, constants.F90:178)

"""
    read_uhu(path; num_bands, num_kpts, nntot) -> uHu

Read a formatted `.uHu` file: `uHu[m, n, b1, b2, q] = <u_{m,q+b1}|H_q|u_{n,q+b2}>`
(the transpose applied on read, as the reference does for pw2wannier90's ordering).
"""
function read_uhu(path::AbstractString; num_bands::Int, num_kpts::Int, nntot::Int)
    open(path, "r") do io
        readline(io)                                    # header
        nb, nk, nn = parse.(Int, split(readline(io)))
        (nb == num_bands && nk == num_kpts && nn == nntot) ||
            error(".uHu dims ($nb,$nk,$nn) don't match model ($num_bands,$num_kpts,$nntot)")
        u = Array{ComplexF64,5}(undef, nb, nb, nntot, nntot, nk)
        for q in 1:nk, b2 in 1:nntot, b1 in 1:nntot
            # file stores Ho(n,m) with n inner; transpose ⇒ row m ↔ bra at q+b1
            for m in 1:nb, n in 1:nb
                t = split(readline(io))
                u[m, n, b1, b2, q] = complex(parse(Float64, t[1]), parse(Float64, t[2]))
            end
        end
        return u
    end
end

"""
    MorbModel

`BerryModel` plus the two H-weighted position operators B(R) (3 components) and C(R)
(3×3 components) needed for the orbital magnetisation.
"""
struct MorbModel
    bm::BerryModel
    Br::Array{ComplexF64,4}          # (nw, nw, nr, 3)
    Cr::Array{ComplexF64,5}          # (nw, nw, nr, 3, 3)
end

Base.show(io::IO, ::MIME"text/plain", m::MorbModel) =
    print(io, "MorbModel: ", num_wann(m.bm), " WF, ", length(m.bm.irvec), " R-vectors (H, A, B, C)")

"""
    MorbModel(seedname) -> MorbModel

Assemble from `seedname.{win,mmn,eig,chk,uHu}` (formatted `.uHu`).
"""
function MorbModel(seedname::AbstractString)
    win = read_win(seedname * ".win")
    chk = isfile(seedname * ".chk") ? read_chk(seedname * ".chk") :
          read_chk_fmt(seedname * ".chk.fmt")
    eig = read_eig(seedname * ".eig")
    M, kpb, gpb, nb, nk, nntot = read_mmn(seedname * ".mmn")
    lattice = Lattice(win.unit_cell)
    kgrid = KGrid(win.kpoints, win.mp_grid)
    bv = build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol=win.kmesh_tol)
    bm = BerryModel(chk, eig, bv, kgrid, lattice)
    uhu = read_uhu(seedname * ".uHu"; num_bands=nb, num_kpts=nk, nntot=nntot)

    nw = num_wann(chk)
    # windowed v(q) = U_opt·U on the window rows, and the window band indices per q
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

    # B_a(q) = i Σ_b w_b b_a V(q)† diag(ε_win(q)) S_o(win_q, win_qb) V(q+b)
    Bq = zeros(ComplexF64, nw, nw, nk, 3)
    for q in 1:nk
        εw = eig[winidx[q], q]
        for b in 1:nntot
            qb = kpb[b, q]
            S = M[winidx[q], winidx[qb], b, q]
            core = vs[q]' * (Diagonal(εw) * S) * vs[qb]
            w = bv.wb[b, q]
            for a in 1:3
                @views Bq[:, :, q, a] .+= (im * w * bv.bvec[a, b, q]) .* core
            end
        end
    end

    # C_ab(q) = Σ_{b1,b2} w_b1 b1_a w_b2 b2_b V(q+b1)† uHu(b1,b2) V(q+b2), a ≤ b then
    # Hermitian completion.
    Cq = zeros(ComplexF64, nw, nw, nk, 3, 3)
    for q in 1:nk, b2 in 1:nntot, b1 in 1:nntot
        qb1, qb2 = kpb[b1, q], kpb[b2, q]
        core = vs[qb1]' * uhu[winidx[qb1], winidx[qb2], b1, b2, q] * vs[qb2]
        for b in 1:3, a in 1:b
            fac = bv.wb[b1, q] * bv.bvec[a, b1, q] * bv.wb[b2, q] * bv.bvec[b, b2, q]
            @views Cq[:, :, q, a, b] .+= fac .* core
        end
    end
    for q in 1:nk, b in 1:3, a in 1:b-1
        @views Cq[:, :, q, b, a] .= (Cq[:, :, q, a, b])'
    end

    # q → R on the Wigner–Seitz set
    nr = length(bm.irvec)
    Br = zeros(ComplexF64, nw, nw, nr, 3)
    Cr = zeros(ComplexF64, nw, nw, nr, 3, 3)
    Threads.@threads for ir in 1:nr
        R = SVector{3,Float64}(bm.irvec[ir]...)
        for q in 1:nk
            fac = cis(-TWOPI * dot(kgrid.frac[q], R)) / nk
            @views for a in 1:3
                Br[:, :, ir, a] .+= fac .* Bq[:, :, q, a]
                for b in 1:3
                    Cr[:, :, ir, a, b] .+= fac .* Cq[:, :, q, a, b]
                end
            end
        end
    end
    return MorbModel(bm, Br, Cr)
end

"""
    orbital_magnetisation(mm; fermi_energy, kmesh) -> SVector{3,Float64}

Orbital magnetisation (μ_B per cell, x/y/z) from the LVTS12 trace formulas on a uniform k-mesh
(no adaptive refinement, matching postw90).
"""
function orbital_magnetisation(mm::MorbModel; fermi_energy::Float64,
                               kmesh::NTuple{3,Int}=(25, 25, 25))
    nktot = prod(kmesh)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    per_k = Vector{SVector{3,Float64}}(undef, nktot)
    # NB: the k-point body lives in its own top-level function (_morb_kdata) rather than inline
    # in this loop — a large inlined @threads body with many captured locals produced
    # nondeterministic results on Julia 1.12 (closure boxing), while the factored form is
    # deterministic and matches the reference. Keep it factored.
    Threads.@threads for idx in 1:nktot
        per_k[idx] = _morb_kdata(mm, kl[idx], fermi_energy)
    end
    fac = -EV_AU / BOHR^2
    return (fac / nktot) .* sum(per_k)
end

"Per-k morb integrand: img + imh − 2E_F·imf (J0+J1+J2 summed), Wannier-gauge traces."
function _morb_kdata(mm::MorbModel, kf::SVector{3,Float64}, ef::Float64)
    bm = mm.bm
    nw = num_wann(bm)
    kd = _berry_kdata(bm, kf)
    E, U = kd.E, kd.U
    H = U * Diagonal(E) * U'

    B = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    C = [zeros(ComplexF64, nw, nw) for _ in 1:3, _ in 1:3]
    for ir in 1:length(bm.irvec)
        fac = cis(TWOPI * dot(kf, SVector{3,Float64}(bm.irvec[ir]...))) / bm.ndegen[ir]
        @views for a in 1:3
            B[a] .+= fac .* mm.Br[:, :, ir, a]
            for b in a:3
                C[a, b] .+= fac .* mm.Cr[:, :, ir, a, b]
            end
        end
    end
    for b in 1:3, a in 1:b-1
        C[b, a] = C[a, b]'
    end

    occ = Float64.(E .< ef)
    f = U * Diagonal(occ) * U'
    JJp = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    JJm = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    for c in 1:3
        for m in 1:nw, n in 1:nw
            if E[n] > ef && E[m] < ef
                JJp[c][n, m] = im * kd.dHh[c][n, m] / (E[m] - E[n])
                JJm[c][m, n] = im * kd.dHh[c][m, n] / (E[n] - E[m])
            end
        end
        JJp[c] = U * JJp[c] * U'
        JJm[c] = U * JJm[c] * U'
    end

    imf = MVector{3,Float64}(0, 0, 0)
    img = MVector{3,Float64}(0, 0, 0)
    imh = MVector{3,Float64}(0, 0, 0)
    for i in 1:3
        α, β = ALPHA_A[i], BETA_A[i]
        Λ = im .* (C[α, β] .- C[α, β]')
        HA = H * kd.A[α]
        s = 2.0 * imag(tr(f * (HA * f * kd.A[β])))
        # J0
        imf[i] += real(tr(f * kd.Ω̄[i]))
        img[i] += real(tr(f * Λ)) - s
        imh[i] += real(tr(f * (H * kd.Ω̄[i]))) + s
        # J1
        imf[i] += -2.0 * (imag(tr(kd.A[α] * JJp[β])) + imag(tr(JJm[α] * kd.A[β])))
        img[i] += -2.0 * (imag(tr(JJm[α] * B[β])) - imag(tr(JJm[β] * B[α])))
        imh[i] += -2.0 * (imag(tr(HA * JJp[β])) + imag(tr((H * JJm[α]) * kd.A[β])))
        # J2
        imf[i] += -2.0 * imag(tr(JJm[α] * JJp[β]))
        img[i] += -2.0 * imag(tr((JJm[α] * H) * JJp[β]))
        imh[i] += -2.0 * imag(tr((H * JJm[α]) * JJp[β]))
    end
    return SVector((img .+ imh .- 2.0 * ef .* imf)...)
end
