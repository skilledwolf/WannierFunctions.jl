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

"""
    read_uhu(path; num_bands, num_kpts, nntot) -> uHu

Read a `.uHu` file: `uHu[m, n, b1, b2, q] = <u_{m,q+b1}|H_q|u_{n,q+b2}>`
(the transpose applied on read, as the reference does for pw2wannier90's ordering).
Formatted and Fortran-unformatted files are distinguished automatically.
"""
function read_uhu(path::AbstractString; num_bands::Int, num_kpts::Int, nntot::Int)
    # Sniff Fortran unformatted (pw2wannier90's default): the file starts with a 4-byte
    # little-endian record marker containing NUL bytes; a formatted file starts with text.
    isbinary = open(io -> any(==(0x00), read(io, 4)), path, "r")
    isbinary && return _read_uhu_unformatted(path, num_bands, num_kpts, nntot)
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

"Read a Fortran-unformatted `.uHu` (sequential records, 4-byte length markers)."
function _read_uhu_unformatted(path::AbstractString, num_bands::Int, num_kpts::Int, nntot::Int)
    open(path, "r") do io
        record() = begin
            len = read(io, Int32)
            data = read(io, Int(len))
            read(io, Int32) == len || error("corrupt Fortran record in $path")
            data
        end
        record()                                        # header string
        dims = reinterpret(Int32, record())
        nb, nk, nn = Int(dims[1]), Int(dims[2]), Int(dims[3])
        (nb == num_bands && nk == num_kpts && nn == nntot) ||
            error(".uHu dims ($nb,$nk,$nn) don't match model ($num_bands,$num_kpts,$nntot)")
        u = Array{ComplexF64,5}(undef, nb, nb, nntot, nntot, nk)
        for q in 1:nk, b2 in 1:nntot, b1 in 1:nntot
            blk = reinterpret(ComplexF64, record())     # Ho(n,m), n inner — one record per block
            length(blk) == nb * nb || error("unexpected .uHu record size in $path")
            idx = 0
            for m in 1:nb, n in 1:nb
                idx += 1
                u[m, n, b1, b2, q] = blk[idx]
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
    MorbModel(seedname; transl_inv_full=false) -> MorbModel

Assemble from `seedname.{win,mmn,eig,chk,uHu}` (formatted `.uHu`). With
`transl_inv_full = true`, A(R), B(R) and C(R) use the one-shell translation-invariant scheme
(e^{ib·r₀} phases, canonical b-slots on the expanded minimal-image R-set with e^{−ib·R̃/2},
plus the H-weighted correction terms of get_oper.F90; C(R) is then NOT Hermitian-paired).
"""
function MorbModel(seedname::AbstractString; transl_inv_full::Bool=false)
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
    uhu = read_uhu(seedname * ".uHu"; num_bands=nb, num_kpts=nk, nntot=nntot)

    nw = num_wann(chk)
    vs, winidx = gauge_v_windows(chk, nb)
    nr = length(bm.irvec)

    if !transl_inv_full
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
        Br = zeros(ComplexF64, nw, nw, nr, 3)
        Cr = zeros(ComplexF64, nw, nw, nr, 3, 3)
        for a in 1:3
            Br[:, :, :, a] = fourier_q_to_R((@view Bq[:, :, :, a]), kgrid, bm.irvec)
            for b in 1:3
                Cr[:, :, :, a, b] = fourier_q_to_R((@view Cq[:, :, :, a, b]), kgrid, bm.irvec)
            end
        end
        return MorbModel(bm, Br, Cr)
    end

    # ---- transl_inv_full: bm lives on the expanded R-set (pre-divided convention) ----
    irvec, ndegen = wigner_seitz(lattice, kgrid.mp_grid)
    ws0 = win.use_ws_distance ?
          ws_translate_dist(irvec, chk.centres, lattice, kgrid.mp_grid) : nothing
    irvec90, idx90 = _pw90_rset(irvec, ws0)
    @assert irvec90 == bm.irvec
    r0 = [(chk.centres[a, i] + chk.centres[a, j]) / 2 for i in 1:nw, j in 1:nw, a in 1:3]
    slot = _bvec_canonical_slots(bv)
    ph1 = (bb, i, j) -> cis(bb[1] * r0[i, j, 1] + bb[2] * r0[i, j, 2] + bb[3] * r0[i, j, 3])
    bk1 = b0 -> SVector(bv.bvec[1, b0, 1], bv.bvec[2, b0, 1], bv.bvec[3, b0, 1])

    # Expanded-set H(R) for the corrections (same convention as bm.Hr)
    Hr90 = bm.Hr

    # B: slot-resolved accumulation with e^{ib·r₀}, per-slot scatter + e^{−ib·R̃/2},
    # then B(R̃) += (r₀_a − R̃_a/2)·H(R̃)
    Bq_b = zeros(ComplexF64, nw, nw, nk, nntot, 3)
    for q in 1:nk
        εw = eig[winidx[q], q]
        for b in 1:nntot
            qb = kpb[b, q]
            bb = SVector(bv.bvec[1, b, q], bv.bvec[2, b, q], bv.bvec[3, b, q])
            b0 = slot[b, q]
            core = vs[q]' * (Diagonal(εw) * M[winidx[q], winidx[qb], b, q]) * vs[qb]
            w = bv.wb[b, q]
            for j in 1:nw, i in 1:nw
                c = im * w * ph1(bb, i, j) * core[i, j]
                for a in 1:3
                    Bq_b[i, j, q, b0, a] += bb[a] * c
                end
            end
        end
    end
    Br = zeros(ComplexF64, nw, nw, nr, 3)
    for a in 1:3, b0 in 1:nntot
        Xb = _scatter_ws(fourier_q_to_R((@view Bq_b[:, :, :, b0, a]), kgrid, irvec),
                         irvec, ndegen, ws0, idx90, nr)
        for ir in 1:nr
            p2 = cis(-dot(bk1(b0), bm.Rcart[ir]) / 2)
            @views Br[:, :, ir, a] .+= p2 .* Xb[:, :, ir]
        end
    end
    for ir in 1:nr, a in 1:3, j in 1:nw, i in 1:nw
        Br[i, j, ir, a] += (r0[i, j, a] - bm.Rcart[ir][a] / 2) * Hr90[i, j, ir]
    end

    # C: slot-pair accumulation (a ≤ b only) with e^{i(b₂−b₁)·r₀}, per-pair scatter with
    # e^{−i(b₁+b₂)·R̃/2}, then the three correction terms (a > b filled only by these).
    Cq_b = zeros(ComplexF64, nw, nw, nk, nntot, nntot, 3, 3)
    for q in 1:nk, b2 in 1:nntot, b1 in 1:nntot
        qb1, qb2 = kpb[b1, q], kpb[b2, q]
        bb1 = SVector(bv.bvec[1, b1, q], bv.bvec[2, b1, q], bv.bvec[3, b1, q])
        bb2 = SVector(bv.bvec[1, b2, q], bv.bvec[2, b2, q], bv.bvec[3, b2, q])
        s1, s2 = slot[b1, q], slot[b2, q]
        core = vs[qb1]' * uhu[winidx[qb1], winidx[qb2], b1, b2, q] * vs[qb2]
        w12 = bv.wb[b1, q] * bv.wb[b2, q]
        for j in 1:nw, i in 1:nw
            c = w12 * ph1(bb2 - bb1, i, j) * core[i, j]
            for b in 1:3, a in 1:b
                Cq_b[i, j, q, s1, s2, a, b] += bb1[a] * bb2[b] * c
            end
        end
    end
    Cr = zeros(ComplexF64, nw, nw, nr, 3, 3)
    for b in 1:3, a in 1:b, s2 in 1:nntot, s1 in 1:nntot
        Xb = _scatter_ws(fourier_q_to_R((@view Cq_b[:, :, :, s1, s2, a, b]), kgrid, irvec),
                         irvec, ndegen, ws0, idx90, nr)
        bsum = bk1(s1) + bk1(s2)
        for ir in 1:nr
            p2 = cis(-dot(bsum, bm.Rcart[ir]) / 2)
            @views Cr[:, :, ir, a, b] .+= p2 .* Xb[:, :, ir]
        end
    end
    irneg = [get(idx90, (-r[1], -r[2], -r[3]), 0) for r in bm.irvec]
    for ir in 1:nr, b in 1:3, a in 1:3
        Rc = bm.Rcart[ir]
        @views begin
            Cr[:, :, ir, a, b] .+= (r0[:, :, a] .+ Rc[a] / 2) .* Br[:, :, ir, b]
            irneg[ir] != 0 &&
                (Cr[:, :, ir, a, b] .+= (Br[:, :, irneg[ir], a])' .* (r0[:, :, b] .- Rc[b] / 2))
            Cr[:, :, ir, a, b] .+= ((r0[:, :, a] .+ Rc[a] / 2) .* Rc[b]) .* Hr90[:, :, ir]
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
    nw = num_wann(mm.bm)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    # NB: the k-point body lives in its own top-level function (_morb_kdata) rather than inline
    # in this loop — a large inlined threaded body with many captured locals produced
    # nondeterministic results on Julia 1.12 (closure boxing), while the factored form is
    # deterministic and matches the reference. Keep it factored.
    states = threaded_ksum(
        (st, idx) -> (st.acc .+= _morb_kdata(mm, kl[idx], fermi_energy, st.work)),
        () -> (work=BerryKWork(nw), acc=MVector{3,Float64}(0, 0, 0)),
        nktot)
    fac = -EV_AU / BOHR^2
    return (fac / nktot) .* SVector(sum(st.acc for st in states))
end

"Per-k morb integrand: img + imh − 2E_F·imf (J0+J1+J2 summed), Wannier-gauge traces."
function _morb_kdata(mm::MorbModel, kf::SVector{3,Float64}, ef::Float64,
                     w::BerryKWork=BerryKWork(num_wann(mm.bm)))
    imf, img, imh = _imfgh_kdata(mm, kf, ef, w)
    return SVector((img .+ imh .- 2.0 * ef .* imf)...)
end

"Per-k LVTS12 trace triple (imf, img, imh), each the J0+J1+J2 sum per axial component."
function _imfgh_kdata(mm::MorbModel, kf::SVector{3,Float64}, ef::Float64,
                      w::BerryKWork=BerryKWork(num_wann(mm.bm)))
    kd, H, B, C = _imfgh_setup(mm, kf, w)
    map!(e -> Float64(e < ef), w.occ, kd.E)
    return _imfgh_occ(kd, H, B, C, w.occ, w)
end

"Interpolate everything the imfgh traces need at one k (shared across occupation choices)."
function _imfgh_setup(mm::MorbModel, kf::SVector{3,Float64},
                      w::BerryKWork=BerryKWork(num_wann(mm.bm)))
    bm = mm.bm
    kd = _berry_kdata!(w, bm, kf)
    w.tmp .= kd.U .* kd.E'                  # U · diag(E)
    H = mul!(w.Hh, w.tmp, kd.U')            # H in the Hamiltonian gauge
    # _ft_op is Wigner–Seitz-distance aware (use_ws_distance), unlike a plain R-sum.
    # All 9 C components are interpolated directly: with transl_inv_full C(R) is not
    # Hermitian-paired, and in the plain case the stored components already are.
    B = [_ft_op(bm, (@view mm.Br[:, :, :, a]), kf) for a in 1:3]
    C = [_ft_op(bm, (@view mm.Cr[:, :, :, a, b]), kf) for a in 1:3, b in 1:3]
    return kd, H, B, C
end

function _imfgh_occ(kd, H, B, C, occ::Vector{Float64},
                    w::BerryKWork=BerryKWork(length(kd.E)))
    E, U = kd.E, kd.U
    nw = length(E)
    w.tmp .= U .* occ'                      # U · diag(occ)
    f = mul!(w.f, w.tmp, U')
    _jj_matrices!(w.JJp, w.JJm, w.tmp, kd, occ)
    JJp, JJm = w.JJp, w.JJm

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
    return SVector(imf), SVector(img), SVector(imh)
end
