# Shift current σ_abc(0; ω, −ω) by Wannier interpolation (postw90's berry_task = sc;
# Ibañez-Azpiroz–Tsirkin–Souza PRB 97, 245143 (2018) [IATS18], with the finite-η correction of
# PRB 103, 247101 (2021)). Exact conventions in docs/reference-notes/shift-current.md.
#
# Per k and ordered band pair (n,m): I_nm(a,bc) = Im[r_mn(b)·r^gen_nm(c;a) + (b↔c)], deposited
# on the frequency grid with Gaussian deltas δ(E_nm ∓ ω); the generalised derivative r^gen is
# the 8-term IATS18 Eq. 34 built from A, ∂A, ∂H, ∂²H, the η-regularised (principal-value) and
# bare D matrices, and band velocities. Result in A/V².

using LinearAlgebra
using StaticArrays
using Printf

const EV_SECONDS = 6.582119e-16          # ħ/e in eV·s — CODATA2006 set (7 digits!)

# (b,c) packing of the 6 symmetric polarisation pairs: xx yy zz xy xz yz
const ALPHA_SC = (1, 2, 3, 1, 1, 2)
const BETA_SC = (1, 2, 3, 2, 3, 3)

"""
    shift_current(bm, centres; fermi_energy, freqs, kmesh, phase_conv=1, sc_eta=0.04,
                  w_thr=5.0, eta_corr=true, adaptive=true, adpt_fac=√2, adpt_max=1.0,
                  smr_width=0.0, eigval_max=Inf) -> (; freqs, sc)

Shift-current tensor on a uniform frequency grid (`freqs` must be uniformly spaced, as the
reference assumes). `centres` is the 3×nw matrix of Wannier centres (Å, from the checkpoint) —
used by the tight-binding phase convention (`phase_conv = 1`); `phase_conv = 2` is the plain
Wannier90 convention. `sc` is 3×6×nfreq: `sc[a, bc, :]` = σ_abc with a the current direction
and bc packing (xx, yy, zz, xy, xz, yz), in A/V².
"""
function shift_current(bm::BerryModel, centres::AbstractMatrix;
                       fermi_energy::Float64, freqs::Vector{Float64},
                       kmesh::NTuple{3,Int}=(25, 25, 25), phase_conv::Int=1,
                       sc_eta::Float64=0.04, w_thr::Float64=5.0, eta_corr::Bool=true,
                       adaptive::Bool=true, adpt_fac::Float64=sqrt(2.0),
                       adpt_max::Float64=1.0, smr_width::Float64=0.0,
                       eigval_max::Float64=Inf)
    phase_conv in (1, 2) || error("shift_current: phase_conv must be 1 or 2")
    size(bm.Ar, 4) == 3 || error("shift_current: the BerryModel needs A(R) (.mmn)")
    ops = _sc_ops(bm, centres, phase_conv)
    nfreq = length(freqs)
    Δk = kmesh_spacing(bm.lattice, kmesh)
    nktot = prod(kmesh)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    per_k = Vector{Array{Float64,3}}(undef, nktot)
    Threads.@threads for idx in 1:nktot
        per_k[idx] = _sc_kpoint(ops, kl[idx], freqs, fermi_energy, Δk, sc_eta, w_thr,
                                eta_corr, adaptive, adpt_fac, adpt_max, smr_width, eigval_max)
    end
    fac = EV_SECONDS * pi * ELEM_CHARGE_SI^3 /
          (4.0 * HBAR_SI^2 * cell_volume(bm.lattice)) / nktot
    return (; freqs, sc=fac .* sum(per_k))
end

"""
    shift_current(seedname; kwargs...) -> (; freqs, sc)

Convenience: assemble the `BerryModel` and centres from `seedname.{win,mmn,eig,chk}`.
"""
function shift_current(seedname::AbstractString; kwargs...)
    chk = isfile(seedname * ".chk") ? read_chk(seedname * ".chk") :
          read_chk_fmt(seedname * ".chk.fmt")
    return shift_current(BerryModel(seedname), chk.centres; kwargs...)
end

# Effective R-space operators for the sc interpolation: ws_distance pre-folded onto the
# expanded minimal-image set (as the reference's operator_wigner_setup does), weights already
# divided out; TB convention (1) removes each WF's own centre from the R = 0 diagonal of A(R)
# and carries τ in the phases, convention (2) sets τ = 0.
struct _ScOps
    lattice::Lattice
    irvec::Vector{NTuple{3,Int}}
    Rcart::Vector{SVector{3,Float64}}
    Hr::Array{ComplexF64,3}
    Ar::Array{ComplexF64,4}
    τc::Matrix{Float64}       # 3×nw Cartesian (zero for phase_conv = 2)
    τf::Matrix{Float64}       # 3×nw fractional
end

function _sc_ops(bm::BerryModel, centres::AbstractMatrix, phase_conv::Int)
    nw = num_wann(bm)
    if bm.wsdist === nothing
        irvec = bm.irvec
        Hr = similar(bm.Hr)
        Ar = similar(bm.Ar)
        for ir in 1:length(irvec)
            @views Hr[:, :, ir] .= bm.Hr[:, :, ir] ./ bm.ndegen[ir]
            @views Ar[:, :, ir, :] .= bm.Ar[:, :, ir, :] ./ bm.ndegen[ir]
        end
    else
        irvec90, idx90 = _pw90_rset(bm.irvec, bm.wsdist)
        irvec = irvec90
        Hr = _scatter_ws(bm.Hr, bm.irvec, bm.ndegen, bm.wsdist, idx90, length(irvec90))
        Ar = zeros(ComplexF64, nw, nw, length(irvec90), 3)
        for c in 1:3
            Ar[:, :, :, c] = _scatter_ws((@view bm.Ar[:, :, :, c]), bm.irvec, bm.ndegen,
                                         bm.wsdist, idx90, length(irvec90))
        end
    end
    Rcart = [bm.lattice.A * SVector{3,Float64}(r...) for r in irvec]
    τc = phase_conv == 1 ? Matrix{Float64}(centres) : zeros(3, nw)
    τf = phase_conv == 1 ? Matrix{Float64}(bm.lattice.A \ centres) : zeros(3, nw)
    if phase_conv == 1
        ir0 = findfirst(==((0, 0, 0)), irvec)
        for c in 1:3, i in 1:nw
            Ar[i, i, ir0, c] -= τc[c, i]
        end
    end
    return _ScOps(bm.lattice, irvec, Rcart, Hr, Ar, τc, τf)
end

# Gather all W-gauge k-matrices with the per-pair TB phases: H, ∂H, ∂²H, A, ∂A.
function _sc_gather(ops::_ScOps, kf::SVector{3,Float64})
    nw = size(ops.Hr, 1)
    HH = zeros(ComplexF64, nw, nw)
    HHda = zeros(ComplexF64, nw, nw, 3)
    HHdadb = zeros(ComplexF64, nw, nw, 3, 3)
    AA = zeros(ComplexF64, nw, nw, 3)
    AAda = zeros(ComplexF64, nw, nw, 3, 3)
    for ir in 1:length(ops.irvec)
        R = SVector{3,Float64}(ops.irvec[ir]...)
        Rc = ops.Rcart[ir]
        for j in 1:nw, i in 1:nw
            ph = cis(TWOPI * (dot(kf, R) +
                              kf[1] * (ops.τf[1, j] - ops.τf[1, i]) +
                              kf[2] * (ops.τf[2, j] - ops.τf[2, i]) +
                              kf[3] * (ops.τf[3, j] - ops.τf[3, i])))
            h = ph * ops.Hr[i, j, ir]
            HH[i, j] += h
            rc = SVector(Rc[1] + ops.τc[1, j] - ops.τc[1, i],
                         Rc[2] + ops.τc[2, j] - ops.τc[2, i],
                         Rc[3] + ops.τc[3, j] - ops.τc[3, i])
            for a in 1:3
                HHda[i, j, a] += im * rc[a] * h
                for c in 1:3
                    HHdadb[i, j, c, a] -= rc[c] * rc[a] * h
                end
            end
            for c in 1:3
                av = ph * ops.Ar[i, j, ir, c]
                AA[i, j, c] += av
                for a in 1:3
                    AAda[i, j, c, a] += im * rc[a] * av
                end
            end
        end
    end
    return HH, HHda, HHdadb, AA, AAda
end

function _sc_kpoint(ops::_ScOps, kf::SVector{3,Float64}, ω::Vector{Float64}, ef::Float64,
                    Δk::Float64, sc_eta::Float64, w_thr::Float64, eta_corr::Bool,
                    adaptive::Bool, adpt_fac::Float64, adpt_max::Float64,
                    smr_width::Float64, eigval_max::Float64)
    nw = size(ops.Hr, 1)
    nfreq = length(ω)
    out = zeros(3, 6, nfreq)
    HH, HHda, HHdadb, AA, AAda = _sc_gather(ops, kf)
    F = eigen(Hermitian(HH))
    E, U = F.values, F.vectors

    Ab = [U' * (@view AA[:, :, c]) * U for c in 1:3]
    Hdb = [U' * (@view HHda[:, :, a]) * U for a in 1:3]
    Adab = [U' * (@view AAda[:, :, c, a]) * U for c in 1:3, a in 1:3]
    Hdadb = [U' * (@view HHdadb[:, :, c, a]) * U for c in 1:3, a in 1:3]
    dE = [real(Hdb[a][n, n]) for a in 1:3, n in 1:nw]
    # D matrices: principal-value (sc_eta) and bare (1e-7 degeneracy guard)
    Dh = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    Dh0 = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    for a in 1:3, m in 1:nw, n in 1:nw
        n == m && continue
        ΔE = E[m] - E[n]
        Dh[a][n, m] = Hdb[a][n, m] * ΔE / (ΔE^2 + sc_eta^2)
        abs(ΔE) >= 1e-7 && (Dh0[a][n, m] = Hdb[a][n, m] / ΔE)
    end
    occ = Float64.(E .< ef)

    wmin, wmax = ω[1], ω[end]
    wstep = length(ω) > 1 ? ω[2] - ω[1] : 1.0
    Inm = zeros(3, 6)
    genr = zeros(ComplexF64, 3)
    for m in 1:nw, n in 1:nw
        n == m && continue
        (E[m] > eigval_max || E[n] > eigval_max) && continue
        occ_fac = occ[n] - occ[m]
        abs(occ_fac) < 1e-10 && continue
        η = adaptive ?
            min(adpt_fac * norm(SVector(dE[1, m] - dE[1, n], dE[2, m] - dE[2, n],
                                        dE[3, m] - dE[3, n])) * Δk, adpt_max) : smr_width
        Enm = E[n] - E[m]
        thr = w_thr * η
        ((Enm + thr < wmin || Enm - thr > wmax) &&
         (-Enm + thr < wmin || -Enm - thr > wmax)) && continue

        for a in 1:3
            for c in 1:3
                # sums over p ∉ {n, m} (IATS18 Eqs. 30/32)
                sAD = -Ab[c][n, n] * Dh[a][n, m] + Dh[a][n, m] * Ab[c][m, m]
                sHD = -Hdb[c][n, n] * Dh[a][n, m] + Dh[a][n, m] * Hdb[c][m, m]
                for p in 1:nw
                    sAD += Ab[c][n, p] * Dh[a][p, m] - Dh[a][n, p] * Ab[c][p, m]
                    sHD += Hdb[c][n, p] * Dh[a][p, m] - Dh[a][n, p] * Hdb[c][p, m]
                end
                g = Adab[c, a][n, m] +
                    (Ab[c][n, n] - Ab[c][m, m]) * Dh0[a][n, m] +
                    (Ab[a][n, n] - Ab[a][m, m]) * Dh0[c][n, m] -
                    im * Ab[c][n, m] * (Ab[a][n, n] - Ab[a][m, m]) +
                    sAD +
                    im * (Hdadb[c, a][n, m] + sHD +
                          Dh0[c][n, m] * (dE[a, n] - dE[a, m]) +
                          Dh0[a][n, m] * (dE[c, n] - dE[c, m])) / (E[m] - E[n])
                if eta_corr
                    for p in 1:nw
                        (p == n || p == m) && continue
                        g -= sc_eta^2 / ((E[p] - E[m])^2 + sc_eta^2) / (E[n] - E[m]) *
                             (Ab[c][n, p] * Hdb[a][p, m] -
                              (Hdb[c][n, p] + im * (E[n] - E[p]) * Ab[c][n, p]) * Ab[a][p, m])
                        g += sc_eta^2 / ((E[n] - E[p])^2 + sc_eta^2) / (E[n] - E[m]) *
                             (Hdb[a][n, p] * Ab[c][p, m] -
                              Ab[a][n, p] * (Hdb[c][p, m] + im * (E[p] - E[m]) * Ab[c][p, m]))
                    end
                end
                genr[c] = g
            end
            for bc in 1:6
                b, c = ALPHA_SC[bc], BETA_SC[bc]
                rb = Ab[b][m, n] + im * Dh0[b][m, n]
                rc = Ab[c][m, n] + im * Dh0[c][m, n]
                Inm[a, bc] = imag(rb * genr[c] + rc * genr[b])
            end
        end

        # deposit on the ω grid: δ(ω − E_nm) and δ(ω − E_mn), int()-truncated windows
        for s in (1.0, -1.0)
            Es = s * Enm
            istart = max(trunc(Int, (Es - thr - wmin) / wstep + 1), 1)
            iend = min(trunc(Int, (Es + thr - wmin) / wstep + 1), nfreq)
            for i in istart:iend
                δ = _w0gauss((ω[i] - Es) / η) / η
                for bc in 1:6, a in 1:3
                    out[a, bc, i] += occ_fac * Inm[a, bc] * δ
                end
            end
        end
    end
    return out
end

"""
    write_shift_current(seedname, freqs, sc) -> seedname

Write the 18 `seedname-sc_<abc>.dat` files (a = current direction, bc the symmetric pair),
`(2E18.8E3)` rows of ω [eV] and σ_abc [A/V²].
"""
function write_shift_current(seedname::AbstractString, freqs::Vector{Float64},
                             sc::Array{Float64,3})
    dirs = ("x", "y", "z")
    for a in 1:3, bc in 1:6
        name = seedname * "-sc_" * dirs[a] * dirs[ALPHA_SC[bc]] * dirs[BETA_SC[bc]] * ".dat"
        open(name, "w") do io
            for i in 1:length(freqs)
                println(io, fortran_e(freqs[i], 18, 8; edigits=3),
                        fortran_e(sc[a, bc, i], 18, 8; edigits=3))
            end
        end
    end
    return seedname
end
