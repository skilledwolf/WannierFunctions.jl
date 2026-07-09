# Gyrotropic responses (postw90's gyrotropic module; Tsirkin–Aguado-Puente–Souza,
# PRB 97, 035158 (2018) [TAS17]): Berry-curvature dipole D and its frequency-dependent
# generalisation tildeD(ω), the C tensor, the kinetic magnetoelectric K tensor (orbital and
# spin parts), natural optical activity γ_abc (orbital and spin), and the Fermi-level DOS —
# all as Fermi-surface integrals over a (possibly reduced) k-box. Exact conventions in
# docs/reference-notes/gyrotropic-module.md.

using LinearAlgebra
using StaticArrays
using Printf

"""
    gyrotropic(m; tasks, fermi_energies, freqs=[0.0], kmesh, smr_width,
               box=I, box_corner=zeros(3), smr_max_arg=5.0, degen_thresh=0.0,
               band_list=nothing, eigval_max=Inf, spin=nothing) -> NamedTuple

Gyrotropic responses on a uniform grid over the fractional box spanned by the rows of `box`
from `box_corner`. `tasks` ⊆ `(:D0, :Dw, :C, :K, :NOA, :dos)`; `:K` needs a `MorbModel`
(`.uHu`), the others a `BerryModel`; a `SpinModel` in `spin` adds the spin parts of K and NOA.
Frequencies are broadened internally as ω + i·`smr_width` (the same η that smears the
Fermi-surface delta — Gaussian, as postw90's default smearing type). Returns
`(; fermi_energies, freqs, D, Dw, C, K_orb, K_spn, NOA_orb, NOA_spn, dos)` in the reference
units (D/tildeD dimensionless, C A/cm, K A, NOA Å, DOS eV⁻¹Å⁻³); tensors are indexed
`[velocity_dir, quantity_dir, ifermi(, ifreq)]`, NOA as `[ab_axial, c, ifermi, ifreq]`.
"""
function gyrotropic(m::Union{BerryModel,MorbModel};
                    tasks=(:D0, :Dw, :C, :K, :NOA, :dos),
                    fermi_energies::Vector{Float64}, freqs::Vector{Float64}=[0.0],
                    kmesh::NTuple{3,Int}=(25, 25, 25), smr_width::Float64,
                    box::AbstractMatrix=Matrix{Float64}(LinearAlgebra.I, 3, 3),
                    box_corner::AbstractVector=[0.0, 0.0, 0.0],
                    smr_max_arg::Float64=5.0, degen_thresh::Float64=0.0,
                    band_list::Union{Nothing,Vector{Int}}=nothing,
                    eigval_max::Float64=Inf,
                    spin::Union{Nothing,SpinModel}=nothing)
    bm = m isa BerryModel ? m : m.bm
    all(t -> t in (:D0, :Dw, :C, :K, :NOA, :dos), tasks) ||
        error("gyrotropic: tasks must be a subset of (:D0, :Dw, :C, :K, :NOA, :dos)")
    :K in tasks && !(m isa MorbModel) && error("gyrotropic: the K tensor needs a MorbModel (.uHu)")
    nw = num_wann(bm)
    bl = band_list === nothing ? collect(1:nw) : band_list
    bx = SMatrix{3,3,Float64}(box)
    corner = SVector{3,Float64}(box_corner...)
    ωc = ComplexF64.(freqs) .+ im * smr_width
    nf, nω = length(fermi_energies), length(freqs)

    nktot = prod(kmesh)
    kweight = det(bx) / nktot
    kl = Vector{SVector{3,Float64}}(undef, nktot)
    idx = 0
    for lx in 0:kmesh[1]-1, ly in 0:kmesh[2]-1, lz in 0:kmesh[3]-1
        f = SVector(lx / kmesh[1], ly / kmesh[2], lz / kmesh[3])
        kl[idx += 1] = corner + bx' * f          # k_j = corner_j + Σ_i f_i box[i,j]
    end

    states = threaded_ksum(
        (st, ik) -> _gyro_kpoint!(st.acc, m, bm, kl[ik], kweight, tasks, fermi_energies, ωc,
                                  smr_width, smr_max_arg, degen_thresh, bl, eigval_max, spin,
                                  st.work, st.occ),
        () -> (acc=_gyro_acc(nf, nω), work=BerryKWork(nw), occ=zeros(nw)),
        nktot)
    tot = states[1].acc
    for ic in 2:length(states), fld in 1:length(tot)
        tot[fld] .+= states[ic].acc[fld]
    end
    D, Dw, C, Korb, Kspn, dos, NOAorb, NOAspn = tot

    V = cell_volume(bm.lattice)
    D ./= V
    Dw ./= V
    C .*= 1.0e8 * ELEM_CHARGE_SI^2 / (TWOPI * HBAR_SI * V)
    Korb .*= ELEM_CHARGE_SI^2 / (2.0 * HBAR_SI * V)
    Kspn .*= -1.0e20 * ELEM_CHARGE_SI * HBAR_SI / (2.0 * ELEC_MASS_SI * V)
    NOAorb .*= 1.0e10 * ELEM_CHARGE_SI / (V * EPS0_SI)
    NOAspn .*= 1.0e30 * HBAR_SI^2 / (V * EPS0_SI * ELEC_MASS_SI)
    dos ./= V
    return (; fermi_energies, freqs, D, Dw, C, K_orb=Korb, K_spn=Kspn,
            NOA_orb=NOAorb, NOA_spn=NOAspn, dos)
end

_gyro_acc(nf::Int, nω::Int) = (zeros(3, 3, nf), zeros(3, 3, nf, nω), zeros(3, 3, nf),
                               zeros(3, 3, nf), zeros(3, 3, nf), zeros(nf),
                               zeros(3, 3, nf, nω), zeros(3, 3, nf, nω))

"Pauli expectation matrices S_h[j] = U†S_j(k)U (full matrices, j = x,y,z)."
function _spin_S3(sp::SpinModel, kf::SVector{3,Float64}, U::AbstractMatrix)
    return [U' * _ft_op(sp.bm, (@view sp.SSr[:, :, :, j]), kf) * U for j in 1:3]
end

# Top-level per-k body (closure-boxing hazard in large @threads bodies on Julia 1.12).
function _gyro_kpoint!(out, m, bm::BerryModel, kf::SVector{3,Float64}, kweight::Float64,
                       tasks, efs::Vector{Float64}, ωc::Vector{ComplexF64}, η::Float64,
                       max_arg::Float64, degen_thresh::Float64, bl::Vector{Int},
                       eigval_max::Float64, spin,
                       w::BerryKWork=BerryKWork(num_wann(bm)),
                       occ::Vector{Float64}=zeros(num_wann(bm)))
    D, Dw, C, Korb, Kspn, dos, NOAorb, NOAspn = out
    # With :K the kdata is (re)built inside _imfgh_setup on the same workspace — the
    # aliased kd fields are identical either way.
    setup = :K in tasks ? _imfgh_setup(m, kf, w) : nothing
    kd = setup === nothing ? _berry_kdata!(w, bm, kf) : setup[1]
    E, U = kd.E, kd.U
    nw = length(E)
    dE = [real(kd.dHh[c][n, n]) for c in 1:3, n in 1:nw]

    AA = nothing
    if :Dw in tasks || :NOA in tasks
        Dh = [zeros(ComplexF64, nw, nw) for _ in 1:3]
        for c in 1:3, mm in 1:nw, n in 1:nw
            (n != mm && abs(E[mm] - E[n]) > 1e-7) &&
                (Dh[c][n, mm] = kd.dHh[c][n, mm] / (E[mm] - E[n]))
        end
        AA = [U' * kd.A[c] * U .+ im .* Dh[c] for c in 1:3]
    end
    S3 = spin === nothing ? nothing : _spin_S3(spin, kf, U)

    for n in bl
        (n > 1 && E[n] - E[n-1] <= degen_thresh) && continue
        (n < nw && E[n+1] - E[n] <= degen_thresh) && continue
        fill!(occ, 0.0)
        occ[n] = 1.0
        curv = orb = nothing
        if :K in tasks
            f, g, h = _imfgh_occ(setup..., occ, w)
            curv = f
            orb = h .- g
        elseif :D0 in tasks || :C in tasks || :dos in tasks
            :D0 in tasks && (curv = _imf_occ!(w, kd, occ))
        end
        curvw = :Dw in tasks ? _gyro_curv_w(E, AA, n, bl, ωc) : nothing     # (nω, 3)

        for (fi, ef) in enumerate(efs)
            arg = (E[n] - ef) / η
            abs(arg) > max_arg && continue
            δ = _w0gauss(arg) / η * kweight
            for i in 1:3
                v = dE[i, n] * δ
                for j in 1:3
                    (:D0 in tasks || :K in tasks) && curv !== nothing &&
                        (D[i, j, fi] += v * curv[j])
                    :K in tasks && (Korb[i, j, fi] += v * orb[j])
                    :K in tasks && S3 !== nothing && (Kspn[i, j, fi] += v * real(S3[j][n, n]))
                    :C in tasks && (C[i, j, fi] += v * dE[j, n])
                    if curvw !== nothing
                        for iw in 1:size(curvw, 1)
                            Dw[i, j, fi, iw] += v * curvw[iw, j]
                        end
                    end
                end
            end
            :dos in tasks && (dos[fi] += δ)
        end
    end

    :NOA in tasks && _gyro_noa!(NOAorb, NOAspn, E, dE, AA, bl, efs, ωc, kweight,
                                eigval_max, S3)
    return out
end

"tildeΩ(ω) for band n: −Σ_{m≠n} 2 Im[A_α(n,m)A_β(m,n)]·Re[w_mn²/(w_mn²−ω²)] (TAS17 Eq. 12)."
function _gyro_curv_w(E, AA, n::Int, bl::Vector{Int}, ωc::Vector{ComplexF64})
    curvw = zeros(length(ωc), 3)
    for m in bl
        m == n && continue
        wmn = E[m] - E[n]
        for i in 1:3
            t = -2.0 * imag(AA[ALPHA_A[i]][n, m] * AA[BETA_A[i]][m, n])
            for (iw, w) in enumerate(ωc)
                curvw[iw, i] += t * real(wmn^2 / (wmn^2 - w^2))
            end
        end
    end
    return curvw
end

# Natural optical activity γ_abc (TAS17 Eq. C12): interband Fermi-sea/surface sum over
# occupied n / unoccupied l pairs with Re-of-complex-frequency regularisation.
function _gyro_noa!(NOAorb, NOAspn, E, dE, AA, bl::Vector{Int}, efs::Vector{Float64},
                    ωc::Vector{ComplexF64}, kweight::Float64, eigval_max::Float64, S3)
    nω = length(ωc)
    for (fi, ef) in enumerate(efs)
        occ = [n for n in bl if E[n] < ef]
        unocc = [n for n in bl if ef <= E[n] < eigval_max]
        (isempty(occ) || isempty(unocc)) && continue
        for n in occ, l in unocc
            wln = E[l] - E[n]
            Borb = zeros(ComplexF64, 3, 3)
            for a in 1:3, c in 1:3
                s = -im * (dE[a, n] + dE[a, l]) * AA[c][n, l]
                for mm in bl
                    s += (E[n] - E[mm]) * AA[a][n, mm] * AA[c][mm, l] -
                         (E[l] - E[mm]) * AA[c][n, mm] * AA[a][mm, l]
                end
                Borb[a, c] = s
            end
            Bspn = S3 === nothing ? nothing : zeros(ComplexF64, 3, 3)
            if S3 !== nothing
                for b in 1:3
                    Bspn[BETA_A[b], ALPHA_A[b]] = -im * S3[b][n, l]
                end
            end
            for (iw, w) in enumerate(ωc)
                mW1 = 1.0 / (wln^2 - w^2)
                mWm = real(mW1) * kweight
                mWe = real(-mW1 * (2.0 * wln^2 * mW1 + 1.0)) * kweight
                for ab in 1:3
                    a, b = ALPHA_A[ab], BETA_A[ab]
                    for c in 1:3
                        NOAorb[ab, c, fi, iw] +=
                            mWm * real(AA[b][l, n] * Borb[a, c] - AA[a][l, n] * Borb[b, c]) +
                            mWe * (dE[c, n] + dE[c, l]) * imag(AA[a][n, l] * AA[b][l, n])
                        Bspn !== nothing && (NOAspn[ab, c, fi, iw] +=
                            mWm * real(AA[b][l, n] * Bspn[a, c] - AA[a][l, n] * Bspn[b, c]))
                    end
                end
            end
        end
    end
    return nothing
end

const _GYRO_SYMHDR = "#                             |                                      symmetric part                                     ||              asymmetric part              |"
const _GYRO_COLHDR = "   # EFERMI(eV)      omega(eV)             xx             yy             zz             xy             xz             yz              x              y              z"
const _GYRO_NOAHDR = "   # EFERMI(eV)      omega(eV)            yzx            zxy            xyz            yzy            yzz            zxz            xyy            xyx            zxx"

"""
    write_gyrotropic(seedname, res; tasks=...) -> seedname

Write the postw90 `seedname-gyrotropic-*.dat` files for the quantities in `res` (from
[`gyrotropic`](@ref)): symmetrised 11-column blocks for D/tildeD/C/K, raw γ_abc columns for
NOA, and the 2-column DOS — exact reference headers and E15.6 formats.
"""
function write_gyrotropic(seedname::AbstractString, res;
                          tasks=(:D0, :Dw, :C, :K, :NOA, :dos), spin::Bool=false)
    efs, ωs = res.fermi_energies, res.freqs
    e15(x) = fortran_e(x, 15, 6)
    sym3(io, T, fi, ω) = begin
        print(io, e15(efs[fi]), e15(ω))
        print(io, e15(T[1, 1]), e15(T[2, 2]), e15(T[3, 3]),
              e15((T[1, 2] + T[2, 1]) / 2), e15((T[1, 3] + T[3, 1]) / 2),
              e15((T[2, 3] + T[3, 2]) / 2), e15((T[2, 3] - T[3, 2]) / 2),
              e15((T[3, 1] - T[1, 3]) / 2), e15((T[1, 2] - T[2, 1]) / 2))
        println(io)
    end
    function write_sym(name, comment, units, blocks)  # blocks: [(ω, T(:,:,fi))]
        open("$seedname-gyrotropic-$name.dat", "w") do io
            println(io, " #", comment)
            println(io, " # in units of [ ", units, " ] ")
            for (ω, Tf) in blocks
                println(io, _GYRO_SYMHDR)
                println(io, _GYRO_COLHDR)
                for fi in 1:length(efs)
                    sym3(io, view(Tf, :, :, fi), fi, ω)
                end
                println(io)
                println(io)
            end
        end
    end
    :D0 in tasks && write_sym("D", "the D tensor -- Eq. 2 of TAS17", "dimensionless",
                              [(0.0, res.D)])
    :Dw in tasks && write_sym("tildeD", "the tildeD tensor -- Eq. 12 of TAS17", "dimensionless",
                              [(ωs[iw], view(res.Dw, :, :, :, iw)) for iw in 1:length(ωs)])
    :C in tasks && write_sym("C", "the C tensor -- Eq. B6 of TAS17", "Ampere/cm", [(0.0, res.C)])
    :K in tasks && write_sym("K_orb", "orbital part of the K tensor -- Eq. 3 of TAS17",
                             "Ampere", [(0.0, res.K_orb)])
    :K in tasks && spin && write_sym("K_spin", "spin part of the K tensor -- Eq. 3 of TAS17",
                                     "Ampere", [(0.0, res.K_spn)])
    noa_order = ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3), (3, 2), (3, 1), (2, 1))
    function write_noa(name, comment, T)
        open("$seedname-gyrotropic-$name.dat", "w") do io
            println(io, " #", comment)
            println(io, " # in units of [ Ang ] ")
            for iw in 1:length(ωs)
                println(io, _GYRO_NOAHDR)
                for fi in 1:length(efs)
                    print(io, e15(efs[fi]), e15(ωs[iw]))
                    for (ab, c) in noa_order
                        print(io, e15(T[ab, c, fi, iw]))
                    end
                    println(io)
                end
                println(io)
                println(io)
            end
        end
    end
    :NOA in tasks && write_noa("NOA_orb", "the tensor \$gamma_{abc}^{orb}\$ (Eq. C12,C14 of TAS17)",
                               res.NOA_orb)
    :NOA in tasks && spin && write_noa("NOA_spin",
                                       "the tensor \$gamma_{abc}^{spin}\$ (Eq. C12,C15 of TAS17)",
                                       res.NOA_spn)
    if :dos in tasks
        open("$seedname-gyrotropic-DOS.dat", "w") do io
            println(io, " #density of states")
            println(io, " # in units of [ eV^{-1}.Ang^{-3} ] ")
            println(io, "  # EFERMI(eV) ")
            for fi in 1:length(efs)
                println(io, e15(efs[fi]), e15(res.dos[fi]))
            end
            println(io)
            println(io)
        end
    end
    return seedname
end
