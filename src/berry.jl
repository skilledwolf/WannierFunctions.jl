# Berry curvature and the intrinsic anomalous Hall conductivity, Wannier-interpolated.
#
# Follows postw90's berry module (Wang–Yates–Souza–Vanderbilt PRB 74, 195118 (2006);
# Lopez–Vanderbilt–Thonhauser–Souza PRB 85, 014435 (2012)); exact conventions in
# docs/reference-notes/berry-ahc.md. The occupied-manifold Berry curvature at each k is the
# three-piece sum (axial packing i: 1=(y,z), 2=(z,x), 3=(x,y); α=alpha_A(i), β=beta_A(i)):
#
#   J0_i = Re Tr[f Ω̄_i]                                   (Ω̄ term)
#   J1_i = −2( Im Tr[A_α J⁺_β] + Im Tr[J⁻_α A_β] )        (D–A term)
#   J2_i = −2  Im Tr[J⁻_α J⁺_β]                            (D–D term)
#
# with f the occupation matrix in the Wannier gauge, A the interpolated Berry connection,
# Ω̄ its curl, and J± the energy-denominator (velocity/(E_m−E_n)) matrices. The AHC is the
# k-average scaled by −10⁸ e²/(ħ V_cell) → S/cm.

using LinearAlgebra
using StaticArrays

# CODATA2006 (matches the reference default build; see constants.jl for the bohr choice).
const ELEM_CHARGE_SI = 1.602176487e-19
const HBAR_SI = 1.054571628e-34

const ALPHA_A = (2, 3, 1)     # i=1 → (y,z), 2 → (z,x), 3 → (x,y)
const BETA_A = (3, 1, 2)

"""
    BerryModel

Everything the Berry-curvature interpolation needs, built once from a checkpoint (or a
completed run): H(R), the Berry-connection matrix A(R) (postw90 convention: linear
finite-difference for ALL elements, Hermitised at each q), the Wigner–Seitz R-set, and the
lattice.
"""
struct BerryModel
    lattice::Lattice
    irvec::Vector{NTuple{3,Int}}
    ndegen::Vector{Int}
    Rcart::Vector{SVector{3,Float64}}
    Hr::Array{ComplexF64,3}            # (nw, nw, nr)
    Ar::Array{ComplexF64,4}            # (nw, nw, nr, 3)
end

num_wann(b::BerryModel) = size(b.Hr, 1)

function Base.show(io::IO, ::MIME"text/plain", b::BerryModel)
    print(io, "BerryModel: ", num_wann(b), " WF, ", length(b.irvec), " R-vectors")
end

"""
    BerryModel(chk::Checkpoint, eig, bv::BVectors, kgrid::KGrid, lattice) -> BerryModel

Assemble H(R) and A(R) from a Wannier90 checkpoint plus the band energies and the b-vector
geometry (from the `.mmn` connectivity, so slot order matches `chk.m_matrix`).

- `H(q) = v† diag(ε_window) v` with `v = U_opt·U` (window rows) for disentangled runs,
  `U† diag(ε) U` otherwise.
- `A_α(q) = Σ_b i w_b b_α M̃(q,b)` for **all** elements (postw90 default, `transl_inv=F` —
  NB this differs from the `_tb.dat` convention which uses −w_b b Im ln M̃_nn on the diagonal),
  then Hermitised: A ← (A + A†)/2.
- Both transformed with `O(R) = (1/N_q) Σ_q e^{−i2πq·R} O(q)` on the Wigner–Seitz R-set.
"""
function BerryModel(chk::Checkpoint, eig::Matrix{Float64}, bv::BVectors, kgrid::KGrid,
                    lattice::Lattice)
    nw = num_wann(chk)
    nk = nkpt(kgrid)
    nb = num_bands(chk)

    # H(q)
    Hq = Array{ComplexF64,3}(undef, nw, nw, nk)
    for q in 1:nk
        if chk.have_disentangled
            nd = chk.ndimwin[q]
            v = chk.u_matrix_opt[1:nd, :, q] * chk.u_matrix[:, :, q]     # ndimwin × nw
            εw = eig[findall(@view chk.lwindow[:, q]), q]
            Hq[:, :, q] = v' * Diagonal(εw) * v
        else
            Uq = @view chk.u_matrix[:, :, q]
            Hq[:, :, q] = Uq' * Diagonal(@view eig[:, q]) * Uq
        end
    end

    # A(q): linear finite difference on the final-gauge overlaps, then Hermitise.
    Aq = zeros(ComplexF64, nw, nw, nk, 3)
    for q in 1:nk, b in 1:bv.nntot
        w = bv.wb[b, q]
        bb = SVector{3,Float64}(bv.bvec[1, b, q], bv.bvec[2, b, q], bv.bvec[3, b, q])
        @views for c in 1:3
            Aq[:, :, q, c] .+= (im * w * bb[c]) .* chk.m_matrix[:, :, b, q]
        end
    end
    for q in 1:nk, c in 1:3
        @views Aq[:, :, q, c] .= (Aq[:, :, q, c] .+ Aq[:, :, q, c]') ./ 2
    end

    irvec, ndegen = wigner_seitz(lattice, kgrid.mp_grid)
    nr = length(irvec)
    Hr = zeros(ComplexF64, nw, nw, nr)
    Ar = zeros(ComplexF64, nw, nw, nr, 3)
    Threads.@threads for ir in 1:nr
        R = SVector{3,Float64}(irvec[ir]...)
        for q in 1:nk
            fac = cis(-TWOPI * dot(kgrid.frac[q], R)) / nk
            @views Hr[:, :, ir] .+= fac .* Hq[:, :, q]
            @views for c in 1:3
                Ar[:, :, ir, c] .+= fac .* Aq[:, :, q, c]
            end
        end
    end
    Rcart = [lattice.A * SVector{3,Float64}(r...) for r in irvec]
    return BerryModel(lattice, irvec, ndegen, Rcart, Hr, Ar)
end

"""
    BerryModel(seedname) -> BerryModel

Convenience: assemble from `seedname.{win,mmn,eig,chk}` (a completed Wannier90 or
WannierFunctions run).
"""
function BerryModel(seedname::AbstractString)
    win = read_win(seedname * ".win")
    chk = isfile(seedname * ".chk") ? read_chk(seedname * ".chk") :
          read_chk_fmt(seedname * ".chk.fmt")
    eig = read_eig(seedname * ".eig")
    _, kpb, gpb, _, nk, _ = read_mmn(seedname * ".mmn")
    lattice = Lattice(win.unit_cell)
    kgrid = KGrid(win.kpoints, win.mp_grid)
    bv = build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol=win.kmesh_tol)
    return BerryModel(chk, eig, bv, kgrid, lattice)
end

"Interpolated k-space data shared by all Fermi levels at one k: E, U, A, Ω̄, rotated velocity."
function _berry_kdata(bm::BerryModel, kfrac::AbstractVector)
    nw = num_wann(bm)
    H = zeros(ComplexF64, nw, nw)
    dH = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    A = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    Ω̄ = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    kf = SVector{3,Float64}(kfrac...)
    for ir in 1:length(bm.irvec)
        fac = cis(TWOPI * dot(kf, SVector{3,Float64}(bm.irvec[ir]...))) / bm.ndegen[ir]
        R = bm.Rcart[ir]
        @views H .+= fac .* bm.Hr[:, :, ir]
        @views for c in 1:3
            dH[c] .+= (fac * im * R[c]) .* bm.Hr[:, :, ir]
            A[c] .+= fac .* bm.Ar[:, :, ir, c]
        end
        # Ω̄_i = Σ_R e^{ikR} i (R_α A_β − R_β A_α)
        @views for i in 1:3
            α, β = ALPHA_A[i], BETA_A[i]
            Ω̄[i] .+= (fac * im) .* (R[α] .* bm.Ar[:, :, ir, β] .- R[β] .* bm.Ar[:, :, ir, α])
        end
    end
    F = eigen(Hermitian((H + H') / 2))
    E, U = F.values, F.vectors
    dHh = [U' * dH[c] * U for c in 1:3]     # velocity in the Hamiltonian gauge
    return (; E, U, A, Ω̄, dHh)
end

"Occupied-manifold curvature (J0+J1+J2, Ų) at one Fermi level from shared k-data."
function _imf_kdata(kd, fermi_energy::Float64)
    E, U = kd.E, kd.U
    nw = length(E)
    occ = Float64.(E .< fermi_energy)
    f = U * Diagonal(occ) * U'
    JJp = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    JJm = [zeros(ComplexF64, nw, nw) for _ in 1:3]
    for c in 1:3
        for m in 1:nw, n in 1:nw
            if E[n] > fermi_energy && E[m] < fermi_energy       # n empty, m occupied
                JJp[c][n, m] = im * kd.dHh[c][n, m] / (E[m] - E[n])
                JJm[c][m, n] = im * kd.dHh[c][m, n] / (E[n] - E[m])
            end
        end
        JJp[c] = U * JJp[c] * U'
        JJm[c] = U * JJm[c] * U'
    end
    imf = MVector{3,Float64}(0, 0, 0)
    for i in 1:3
        α, β = ALPHA_A[i], BETA_A[i]
        J0 = real(tr(f * kd.Ω̄[i]))
        J1 = -2.0 * (imag(tr(kd.A[α] * JJp[β])) + imag(tr(JJm[α] * kd.A[β])))
        J2 = -2.0 * imag(tr(JJm[α] * JJp[β]))
        imf[i] = J0 + J1 + J2
    end
    return SVector(imf)
end

"""
    berry_curvature_k(bm, kfrac, fermi_energy) -> SVector{3,Float64}

Occupied-manifold Berry curvature −2 Im f(k) (Ų, axial components yz/zx/xy) at one
fractional k-point, as the J0+J1+J2 sum.
"""
berry_curvature_k(bm::BerryModel, kfrac::AbstractVector, fermi_energy::Float64) =
    _imf_kdata(_berry_kdata(bm, kfrac), fermi_energy)

"""
    ahc_fermiscan(bm; fermi_energies, kmesh, adpt_kmesh=1, adpt_thresh=100.0) -> Matrix

AHC (S/cm) for a list of Fermi energies (3 × nf matrix), with postw90's adaptive k-mesh
refinement: any coarse point whose curvature norm exceeds `adpt_thresh` (Ų) at some Fermi
level is re-evaluated on an `adpt_kmesh`³ sub-mesh centred on it (replacing the coarse value
for the triggered levels only). `adpt_kmesh = 1` disables refinement.
"""
function ahc_fermiscan(bm::BerryModel; fermi_energies::Vector{Float64},
                       kmesh::NTuple{3,Int}=(25, 25, 25),
                       adpt_kmesh::Int=1, adpt_thresh::Float64=100.0)
    nf = length(fermi_energies)
    nktot = prod(kmesh)
    db = SVector(1.0 / kmesh[1], 1.0 / kmesh[2], 1.0 / kmesh[3])
    na = adpt_kmesh
    # sub-mesh offsets centred on the coarse point (berry.F90:617-626)
    adkpt = [SVector(db[1] * ((i + 0.5) / na - 0.5),
                     db[2] * ((j + 0.5) / na - 0.5),
                     db[3] * ((k + 0.5) / na - 0.5))
             for i in 0:na-1 for j in 0:na-1 for k in 0:na-1]

    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    kw = 1.0 / nktot
    per_k = Vector{Matrix{Float64}}(undef, nktot)
    Threads.@threads for idx in 1:nktot
        acc = zeros(3, nf)
        kd = _berry_kdata(bm, kl[idx])
        ladpt = falses(nf)
        for (ife, ef) in enumerate(fermi_energies)
            imf = _imf_kdata(kd, ef)
            if na > 1 && norm(imf) > adpt_thresh
                ladpt[ife] = true                       # coarse value discarded
            else
                acc[:, ife] .+= imf .* kw
            end
        end
        if any(ladpt)
            kwa = kw / na^3
            for off in adkpt
                kda = _berry_kdata(bm, kl[idx] + off)
                for (ife, ef) in enumerate(fermi_energies)
                    ladpt[ife] || continue
                    acc[:, ife] .+= _imf_kdata(kda, ef) .* kwa
                end
            end
        end
        per_k[idx] = acc
    end
    imf_tot = sum(per_k)
    fac = -1.0e8 * ELEM_CHARGE_SI^2 / (HBAR_SI * cell_volume(bm.lattice))
    return fac .* imf_tot
end

"""
    anomalous_hall(bm; fermi_energy, kmesh) -> SVector{3,Float64}

Intrinsic anomalous Hall conductivity (S/cm, axial components x/y/z) as the Berry-curvature
average over a uniform `kmesh` (postw90's `berry_task = ahc` with no adaptive refinement).
"""
function anomalous_hall(bm::BerryModel; fermi_energy::Float64,
                        kmesh::NTuple{3,Int}=(25, 25, 25))
    out = ahc_fermiscan(bm; fermi_energies=[fermi_energy], kmesh=kmesh)
    return SVector(out[1, 1], out[2, 1], out[3, 1])
end
