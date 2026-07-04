# Spin expectation values on Wannier-interpolated bands (postw90's spin module): the total
# spin magnetic moment and per-band ⟨σ·n̂⟩_nk (used for the spin-decomposed DOS and for
# band-structure colouring). Exact conventions in docs/reference-notes/spin-module.md.
#
# The operator is the Pauli vector σ (eigenvalues in [−1, 1], no ħ/2 and no g-factor);
# the moment is m = −(1/N_k) Σ_k Σ_{n: ε_nk < E_F} ⟨σ⟩_nk in Bohr magnetons per cell.

using LinearAlgebra
using StaticArrays

"""
    SpinModel

`BerryModel` plus the Wannier-gauge Pauli matrices S(R) (3 components) from a `.spn` file.
"""
struct SpinModel
    bm::BerryModel
    SSr::Array{ComplexF64,4}       # (nw, nw, nr, 3)
end

Base.show(io::IO, ::MIME"text/plain", m::SpinModel) =
    print(io, "SpinModel: ", num_wann(m.bm), " WF, ", length(m.bm.irvec), " R-vectors (H, S)")

"""
    SpinModel(bm, chk, spn, kgrid) -> SpinModel

Assemble S(R) from raw `.spn` matrices: `S(q) = v† σ(win_q, win_q) v` on the disentanglement
window, then Fourier to the Wigner–Seitz R-set of `bm`.
"""
function SpinModel(bm::BerryModel, chk::Checkpoint, spn::Array{ComplexF64,4}, kgrid::KGrid)
    nb, nk = size(spn, 1), size(spn, 3)
    nw = num_wann(chk)
    vs, winidx = gauge_v_windows(chk, nb)
    SSr = zeros(ComplexF64, nw, nw, length(bm.irvec), 3)
    for s in 1:3
        SSq = Array{ComplexF64,3}(undef, nw, nw, nk)
        for q in 1:nk
            SSq[:, :, q] = vs[q]' * spn[winidx[q], winidx[q], q, s] * vs[q]
        end
        SSr[:, :, :, s] = fourier_q_to_R(SSq, kgrid, bm.irvec)
    end
    return SpinModel(bm, SSr)
end

"""
    SpinModel(seedname) -> SpinModel

Assemble from `seedname.{win,eig,chk,spn}` (plus `.mmn` for the Berry connection when
present — not needed for the spin quantities themselves).
"""
function SpinModel(seedname::AbstractString)
    win = read_win(seedname * ".win")
    chk = isfile(seedname * ".chk") ? read_chk(seedname * ".chk") :
          read_chk_fmt(seedname * ".chk.fmt")
    eig = read_eig(seedname * ".eig")
    lattice = Lattice(win.unit_cell)
    kgrid = KGrid(win.kpoints, win.mp_grid)
    bm = BerryModel(chk, eig, nothing, kgrid, lattice; use_ws_distance=win.use_ws_distance)
    spn = read_spn(seedname * ".spn"; num_bands=num_bands(chk), num_kpts=nkpt(kgrid))
    return SpinModel(bm, chk, spn, kgrid)
end

"Diagonal ⟨σ·n̂⟩ in the eigenbasis U at k (spin axis from polar/azimuth angles in degrees)."
function _spin_diag(sm::SpinModel, kf::SVector{3,Float64}, U::AbstractMatrix;
                    polar::Float64=0.0, azimuth::Float64=0.0)
    n̂ = SVector(sind(polar) * cosd(azimuth), sind(polar) * sind(azimuth), cosd(polar))
    nw = num_wann(sm.bm)
    Sn = zeros(ComplexF64, nw, nw)
    for s in 1:3
        n̂[s] == 0.0 && continue
        Sn .+= n̂[s] .* _ft_op(sm.bm, (@view sm.SSr[:, :, :, s]), kf)
    end
    Sh = U' * Sn * U
    return [real(Sh[n, n]) for n in 1:nw]
end

"""
    spin_expectation(sm, kfrac; polar=0.0, azimuth=0.0) -> (E, spn)

Band energies and per-band ⟨σ·n̂⟩_nk ∈ [−1, 1] at one k-point.
"""
function spin_expectation(sm::SpinModel, kfrac::AbstractVector;
                          polar::Float64=0.0, azimuth::Float64=0.0)
    kf = SVector{3,Float64}(kfrac...)
    E, _, U = eig_deleig_vec(sm.bm, kf; deriv=false)
    return E, _spin_diag(sm, kf, U; polar=polar, azimuth=azimuth)
end

"""
    spin_moment(sm; fermi_energy, kmesh=(25,25,25)) -> (; moment, theta, phi)

Total spin magnetic moment (Bohr magnetons per cell) and its polar/azimuthal angles in
degrees (postw90 convention: θ = acos(m_z/|m|), φ = atan(m_y/m_x) — plain atan, as in pwscf).
Occupations are the T = 0 step with strict ε < E_F.
"""
function spin_moment(sm::SpinModel; fermi_energy::Float64, kmesh::NTuple{3,Int}=(25, 25, 25))
    nktot = prod(kmesh)
    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    per_k = Vector{SVector{3,Float64}}(undef, nktot)
    Threads.@threads for idx in 1:nktot
        per_k[idx] = _spin_moment_k(sm, kl[idx], fermi_energy)
    end
    m = -sum(per_k) / nktot
    theta = acosd(m[3] / norm(m))
    phi = atand(m[2] / m[1])
    return (; moment=m, theta, phi)
end

function _spin_moment_k(sm::SpinModel, kf::SVector{3,Float64}, ef::Float64)
    E, _, U = eig_deleig_vec(sm.bm, kf; deriv=false)
    occ = findall(<(ef), E)
    acc = MVector{3,Float64}(0.0, 0.0, 0.0)
    for s in 1:3
        Sh = U' * _ft_op(sm.bm, (@view sm.SSr[:, :, :, s]), kf) * U
        for n in occ
            acc[s] += real(Sh[n, n])
        end
    end
    return SVector(acc)
end
