# Kubo–Greenwood optical conductivity and joint density of states, Wannier-interpolated
# (postw90's `berry_task = kubo`; WYSV PRB 74 195118 (2006) Sec. VI, YWVS PRB 75 195121 (2007)
# adaptive smearing). Exact conventions in docs/reference-notes/kubo-morb-geninterp.md.
#
# Interband, per k and band pair n≠m (T = 0 occupations f at the Fermi energy):
#   σ^H_ij(ω)  += −π (f_m−f_n)(ε_m−ε_n) δ_η(ε_m−ε_n−ω) A_nm,i A_mn,j     (Hermitian part)
#   σ^AH_ij(ω) += i (f_m−f_n)(ε_m−ε_n)/(ε_m−ε_n−ω_c) A_nm,i A_mn,j       (anti-Hermitian)
#   jdos(ω)    += f_n (1−f_m) δ_η(ε_m−ε_n−ω)
# with A = U†A_W U + i D_h (WYSV Eq. 25), D_h(n,m) = (U†∂HU)_nm/(ε_m−ε_n), Gaussian δ, and
# per-pair adaptive width η = min(fac·|∇ε_m−∇ε_n|·Δk, η_max). σ scaled by 10⁸e²/(ħV) → S/cm.

using LinearAlgebra
using StaticArrays

"Gaussian delta representation: w0gauss(x, 0) = e^{−x²}/√π."
_w0gauss(x::Float64) = abs(x) > 6.0 ? 0.0 : exp(-x^2) / sqrt(π)

"Largest reciprocal-mesh spacing max_i |b_i|/N_i (postw90's kmesh_spacing)."
function kmesh_spacing(lattice::Lattice, mesh::NTuple{3,Int})
    return maximum(norm(lattice.B[:, i]) / mesh[i] for i in 1:3)
end

"""
    KuboResult

Frequency grid (eV), Hermitian and anti-Hermitian conductivity tensors (3×3×nfreq, S/cm), and
the joint density of states (states/eV, unscaled).
"""
struct KuboResult
    freqs::Vector{Float64}
    H::Array{ComplexF64,3}
    AH::Array{ComplexF64,3}
    jdos::Vector{Float64}
end

"σ_S(i,j,ω) = Re[(H_ij+H_ji)/2] + i·Im[(AH_ij+AH_ji)/2] — the symmetric output combination."
kubo_S(r::KuboResult, i, j, f) = complex(real(r.H[i, j, f] + r.H[j, i, f]) / 2,
                                         imag(r.AH[i, j, f] + r.AH[j, i, f]) / 2)
"σ_A(i,j,ω) = Re[(AH_ij−AH_ji)/2] + i·Im[(H_ij−H_ji)/2] — the antisymmetric combination."
kubo_A(r::KuboResult, i, j, f) = complex(real(r.AH[i, j, f] - r.AH[j, i, f]) / 2,
                                         imag(r.H[i, j, f] - r.H[j, i, f]) / 2)

"""
    optical_conductivity(bm; fermi_energy, kmesh, freqs=0.0:0.01:..., eigval_max=Inf,
                         adaptive=true, adpt_fac=√2, adpt_max=1.0, smr_width=0.0,
                         ) -> KuboResult

Interband optical conductivity σ(ħω) and JDOS on a uniform k-mesh, with postw90's defaults
(Gaussian broadening; per-pair adaptive width from the band-velocity difference).
"""
function optical_conductivity(bm::BerryModel; fermi_energy::Float64,
                              kmesh::NTuple{3,Int}=(25, 25, 25),
                              freqs::AbstractVector{<:Real}=0.0:0.01:7.0,
                              eigval_max::Float64=Inf,
                              adaptive::Bool=true, adpt_fac::Float64=sqrt(2.0),
                              adpt_max::Float64=1.0, smr_width::Float64=0.0)
    nf = length(freqs)
    ωs = collect(Float64, freqs)
    nktot = prod(kmesh)
    Δk = kmesh_spacing(bm.lattice, kmesh)
    nw = num_wann(bm)

    kl = [SVector(i / kmesh[1], j / kmesh[2], k / kmesh[3])
          for i in 0:kmesh[1]-1 for j in 0:kmesh[2]-1 for k in 0:kmesh[3]-1]
    states = threaded_ksum(
        (st, idx) -> _kubo_kpoint!(st, bm, kl[idx], fermi_energy, ωs, Δk, eigval_max,
                                   adaptive, adpt_fac, adpt_max, smr_width),
        () -> (work=BerryKWork(nw), H=zeros(ComplexF64, 3, 3, nf),
               AH=zeros(ComplexF64, 3, 3, nf), J=zeros(nf),
               A=[zeros(ComplexF64, nw, nw) for _ in 1:3]),
        nktot)
    fac = 1.0e8 * ELEM_CHARGE_SI^2 / (HBAR_SI * cell_volume(bm.lattice))
    Hs = sum(st.H for st in states) .* (fac / nktot)
    AHs = sum(st.AH for st in states) .* (fac / nktot)
    Js = sum(st.J for st in states) ./ nktot
    return KuboResult(ωs, Hs, AHs, Js)
end

# Per-k body of optical_conductivity (top-level: closure-boxing hazard in threaded bodies).
function _kubo_kpoint!(st, bm::BerryModel, kf::SVector{3,Float64}, fermi_energy::Float64,
                       ωs::Vector{Float64}, Δk::Float64, eigval_max::Float64, adaptive::Bool,
                       adpt_fac::Float64, adpt_max::Float64, smr_width::Float64)
    w = st.work
    nw = num_wann(bm)
    nf = length(ωs)
    kd = _berry_kdata!(w, bm, kf)
    E, U = kd.E, kd.U
    # band velocities from the diagonal of the rotated velocity matrices
    dE = [real(kd.dHh[c][n, n]) for c in 1:3, n in 1:nw]
    # A(k) in the Hamiltonian gauge: U†A_W U + i D_h, with D_h(n,m) = (U†∂HU)_nm/(ε_m−ε_n)
    A = st.A
    for c in 1:3
        mul!(w.tmp, U', kd.A[c])
        mul!(A[c], w.tmp, U)
        for m in 1:nw, n in 1:nw
            (n != m && abs(E[m] - E[n]) > 1e-7) &&
                (A[c][n, m] += im * kd.dHh[c][n, m] / (E[m] - E[n]))
        end
    end
    map!(e -> Float64(e < fermi_energy), w.occ, E)
    occ = w.occ

    H = st.H; AH = st.AH; J = st.J
    for m in 1:nw, n in 1:nw
        n == m && continue
        (E[m] > eigval_max || E[n] > eigval_max) && continue
        η = adaptive ?
            min(adpt_fac * norm(SVector(dE[1, m] - dE[1, n], dE[2, m] - dE[2, n],
                                        dE[3, m] - dE[3, n])) * Δk, adpt_max) : smr_width
        η <= 0 && (η = 1e-10)
        rfac1 = (occ[m] - occ[n]) * (E[m] - E[n])
        occ_prod = occ[n] * (1.0 - occ[m])
        for f in 1:nf
            ωc = complex(ωs[f], adaptive ? η : smr_width)
            δ = _w0gauss((E[m] - E[n] - ωs[f]) / η) / η
            J[f] += occ_prod * δ
            rfac2 = -π * rfac1 * δ
            cfac = im * rfac1 / (E[m] - E[n] - ωc)
            for j in 1:3, i in 1:3
                aa = A[i][n, m] * A[j][m, n]
                H[i, j, f] += rfac2 * aa
                AH[i, j, f] += cfac * aa
            end
        end
    end
    return nothing
end
