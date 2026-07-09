# Tetrahedron-method spin Hall conductivity (Ghim–Park PRB 106, 075126 (2022) with the
# Kawamura PRB 89, 094515 (2014) optimized-tetrahedron correction), as in postw90's
# berry_get_shc_tetrahedron. Replaces the Gaussian energy-smearing of shc_fermiscan with an
# analytic per-tetrahedron integration of the band-resolved integrand imjv (see shc_imjv),
# needing no smearing width. Exact conventions in docs/reference-notes/tetrahedron-shc.md.
#
# The pure-numerical integration kernels (P-matrix, stencil, tet_spinhall, tet_integral, …) are
# a verbatim port of tetrahedron.F90 in tetrahedron_kernels.jl.

using LinearAlgebra
using StaticArrays

const MESH_SHIFT = 0.5      # half-cell grid shift (avoids Γ)

# imjv exactly as berry_get_shc_tetrahedron builds it (Qiao method): the bare-velocity
# anticommutator spin-velocity ½{v_α, σ_γ} with the SR_R/SHR_R corrections, contracted with
# the covariant velocity v_β = VV0_β − i A_β (ε_m−ε_n). This differs on the discrete grid from
# the QZYZ18 D_h form of shc_imjv (equal only at exact interpolation), so the tetrahedron path
# must use this construction to reproduce the reference.
function _tetra_imjv(sm::ShcModel, kf::SVector{3,Float64}, γ::Int, α::Int, β::Int,
                     kw::BerryKWork=BerryKWork(num_wann(sm.bm)))
    bm = sm.bm
    nw = num_wann(bm)
    kd = _berry_kdata!(kw, bm, kf)
    E, U = kd.E, kd.U
    VV0α = kd.dHh[α]                                    # U†∂_αH U (bare velocity)
    VV0β = kd.dHh[β]
    SS = U' * _ft_op(bm, (@view sm.SSr[:, :, :, γ]), kf) * U
    SAA = U' * _ft_op(bm, (@view sm.SRr[:, :, :, γ, α]), kf) * U
    SBB = U' * _ft_op(bm, (@view sm.SHRr[:, :, :, γ, α]), kf) * U
    AAβ = U' * kd.A[β] * U                              # rotated Berry connection (no D_h)
    spinVel0 = VV0α * SS .+ SS * VV0α
    spinVel = similar(spinVel0)
    for m in 1:nw, n in 1:nw
        spinVel[n, m] = (spinVel0[n, m]
                         - im * (E[m] * SAA[n, m] - SBB[n, m])
                         + im * (E[n] * conj(SAA[m, n]) - conj(SBB[m, n]))) / 2
    end
    imjv = zeros(nw, nw)
    for m in 1:nw, n in 1:nw
        vβ = VV0β[m, n] - im * AAβ[m, n] * (E[n] - E[m])   # VV(m,n,β), reference berry.F90:3232
        imjv[n, m] = imag(spinVel[n, m] * vβ)
    end
    return E, imjv
end

# Build imjv[n,m] and eig on the padded cell-offset grid cx,cy,cz ∈ -1 .. mesh+1
# (k = (offset+0.5)/mesh, wrapped). Flat arrays (nw[, nw], nx, ny, nz): one contiguous
# allocation instead of O(nx·ny·nz) heap objects, and cache-friendly cube gathers.
function _tetra_grid(sm, kmesh::NTuple{3,Int}, γ::Int, α::Int, β::Int)
    m1, m2, m3 = kmesh
    ax = collect(-1:m1+1); ay = collect(-1:m2+1); az = collect(-1:m3+1)
    nx, ny, nz = length(ax), length(ay), length(az)
    nw = num_wann(sm.bm)
    imjv = Array{Float64,5}(undef, nw, nw, nx, ny, nz)
    eig = Array{Float64,4}(undef, nw, nx, ny, nz)
    pts = [(ix, iy, iz) for ix in 1:nx for iy in 1:ny for iz in 1:nz]
    threaded_ksum(
        (st, p) -> begin
            ix, iy, iz = pts[p]
            kf = SVector((ax[ix] + MESH_SHIFT) / m1, (ay[iy] + MESH_SHIFT) / m2,
                         (az[iz] + MESH_SHIFT) / m3)
            E, jv = _tetra_imjv(sm, kf, γ, α, β, st.work)
            eig[:, ix, iy, iz] = E
            imjv[:, :, ix, iy, iz] = jv
        end,
        () -> (work=BerryKWork(nw),),
        length(pts))
    return imjv, eig, ax, ay, az
end

# cell-offset → padded array index (offset o ∈ -1..mesh+1 → index o+2)
@inline _cidx(o::Int) = o + 2

"""
    shc_tetra(sm; kmesh, fermi_energies=nothing, freqs=nothing, fermi_energy=nothing,
              γ=3, α=1, β=2, cutoff=1e-4, avoid_deg=3e-4) -> Vector

Tetrahedron-method spin Hall conductivity in (ħ/e)·S/cm. Pass `fermi_energies` for a Fermi
scan (ω = 0, real output) or `freqs` + `fermi_energy` for a frequency scan (complex output).
`cutoff` = `tetrahedron_cutoff`, `avoid_deg` = `tetrahedron_avoid_degeneracy`.
"""
function shc_tetra(sm::Union{ShcModel,ShcRyooModel}; kmesh::NTuple{3,Int},
                   fermi_energies::Union{Nothing,Vector{Float64}}=nothing,
                   freqs::Union{Nothing,Vector{Float64}}=nothing,
                   fermi_energy::Union{Nothing,Float64}=nothing,
                   γ::Int=3, α::Int=1, β::Int=2,
                   cutoff::Float64=1e-4, avoid_deg::Float64=3e-4)
    freq_scan = freqs !== nothing
    if freq_scan
        fermi_energy !== nothing || error("shc_tetra freq scan needs a single fermi_energy")
        scan = freqs
    else
        fermi_energies !== nothing || error("shc_tetra needs fermi_energies or freqs")
        scan = fermi_energies
    end
    nscan = length(scan)
    bm = sm.bm
    nw = num_wann(bm)
    m1, m2, m3 = kmesh

    imjv, eig, ax, ay, az = _tetra_grid(sm, kmesh, γ, α, β)
    P = tet_p_matrix()
    TA = tet_array()

    # accumulate per interior cube (chunked over threads, per-chunk buffers reused)
    cubes = [(cx, cy, cz) for cx in 0:m1-1 for cy in 0:m2-1 for cz in 0:m3-1]
    ef0 = freq_scan ? fermi_energy : 0.0
    states = threaded_ksum(
        (st, ci) -> _tetra_cube!(st, cubes[ci], imjv, eig, P, TA, kmesh, scan, freq_scan,
                                 ef0, nw, cutoff, avoid_deg),
        () -> (acc=zeros(ComplexF64, nscan), gidx=Matrix{Int}(undef, 3, 64),
               kc=Matrix{Float64}(undef, 3, 64), Fopt=zero(MVector{20,Float64}),
               E1opt=zero(MVector{20,Float64}), E2opt=zero(MVector{20,Float64})),
        length(cubes))
    total = sum(st.acc for st in states)
    fac = 1.0e8 * ELEM_CHARGE_SI^2 / (HBAR_SI * cell_volume(bm.lattice)) / 2.0
    total .*= fac
    return freq_scan ? total : real.(total)
end

# Per-cube body: gather the 64-point super-block (idx = 16*l + 4*k + i + 1, i,k,l ∈ 0..3),
# contract with the Kawamura P-matrix, and integrate the 6 tetrahedra for every band pair.
function _tetra_cube!(st, cube::NTuple{3,Int}, imjv::Array{Float64,5}, eig::Array{Float64,4},
                      P::SMatrix{4,20,Float64,80}, TA::SMatrix{6,20,Int,120},
                      kmesh::NTuple{3,Int}, scan::Vector{Float64}, freq_scan::Bool,
                      ef0::Float64, nw::Int, cutoff::Float64, avoid_deg::Float64)
    m1, m2, m3 = kmesh
    cx, cy, cz = cube
    acc = st.acc
    gidx = st.gidx
    kc = st.kc
    for l in 0:3, k in 0:3, i in 0:3
        idx = 16l + 4k + i + 1
        gidx[1, idx] = _cidx(cx + i - 1)
        gidx[2, idx] = _cidx(cy + k - 1)
        gidx[3, idx] = _cidx(cz + l - 1)
        kc[1, idx] = (cx + i - 1 + MESH_SHIFT) / m1
        kc[2, idx] = (cy + k - 1 + MESH_SHIFT) / m2
        kc[3, idx] = (cz + l - 1 + MESH_SHIFT) / m3
    end
    Fopt = st.Fopt
    E1opt = st.E1opt
    E2opt = st.E2opt
    nscan = length(scan)
    @inbounds for itet in 1:6
        v1, v2, v3, v4 = TA[itet, 1], TA[itet, 2], TA[itet, 3], TA[itet, 4]
        # column kk = kc[:, verts[kk+1]] − kc[:, verts[1]] (column-major constructor)
        tt = SMatrix{3,3,Float64,9}(
            kc[1, v2] - kc[1, v1], kc[2, v2] - kc[2, v1], kc[3, v2] - kc[3, v1],
            kc[1, v3] - kc[1, v1], kc[2, v3] - kc[2, v1], kc[3, v3] - kc[3, v1],
            kc[1, v4] - kc[1, v1], kc[2, v4] - kc[2, v1], kc[3, v4] - kc[3, v1])
        for n in 1:nw, m in 1:nw
            n == m && continue
            for pt in 1:20
                q = TA[itet, pt]
                gx, gy, gz = gidx[1, q], gidx[2, q], gidx[3, q]
                Fopt[pt] = imjv[n, m, gx, gy, gz]
                E1opt[pt] = eig[n, gx, gy, gz]
                E2opt[pt] = eig[m, gx, gy, gz]
            end
            # P-matrix contraction to the 4 corrected corner values (one static matvec)
            Ftet = P * SVector(Fopt)
            E1tet = P * SVector(E1opt)
            E2tet = P * SVector(E2opt)
            for s in 1:nscan
                ω = freq_scan ? scan[s] : 0.0
                Ef = freq_scan ? ef0 : scan[s]
                if ω == 0.0
                    v = tet_spinhall(Ftet, E1tet, E2tet, tt, 0.0, Ef, 3, cutoff, avoid_deg)
                    acc[s] -= v
                else
                    re = (tet_spinhall(Ftet, E1tet, E2tet, tt, -ω, Ef, 1, cutoff, avoid_deg) -
                          tet_spinhall(Ftet, E1tet, E2tet, tt, ω, Ef, 1, cutoff, avoid_deg)) / (2ω)
                    im_ = pi * (tet_spinhall(Ftet, E1tet, E2tet, tt, -ω, Ef, 2, cutoff, avoid_deg) +
                                tet_spinhall(Ftet, E1tet, E2tet, tt, ω, Ef, 2, cutoff, avoid_deg)) / (2ω)
                    acc[s] -= complex(re, im_)
                end
            end
        end
    end
    return nothing
end
