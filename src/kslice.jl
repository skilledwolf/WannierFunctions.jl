# kslice — quantities on a 2-D slice through the BZ (postw90's kslice module): band energies
# (for Fermi-line plots) and the occupied-manifold Berry curvature at each point of an
# (N1+1)×(N2+1) inclusive grid spanned by two fractional vectors from a corner.

using Printf
using StaticArrays

"""
    kslice(bm; corner, b1, b2, mesh=(50,50), fermi_energy=nothing, curvature=false)
        -> (kpts, coords, bands[, curv])

Grid k(i,j) = corner + (i/N1)·b1 + (j/N2)·b2 (fractional; i = 0..N1 fastest). `coords` are the
2-D Cartesian coordinates in the slice plane (Å⁻¹). With `curvature = true` (needs
`fermi_energy` and the Berry connection), also returns **−**(J0+J1+J2) per point — postw90's
sign convention for `kslice_task = curv`.
"""
function kslice(bm::BerryModel; corner::AbstractVector=[0.0, 0.0, 0.0],
                b1::AbstractVector, b2::AbstractVector, mesh::NTuple{2,Int}=(50, 50),
                fermi_energy::Union{Nothing,Float64}=nothing, curvature::Bool=false)
    c = SVector{3,Float64}(corner...)
    v1 = SVector{3,Float64}(b1...)
    v2 = SVector{3,Float64}(b2...)
    N1, N2 = mesh
    kpts = [c + (i / N1) * v1 + (j / N2) * v2 for j in 0:N2 for i in 0:N1]   # i fastest

    # orthonormal in-plane basis for the 2-D coordinates
    e1c = bm.lattice.B * v1
    e2c = bm.lattice.B * v2
    ê1 = e1c / norm(e1c)
    e2p = e2c - dot(e2c, ê1) * ê1
    ê2 = e2p / norm(e2p)
    coords = [(dot(bm.lattice.B * (k - c), ê1), dot(bm.lattice.B * (k - c), ê2)) for k in kpts]

    np = length(kpts)
    nw = num_wann(bm)
    bands = Matrix{Float64}(undef, nw, np)
    curv = curvature ? Matrix{Float64}(undef, 3, np) : nothing
    Threads.@threads for p in 1:np
        if curvature
            kd = _berry_kdata(bm, kpts[p])
            bands[:, p] = kd.E
            curv[:, p] = -_imf_kdata(kd, fermi_energy)
        else
            E, _ = eig_deleig(bm, kpts[p]; deriv=false)
            bands[:, p] = E
        end
    end
    return curvature ? (kpts, coords, bands, curv) : (kpts, coords, bands)
end

"""
    write_kslice(seedname, coords, bands; curv=nothing)

Write the postw90 kslice output files: `-kslice-coord.dat`, `-kslice-bands.dat`, and (when
given) `-kslice-curv.dat`.
"""
function write_kslice(seedname::AbstractString, coords, bands; curv=nothing)
    open(seedname * "-kslice-coord.dat", "w") do io
        for (x, y) in coords
            println(io, fortran_e(x, 16, 8), fortran_e(y, 16, 8))
        end
    end
    open(seedname * "-kslice-bands.dat", "w") do io
        for p in 1:size(bands, 2), n in 1:size(bands, 1)
            println(io, fortran_e(bands[n, p], 16, 8))
        end
    end
    if curv !== nothing
        open(seedname * "-kslice-curv.dat", "w") do io
            for p in 1:size(curv, 2)
                println(io, fortran_e(curv[1, p], 16, 8), fortran_e(curv[2, p], 16, 8),
                        fortran_e(curv[3, p], 16, 8))
            end
        end
    end
    return seedname
end
