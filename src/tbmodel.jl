# Tight-binding model input: read a `_hr.dat` or `_tb.dat` file back as an interpolation
# source, so the whole postw90 stack (bands, DOS, Berry curvature, AHC, …) runs on a model
# Hamiltonian without a DFT run — the WannierBerri `System_tb` use case. The readers are the
# exact inverse of `write_hr` / `write_tb` (src/output.jl); a `_tb.dat` additionally supplies
# the position operator r(R), enabling the Berry-connection quantities.

using LinearAlgebra
using StaticArrays

"""
    read_tb(path) -> (lattice, num_wann, irvec, ndegen, Hr, pos)

Read a `seedname_tb.dat`: the lattice vectors (Å), the Wannier Hamiltonian `Hr[j,i,ir]` =
⟨0j|H|R i⟩, and the position operator `pos[j,i,ir,α]` = ⟨0j|r_α|R i⟩ (Å). Inverse of
[`write_tb`](@ref).
"""
function read_tb(path::AbstractString)
    lines = readlines(path)
    idx = 2                                        # skip header comment
    A = Matrix{Float64}(undef, 3, 3)
    for k in 1:3
        A[:, k] = parse.(Float64, split(strip(lines[idx])))   # column k = a_k
        idx += 1
    end
    lattice = Lattice(SMatrix{3,3,Float64}(A))
    num_wann = parse(Int, strip(lines[idx])); idx += 1
    nrpts = parse(Int, strip(lines[idx])); idx += 1
    ndegen = Int[]
    while length(ndegen) < nrpts
        append!(ndegen, parse.(Int, split(strip(lines[idx]))))
        idx += 1
    end
    length(ndegen) == nrpts || error("read_tb: parsed $(length(ndegen)) ndegen, expected $nrpts")

    irvec = Vector{NTuple{3,Int}}(undef, nrpts)
    Hr = zeros(ComplexF64, num_wann, num_wann, nrpts)
    for ir in 1:nrpts
        isempty(strip(lines[idx])) && (idx += 1)   # blank line before each R block
        t = split(strip(lines[idx])); idx += 1
        irvec[ir] = (parse(Int, t[1]), parse(Int, t[2]), parse(Int, t[3]))
        for i in 1:num_wann, j in 1:num_wann
            tt = split(strip(lines[idx])); idx += 1
            jj, ii = parse(Int, tt[1]), parse(Int, tt[2])
            Hr[jj, ii, ir] = complex(parse(Float64, tt[3]), parse(Float64, tt[4]))
        end
    end
    pos = zeros(ComplexF64, num_wann, num_wann, nrpts, 3)
    for ir in 1:nrpts
        isempty(strip(lines[idx])) && (idx += 1)
        idx += 1                                   # R line (same order as H block)
        for i in 1:num_wann, j in 1:num_wann
            tt = split(strip(lines[idx])); idx += 1
            jj, ii = parse(Int, tt[1]), parse(Int, tt[2])
            for α in 1:3
                pos[jj, ii, ir, α] = complex(parse(Float64, tt[2+2α-1]), parse(Float64, tt[2+2α]))
            end
        end
    end
    return lattice, num_wann, irvec, ndegen, Hr, pos
end

"""
    BerryModel(lattice, irvec, ndegen, Hr; pos=nothing) -> BerryModel

Build an interpolation model directly from a Wannier Hamiltonian H(R) (and, optionally, the
position operator r(R) as `pos[j,i,ir,α]`). With `pos`, the full Berry-connection stack (AHC,
Kubo, curvature, kpath curv/morb) is available; without it, only H(R)-derived quantities
(bands, DOS, velocities, BoltzWann, geninterp) work. This is the TB-model entry point — the R
list carries its own Wigner–Seitz degeneracies, so no k-mesh or checkpoint is needed.
"""
function BerryModel(lattice::Lattice, irvec::Vector{NTuple{3,Int}}, ndegen::Vector{Int},
                    Hr::Array{ComplexF64,3};
                    pos::Union{Nothing,Array{ComplexF64,4}}=nothing)
    nw = size(Hr, 1)
    Rcart = [lattice.A * SVector{3,Float64}(r...) for r in irvec]
    Ar = pos === nothing ? zeros(ComplexF64, nw, nw, length(irvec), 0) :
         Array{ComplexF64,4}(pos)
    return BerryModel(lattice, irvec, ndegen, Rcart, Hr, Ar, nothing)
end

"""
    tb_model(path; lattice=nothing) -> BerryModel

Load a tight-binding model as a `BerryModel`. A `_tb.dat` supplies the lattice and r(R)
automatically; a `_hr.dat` needs the `lattice` (a `Lattice`) passed explicitly and yields an
H(R)-only model.
"""
function tb_model(path::AbstractString; lattice::Union{Nothing,Lattice}=nothing)
    if endswith(path, "_tb.dat") || endswith(path, ".tb.dat")
        lat, _, irvec, ndegen, Hr, pos = read_tb(path)
        return BerryModel(lat, irvec, ndegen, Hr; pos=pos)
    else
        lattice === nothing &&
            error("tb_model: reading a _hr.dat needs `lattice` (a _hr.dat carries no cell)")
        _, irvec, ndegen, Hr = read_hr(path)
        return BerryModel(lattice, irvec, ndegen, Hr)
    end
end
