# In-memory bridge: build a `Model` directly from arrays a DFT code holds in memory (overlaps,
# projections, eigenvalues, neighbour geometry) with no `.amn`/`.mmn`/`.eig` file round-trip.
# This is the entry point a Julia DFT code (e.g. DFTK.jl) uses to hand off to the wannieriser;
# the DFTK-specific glue lives in the package extension `ext/WannierFunctionsDFTKExt.jl`.

using StaticArrays

"""
    wannier_model(; unit_cell, kpoints, mp_grid, num_wann, M, A, kpb, gpb, eig=nothing,
                  kmesh_tol=1e-6, seedname="model") -> Model

Construct a wannierisation `Model` from in-memory arrays (the same content as the
`.win`/`.mmn`/`.amn`/`.eig` files, but passed directly):

- `unit_cell` — 3×3 real lattice, columns = a₁,a₂,a₃ (Å).
- `kpoints` — vector of fractional k-points; `mp_grid` — the Monkhorst–Pack dimensions.
- `M[b,b′,nn,k]` — the `.mmn` overlaps `⟨u_{b,k}|u_{b′,k+nn}⟩`; `kpb[nn,k]` the neighbour
  k-index and `gpb[:,nn,k]` the reciprocal-lattice shift folding k+b back into the list.
- `A[b,w,k]` — the `.amn` projections; `num_wann` the target count.
- `eig[b,k]` — band energies (eV); `nothing` for the isolated case with no interpolation.

Feeds `run_wannier(model)` exactly as a file-read model would.
"""
function wannier_model(; unit_cell::AbstractMatrix, kpoints::Vector, mp_grid::NTuple{3,Int},
                       num_wann::Int, M::Array{ComplexF64,4}, A::Array{ComplexF64,3},
                       kpb::Matrix{Int}, gpb::Array{Int,3},
                       eig::Union{Nothing,Matrix{Float64}}=nothing,
                       kmesh_tol::Float64=1e-6, seedname::AbstractString="model")
    lattice = Lattice(SMatrix{3,3,Float64}(unit_cell))
    kfrac = [SVector{3,Float64}(k...) for k in kpoints]
    kgrid = KGrid(kfrac, mp_grid)
    bvectors = build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol=kmesh_tol)
    num_bands = size(M, 1)
    return Model(lattice, kgrid, bvectors, num_bands, num_wann, M, A, eig, String(seedname))
end
