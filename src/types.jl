# Core data structures for the Wannierisation pipeline.
#
# Index conventions (1-based, matching the physics literature and the reference code):
#   m, n  band / Wannier indices
#   k     k-point index (1..nkpt)
#   b     neighbour-shell b-vector slot (1..nntot)
# Arrays put the fast-varying physical index first (column-major friendly).

using StaticArrays
using LinearAlgebra

"""
    Lattice(A)

Real- and reciprocal-space lattice. `A` has the three real-space lattice vectors as its
**columns** (Ångström). The reciprocal lattice `B` (columns b₁,b₂,b₃, units Å⁻¹) satisfies the
crystallographic convention aᵢ·bⱼ = 2π δᵢⱼ, i.e. `B = 2π (A⁻¹)ᵀ`.
"""
struct Lattice
    A::SMatrix{3,3,Float64,9}
    B::SMatrix{3,3,Float64,9}
end

function Lattice(A::AbstractMatrix)
    Am = SMatrix{3,3,Float64}(A)
    # b_i = 2π (a_j × a_k) / V — the adjugate/cross-product form, matching the reference's
    # utility_recip_lattice to the last ULP (an LU-based inv() can differ in the final digit,
    # which shows up in fixed-format file output).
    a1, a2, a3 = Am[:, 1], Am[:, 2], Am[:, 3]
    V = dot(a1, cross(a2, a3))
    B = SMatrix{3,3,Float64}(hcat(TWOPI / V * cross(a2, a3),
                                  TWOPI / V * cross(a3, a1),
                                  TWOPI / V * cross(a1, a2)))
    return Lattice(Am, B)
end

"Unit-cell volume (Ų·Å = ų)."
cell_volume(l::Lattice) = abs(det(l.A))

"""
    KGrid

The Monkhorst–Pack k-point grid: `frac` holds fractional coordinates (length `nkpt`),
`mp_grid` the grid subdivisions.
"""
struct KGrid
    frac::Vector{SVector{3,Float64}}
    mp_grid::NTuple{3,Int}
end

nkpt(k::KGrid) = length(k.frac)

"""
    BVectors

Finite-difference neighbour geometry. For each k-point, `nntot` neighbours are recorded:

- `kpb[b,k]`     : index of the neighbour k-point k+b (in the k list)
- `gpb[:,b,k]`   : integer reciprocal-lattice shift G folding k+b back into the list
- `bvec[:,b,k]`  : Cartesian b-vector (Å⁻¹), = B·(kfrac[kpb] + G − kfrac[k])
- `wb[b,k]`      : finite-difference weight for that b (Ų)

The weights satisfy the B1 completeness relation Σ_b w_b b_α b_β = δ_αβ (PRB 56, 12847).
"""
struct BVectors
    nntot::Int
    kpb::Matrix{Int}                 # (nntot, nkpt)
    gpb::Array{Int,3}                # (3, nntot, nkpt)
    bvec::Array{Float64,3}           # (3, nntot, nkpt)
    wb::Matrix{Float64}              # (nntot, nkpt)
    shells::Vector{Float64}          # distinct shell radii used (Å⁻¹)
    shell_weight::Vector{Float64}    # weight per shell (Ų)
end

"""
    Model

Everything needed to Wannierise a set of bands: geometry, k-mesh, neighbour weights, the
overlap matrices M, the projection matrices A, and (optionally) the band energies.

- `M[m,n,b,k]` = ⟨u_{m,k} | u_{n,k+b}⟩   (num_bands × num_bands × nntot × nkpt)
- `A[m,n,k]`   = ⟨ψ_{m,k} | g_n⟩         (num_bands × num_wann × nkpt), trial-projection overlaps
- `eig[m,k]`   = ε_{m,k} in eV           (num_bands × nkpt), optional (disentanglement/interp)
"""
mutable struct Model
    lattice::Lattice
    kgrid::KGrid
    bvectors::BVectors
    num_bands::Int
    num_wann::Int
    M::Array{ComplexF64,4}
    A::Array{ComplexF64,3}
    eig::Union{Nothing,Matrix{Float64}}
    seedname::String
end
