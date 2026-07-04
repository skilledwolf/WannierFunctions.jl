# Finite-difference b-vectors and their B1 completeness weights.
#
# The neighbour connectivity (which k-point is k+b, and the reciprocal fold G) is taken from the
# .mmn file so that the b-vector ordering aligns exactly with the stored overlap matrices. The
# *weights* are solved here from the B1 completeness relation Σ_b w_b b_α b_β = δ_αβ
# (Marzari–Vanderbilt, PRB 56, 12847 (1997), Eq. B1), which is the gauge-independent "science"
# of the k-mesh and must reproduce the reference shell weights.

using LinearAlgebra
using StaticArrays

"""
    build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol) -> BVectors

Compute Cartesian b-vectors for the given neighbour connectivity and solve the B1 relation for
the per-shell finite-difference weights. Verifies completeness and errors if it is not met.
"""
function build_bvectors(kgrid::KGrid, lattice::Lattice,
                        kpb::Matrix{Int}, gpb::Array{Int,3};
                        kmesh_tol::Float64=KMESH_TOL_DEFAULT)
    nk = nkpt(kgrid)
    nntot = size(kpb, 1)
    bvec = Array{Float64,3}(undef, 3, nntot, nk)
    for k in 1:nk, b in 1:nntot
        dfrac = kgrid.frac[kpb[b, k]] + SVector{3,Float64}(gpb[1, b, k], gpb[2, b, k], gpb[3, b, k]) - kgrid.frac[k]
        bvec[:, b, k] = lattice.B * dfrac
    end

    # Distinct shell radii from the k=1 neighbour set (identical across k for a uniform mesh).
    b1(b) = SVector{3,Float64}(bvec[1, b, 1], bvec[2, b, 1], bvec[3, b, 1])
    radii = Float64[]
    for b in 1:nntot
        r = norm(b1(b))
        any(x -> abs(x - r) < kmesh_tol, radii) || push!(radii, r)
    end
    sort!(radii)
    shell_of(r) = findfirst(x -> abs(x - r) < kmesh_tol, radii)::Int
    nsh = length(radii)

    # Least-squares solve of the six independent components of the B1 tensor relation.
    Amat = zeros(6, nsh)
    for b in 1:nntot
        v = b1(b)
        s = shell_of(norm(v))
        Amat[1, s] += v[1]^2; Amat[2, s] += v[2]^2; Amat[3, s] += v[3]^2
        Amat[4, s] += v[1]*v[2]; Amat[5, s] += v[1]*v[3]; Amat[6, s] += v[2]*v[3]
    end
    target = SVector{6,Float64}(1, 1, 1, 0, 0, 0)
    shell_weight = Amat \ Vector(target)             # SVD/QR least squares

    # Verify B1 completeness at k=1.
    T = zeros(3, 3)
    for b in 1:nntot
        v = b1(b); w = shell_weight[shell_of(norm(v))]
        T .+= w .* (v * v')
    end
    resid = norm(T - I)
    resid < 1e-6 || error("B1 completeness not satisfied (‖Σ w b⊗b − I‖ = $resid). " *
                          "Neighbour shells may be insufficient.")

    wb = Matrix{Float64}(undef, nntot, nk)
    for k in 1:nk, b in 1:nntot
        r = norm(SVector{3,Float64}(bvec[1, b, k], bvec[2, b, k], bvec[3, b, k]))
        wb[b, k] = shell_weight[shell_of(r)]
    end

    return BVectors(nntot, kpb, gpb, bvec, wb, radii, shell_weight)
end
