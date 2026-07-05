# Finite-difference b-vectors and their B1 completeness weights.
#
# The neighbour connectivity (which k-point is k+b, and the reciprocal fold G) is taken from the
# .mmn file so that the b-vector ordering aligns exactly with the stored overlap matrices. The
# *weights* are solved here from the B1 completeness relation Σ_b w_b b_α b_β = δ_αβ
# (Marzari–Vanderbilt, PRB 56, 12847 (1997), Eq. B1), which is the gauge-independent "science"
# of the k-mesh and must reproduce the reference shell weights.
#
# Higher-order finite differences (`higher_order_n = N`, Lihm's scheme in the reference's
# kmesh.F90): the b-vector list consists of N consecutive blocks, block m holding m·b for every
# first-order b (this is the order pw2wannier90 writes the .mmn in). B1 alone does not fix the
# weights of parallel shells (w_{2b} = 0 also satisfies it), so the first-order weights are
# solved from B1 on block 1 and the multiples take the 1D central-difference factors
# w_{mb} = w_b · (1/m²) Π_{j≠m} j²/(j²−m²), which cancel the O(b²)…O(b^{2N-2}) error terms and
# preserve the B1 sum (Σ_m m²·fact_m = 1).

using LinearAlgebra
using StaticArrays

"Central-difference weight factor for the m-th multiple in an order-N scheme."
function _higher_order_factor(m::Int, N::Int)
    fact = 1.0 / m^2
    for j in 1:N
        j == m && continue
        fact *= j^2 / (j^2 - m^2)
    end
    return fact
end

"""
Detect the higher-order structure of the k=1 b-vector list, independent of ordering: each b is
assigned the integer multiple `m` of the shortest vector parallel to it in the list. Returns
`(N, mult)` where `mult[b]` is that multiple and `N = maximum(mult)`; a standard mesh gives
`N = 1`. Errors if the multiples are not a complete 1..N per direction (not a Lihm mesh).
"""
function _detect_higher_order(bvec::Array{Float64,3}, nntot::Int, tol::Float64)
    b1(b) = SVector{3,Float64}(bvec[1, b, 1], bvec[2, b, 1], bvec[3, b, 1])
    # shortest same-direction partner of each b (b itself if none shorter)
    mult = ones(Int, nntot)
    prim = collect(1:nntot)
    for b in 1:nntot
        v = b1(b)
        for c in 1:nntot
            u = b1(c)
            norm(u) < norm(v) - tol || continue
            norm(cross(u, v)) < tol * norm(v) && dot(u, v) > 0 || continue
            if prim[b] == b || norm(u) < norm(b1(prim[b]))
                prim[b] = c
            end
        end
        if prim[b] != b
            r = norm(v) / norm(b1(prim[b]))
            m = round(Int, r)
            abs(r - m) < tol && norm(v - m * b1(prim[b])) < tol * m ||
                error("b-vector $b is parallel to a shorter one but not an integer multiple " *
                      "(ratio $r) — unsupported k-mesh")
            mult[b] = m
        end
    end
    N = maximum(mult)
    if N > 1
        # every direction family must carry the complete multiples 1..N
        for b in 1:nntot
            fam = sort([mult[c] for c in 1:nntot if prim[c] == prim[b] || c == prim[b]])
            fam == collect(1:N) ||
                error("higher-order k-mesh: direction family of b-vector $b has multiples " *
                      "$fam, expected 1:$N")
        end
    end
    return N, mult
end

"""
    build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol) -> BVectors

Compute Cartesian b-vectors for the given neighbour connectivity and solve the B1 relation for
the per-shell finite-difference weights. Verifies completeness and errors if it is not met.
Higher-order (`higher_order_n`) meshes are detected from the block structure of the b-list.
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
    b1(b) = SVector{3,Float64}(bvec[1, b, 1], bvec[2, b, 1], bvec[3, b, 1])

    # Higher-order structure: only the first-order (multiple = 1) weights are free parameters.
    horder, mult = _detect_higher_order(bvec, nntot, kmesh_tol)

    # Distinct shell radii of the FIRST-ORDER b's (identical across k for a uniform mesh).
    radii = Float64[]
    for b in 1:nntot
        mult[b] == 1 || continue
        r = norm(b1(b))
        any(x -> abs(x - r) < kmesh_tol, radii) || push!(radii, r)
    end
    sort!(radii)
    shell_of(r) = findfirst(x -> abs(x - r) < kmesh_tol, radii)::Int
    nsh1 = length(radii)

    # Least-squares solve of the six independent components of the B1 tensor relation on the
    # first-order b's; the multiples carry the fixed central-difference factors, so their
    # contribution to B1 is folded in through Σ_m m²·fact_m = 1.
    Amat = zeros(6, nsh1)
    for b in 1:nntot
        mult[b] == 1 || continue
        v = b1(b)
        s = shell_of(norm(v))
        Amat[1, s] += v[1]^2; Amat[2, s] += v[2]^2; Amat[3, s] += v[3]^2
        Amat[4, s] += v[1]*v[2]; Amat[5, s] += v[1]*v[3]; Amat[6, s] += v[2]*v[3]
    end
    target = SVector{6,Float64}(1, 1, 1, 0, 0, 0)
    w1 = Amat \ Vector(target)                       # SVD/QR least squares

    # Extended per-shell weights and per-position shell index [1× shells, 2× shells, …].
    shell_weight = vcat((w1 .* _higher_order_factor(m, horder) for m in 1:horder)...)
    shellidx = Vector{Int}(undef, nntot)
    for b in 1:nntot
        m = mult[b]
        shellidx[b] = (m - 1) * nsh1 + shell_of(norm(b1(b)) / m)
    end
    radii_ext = vcat((radii .* m for m in 1:horder)...)

    # Verify B1 completeness over the full (possibly extended) set at k=1.
    T = zeros(3, 3)
    for b in 1:nntot
        v = b1(b)
        T .+= shell_weight[shellidx[b]] .* (v * v')
    end
    resid = norm(T - I)
    resid < 1e-6 || error("B1 completeness not satisfied (‖Σ w b⊗b − I‖ = $resid). " *
                          "Neighbour shells may be insufficient.")

    wb = Matrix{Float64}(undef, nntot, nk)
    for k in 1:nk, b in 1:nntot
        wb[b, k] = shell_weight[shellidx[b]]
    end

    return BVectors(nntot, kpb, gpb, bvec, wb, radii_ext, shell_weight)
end
