# Crystal symmetry operations and Brillouin-zone reduction. Reads the `.sym` file written by
# pw2wannier90 (the space-group operations {R, t}) and reduces a uniform k-mesh to its
# irreducible wedge — used for symmetry-adapted Wannier functions and for symmetrised
# Brillouin-zone integration (e.g. anomalous_hall_sym).

using LinearAlgebra
using StaticArrays

"""
    SymmetryOps

Space-group operations from a `.sym` file: `rot[s]` is the 3×3 rotation of symmetry `s`
(fractional/crystal basis, acting on fractional coordinates) and `trans[s]` its fractional
translation. `nsym` operations in total.
"""
struct SymmetryOps
    rot::Vector{SMatrix{3,3,Float64,9}}
    trans::Vector{SVector{3,Float64}}
end
nsym(s::SymmetryOps) = length(s.rot)

Base.show(io::IO, ::MIME"text/plain", s::SymmetryOps) =
    print(io, "SymmetryOps: ", nsym(s), " operations")

"""
    read_sym(path) -> SymmetryOps

Read a `seedname.sym`: line 1 = number of operations; then, per operation, three rows of the
3×3 rotation matrix followed by one row of the 3-vector fractional translation (blank lines
between blocks are ignored).
"""
function read_sym(path::AbstractString)
    toks = String[]
    for ln in eachline(path)
        s = strip(ln)
        isempty(s) && continue
        append!(toks, split(s))
    end
    n = parse(Int, toks[1])
    rot = Vector{SMatrix{3,3,Float64,9}}(undef, n)
    trans = Vector{SVector{3,Float64}}(undef, n)
    idx = 2
    for s in 1:n
        M = Matrix{Float64}(undef, 3, 3)
        for i in 1:3, j in 1:3           # row i, column j — file rows are matrix rows
            M[i, j] = parse(Float64, toks[idx]); idx += 1
        end
        rot[s] = SMatrix{3,3,Float64}(M)
        trans[s] = SVector{3,Float64}(parse(Float64, toks[idx]), parse(Float64, toks[idx+1]),
                                      parse(Float64, toks[idx+2]))
        idx += 3
    end
    return SymmetryOps(rot, trans)
end

"""
    cubic_point_group() -> SymmetryOps

The 48 operations of the cubic point group Oₕ as Cartesian 3×3 matrices (all signed 3×3
permutation matrices), with zero translation. Useful as a symmetry source for cubic crystals
when no `.sym` file is available; the magnetic subgroup of a given quantity can be filtered
from it numerically.
"""
function cubic_point_group()
    perms = ([1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1])
    rot = SMatrix{3,3,Float64,9}[]
    for p in perms, s1 in (1, -1), s2 in (1, -1), s3 in (1, -1)
        M = zeros(3, 3)
        M[1, p[1]] = s1; M[2, p[2]] = s2; M[3, p[3]] = s3
        push!(rot, SMatrix{3,3,Float64}(M))
    end
    return SymmetryOps(rot, [SVector(0.0, 0.0, 0.0) for _ in rot])
end

"Subset of `sym` whose Cartesian rotations `S` satisfy Ω(Sk)=det(S)·S·Ω(k) for the axial field."
function _pseudovector_subgroup(bm::BerryModel, sym::SymmetryOps, kmesh::NTuple{3,Int},
                                fermi_energy::Float64; tol::Float64=1e-4)
    n1, n2, n3 = kmesh
    ks = [SVector(i / n1, j / n2, k / n3) for i in 0:n1-1 for j in 0:n2-1 for k in 0:n3-1]
    Ω = Dict{NTuple{3,Int},SVector{3,Float64}}()
    key(k) = (round(Int, mod(k[1], 1) * n1) % n1, round(Int, mod(k[2], 1) * n2) % n2,
              round(Int, mod(k[3], 1) * n3) % n3)
    for k in ks
        Ω[key(k)] = _imf_kdata(_berry_kdata(bm, k), fermi_energy)
    end
    keep = Int[]
    B = bm.lattice.B
    for s in 1:nsym(sym)
        R = sym.rot[s]
        d = det(R)
        ok = true
        for k in ks
            Sk_frac = B \ (R * (B * k))
            got = Ω[key(Sk_frac)]
            want = d .* (R * Ω[key(k)])
            if norm(got - want) > tol * (1 + norm(want))
                ok = false; break
            end
        end
        ok && push!(keep, s)
    end
    return SymmetryOps(sym.rot[keep], sym.trans[keep])
end

"""
    anomalous_hall_sym(bm, sym; fermi_energy, kmesh) -> (ahc, ninfo)

Anomalous Hall conductivity (S/cm, axial x/y/z) computed on the irreducible wedge of `kmesh`
under `sym` and symmetrised — the Berry curvature is a pseudovector, so each irreducible point
contributes `w · (1/|G|) Σ_s det(R_s)·R_s·Ω(k)`. Equal to the full-BZ [`anomalous_hall`] but
evaluates only the irreducible k-points. `ninfo = (n_irreducible, n_full)`.
"""
function anomalous_hall_sym(bm::BerryModel, sym::SymmetryOps; fermi_energy::Float64,
                            kmesh::NTuple{3,Int}=(25, 25, 25))
    reps, wts, _ = irreducible_kmesh(kmesh, sym; kaction=:cart, lattice=bm.lattice)
    ng = nsym(sym)
    nktot = prod(kmesh)
    acc = zeros(SVector{3,Float64})
    per = Vector{SVector{3,Float64}}(undef, length(reps))
    Threads.@threads for i in 1:length(reps)
        Ω = _imf_kdata(_berry_kdata(bm, reps[i]), fermi_energy)
        s = zeros(SVector{3,Float64})
        for r in 1:ng
            R = sym.rot[r]
            s += det(R) .* (R * Ω)
        end
        per[i] = (wts[i] / ng) .* s
    end
    acc = reduce(+, per)
    fac = -1.0e8 * ELEM_CHARGE_SI^2 / (HBAR_SI * cell_volume(bm.lattice))
    return (fac / nktot) .* acc, (length(reps), nktot)
end

"Fold a fractional k-point into [0,1)³ (numerical tolerance `tol`)."
_wrapk(k::SVector{3,Float64}; tol::Float64=1e-7) =
    SVector(mod(k[1] + tol, 1.0) - tol, mod(k[2] + tol, 1.0) - tol, mod(k[3] + tol, 1.0) - tol)

"""
    irreducible_kmesh(kmesh, sym; kaction=:frac) -> (kpts, weights, reps)

Reduce the uniform Γ-centred mesh `k = (i/n1, j/n2, k/n3)` to its irreducible wedge under the
point-group parts of `sym`. Returns the irreducible fractional k-points, their integer
multiplicities (orbit sizes, summing to `prod(kmesh)`), and `reps` mapping each irreducible
index to the list of symmetry indices generating distinct star members (for tensor
symmetrisation). `kaction`: how a rotation acts on a fractional k — `:frac` uses `R·k`
directly (the pw2wannier90 convention), `:fracT` uses `Rᵀ·k`.
"""
# The fractional k-action of symmetry `s`. `.sym` rotations are Cartesian, so the action on a
# fractional k is B⁻¹·R·B (`:cart`); `:frac`/`:fracT` treat R as already crystal-basis.
function _kact(sym::SymmetryOps, s::Int, k::SVector{3,Float64}, kaction::Symbol, Bmat)
    R = sym.rot[s]
    kaction === :cart && return Bmat \ (R * (Bmat * k))
    kaction === :fracT && return R' * k
    return R * k
end

function irreducible_kmesh(kmesh::NTuple{3,Int}, sym::SymmetryOps;
                           kaction::Symbol=:cart, lattice=nothing)
    Bmat = kaction === :cart ?
           (lattice === nothing ? error("irreducible_kmesh(:cart) needs `lattice`") : lattice.B) :
           nothing
    n1, n2, n3 = kmesh
    all_k = [SVector(i / n1, j / n2, k / n3) for i in 0:n1-1 for j in 0:n2-1 for k in 0:n3-1]
    key(k) = (round(Int, k[1] * n1), round(Int, k[2] * n2), round(Int, k[3] * n3))
    seen = Set{NTuple{3,Int}}()
    reps = Vector{SVector{3,Float64}}()
    weights = Int[]
    for k in all_k
        kk = key(_wrapk(k))
        kk in seen && continue
        star = Set{NTuple{3,Int}}()          # distinct images (mod 1)
        for s in 1:nsym(sym)
            push!(star, key(_wrapk(_kact(sym, s, k, kaction, Bmat))))
        end
        for m in star
            push!(seen, m)
        end
        push!(reps, k)
        push!(weights, length(star))
    end
    return reps, weights, sym
end
