# The operator-centric core of Wannier interpolation.
#
# Every Wannier-interpolable quantity — the Hamiltonian, the position operator, and (later) spin
# or Berry-connection blocks — is the same object: a matrix-valued lattice operator O(R) on the
# Wigner–Seitz R-vector set, Fourier-interpolable to any k. `TBOperator` captures that once;
# the Hamiltonian and r(R) are just instances.
#
#   O(R)  = (1/N_k) Σ_k e^{-i 2π k·R} O(k)         (fourier_to_R)
#   O(k') = Σ_R e^{+i 2π k'·R} O(R) / ndegen(R)     (evaluation: `op(k')`)

using LinearAlgebra
using StaticArrays

"""
    TBOperator

A tight-binding (Wannier-gauge) operator: `data[m, n, ir, c]` is the (m,n) matrix element of
component `c` at lattice vector `irvec[ir]`, stored **undivided** by `ndegen` (the Wannier90
convention). Component count is 1 for the Hamiltonian, 3 (x,y,z) for the position operator.

Evaluate at a fractional k-point by calling it: `op(k)` returns a `num_wann × num_wann` matrix
(1 component) or a 3-vector of matrices.
"""
struct TBOperator
    name::Symbol
    lattice::Lattice
    irvec::Vector{NTuple{3,Int}}
    ndegen::Vector{Int}
    data::Array{ComplexF64,4}          # (nw, nw, nrpts, ncomp)
end

num_wann(op::TBOperator) = size(op.data, 1)
ncomponents(op::TBOperator) = size(op.data, 4)

function Base.show(io::IO, ::MIME"text/plain", op::TBOperator)
    print(io, "TBOperator :", op.name, "  (", num_wann(op), "×", num_wann(op), ", ",
          length(op.irvec), " R-vectors, ", ncomponents(op), " component",
          ncomponents(op) == 1 ? "" : "s", ")")
end

"""
    fourier_to_R(Ok, kgrid, irvec) -> OR

Generic k→R transform of a matrix-valued operator: `OR[:,:,ir,c] = (1/N_k) Σ_k e^{-i2πk·R} Ok[:,:,k,c]`.
`Ok` is (nw × nw × nkpt × ncomp).
"""
function fourier_to_R(Ok::Array{ComplexF64,4}, kgrid::KGrid, irvec::Vector{NTuple{3,Int}})
    nw, _, nk, nc = size(Ok)
    OR = zeros(ComplexF64, nw, nw, length(irvec), nc)
    @maybe_threads (length(irvec) >= THREAD_MIN) for ir in 1:length(irvec)
        R = SVector{3,Float64}(irvec[ir]...)
        for k in 1:nk
            fac = cis(-TWOPI * dot(kgrid.frac[k], R)) / nk
            @views for c in 1:nc
                OR[:, :, ir, c] .+= fac .* Ok[:, :, k, c]
            end
        end
    end
    return OR
end

"Evaluate the operator at fractional k: Σ_R e^{+i2πk·R} O(R)/ndegen(R)."
function (op::TBOperator)(kfrac::AbstractVector)
    kf = SVector{3,Float64}(kfrac...)
    nw = num_wann(op); nc = ncomponents(op)
    out = zeros(ComplexF64, nw, nw, nc)
    for ir in 1:length(op.irvec)
        R = SVector{3,Float64}(op.irvec[ir]...)
        fac = cis(TWOPI * dot(kf, R)) / op.ndegen[ir]
        @views for c in 1:nc
            out[:, :, c] .+= fac .* op.data[:, :, ir, c]
        end
    end
    return nc == 1 ? out[:, :, 1] : [out[:, :, c] for c in 1:nc]
end

"""
    hamiltonian_operator(model, result) -> TBOperator

The interpolable Hamiltonian H(R) from a wannierisation result (isolated or disentangled).
"""
function hamiltonian_operator(model::Model, result::WannierResult)
    result.eig_interp !== nothing ||
        error("no band energies (.eig) available to build the Hamiltonian operator")
    irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
    Hr, _ = build_hr(result.U, result.eig_interp, model.kgrid, irvec)
    data = reshape(Hr, size(Hr)..., 1)
    return TBOperator(:hamiltonian, model.lattice, irvec, ndegen, data)
end

"""
    position_operator(model, result) -> TBOperator

The position operator ⟨0m|r|Rn⟩ (3 Cartesian components, Å) from the final-gauge overlap
matrices, using the finite-difference Berry-connection formula the reference uses in `_tb.dat`:
off-diagonal `i w_b b M̃_mn` (Wang–Yates–Souza–Vanderbilt PRB 74, 195118 (2006), Eq. 44) and
diagonal `−w_b b Im ln M̃_nn` (Marzari–Vanderbilt PRB 56, 12847 (1997), Eq. 32), then Fourier
k→R.
"""
function position_operator(model::Model, result::WannierResult)
    bv = model.bvectors
    Mrot = result.Mrot
    nw = size(Mrot, 1)
    nk = nkpt(model.kgrid)
    # Berry-connection matrix A(k) per component, then generic k→R.
    Ak = zeros(ComplexF64, nw, nw, nk, 3)
    for k in 1:nk, nn in 1:bv.nntot
        w = bv.wb[nn, k]
        b = SVector{3,Float64}(bv.bvec[1, nn, k], bv.bvec[2, nn, k], bv.bvec[3, nn, k])
        @inbounds for i in 1:nw, j in 1:nw
            if i == j
                v = -w * imag(log(Mrot[i, i, nn, k]))
                for c in 1:3
                    Ak[i, i, k, c] += v * b[c]
                end
            else
                v = im * w * Mrot[j, i, nn, k]      # element (j,i), matching the reference loop
                for c in 1:3
                    Ak[j, i, k, c] += v * b[c]
                end
            end
        end
    end
    irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
    data = fourier_to_R(Ak, model.kgrid, irvec)
    return TBOperator(:position, model.lattice, irvec, ndegen, data)
end

"""
    bands(H::TBOperator, kpts) -> energies

Interpolated eigenvalues (num_wann × npts, ascending per column) of a 1-component operator.
"""
function bands(H::TBOperator, kpts::AbstractVector)
    ncomponents(H) == 1 || error("bands: need a scalar (1-component) operator, got $(H.name)")
    nw = num_wann(H)
    E = Matrix{Float64}(undef, nw, length(kpts))
    @maybe_threads (length(kpts) >= 32) for ik in 1:length(kpts)
        Hk = H(kpts[ik])
        E[:, ik] = eigvals(Hermitian((Hk + Hk') / 2))
    end
    return E
end
