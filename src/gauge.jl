# Gauge construction and rotation.
#
# The initial Wannier gauge comes from the trial projections A via Löwdin orthonormalisation
# U_k = A_k (A_k† A_k)^{-1/2}, computed stably as the unitary/isometric polar factor of A_k from
# its SVD (A = VΣW† ⇒ U = VW†). Overlaps are then rotated into the current gauge,
# M̃_k,b = U_k† M_k,b U_{k+b}.

using LinearAlgebra

"""
    lowdin(A) -> U

Löwdin-orthonormalised gauge from an (num_bands × num_wann) projection block: the closest
isometry to `A` in Frobenius norm, `U = V W†` where `A = V Σ W†`.
"""
function lowdin(A::AbstractMatrix{ComplexF64})
    F = svd(A)
    return F.U * F.Vt
end

"""
    initial_gauge(A) -> U

Per-k Löwdin gauge for the whole model. `A` is (num_bands × num_wann × num_kpts); returns
`U` of the same shape.
"""
function initial_gauge(A::Array{ComplexF64,3})
    nb, nw, nk = size(A)
    U = Array{ComplexF64,3}(undef, nb, nw, nk)
    for k in 1:nk
        U[:, :, k] = lowdin(@view A[:, :, k])
    end
    return U
end

"""
    rotate_overlaps(M, U, kpb) -> Mrot

Rotate the band-space overlaps into the Wannier gauge:
`Mrot[:,:,b,k] = U_k† · M[:,:,b,k] · U_{k+b}`, shape (num_wann × num_wann × nntot × num_kpts).
"""
function rotate_overlaps(M::Array{ComplexF64,4}, U::Array{ComplexF64,3}, kpb::Matrix{Int})
    nb, nw, nk = size(U)
    nntot = size(M, 3)
    Mrot = Array{ComplexF64,4}(undef, nw, nw, nntot, nk)
    for k in 1:nk, b in 1:nntot
        kb = kpb[b, k]
        Mrot[:, :, b, k] = (@view U[:, :, k])' * (@view M[:, :, b, k]) * (@view U[:, :, kb])
    end
    return Mrot
end
