# Symmetry-adapted Wannier functions (SAWF / site_symmetry; Sakuma PRB 87, 235109 (2013)).
# Reads the `.dmn` file written by pw2wannier90 — the k-point↔irreducible map and the band- and
# Wannier-representation matrices — and symmetrises the gauge so the Wannier functions transform
# according to the chosen site-symmetry representation. The symmetrisation projects U(k) onto the
# symmetric subspace at the irreducible k and reconstructs the star; the gradient is projected the
# same way. Exact conventions in docs/reference-notes/sawf-sitesym.md. Validation is
# gauge-invariant only (Ω, centres, bands) — the symmetric U differs from the reference's.

using LinearAlgebra
using StaticArrays

"""
    Sitesym

Site-symmetry data from a `.dmn`: `ik2ir[ik]` (irreducible index of k), `ir2ik[ir]`
(representative k of irreducible index ir), `kptsym[isym, ir]` (the k reached by applying
`isym` to irreducible point `ir`), and the band/Wannier representation matrices
`d_band[:,:,isym,ir]` (num_bands²) and `d_wann[:,:,isym,ir]` (num_wann²).
"""
struct Sitesym
    nsym::Int
    nkptirr::Int
    ik2ir::Vector{Int}
    ir2ik::Vector{Int}
    kptsym::Matrix{Int}                  # (nsym, nkptirr)
    d_wann::Array{ComplexF64,4}          # (nw, nw, nsym, nkptirr)
    d_band::Array{ComplexF64,4}          # (nb, nb, nsym, nkptirr)
    eps::Float64                         # symmetrize_eps convergence threshold
end

Sitesym(nsym, nkptirr, ik2ir, ir2ik, kptsym, d_wann, d_band) =
    Sitesym(nsym, nkptirr, ik2ir, ir2ik, kptsym, d_wann, d_band, 1e-3)

"After disentanglement the gauge is square (num_wann): the band representation is replaced by
the Wannier one for the localisation phase (the reference's `sitesym_replace_d_matrix_band`)."
replace_d_matrix_band(s::Sitesym) =
    Sitesym(s.nsym, s.nkptirr, s.ik2ir, s.ir2ik, s.kptsym, s.d_wann, s.d_wann, s.eps)

Base.show(io::IO, ::MIME"text/plain", s::Sitesym) =
    print(io, "Sitesym: ", s.nsym, " symmetries, ", s.nkptirr, " irreducible k")

"""
    read_dmn(path, num_bands, num_wann) -> Sitesym

Read a `seedname.dmn` (formatted, list-directed). Line 2 is `num_bands nsym nkptirr num_kpts`;
then `ik2ir` (num_kpts), `ir2ik` (nkptirr), `kptsym` (nsym×nkptirr, isym fastest), `d_wann`
(nw²×nsym×nkptirr, complex, i fastest), `d_band` (nb²×nsym×nkptirr, complex).
"""
function read_dmn(path::AbstractString, num_bands::Int, num_wann::Int; eps::Float64=1e-3)
    lines = readlines(path)
    # Line 1 is the header comment; everything after is list-directed data. Strip () and commas.
    raw = replace(join(lines[2:end], " "), '(' => ' ', ')' => ' ', ',' => ' ')
    f = split(raw)
    p = 1
    nextint() = (v = parse(Int, f[p]); p += 1; v)
    nb = nextint(); nsym = nextint(); nkptirr = nextint(); nk = nextint()
    nb == num_bands || error("read_dmn: num_bands $nb ≠ $num_bands")
    ik2ir = [nextint() for _ in 1:nk]
    ir2ik = [nextint() for _ in 1:nkptirr]
    kptsym = Matrix{Int}(undef, nsym, nkptirr)
    for ir in 1:nkptirr, is in 1:nsym          # isym fastest
        kptsym[is, ir] = nextint()
    end
    nextcplx() = (re = parse(Float64, f[p]); im = parse(Float64, f[p+1]); p += 2; complex(re, im))
    d_wann = Array{ComplexF64,4}(undef, num_wann, num_wann, nsym, nkptirr)
    for ir in 1:nkptirr, is in 1:nsym, j in 1:num_wann, i in 1:num_wann   # i fastest
        d_wann[i, j, is, ir] = nextcplx()
    end
    d_band = Array{ComplexF64,4}(undef, num_bands, num_bands, nsym, nkptirr)
    for ir in 1:nkptirr, is in 1:nsym, j in 1:num_bands, i in 1:num_bands
        d_band[i, j, is, ir] = nextcplx()
    end
    return Sitesym(nsym, nkptirr, ik2ir, ir2ik, kptsym, d_wann, d_band, eps)
end

"Löwdin orthonormalisation of a square matrix: U = W·V† from the SVD A = W·Σ·V†."
function _ortho_lowdin(A::AbstractMatrix{ComplexF64})
    F = svd(A)
    return F.U * F.Vt
end

"Little-group symmetry indices of irreducible point `ir` (those fixing its representative k)."
_little_group(s::Sitesym, ir::Int) = findall(is -> s.kptsym[is, ir] == s.ir2ik[ir], 1:s.nsym)

"""
    symmetrize_ukirr(U0, s, ir, dband, dwann, n) -> U

Project the (leading n×n block of the) gauge `U0` at irreducible point `ir` onto the symmetric
subspace: iterate `ũ = (1/|G_k|) Σ_{R'∈G_k} d(R')† · U · D(R')` then Löwdin-orthonormalise,
until convergence. `dband`/`dwann` select the representation (num_bands / num_wann).
"""
function symmetrize_ukirr(U0::AbstractMatrix{ComplexF64}, s::Sitesym, ir::Int,
                          dband::Array{ComplexF64,4}, dwann::Array{ComplexF64,4}, n::Int;
                          eps::Float64=s.eps, maxiter::Int=100)
    lg = _little_group(s, ir)
    nw = size(U0, 2)
    U = copy(U0[1:n, :])
    ngk = length(lg)
    ngk <= 1 && return _ortho_lowdin(U)
    conv = false
    for _ in 1:maxiter
        usum = zeros(ComplexF64, n, nw)
        cmat2 = ComplexF64.(ngk * Matrix{ComplexF64}(I, nw, nw))
        for is in lg
            cmat = (@view dband[1:n, 1:n, is, ir])' * U * (@view dwann[:, :, is, ir])  # d†·U·D
            usum .+= cmat
            cmat2 .-= U' * cmat                       # ngk·I − Σ U†·(d†UD), on the CURRENT U
        end
        if sum(abs.(cmat2)) < eps                     # converged (checked before the update)
            conv = true
            break
        end
        usum ./= ngk
        U = _ortho_lowdin(usum)
    end
    return U
end

"""
    symmetrize_u!(U, s, dband, dwann; n=size(U,1))

Symmetrise a full gauge stack in place: at each irreducible representative project onto the
symmetric subspace, then reconstruct the star via `U(Rk) = d(R)·U(k)·D(R)†`. `U` is
`ndim × num_wann × num_kpts` (ndim = n).
"""
function symmetrize_u!(U::Array{ComplexF64,3}, s::Sitesym,
                       dband::Array{ComplexF64,4}, dwann::Array{ComplexF64,4}; n::Int=size(U, 1))
    nw = size(U, 2)
    done = falses(size(U, 3))
    for ir in 1:s.nkptirr
        ik = s.ir2ik[ir]
        Uk = symmetrize_ukirr((@view U[:, :, ik]), s, ir, dband, dwann, n)
        U[1:n, :, ik] = Uk
        done[ik] = true
        for is in 1:s.nsym
            irk = s.kptsym[is, ir]
            done[irk] && continue
            U[1:n, :, irk] = (@view dband[1:n, 1:n, is, ir]) * Uk * (@view dwann[:, :, is, ir])'
            done[irk] = true
        end
    end
    return U
end

"""
    symmetrize_rotation!(d, s)

Propagate the CG rotation generator from each irreducible representative to its star,
`d(Rk) = D(R)·d(k)·D(R)†`, so applying `exp(d)` keeps U symmetry-adapted (the representative-k
gradient is the only nonzero one after mode-1 gradient symmetrisation).
"""
function symmetrize_rotation!(d::Array{ComplexF64,3}, s::Sitesym)
    D = s.d_wann
    for ir in 1:s.nkptirr
        ik = s.ir2ik[ir]
        done = Set{Int}([ik])
        for is in 1:s.nsym
            irk = s.kptsym[is, ir]
            irk in done && continue
            push!(done, irk)
            d[:, :, irk] = (@view D[:, :, is, ir]) * (@view d[:, :, ik]) * (@view D[:, :, is, ir])'
        end
    end
    return d
end

"""
    symmetrize_zmatrix!(Z, s)

Symmetrise the disentanglement Z matrices across each k-star and little group,
`Z(k) ← (1/|G_k|) Σ_{R'∈G_k∪1} d†(R') [Σ_star d†(R) Z(Rk) d(R)] d(R')`, updating only the
irreducible representatives (non-representative entries are never used afterwards). Each
distinct star member contributes once, through the first symmetry that reaches it.
"""
function symmetrize_zmatrix!(Z::Vector{Matrix{ComplexF64}}, s::Sitesym)
    lfound = falses(length(Z))
    for ir in 1:s.nkptirr
        ik = s.ir2ik[ir]
        nd = size(Z[ik], 1)
        lfound[ik] = true
        acc = copy(Z[ik])
        for is in 2:s.nsym
            irk = s.kptsym[is, ir]
            lfound[irk] && continue
            lfound[irk] = true
            d = @view s.d_band[1:nd, 1:nd, is, ir]
            acc .+= d' * Z[irk] * d
        end
        tmp = copy(acc)
        ngk = 1
        for is in 2:s.nsym
            s.kptsym[is, ir] == ik || continue
            ngk += 1
            d = @view s.d_band[1:nd, 1:nd, is, ir]
            acc .+= d' * tmp * d
        end
        Z[ik] = acc ./ ngk
    end
    return Z
end

"""
    dis_extract_symmetry!(U, Z, s, ik, n) -> λ

Constrained Ω_I step at irreducible representative `ik` (the reference's
`sitesym_dis_extract_symmetry`): steepest-descent on the subspace embedding `U` (n × num_wann)
along `ΔU = Z·U − U·λ` with `λ = U†ZU`, maximising the Rayleigh quotient band-by-band in the
2-dimensional span {u_i, Δu_i} (generalized 2×2 eigenproblem, larger eigenvalue), then
re-projecting onto the symmetric manifold each sweep. Returns the final `λ` — its real trace
is the k-point's Z-eigenvalue-sum surrogate for the Ω_I bookkeeping.
"""
function dis_extract_symmetry!(U::AbstractMatrix{ComplexF64}, Z::AbstractMatrix{ComplexF64},
                               s::Sitesym, ik::Int, n::Int)
    nw = size(U, 2)
    λ = zeros(ComplexF64, nw, nw)
    ir = s.ik2ir[ik]
    Unew = similar(U, n, nw)
    for _ in 1:50
        ZU = Z * U                               # n × nw
        λ = U' * ZU
        ΔU = ZU - U * λ
        sum(abs, ΔU) < 1e-10 && return λ
        for i in 1:nw
            u = @view U[:, i]
            du = @view ΔU[:, i]
            zu = @view ZU[:, i]
            s22 = real(dot(du, du))
            if abs(s22) < 1e-10
                Unew[:, i] = u
                continue
            end
            h12 = dot(zu, du)                    # ⟨Zu|Δu⟩ (Z hermitian)
            s12 = dot(u, du)
            H2 = ComplexF64[real(dot(u, zu)) h12; conj(h12) real(dot(du, Z * du))]
            S2 = ComplexF64[real(dot(u, u)) s12; conj(s12) s22]
            F = eigen(Hermitian(H2), Hermitian(S2))          # ascending; v†S v = 1
            v = @view F.vectors[:, 2]                        # larger eigenvalue
            @. Unew[:, i] = v[1] * u + v[2] * du
        end
        U .= symmetrize_ukirr(Unew, s, ir, s.d_band, s.d_wann, n)
    end
    return λ
end

"""
    symmetrize_gradient!(G, s)

Project the num_wann gradient field onto the symmetric-gauge tangent space (localisation
phase, d = d_wann): mode 1 accumulates the star into each representative and zeroes the rest;
mode 2 averages over the little group. Applied mode 1 then mode 2, matching the reference.
"""
function symmetrize_gradient!(G::Array{ComplexF64,3}, s::Sitesym)
    D = s.d_wann
    # mode 1: G[ik] += Σ_{isym} D(isym,ir)† G[irk] D(isym,ir); then zero non-representatives
    for ir in 1:s.nkptirr
        ik = s.ir2ik[ir]
        done = Set{Int}([ik])
        for is in 2:s.nsym
            irk = s.kptsym[is, ir]
            irk in done && continue
            push!(done, irk)
            G[:, :, ik] .+= (@view D[:, :, is, ir])' * (@view G[:, :, irk]) * (@view D[:, :, is, ir])
        end
    end
    for k in 1:size(G, 3)
        s.ir2ik[s.ik2ir[k]] != k && (G[:, :, k] .= 0)
    end
    # mode 2: little-group average at each representative
    for ir in 1:s.nkptirr
        ik = s.ir2ik[ir]
        lg = _little_group(s, ir)
        length(lg) <= 1 && continue
        acc = copy(@view G[:, :, ik])                 # the isym = 1 (identity) term
        for is in lg
            is == 1 && continue
            acc .+= (@view D[:, :, is, ir])' * (@view G[:, :, ik]) * (@view D[:, :, is, ir])
        end
        G[:, :, ik] = acc ./ length(lg)
    end
    return G
end
