# Wannier interpolation: real-space Hamiltonian H(R) from the wannierised gauge, the
# Wigner–Seitz R-vector set with degeneracy weights, and band interpolation on a k-path.
# Conventions follow the reference (hamiltonian.F90); see docs/reference-notes/interpolation.md.
#
#   H(k) = U†(k) diag(ε_k) U(k)                    (Wannier-gauge Hamiltonian)
#   H(R) = (1/N_k) Σ_k e^{-i 2π k·R} H(k)          (stored undivided by ndegen)
#   H(k') = Σ_R e^{+i 2π k'·R} H(R) / ndegen(R)    (then diagonalise, ascending)

using LinearAlgebra
using StaticArrays

"""
    wigner_seitz(lattice, mp_grid; ws_search_size=2, tol=1e-5) -> (irvec, ndegen)

Wigner–Seitz lattice vectors (columns of `irvec`, integer) inside the superlattice cell defined
by `mp_grid`, with degeneracy `ndegen` for boundary points. Enforces the sum rule
Σ_R 1/ndegen(R) = ∏ mp_grid.
"""
function wigner_seitz(lattice::Lattice, mp_grid::NTuple{3,Int};
                     ws_search_size::Int=2, tol::Float64=1.0e-5)
    metric = Matrix(transpose(lattice.A) * lattice.A)   # a_i·a_j
    t2 = tol^2
    dist(n) = (v = SVector{3,Float64}(n[1], n[2], n[3]); dot(v, metric * v))

    irvec = NTuple{3,Int}[]
    ndegen = Int[]
    rng = -(ws_search_size + 1):(ws_search_size + 1)
    for n1 in -ws_search_size*mp_grid[1]:ws_search_size*mp_grid[1],
        n2 in -ws_search_size*mp_grid[2]:ws_search_size*mp_grid[2],
        n3 in -ws_search_size*mp_grid[3]:ws_search_size*mp_grid[3]

        dmin = Inf
        d0 = 0.0
        deg = 0
        dists = Float64[]
        for i1 in rng, i2 in rng, i3 in rng
            d = dist((n1 - i1*mp_grid[1], n2 - i2*mp_grid[2], n3 - i3*mp_grid[3]))
            push!(dists, d)
            (i1 == 0 && i2 == 0 && i3 == 0) && (d0 = d)
            d < dmin && (dmin = d)
        end
        if abs(d0 - dmin) < t2
            deg = count(d -> abs(d - dmin) < t2, dists)
            push!(irvec, (n1, n2, n3))
            push!(ndegen, deg)
        end
    end

    sr = sum(1.0 / d for d in ndegen)
    abs(sr - prod(mp_grid)) < 1e-8 ||
        error("Wigner–Seitz sum rule failed: Σ 1/ndegen = $sr ≠ ∏ mp_grid = $(prod(mp_grid))")
    return irvec, ndegen
end

"""
    build_hr(U, eig, kgrid, irvec) -> (Hr, Hk)

Wannier-gauge Hamiltonian in reciprocal space `Hk[:,:,k] = U_k† diag(ε_k) U_k` and its Fourier
transform to real space `Hr[:,:,R]` (undivided by ndegen). `U` is (num_wann × num_wann × nkpt),
`eig` is (num_bands × nkpt) with num_bands == num_wann for the isolated case.
"""
function build_hr(U::Array{ComplexF64,3}, eig::Matrix{Float64}, kgrid::KGrid,
                 irvec::Vector{NTuple{3,Int}})
    nw = size(U, 1); nk = size(U, 3); nr = length(irvec)
    Hk = Array{ComplexF64,3}(undef, nw, nw, nk)
    for k in 1:nk
        Uk = @view U[:, :, k]
        Hk[:, :, k] = Uk' * Diagonal(@view eig[:, k]) * Uk
    end
    Hr = zeros(ComplexF64, nw, nw, nr)
    for ir in 1:nr
        R = SVector{3,Float64}(irvec[ir]...)
        for k in 1:nk
            fac = cis(-TWOPI * dot(kgrid.frac[k], R)) / nk
            @views Hr[:, :, ir] .+= fac .* Hk[:, :, k]
        end
    end
    return Hr, Hk
end

"""
    interpolate_hk(Hr, irvec, ndegen, kfrac) -> H

Interpolated Wannier-gauge Hamiltonian at a single fractional k-point `kfrac`
(num_wann × num_wann, Hermitian).
"""
function interpolate_hk(Hr::Array{ComplexF64,3}, irvec::Vector{NTuple{3,Int}},
                       ndegen::Vector{Int}, kfrac::SVector{3,Float64})
    nw = size(Hr, 1)
    H = zeros(ComplexF64, nw, nw)
    for ir in 1:length(irvec)
        R = SVector{3,Float64}(irvec[ir]...)
        fac = cis(TWOPI * dot(kfrac, R)) / ndegen[ir]
        @views H .+= fac .* Hr[:, :, ir]
    end
    return H
end

"""
    generate_kpath(win, lattice; bands_num_points=100) -> (kpts, xvals, labels, label_idx)

Sample the `kpoint_path` block of a `.win` in segment mode: each segment gets a point count
proportional to its length (the first segment gets `bands_num_points`), linearly interpolated in
fractional coordinates. `xvals` is the cumulative Cartesian path length (Å⁻¹). `labels`/`label_idx`
mark the special points for axis ticks.
"""
function generate_kpath(win::WinInput, lattice::Lattice; bands_num_points::Int=100)
    haskey(win.blocks, "kpoint_path") || return (SVector{3,Float64}[], Float64[], String[], Int[])
    segs = Tuple{String,SVector{3,Float64},String,SVector{3,Float64}}[]
    for ln in win.blocks["kpoint_path"]
        t = split(ln)
        length(t) >= 8 || continue
        push!(segs, (String(t[1]), SVector{3,Float64}(parse_f64.(t[2:4])...),
                     String(t[5]), SVector{3,Float64}(parse_f64.(t[6:8])...)))
    end
    isempty(segs) && return (SVector{3,Float64}[], Float64[], String[], Int[])
    seglen = [norm(lattice.B * (s[4] - s[2])) for s in segs]
    base = seglen[1] > 0 ? seglen[1] : 1.0
    # Fortran `nint` rounds halves away from zero (Julia's default rounds to even).
    npts = [max(1, round(Int, bands_num_points * l / base, RoundNearestTiesAway)) for l in seglen]

    kpts = SVector{3,Float64}[]; xvals = Float64[]; labels = String[]; lidx = Int[]
    x = 0.0
    for (i, s) in enumerate(segs)
        push!(labels, s[1]); push!(lidx, length(kpts) + 1)     # segment-start special point
        for j in 0:npts[i]-1
            f = j / npts[i]
            push!(kpts, s[2] + (s[4] - s[2]) * f)
            push!(xvals, x + seglen[i] * f)
        end
        x += seglen[i]
        # Emit the segment endpoint here only if the path does not continue smoothly into the next
        # segment (last segment, or a discontinuity where the next segment starts elsewhere). For a
        # continuous join, the shared point is emitted as the next segment's first sample.
        is_last = i == length(segs)
        continuous = !is_last && norm(segs[i+1][2] - s[4]) < 1e-6
        if is_last || !continuous
            push!(kpts, s[4]); push!(xvals, x)
            push!(labels, s[3]); push!(lidx, length(kpts))
        end
    end
    return kpts, xvals, labels, lidx
end

"""
    interpolate_bands(Hr, irvec, ndegen, kpts) -> eigs

Interpolated band energies (num_wann × npts, ascending per column) at fractional k-points `kpts`.
"""
function interpolate_bands(Hr::Array{ComplexF64,3}, irvec::Vector{NTuple{3,Int}},
                          ndegen::Vector{Int}, kpts::Vector{SVector{3,Float64}})
    nw = size(Hr, 1)
    eigs = Matrix{Float64}(undef, nw, length(kpts))
    for (ik, kf) in enumerate(kpts)
        H = interpolate_hk(Hr, irvec, ndegen, kf)
        eigs[:, ik] = eigvals(Hermitian((H + H') / 2))
    end
    return eigs
end

# ---------------------------------------------------------------------------
# use_ws_distance: per-Wannier-pair minimal-image improvement (default in the reference).
# For each pair (i,j) and each R, translate WF j by supercell vectors so it lands in the
# Wigner–Seitz cell around WF i, giving a pair-specific set of equivalent R-vectors with
# degeneracy ndeg(i,j,R). Mirrors ws_distance.F90 (R_wz_sc) and plot.F90:786-797.
# ---------------------------------------------------------------------------

"Place `Rin` (Cartesian) in the Wigner–Seitz cell around the origin; return the equivalent integer
lattice shifts (multiples of mp_grid) and their degeneracy."
function _R_wz_sc(Rin::SVector{3,Float64}, lat::Lattice, mp_grid::NTuple{3,Int},
                 ws_search_size::Int, tol::Float64)
    Rin_f = lat.A \ Rin
    mod2 = sum(abs2, Rin)
    best = (0, 0, 0)
    rng = -(ws_search_size + 1):(ws_search_size + 1)
    for i in rng, j in rng, k in rng
        s = SVector{3,Float64}(i * mp_grid[1], j * mp_grid[2], k * mp_grid[3])
        R = lat.A * (Rin_f + s)
        d2 = sum(abs2, R)
        if d2 < mod2 - 1e-13
            mod2 = d2
            best = (i * mp_grid[1], j * mp_grid[2], k * mp_grid[3])
        end
    end
    if mod2 < tol^2
        return NTuple{3,Int}[best]
    end
    Rbz_f = lat.A \ (lat.A * (Rin_f + SVector{3,Float64}(best...)))   # = fractional of R_bz
    dref = sqrt(mod2)
    shifts = NTuple{3,Int}[]
    for i in rng, j in rng, k in rng
        s = SVector{3,Float64}(i * mp_grid[1], j * mp_grid[2], k * mp_grid[3])
        R = lat.A * (Rbz_f + s)
        if abs(sqrt(sum(abs2, R)) - dref) < tol
            push!(shifts, (best[1] + i * mp_grid[1], best[2] + j * mp_grid[2], best[3] + k * mp_grid[3]))
        end
    end
    return shifts
end

"""
    ws_translate_dist(irvec, centres, lattice, mp_grid; ws_search_size=2, tol=1e-5) -> irdist

For each Wannier pair (i,j) and each R (column of `irvec`), the minimal-image-equivalent integer
R-vectors. `irdist[i,j,ir]` is a `Vector{NTuple{3,Int}}` of length ndeg(i,j,ir). `centres` is the
(3 × num_wann) Cartesian Wannier centres.
"""
function ws_translate_dist(irvec::Vector{NTuple{3,Int}}, centres::Matrix{Float64}, lat::Lattice,
                          mp_grid::NTuple{3,Int}; ws_search_size::Int=2, tol::Float64=1.0e-5)
    nw = size(centres, 2); nr = length(irvec)
    irdist = Array{Vector{NTuple{3,Int}}}(undef, nw, nw, nr)
    for ir in 1:nr
        Rcart = lat.A * SVector{3,Float64}(irvec[ir]...)
        for jw in 1:nw, iw in 1:nw
            v = Rcart + SVector{3,Float64}(centres[:, jw]) - SVector{3,Float64}(centres[:, iw])
            shifts = _R_wz_sc(v, lat, mp_grid, ws_search_size, tol)
            irdist[iw, jw, ir] = [(irvec[ir][1] + s[1], irvec[ir][2] + s[2], irvec[ir][3] + s[3])
                                  for s in shifts]
        end
    end
    return irdist
end

"Interpolated H(k') with the per-pair minimal-image (ws_distance) weighting."
function interpolate_hk_ws(Hr::Array{ComplexF64,3}, ndegen::Vector{Int},
                          irdist::Array{Vector{NTuple{3,Int}}}, kfrac::SVector{3,Float64})
    nw = size(Hr, 1); nr = size(irdist, 3)
    H = zeros(ComplexF64, nw, nw)
    for ir in 1:nr, jw in 1:nw, iw in 1:nw
        dl = irdist[iw, jw, ir]
        nd = length(dl)
        acc = zero(ComplexF64)
        for R in dl
            acc += cis(TWOPI * dot(kfrac, SVector{3,Float64}(R...)))
        end
        H[iw, jw] += acc / (ndegen[ir] * nd) * Hr[iw, jw, ir]
    end
    return H
end

"""
    interpolate_bands_ws(Hr, irvec, ndegen, centres, lattice, mp_grid, kpts) -> eigs

Band interpolation with the `use_ws_distance=.true.` minimal-image improvement (the reference
default). `centres` are the (3 × num_wann) Cartesian Wannier centres.
"""
function interpolate_bands_ws(Hr::Array{ComplexF64,3}, irvec::Vector{NTuple{3,Int}},
                             ndegen::Vector{Int}, centres::Matrix{Float64}, lattice::Lattice,
                             mp_grid::NTuple{3,Int}, kpts::Vector{SVector{3,Float64}})
    irdist = ws_translate_dist(irvec, centres, lattice, mp_grid)
    nw = size(Hr, 1)
    eigs = Matrix{Float64}(undef, nw, length(kpts))
    for (ik, kf) in enumerate(kpts)
        H = interpolate_hk_ws(Hr, ndegen, irdist, kf)
        eigs[:, ik] = eigvals(Hermitian((H + H') / 2))
    end
    return eigs
end
