# From-scratch k-mesh generation and the `.nnkp` file — the `-pp` (post-processing setup) mode
# that lets a DFT interface (pw2wannier90, ...) know which overlaps and projections to compute.
#
# This is the *generation* side of the k-mesh: neighbour shells are searched from the lattice and
# mp_grid alone (no .mmn needed), shells are selected to satisfy the B1 completeness relation, and
# the per-k neighbour list is enumerated in the reference's exact deterministic order so the
# resulting .nnkp matches `wannier90.x -pp` byte-for-byte (modulo the date header).
# Port of kmesh.F90 (kmesh_get / kmesh_supercell_sort / kmesh_shell_automatic / kmesh_write).

using LinearAlgebra
using StaticArrays

const NSUPCELL = 5              # kmesh.F90:68
const EPS5, EPS6, EPS8 = 1e-5, 1e-6, 1e-8

"""
    supercell_cells(B) -> Vector{NTuple{3,Int}}

The (2·NSUPCELL+1)³ reciprocal-superlattice cells, ordered exactly as the reference: (0,0,0)
first, then ascending by |l·b₁+m·b₂+n·b₃|, ties resolved by the reference's extract-max /
lowest-index rule (kmesh_supercell_sort + internal_maxloc).
"""
function supercell_cells(lattice::Lattice)
    n = (2 * NSUPCELL + 1)^3
    cells = Vector{NTuple{3,Int}}(undef, n)
    dist = Vector{Float64}(undef, n)
    cells[1] = (0, 0, 0); dist[1] = 0.0
    c = 1
    for l in -NSUPCELL:NSUPCELL, m in -NSUPCELL:NSUPCELL, nn in -NSUPCELL:NSUPCELL
        (l == 0 && m == 0 && nn == 0) && continue
        c += 1
        cells[c] = (l, m, nn)
        dist[c] = norm(lattice.B * SVector{3,Float64}(l, m, nn))
    end
    # Literal port of the extraction sort: repeatedly take the max (lowest index among eps8
    # ties) and place it at the end. Within a tie group this reverses enumeration order.
    d = copy(dist)
    out = Vector{NTuple{3,Int}}(undef, n)
    for pos in n:-1:1
        mx = maximum(d)
        idx = findfirst(x -> abs(x - mx) < EPS8, d)::Int
        out[pos] = cells[idx]
        d[idx] = -1.0
    end
    return out
end

"""
    shell_distances(kcart, lattice, cells; search_shells, tol) -> (dnn, multi)

Distances of the first `search_shells` neighbour shells of k-point 1 (over all k+G images) and
their multiplicities. Port of kmesh.F90:172-199.
"""
function shell_distances(kcart::Vector{SVector{3,Float64}}, lattice::Lattice,
                        cells::Vector{NTuple{3,Int}}; search_shells::Int, tol::Float64)
    dnn = Float64[]; multi = Int[]
    dnn0 = 0.0
    for _ in 1:search_shells
        dnn1 = Inf; counter = 0
        for k2 in kcart, cell in cells
            v = k2 + lattice.B * SVector{3,Float64}(cell...)
            d = norm(kcart[1] - v)
            if d > tol && d > dnn0 + tol
                if d < dnn1 - tol
                    dnn1 = d; counter = 0
                end
                (dnn1 - tol < d < dnn1 + tol) && (counter += 1)
            end
        end
        isfinite(dnn1) || break
        push!(dnn, dnn1); push!(multi, counter)
        dnn0 = dnn1
    end
    return dnn, multi
end

"b-vectors of one shell around k-point 1, in reference discovery order (kmesh_get_bvectors)."
function shell_bvectors(kcart, lattice, cells, shell_dist, mult; tol)
    bs = SVector{3,Float64}[]
    for cell in cells
        v2 = lattice.B * SVector{3,Float64}(cell...)
        for k2 in kcart
            v = v2 + k2
            d = norm(kcart[1] - v)
            (shell_dist - tol <= d <= shell_dist + tol) && push!(bs, v - kcart[1])
            length(bs) == mult && return bs
        end
    end
    length(bs) == mult || error("shell_bvectors: found $(length(bs)) of $mult b-vectors")
    return bs
end

"""
    select_shells(kcart, lattice, cells, dnn, multi; kmesh_tol) -> (shell_list, weights)

Automatic shell selection satisfying the B1 relation (kmesh_shell_automatic): add shells in
distance order, skipping shells parallel to accepted ones (|cosine| within 1e-6 of 1) and
rejecting additions that make the design matrix singular (singular value < 1e-5), until
Σ_s w_s Σ_b b⊗b = 1 holds within `kmesh_tol`.
"""
function select_shells(kcart, lattice, cells, dnn::Vector{Float64}, multi::Vector{Int};
                      kmesh_tol::Float64)
    target = [1.0, 0.0, 1.0, 0.0, 0.0, 1.0]            # (xx,xy,yy,xz,yz,zz)
    shell_list = Int[]
    bvecs = Vector{Vector{SVector{3,Float64}}}()        # accepted shells' b-vectors
    weights = Float64[]
    for shell in 1:length(dnn)
        bnew = shell_bvectors(kcart, lattice, cells, dnn[shell], multi[shell]; tol=kmesh_tol)
        # parallel-shell rejection
        lpar = false
        for bn in bnew, bset in bvecs, b in bset
            cosang = dot(bn, b) / (norm(bn) * norm(b))
            abs(abs(cosang) - 1.0) < EPS6 && (lpar = true)
        end
        lpar && continue
        push!(shell_list, shell); push!(bvecs, bnew)
        # design matrix rows (xx,xy,yy,xz,yz,zz) per accepted shell
        A = zeros(6, length(shell_list))
        for (s, bset) in enumerate(bvecs), b in bset
            A[1, s] += b[1]^2;    A[2, s] += b[1] * b[2]; A[3, s] += b[2]^2
            A[4, s] += b[1] * b[3]; A[5, s] += b[2] * b[3]; A[6, s] += b[3]^2
        end
        F = svd(A)
        if any(s -> abs(s) < EPS5, F.S)
            length(shell_list) == 1 &&
                error("select_shells: singular design matrix on the first shell")
            pop!(shell_list); pop!(bvecs)               # reject this shell, try the next
            continue
        end
        w = F.V * (Diagonal(1.0 ./ F.S) * (F.U' * target))
        if maximum(abs.(A * w - target)) < kmesh_tol    # B1 satisfied
            weights = w
            break
        end
    end
    isempty(weights) && error("select_shells: B1 completeness not satisfied within the searched shells")
    return shell_list, weights
end

"""
    build_nnlist(kfrac, lattice, cells, dnn, shell_list, multi; tol)
        -> (nnlist, nncell, nntot)

Per-k neighbour list in the reference's exact deterministic order (shell-major, then
distance-sorted supercell cells, then k₂ ascending; early exit per filled shell).
Port of kmesh.F90:423-450. `nnlist[k, nn]` is the neighbour k-index, `nncell[:, k, nn]` the
reciprocal-lattice fold G.
"""
function build_nnlist(kfrac::Vector{SVector{3,Float64}}, lattice::Lattice,
                     cells::Vector{NTuple{3,Int}}, dnn::Vector{Float64},
                     shell_list::Vector{Int}, multi::Vector{Int}; tol::Float64)
    kcart = [lattice.B * k for k in kfrac]
    nk = length(kfrac)
    nntot = sum(multi[s] for s in shell_list)
    nnlist = zeros(Int, nk, nntot)
    nncell = zeros(Int, 3, nk, nntot)
    for k in 1:nk
        nnx = 0
        for s in shell_list
            found = 0
            for cell in cells
                v2 = lattice.B * SVector{3,Float64}(cell...)
                for k2 in 1:nk
                    v = v2 + kcart[k2]
                    d = norm(kcart[k] - v)
                    if dnn[s] - tol <= d <= dnn[s] + tol
                        nnx += 1; found += 1
                        nnlist[k, nnx] = k2
                        nncell[:, k, nnx] .= cell
                    end
                    found == multi[s] && break
                end
                found == multi[s] && break
            end
            found == multi[s] || error("build_nnlist: k=$k shell=$s found $found of $(multi[s])")
        end
    end
    return nnlist, nncell, nntot
end

# ---------------------------------------------------------------------------
# Projections: parse the .win `projections` block into per-orbital entries.
# ---------------------------------------------------------------------------

"""
One projection: site (fractional), angular character (l, mr), radial node count, axes; for
spinor projections `s = ±1` selects the spin component along `s_qaxis` (`s = 0` ⇒ spatial).
"""
struct Projection
    site::SVector{3,Float64}
    l::Int
    mr::Int
    radial::Int
    z::SVector{3,Float64}
    x::SVector{3,Float64}
    zona::Float64
    s::Int
    s_qaxis::SVector{3,Float64}
end
Projection(site, l, mr, radial, z, x, zona) =
    Projection(site, l, mr, radial, z, x, zona, 0, SVector(0.0, 0.0, 1.0))

# orbital label → (l, [mr...]); bare "p"/"d"/"f" expand to all mr. Wannier90 user-guide table 3.1.
const ORBITALS = Dict{String,Tuple{Int,Vector{Int}}}(
    "s" => (0, [1]),
    "p" => (1, [1, 2, 3]), "pz" => (1, [1]), "px" => (1, [2]), "py" => (1, [3]),
    "d" => (2, [1, 2, 3, 4, 5]), "dz2" => (2, [1]), "dxz" => (2, [2]), "dyz" => (2, [3]),
    "dx2-y2" => (2, [4]), "dxy" => (2, [5]),
    "f" => (3, [1, 2, 3, 4, 5, 6, 7]),
    "sp" => (-1, [1, 2]), "sp2" => (-2, [1, 2, 3]), "sp3" => (-3, [1, 2, 3, 4]),
    "sp3d" => (-4, [1, 2, 3, 4, 5]), "sp3d2" => (-5, [1, 2, 3, 4, 5, 6]),
)

"Parse the atoms_frac / atoms_cart block into (species, fractional site) pairs."
function parse_atoms(win::WinInput)
    atoms = Tuple{String,SVector{3,Float64}}[]
    if haskey(win.blocks, "atoms_frac")
        for ln in win.blocks["atoms_frac"]
            t = split(ln)
            length(t) >= 4 || continue
            push!(atoms, (String(t[1]), SVector{3,Float64}(parse_f64.(t[2:4])...)))
        end
    elseif haskey(win.blocks, "atoms_cart")
        scale = 1.0
        Ainv = inv(win.unit_cell)
        for ln in win.blocks["atoms_cart"]
            t = split(ln)
            if length(t) == 1
                u = lowercase(t[1])
                (u == "bohr" || u == "b") && (scale = BOHR)
                continue
            end
            length(t) >= 4 || continue
            cart = scale .* SVector{3,Float64}(parse_f64.(t[2:4])...)
            push!(atoms, (String(t[1]), Ainv * cart))
        end
    end
    return atoms
end

"""
    parse_projections(win) -> Vector{Projection}

Parse the `projections` block: `Species:orb[;orb2...]`, `f=x,y,z:orb`, or `c=x,y,z:orb` sites
(a species expands to every atom of that species, in atoms-block order). Default z-axis (0,0,1),
x-axis (1,0,0), zona 1.0, radial 1 — the reference defaults.
"""
function parse_projections(win::WinInput; spinors::Bool=_getbool(win.raw, "spinors", false))
    haskey(win.blocks, "projections") || return Projection[]
    atoms = parse_atoms(win)
    projs = Projection[]
    for ln in win.blocks["projections"]
        lc = strip(ln)
        (isempty(lc) || lowercase(lc) in ("random", "bohr", "ang")) && continue
        parts = split(lc, ':')
        length(parts) >= 2 || error("projections: cannot parse line `$ln`")
        sitespec = strip(parts[1])
        orbspec = strip(parts[2])            # axis/zona modifiers (parts 3+) not yet supported
        length(parts) > 2 && @warn "projections: axis/zona modifiers not supported, using defaults" line = ln maxlog = 1
        # resolve site(s)
        sites = SVector{3,Float64}[]
        low = lowercase(sitespec)
        if startswith(low, "f=")
            push!(sites, SVector{3,Float64}(parse_f64.(split(sitespec[3:end], ','))...))
        elseif startswith(low, "c=")
            cart = SVector{3,Float64}(parse_f64.(split(sitespec[3:end], ','))...)
            push!(sites, inv(win.unit_cell) * cart)
        else
            for (sp, frac) in atoms
                sp == sitespec && push!(sites, frac)
            end
            isempty(sites) && error("projections: no atoms of species `$sitespec` " *
                                    "(have: $(unique(first.(atoms))))")
        end
        # resolve orbitals (';'- or ','-separated, as the reference accepts both:
        # "O:s,p" ≡ "O:s;p"). The reference emits a line's orbitals in ascending-l order
        # regardless of spec order ("Pt: d;s;p" → s, p, d).
        orbs = Tuple{Int,Vector{Int}}[]
        for orb in split(orbspec, [';', ','])
            key = lowercase(strip(orb))
            haskey(ORBITALS, key) || error("projections: unknown orbital `$orb` " *
                                           "(supported: $(join(sort(collect(keys(ORBITALS))), ", ")))")
            push!(orbs, ORBITALS[key])
        end
        sort!(orbs; by=first)
        for site in sites, (l, mrs) in orbs, mr in mrs
            if spinors
                for spin in (1, -1)     # up/down interleaved, as the reference writes them
                    push!(projs, Projection(site, l, mr, 1, SVector(0.0, 0.0, 1.0),
                                            SVector(1.0, 0.0, 0.0), 1.0, spin,
                                            SVector(0.0, 0.0, 1.0)))
                end
            else
                push!(projs, Projection(site, l, mr, 1,
                                        SVector(0.0, 0.0, 1.0), SVector(1.0, 0.0, 0.0), 1.0))
            end
        end
    end
    return projs
end

# ---------------------------------------------------------------------------
# .nnkp writer (kmesh_write, kmesh.F90:962-1128) — exact Fortran formats.
# ---------------------------------------------------------------------------

_f(x, w, d) = lpad(@sprintf("%.*f", d, x), w)

"""
    write_nnkp(path, lattice, kfrac, projs, nnlist, nncell; exclude_bands=Int[], calc_only_A=false)

Write a `.nnkp` file in the reference format (drop-in for `wannier90.x -pp` output).
"""
function write_nnkp(path::AbstractString, lattice::Lattice, kfrac::Vector{SVector{3,Float64}},
                   projs::Vector{Projection}, nnlist::Matrix{Int}, nncell::Array{Int,3};
                   exclude_bands::Vector{Int}=Int[], calc_only_A::Bool=false)
    nk = length(kfrac)
    nntot = size(nnlist, 2)
    open(path, "w") do io
        println(io, "# File written by WannierFunctions.jl\n")
        println(io, "calc_only_A  :  ", calc_only_A ? "T" : "F", "\n")
        println(io, "begin real_lattice")
        for i in 1:3   # row i = lattice vector a_i (columns of lattice.A)
            println(io, _f(lattice.A[1, i], 12, 7), _f(lattice.A[2, i], 12, 7), _f(lattice.A[3, i], 12, 7))
        end
        println(io, "end real_lattice\n")
        println(io, "begin recip_lattice")
        for i in 1:3   # row i = reciprocal vector b_i
            println(io, _f(lattice.B[1, i], 12, 7), _f(lattice.B[2, i], 12, 7), _f(lattice.B[3, i], 12, 7))
        end
        println(io, "end recip_lattice\n")
        println(io, "begin kpoints")
        println(io, lpad(nk, 6))
        for k in kfrac
            println(io, _f(k[1], 14, 8), _f(k[2], 14, 8), _f(k[3], 14, 8))
        end
        println(io, "end kpoints\n")
        spinor = !isempty(projs) && projs[1].s != 0
        blockname = spinor ? "spinor_projections" : "projections"
        println(io, "begin ", blockname)
        println(io, lpad(length(projs), 6))
        for p in projs
            println(io, _f(p.site[1], 10, 5), " ", _f(p.site[2], 10, 5), " ", _f(p.site[3], 10, 5),
                    "   ", lpad(p.l, 3), lpad(p.mr, 3), lpad(p.radial, 3))
            println(io, "  ", _f(p.z[1], 11, 7), _f(p.z[2], 11, 7), _f(p.z[3], 11, 7), " ",
                    _f(p.x[1], 11, 7), _f(p.x[2], 11, 7), _f(p.x[3], 11, 7), " ", _f(p.zona, 7, 2))
            spinor && println(io, "  ", lpad(p.s, 3), " ", _f(p.s_qaxis[1], 11, 7),
                              _f(p.s_qaxis[2], 11, 7), _f(p.s_qaxis[3], 11, 7))
        end
        println(io, "end ", blockname, "\n")
        println(io, "begin nnkpts")
        println(io, lpad(nntot, 4))
        for k in 1:nk, nn in 1:nntot
            println(io, lpad(k, 8), lpad(nnlist[k, nn], 8), "   ",
                    lpad(nncell[1, k, nn], 4), lpad(nncell[2, k, nn], 4), lpad(nncell[3, k, nn], 4))
        end
        println(io, "end nnkpts\n")
        println(io, "begin exclude_bands")
        println(io, lpad(length(exclude_bands), 4))
        for b in exclude_bands
            println(io, lpad(b, 4))
        end
        println(io, "end exclude_bands")
    end
    return path
end

"""
    parse_exclude_bands(win) -> Vector{Int}

Parse the `exclude_bands` range list (`exclude_bands = 1-4,10`) into a sorted index list.
"""
function parse_exclude_bands(win::WinInput)
    haskey(win.raw, "exclude_bands") || return Int[]
    out = Int[]
    for tok in split(win.raw["exclude_bands"], ',')
        t = strip(tok)
        isempty(t) && continue
        m = match(r"^(\d+)\s*-\s*(\d+)$", t)
        if m !== nothing
            append!(out, parse(Int, m.captures[1]):parse(Int, m.captures[2]))
        else
            push!(out, parse(Int, t))
        end
    end
    return sort!(unique!(out))
end

"""
    write_nnkp(seedname_or_win; out=...) — the `-pp` mode

Generate the k-mesh from `seedname.win` alone (shell search + B1 weights + neighbour list, no
`.mmn` required) and write `seedname.nnkp` for the DFT interface. Returns the output path.
"""
function generate_nnkp(seedname::AbstractString; out::AbstractString=seedname * ".nnkp")
    win = read_win(seedname * ".win")
    lattice = Lattice(win.unit_cell)
    cells = supercell_cells(lattice)
    kcart = [lattice.B * k for k in win.kpoints]
    dnn, multi = shell_distances(kcart, lattice, cells;
                                 search_shells=win.search_shells, tol=win.kmesh_tol)
    shell_list, weights = select_shells(kcart, lattice, cells, dnn, multi; kmesh_tol=win.kmesh_tol)
    nnlist, nncell, nntot = build_nnlist(win.kpoints, lattice, cells, dnn, shell_list, multi;
                                         tol=win.kmesh_tol)
    horder = _getint(win.raw, "higher_order_n", 1)
    if horder > 1
        # Higher-order finite differences (Lihm): append the 2b … Nb multiples of every
        # first-order b, in the reference's canonical order [block1, 2·block1, …]. Each
        # multiple is located by exact match of k + m·b folded into the mesh (kmesh.F90:491).
        win.gamma_only && error("higher_order_n with gamma_only is not supported")
        nk = length(win.kpoints)
        nnlist2 = zeros(Int, nk, nntot * horder)
        nncell2 = zeros(Int, 3, nk, nntot * horder)
        nnlist2[:, 1:nntot] = nnlist
        nncell2[:, :, 1:nntot] = nncell
        for k in 1:nk, m in 2:horder, nn in 1:nntot
            bfrac = m .* (win.kpoints[nnlist[k, nn]] +
                          SVector{3,Float64}(nncell[:, k, nn]...) - win.kpoints[k])
            lmn = floor.(Int, win.kpoints[k] + bfrac .+ 1e-6)
            tgt = win.kpoints[k] + bfrac - lmn
            k2 = findfirst(kf -> norm(kf - tgt) < win.kmesh_tol, win.kpoints)
            k2 === nothing && error("generate_nnkp: could not find k+$(m)b neighbour on the mesh")
            nnx = (m - 1) * nntot + nn
            nnlist2[k, nnx] = k2
            nncell2[:, k, nnx] .= lmn
        end
        nnlist, nncell, nntot = nnlist2, nncell2, nntot * horder
    end
    if win.gamma_only
        # Keep the first of each ±b pair, in discovery order (kmesh.F90:833-861); the DFT
        # interface then computes overlaps for half the b-vectors only.
        length(win.kpoints) == 1 || error("gamma_only requires a single k-point")
        bcart(nn) = lattice.B * (win.kpoints[nnlist[1, nn]] +
                                 SVector{3,Float64}(nncell[:, 1, nn]...) - win.kpoints[1])
        keep = Int[]
        for nn in 1:nntot
            any(k -> norm(bcart(k) + bcart(nn)) < 1e-8, keep) || push!(keep, nn)
        end
        length(keep) == nntot ÷ 2 || error("gamma_only: could not pair b-vectors")
        nnlist = nnlist[:, keep]
        nncell = nncell[:, :, keep]
        nntot = length(keep)
    end
    projs = parse_projections(win)
    excl = parse_exclude_bands(win)
    write_nnkp(out, lattice, win.kpoints, projs, nnlist, nncell; exclude_bands=excl)
    return out, (; shells=shell_list, weights, nntot, dnn=dnn[shell_list], multi=multi[shell_list])
end
