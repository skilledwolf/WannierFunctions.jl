# Real-space Wannier-function plotting (`wannier_plot = .true.`): read the periodic parts
# u_{bk}(r) from UNK files, assemble w_n(r) on a plotting supercell in the converged gauge, and
# write XCrySDen .xsf volumetric files. Follows plot.F90 (plot_wannier); see
# docs/reference-notes/wannier-plot.md for the exact conventions.
#
#   w_n(r) = (1/N_k) Σ_k e^{i k·r} Σ_b [U_opt(k) U(k)]_{b n} u_{b k}(r mod cell)
#
# Grid index n ↔ fractional coordinate (n−1)/ng (index 1 is the origin). The global phase of
# each WF is fixed so the point of maximum |w| is real and positive.

using Printf
using StaticArrays

"""
    read_unk(path) -> (ng, ik, u)

Read a formatted UNK file: `ng = (ngx, ngy, ngz)`, the k-index `ik`, and
`u[npoint, band]` with the x-fastest linear index `npoint = nx + (ny−1)ngx + (nz−1)ngx·ngy`.
"""
function read_unk(path::AbstractString)
    open(path, "r") do io
        hdr = parse.(Int, split(readline(io)))
        ngx, ngy, ngz, ik, nbnd = hdr
        npts = ngx * ngy * ngz
        u = Matrix{ComplexF64}(undef, npts, nbnd)
        for b in 1:nbnd, p in 1:npts
            t = split(readline(io))
            u[p, b] = complex(parse(Float64, t[1]), parse(Float64, t[2]))
        end
        return (ngx, ngy, ngz), ik, u
    end
end

"Parse a Wannier90 range-vector string (`\"1-4,6\"`) into a sorted index list."
function parse_range_list(str::AbstractString)
    out = Int[]
    for tok in split(str, ',')
        t = strip(tok)
        isempty(t) && continue
        m = match(r"^(\d+)\s*-\s*(\d+)$", t)
        m === nothing ? push!(out, parse(Int, t)) :
                        append!(out, parse(Int, m.captures[1]):parse(Int, m.captures[2]))
    end
    return sort!(unique!(out))
end

"""
    wannier_function_grid(model, win, result; list, supercell=(2,2,2), dir=".")
        -> (w, ng, los)

Assemble the plotted Wannier functions on the supercell grid. Returns the complex 4-D array
`w[ix, iy, iz, n]` (axes are the supercell grid offsets `los[d]:los[d]+ngs[d]*ng[d]-1`), the
home-cell grid dims `ng`, and the lower bounds `los`. UNK files are read from `dir`.
"""
function wannier_function_grid(model::Model, win::WinInput, result::WannierResult;
                               list::Vector{Int}=collect(1:model.num_wann),
                               supercell::NTuple{3,Int}=(2, 2, 2), dir::AbstractString=".")
    nk = nkpt(model.kgrid)
    nw = model.num_wann
    all(1 .<= list .<= nw) || error("wannier_plot_list entries must be in 1..$nw")

    # Window bookkeeping for the disentangled case (UNK holds all bands; keep lwindow ones).
    wd = nothing
    if result.disentangled
        wd = dis_windows(model.eig, nw;
                         win_min=win.dis_win_min, win_max=win.dis_win_max,
                         froz_min=win.dis_froz_min,
                         froz_max=(win.dis_froz_max == -Inf ? nothing : win.dis_froz_max))
    end

    # Grid dims from the first UNK.
    ng, _, _ = read_unk(joinpath(dir, @sprintf("UNK%05d.%1d", 1, 1)))
    ngx, ngy, ngz = ng
    los = (-(supercell[1] ÷ 2) * ngx, -(supercell[2] ÷ 2) * ngy, -(supercell[3] ÷ 2) * ngz)
    his = (((supercell[1] + 1) ÷ 2) * ngx - 1, ((supercell[2] + 1) ÷ 2) * ngy - 1,
           ((supercell[3] + 1) ÷ 2) * ngz - 1)
    dims = his .- los .+ 1
    w = zeros(ComplexF64, dims[1], dims[2], dims[3], length(list))

    for k in 1:nk
        ngk, ik, u = read_unk(joinpath(dir, @sprintf("UNK%05d.%1d", k, 1)))
        (ngk == ng && ik == k) || error("UNK file $k: header mismatch")

        # Gauge contraction c[:, w] = Σ_b [U_opt·U]_{b, list[w]} u[:, b].
        if result.disentangled
            keep = findall(i -> wd.nfirstwin[k] <= i < wd.nfirstwin[k] + wd.ndimwin[k],
                           1:model.num_bands)
            V = result.dis.Uopt[k] * result.U[:, :, k]     # ndimwin × num_wann
            c = u[:, keep] * V[:, list]
        else
            c = u * result.U[:, list, k]
        end

        kf = model.kgrid.frac[k]
        phx = [cis(TWOPI * kf[1] * (nxx - 1) / ngx) for nxx in los[1]:his[1]]
        phy = [cis(TWOPI * kf[2] * (nyy - 1) / ngy) for nyy in los[2]:his[2]]
        phz = [cis(TWOPI * kf[3] * (nzz - 1) / ngz) for nzz in los[3]:his[3]]

        Threads.@threads for izz in 1:dims[3]
            nz = mod1(los[3] + izz - 1, ngz)
            for iyy in 1:dims[2]
                ny = mod1(los[2] + iyy - 1, ngy)
                base = (ny - 1) * ngx + (nz - 1) * ngx * ngy
                for ixx in 1:dims[1]
                    nx = mod1(los[1] + ixx - 1, ngx)
                    ph = phx[ixx] * phy[iyy] * phz[izz]
                    p = nx + base
                    @inbounds for n in 1:length(list)
                        w[ixx, iyy, izz, n] += c[p, n] * ph
                    end
                end
            end
        end
    end
    w ./= nk

    # Fix each WF's global phase: max-|w| point real positive.
    for n in 1:length(list)
        v = @view w[:, :, :, n]
        wmod = v[argmax(abs2.(v))]
        v ./= (wmod / abs(wmod))
    end
    return w, ng, los
end

"""
    write_xsf(path, lattice, atoms, w, ng, los, supercell)

Write one WF's real part as an XCrySDen .xsf (crystal mode, Å): CRYSTAL/PRIMVEC/CONVVEC/
PRIMCOORD blocks, then a general 3-D datagrid over the whole plotting supercell.
"""
function write_xsf(path::AbstractString, lattice::Lattice,
                   atoms::Vector{Tuple{String,SVector{3,Float64}}},
                   w::AbstractArray{ComplexF64,3}, ng::NTuple{3,Int}, los::NTuple{3,Int},
                   supercell::NTuple{3,Int})
    A = lattice.A
    open(path, "w") do io
        println(io, "      # Generated by WannierFunctions.jl")
        println(io, "      #")
        println(io, "      #")
        println(io, "      #")
        println(io, "CRYSTAL")
        println(io, "PRIMVEC")
        for i in 1:3
            @printf(io, "%12.7f%12.7f%12.7f\n", A[1, i], A[2, i], A[3, i])
        end
        println(io, "CONVVEC")
        for i in 1:3
            @printf(io, "%12.7f%12.7f%12.7f\n", A[1, i], A[2, i], A[3, i])
        end
        println(io, "PRIMCOORD")
        @printf(io, "%6d  1\n", length(atoms))
        for (sp, frac) in atoms
            cart = A * frac
            @printf(io, "%-2s   %12.7f%12.7f%12.7f\n", sp, cart[1], cart[2], cart[3])
        end
        println(io)
        println(io, "BEGIN_BLOCK_DATAGRID_3D")
        println(io, "3D_field")
        println(io, "BEGIN_DATAGRID_3D_UNKNOWN")
        N = supercell .* ng
        @printf(io, "%6d%6d%6d\n", N[1], N[2], N[3])
        # origin = grid point (los .- 1)/ng in fractional coordinates
        orig = A * SVector((los[1] - 1) / ng[1], (los[2] - 1) / ng[2], (los[3] - 1) / ng[3])
        @printf(io, "%12.6f%12.6f%12.6f\n", orig[1], orig[2], orig[3])
        for j in 1:3
            span = A[:, j] * ((N[j] - 1) / ng[j])
            @printf(io, "%12.7f%12.7f%12.7f\n", span[1], span[2], span[3])
        end
        vals = vec(real.(w))                              # x fastest, matching array order
        for i in 1:6:length(vals)
            println(io, join((fortran_e(v, 13, 5) for v in vals[i:min(i + 5, end)]), ""))
        end
        println(io, "END_DATAGRID_3D")
        println(io, "END_BLOCK_DATAGRID_3D")
    end
    return path
end

"""
    plot_wannier_functions(model, win, result; seedname, dir) -> paths

Assemble and write `seedname_0000n.xsf` for each WF in `wannier_plot_list` (default: all).
Honours `wannier_plot_supercell`. Only the formatted-UNK, xsf, crystal-mode path is
implemented (`wannier_plot_format = cube` is not yet supported).
"""
function plot_wannier_functions(model::Model, win::WinInput, result::WannierResult;
                                seedname::AbstractString=model.seedname,
                                dir::AbstractString=dirname(abspath(seedname * ".win")))
    fmt = lowercase(get(win.raw, "wannier_plot_format", "xcrysden"))
    occursin("xcrys", fmt) || error("wannier_plot_format=$fmt not supported (xsf only)")
    _getbool(win.raw, "wvfn_formatted", false) ||
        @warn "wvfn_formatted=.false.: attempting formatted UNK read anyway" maxlog = 1
    list = haskey(win.raw, "wannier_plot_list") ?
           parse_range_list(win.raw["wannier_plot_list"]) : collect(1:model.num_wann)
    sc = if haskey(win.raw, "wannier_plot_supercell")
        v = parse.(Int, split(win.raw["wannier_plot_supercell"]))
        length(v) == 1 ? (v[1], v[1], v[1]) : (v[1], v[2], v[3])
    else
        (2, 2, 2)
    end
    w, ng, los = wannier_function_grid(model, win, result; list=list, supercell=sc, dir=dir)
    atoms = parse_atoms(win)
    paths = String[]
    for (i, n) in enumerate(list)
        p = @sprintf("%s_%05d.xsf", seedname, n)
        write_xsf(p, model.lattice, atoms, (@view w[:, :, :, i]), ng, los, sc)
        push!(paths, p)
    end
    return paths
end
