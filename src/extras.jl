# Small wannier90.x output writers: Gaussian .cube Wannier-function plots, the `_r.dat`
# position matrix elements (write_rmn), and the XCrySDen .bxsf Fermi-surface grid.
# Exact conventions in docs/reference-notes/w90x-outputs.md.

using Printf
using StaticArrays

const _PERIODIC_TABLE = split(
    "h he li be b c n o f ne na mg al si p s cl ar k ca sc ti v cr mn fe co ni cu zn ga ge " *
    "as se br kr rb sr y zr nb mo tc ru rh pd ag cd in sn sb te i xe cs ba la ce pr nd pm sm " *
    "eu gd tb dy ho er tm yb lu hf ta w re os ir pt au hg tl pb bi po at rn fr ra ac th pa u " *
    "np pu am cm bk cf es fm md no lr rf db sg bh hs mt ds rg cn")

_atomic_number(sym::AbstractString) =
    something(findfirst(==(lowercase(strip(sym))), _PERIODIC_TABLE), 0)

"gfortran-style date/time pair, e.g. (\" 4Jul2026\", \"19:40:52\")."
function _w90_datetime()
    t = _NOW[]
    months = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    cdate = @sprintf("%2d%s%4d", t.day, months[t.month], t.year)
    ctime = @sprintf("%02d:%02d:%02d", t.hour, t.minute, t.second)
    return cdate, ctime
end

# injectable clock (keeps the writers testable without wall-time)
struct _DT
    year::Int
    month::Int
    day::Int
    hour::Int
    minute::Int
    second::Int
end
const _NOW = Ref(_DT(2026, 1, 1, 0, 0, 0))
function _set_now!()
    c = round(Int, time())
    days = c ÷ 86400
    secs = c % 86400
    # civil-from-days (Howard Hinnant's algorithm)
    z = days + 719468
    era = fld(z, 146097)
    doe = z - era * 146097
    yoe = (doe - doe ÷ 1460 + doe ÷ 36524 - doe ÷ 146096) ÷ 365
    y = yoe + era * 400
    doy = doe - (365 * yoe + yoe ÷ 4 - yoe ÷ 100)
    mp = (5 * doy + 2) ÷ 153
    d = doy - (153 * mp + 2) ÷ 5 + 1
    m = mp < 10 ? mp + 3 : mp - 9
    y = m <= 2 ? y + 1 : y
    _NOW[] = _DT(y, m, d, secs ÷ 3600, (secs % 3600) ÷ 60, secs % 60)
    return nothing
end

"""
    write_rmn(seedname, model, Mrot) -> path

Write `seedname_r.dat`: the position matrix elements ⟨0n|r|Rm⟩ (Å) on the Wigner–Seitz R-set,
from the final-gauge overlaps `Mrot` — linear WYSV06 form off the diagonal, −Σ w_b b Im ln M̃_nn
on the diagonal (all R). Reference record layout `(5I5,6F12.6)`, n varying fastest.
"""
function write_rmn(seedname::AbstractString, model::Model, Mrot::Array{ComplexF64,4})
    kgrid, bv = model.kgrid, model.bvectors
    nw = model.num_wann
    nk = nkpt(kgrid)
    irvec, _ = wigner_seitz(model.lattice, kgrid.mp_grid)
    Aq = zeros(ComplexF64, nw, nw, nk, 3)
    for k in 1:nk, b in 1:bv.nntot
        w = bv.wb[b, k]
        for α in 1:3
            f = w * bv.bvec[α, b, k]
            for m in 1:nw, n in 1:nw
                if n == m
                    Aq[n, n, k, α] -= f * imag(log(Mrot[n, n, b, k]))
                else
                    Aq[n, m, k, α] += im * f * Mrot[n, m, b, k]
                end
            end
        end
    end
    pos = [fourier_q_to_R((@view Aq[:, :, :, α]), kgrid, irvec) for α in 1:3]
    _set_now!()
    cdate, ctime = _w90_datetime()
    open(seedname * "_r.dat", "w") do io
        println(io, " written on ", cdate, " at ", ctime, " ")
        @printf(io, "%12d\n", nw)
        @printf(io, "%12d\n", length(irvec))
        for ir in 1:length(irvec), m in 1:nw, n in 1:nw
            @printf(io, "%5d%5d%5d%5d%5d", irvec[ir][1], irvec[ir][2], irvec[ir][3], n, m)
            for α in 1:3
                @printf(io, "%12.6f%12.6f", real(pos[α][n, m, ir]), imag(pos[α][n, m, ir]))
            end
            println(io)
        end
    end
    return seedname * "_r.dat"
end

"""
    hr_diagonal(Hr, irvec, ndegen) -> Vector{Float64}

The on-site Hamiltonian matrix elements ⟨0n|H|0n⟩ (eV) = real diagonal of H(R = 0) / ndegen.
"""
function hr_diagonal(Hr::Array{ComplexF64,3}, irvec::Vector{NTuple{3,Int}}, ndegen::Vector{Int})
    ir0 = findfirst(==((0, 0, 0)), irvec)
    ir0 === nothing && error("hr_diagonal: no R = 0 vector")
    nw = size(Hr, 1)
    return [real(Hr[i, i, ir0]) / ndegen[ir0] for i in 1:nw]
end

"""
    write_hr_diag(io, Hr, irvec, ndegen)

Write the postw90 `write_hr_diag` on-site Hamiltonian table (`⟨0n|H|0n⟩` in eV) to `io`
(default `stdout`), matching the reference stdout format.
"""
function write_hr_diag(io::IO, Hr::Array{ComplexF64,3}, irvec::Vector{NTuple{3,Int}},
                       ndegen::Vector{Int})
    d = hr_diagonal(Hr, irvec, ndegen)
    println(io)
    println(io, " On-site Hamiltonian matrix elements")
    println(io, "     n        <0n|H|0n> (eV)")
    println(io, "   -------------------------")
    for (i, v) in enumerate(d)
        @printf(io, "   %3d     %12.6f\n", i, v)
    end
    println(io)
    return d
end
write_hr_diag(Hr, irvec, ndegen) = write_hr_diag(stdout, Hr, irvec, ndegen)

"Translate a Cartesian vector into the home unit cell [0,1)³ (reference utility_translate_home)."
function translate_home(vec::AbstractVector, lattice::Lattice)
    f = lattice.B * SVector{3,Float64}(vec...) ./ TWOPI      # Cartesian → fractional
    fr = MVector{3,Float64}(f...)
    for i in 1:3
        fr[i] < 0.0 && (fr[i] += ceil(abs(fr[i])))
        fr[i] > 1.0 && (fr[i] -= trunc(fr[i]))
    end
    return lattice.A * SVector(fr)
end

"""
    write_xyz(path, centres, atoms_cart; translate_home_cell=false, lattice=nothing) -> path

Write `seedname_centres.xyz`: Wannier centres as pseudo-atoms `X`, then the real atoms.
`centres` is a 3×nw matrix (Å), `atoms_cart` a vector of `(symbol, position)` pairs (Å). With
`translate_home_cell`, centres are rationalised into the home cell (needs `lattice`).
"""
function write_xyz(path::AbstractString, centres::AbstractMatrix,
                   atoms_cart::Vector{<:Tuple{<:AbstractString,<:AbstractVector}};
                   translate_home_cell::Bool=false, lattice::Union{Nothing,Lattice}=nothing)
    nw = size(centres, 2)
    wc = translate_home_cell ?
         hcat([translate_home(centres[:, i], lattice) for i in 1:nw]...) : centres
    _set_now!()
    cdate, ctime = _w90_datetime()
    open(path, "w") do io
        @printf(io, "%6d\n", nw + length(atoms_cart))
        println(io, " Wannier centres, written by WannierFunctions.jl on", cdate, " at ", ctime)
        for i in 1:nw
            @printf(io, "X      %14.8f   %14.8f   %14.8f\n", wc[1, i], wc[2, i], wc[3, i])
        end
        for (sym, p) in atoms_cart
            @printf(io, "%-2s     %14.8f   %14.8f   %14.8f\n", sym, p[1], p[2], p[3])
        end
    end
    return path
end

"""
    write_bxsf(seedname, lattice, Hr, irvec, ndegen; fermi_energy=0.0, num_points=50) -> path

Write `seedname.bxsf` (XCrySDen Fermi-surface grid): all interpolated bands on the inclusive
(num_points+1)³ fractional grid (x outer, z fastest), reference layout.
"""
function write_bxsf(seedname::AbstractString, lattice::Lattice, Hr::Array{ComplexF64,3},
                    irvec::Vector{NTuple{3,Int}}, ndegen::Vector{Int};
                    fermi_energy::Float64=0.0, num_points::Int=50)
    nw = size(Hr, 1)
    np = num_points
    _set_now!()
    cdate, ctime = _w90_datetime()
    eigs = Array{Float64,2}(undef, nw, (np + 1)^3)
    kpts = [SVector((lx - 1) / np, (ly - 1) / np, (lz - 1) / np)
            for lx in 1:np+1 for ly in 1:np+1 for lz in 1:np+1]     # z fastest
    Threads.@threads for ik in 1:length(kpts)
        H = zeros(ComplexF64, nw, nw)
        for ir in 1:length(irvec)
            fac = cis(TWOPI * dot(kpts[ik], SVector{3,Float64}(irvec[ir]...))) / ndegen[ir]
            @views H .+= fac .* Hr[:, :, ir]
        end
        eigs[:, ik] = eigvals(Hermitian(H))
    end
    B = lattice.B
    open(seedname * ".bxsf", "w") do io
        println(io, "  BEGIN_INFO")
        println(io, "       #")
        println(io, "       # this is a Band-XCRYSDEN-Structure-File")
        println(io, "       # for Fermi Surface Visualisation")
        println(io, "       #")
        println(io, "       # Generated by the Wannier90 code http://www.wannier.org")
        println(io, "       # On ", cdate, "  at ", ctime)
        println(io, "       #")
        @printf(io, "       Fermi Energy:   %.15f     \n", fermi_energy)
        println(io, "  END_INFO")
        println(io)
        println(io, "  BEGIN_BLOCK_BANDGRID_3D")
        println(io, " from_wannier_code")
        println(io, "  BEGIN_BANDGRID_3D_fermi")
        @printf(io, "%12d\n", nw)
        @printf(io, "%12d%12d%12d\n", np + 1, np + 1, np + 1)
        println(io, " 0.0 0.0 0.0")
        for i in 1:3
            @printf(io, "%21.9f%21.9f%21.9f\n", B[1, i], B[2, i], B[3, i])
        end
        for n in 1:nw
            @printf(io, " BAND: %12d\n", n)
            for ik in 1:(np+1)^3
                println(io, fortran_e(eigs[n, ik], 16, 8))
            end
        end
        println(io, " END_BANDGRID_3D")
        println(io, "  END_BLOCK_BANDGRID_3D")
    end
    return seedname * ".bxsf"
end

"""
    tabulate_3d(bm; mesh, colour=nothing) -> (energies, colours)

Tabulate interpolated band energies on a regular (Γ-inclusive, no boundary doubling) 3-D grid
`k = ((i1-1)/n1, (i2-1)/n2, (i3-1)/n3)`. `energies` is `nband × n1 × n2 × n3`. If `colour` is a
function `(bm, kf, E, U) -> Vector{Float64}` (one scalar per band), its values are tabulated in
the same layout and returned as `colours` (else `nothing`). Useful for FermiSurfer surfaces
coloured by velocity, spin, or Berry curvature.
"""
function tabulate_3d(bm::BerryModel; mesh::NTuple{3,Int}, colour=nothing)
    nw = num_wann(bm)
    n1, n2, n3 = mesh
    E = Array{Float64,4}(undef, nw, n1, n2, n3)
    C = colour === nothing ? nothing : Array{Float64,4}(undef, nw, n1, n2, n3)
    kl = [(i1, i2, i3) for i1 in 1:n1 for i2 in 1:n2 for i3 in 1:n3]
    Threads.@threads for idx in 1:length(kl)
        i1, i2, i3 = kl[idx]
        kf = SVector((i1 - 1) / n1, (i2 - 1) / n2, (i3 - 1) / n3)
        Ek, _, U = eig_deleig_vec(bm, kf; deriv=false)
        E[:, i1, i2, i3] = Ek
        colour !== nothing && (C[:, i1, i2, i3] = colour(bm, kf, Ek, U))
    end
    return E, C
end

"""
    write_frmsf(path, lattice, energies; fermi_energy=0.0, colours=nothing) -> path

Write a FermiSurfer `.frmsf` file: grid dims, shift flag 1, band count, the three reciprocal
vectors (Å⁻¹), then band-outer / i3-fastest energies (shifted by `fermi_energy` so the surface
lands at 0), and optionally a colour block in the same order. `energies`/`colours` are the
`nband × n1 × n2 × n3` arrays from [`tabulate_3d`](@ref).
"""
function write_frmsf(path::AbstractString, lattice::Lattice, energies::Array{Float64,4};
                     fermi_energy::Float64=0.0, colours::Union{Nothing,Array{Float64,4}}=nothing)
    nb, n1, n2, n3 = size(energies)
    B = lattice.B
    open(path, "w") do io
        @printf(io, "%d %d %d\n", n1, n2, n3)
        println(io, "1")
        @printf(io, "%d\n", nb)
        for i in 1:3
            @printf(io, "  %14.8f  %14.8f  %14.8f\n", B[1, i], B[2, i], B[3, i])
        end
        # band outer, i3 fastest
        for b in 1:nb, i1 in 1:n1, i2 in 1:n2, i3 in 1:n3
            @printf(io, "%15.8e\n", energies[b, i1, i2, i3] - fermi_energy)
        end
        if colours !== nothing
            for b in 1:nb, i1 in 1:n1, i2 in 1:n2, i3 in 1:n3
                @printf(io, "%15.8e\n", colours[b, i1, i2, i3])
            end
        end
    end
    return path
end

"""
    write_cube(seedname, index, lattice, ng, w, los, centre, atoms_cart;
               radius=3.5, scale=1.0, mode=:crystal, supercell=(2,2,2)) -> path

Write `seedname_<index>.cube` for one Wannier function: `w` is its complex supercell grid
(from [`wannier_function_grid`](@ref), offsets `los`, UNK grid `ng`), `centre` its Cartesian
centre (Å), `atoms_cart` a vector of `(symbol, position)` pairs (Å). Crystal mode keeps every
periodic atom image within `scale·radius` of the WF centre; the data box spans `radius` around
the centre along each lattice direction (reference geometry, everything written in bohr).
"""
function write_cube(seedname::AbstractString, index::Int, lattice::Lattice,
                    ng::NTuple{3,Int}, w::AbstractArray{ComplexF64,3},
                    los::NTuple{3,Int}, centre::AbstractVector,
                    atoms_cart::Vector{<:Tuple{<:AbstractString,<:AbstractVector}};
                    radius::Float64=3.5, scale::Float64=1.0, mode::Symbol=:crystal,
                    supercell::NTuple{3,Int}=(2, 2, 2))
    A, B = lattice.A, lattice.B
    moda = SVector(norm(A[:, 1]), norm(A[:, 2]), norm(A[:, 3]))
    modb = SVector(norm(B[:, 1]), norm(B[:, 2]), norm(B[:, 3]))
    dgrid = SVector(moda[1] / ng[1], moda[2] / ng[2], moda[3] / ng[3])

    proj = SVector((dot(centre, B[:, i]) * moda[i] / TWOPI for i in 1:3)...)
    rstart = SVector((proj[i] - TWOPI * radius / (moda[i] * modb[i]) for i in 1:3)...)
    rend = SVector((proj[i] + TWOPI * radius / (moda[i] * modb[i]) for i in 1:3)...)
    ilength = SVector{3,Int}((ceil(Int, (rend[i] - rstart[i]) / dgrid[i]) for i in 1:3)...)
    istart = SVector{3,Int}((floor(Int, rstart[i] / dgrid[i]) + 1 for i in 1:3)...)
    orig = A * SVector(((istart[i] - 1) / ng[i] for i in 1:3)...)

    his = (los[1] + size(w, 1) - 1, los[2] + size(w, 2) - 1, los[3] + size(w, 3) - 1)
    # fold-up-from-below (reference behaviour; above the supercell is an error)
    function q(nn::Int, d::Int)
        qq = nn + istart[d] - 1
        if qq < los[d]
            qq += (div(abs(qq) - 1, ng[d])) * ng[d]
        end
        (qq < los[d] || qq > his[d]) &&
            error("write_cube: box outside the plot supercell — increase wannier_plot_supercell or decrease radius")
        return qq - los[d] + 1
    end

    kept = Tuple{Int,SVector{3,Float64}}[]
    if mode == :crystal
        for (sym, pc) in atoms_cart
            Z = _atomic_number(sym)
            for nx in -(supercell[1] ÷ 2):((supercell[1] + 1) ÷ 2),
                ny in -(supercell[2] ÷ 2):((supercell[2] + 1) ÷ 2),
                nz in -(supercell[3] ÷ 2):((supercell[3] + 1) ÷ 2)

                pos = SVector{3,Float64}(pc...) + nx * A[:, 1] + ny * A[:, 2] + nz * A[:, 3]
                norm(pos - centre) <= scale * radius && push!(kept, (Z, pos))
            end
        end
    else
        for (sym, pc) in atoms_cart
            push!(kept, (_atomic_number(sym), SVector{3,Float64}(pc...)))
        end
    end

    _set_now!()
    cdate, ctime = _w90_datetime()
    path = @sprintf("%s_%05d.cube", seedname, index)
    open(path, "w") do io
        println(io, "      Generated by Wannier90 code http://www.wannier.org")
        println(io, "      On ", cdate, " at ", ctime, " ")
        @printf(io, "%4d%13.5f%13.5f%13.5f\n", length(kept), orig[1] / BOHR, orig[2] / BOHR,
                orig[3] / BOHR)
        for i in 1:3
            @printf(io, "%4d%13.5f%13.5f%13.5f\n", ilength[i],
                    A[1, i] / (ng[i] * BOHR), A[2, i] / (ng[i] * BOHR), A[3, i] / (ng[i] * BOHR))
        end
        for (Z, pos) in kept
            @printf(io, "%4d%13.5f%13.5f%13.5f%13.5f\n", Z, 1.0,
                    pos[1] / BOHR, pos[2] / BOHR, pos[3] / BOHR)
        end
        for nx in 1:ilength[1], ny in 1:ilength[2]
            iz = 1
            while iz <= ilength[3]
                hi = min(iz + 5, ilength[3])
                for nz in iz:hi
                    print(io, fortran_e(real(w[q(nx, 1), q(ny, 2), q(nz, 3)]), 13, 5))
                end
                println(io)
                iz = hi + 1
            end
        end
    end
    return path
end
