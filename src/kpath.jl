# kpath — quantities along a high-symmetry path (postw90's kpath module): interpolated bands
# (optionally coloured by ⟨σ·n̂⟩ or by the band-resolved spin Hall term), the occupied-manifold
# Berry curvature, the orbital-magnetisation integrand, and the k-resolved spin Hall
# conductivity. Exact conventions in docs/reference-notes/kpath-module.md.

using Printf
using StaticArrays

"""
    kpath_points(segments, lattice; num_points=100) -> (kpts, xvals)

Reference path construction: segment 1 gets `num_points` samples, segment p gets
`nint(num_points·len_p/len_1)`; each segment is sampled at fractions (i−1)/n_p (interior
vertices appear once, as the next segment's first point) and one final point closes the path.
`xvals` reproduce postw90's accumulation: the increment is always the *current* segment's step,
and only the final value is forced to the exact total length (so interior-vertex x-coordinates
differ slightly from the cumulative segment lengths — this matches the reference files).
"""
function kpath_points(segments::Vector{<:Tuple{<:AbstractVector,<:AbstractVector}},
                      lattice::Lattice; num_points::Int=100)
    segs = [(SVector{3,Float64}(a...), SVector{3,Float64}(b...)) for (a, b) in segments]
    seglen = [norm(lattice.B * (b - a)) for (a, b) in segs]
    npts = [p == 1 ? num_points :
            round(Int, num_points * seglen[p] / seglen[1], RoundNearestTiesAway)
            for p in 1:length(segs)]
    kpts = SVector{3,Float64}[]
    xvals = Float64[]
    x = 0.0
    for (p, (a, b)) in enumerate(segs)
        for i in 1:npts[p]
            isempty(xvals) || (x += seglen[p] / npts[p])   # crossing uses the NEW segment's step
            f = (i - 1) / npts[p]
            # fma mirrors the reference binary's contracted a+(b−a)·f (sign-of-zero in -path.kpt)
            push!(kpts, SVector(fma(b[1] - a[1], f, a[1]), fma(b[2] - a[2], f, a[2]),
                                fma(b[3] - a[3], f, a[3])))
            push!(xvals, x)
        end
    end
    push!(kpts, segs[end][2])                 # forced final point at the exact total length
    push!(xvals, sum(seglen))
    return kpts, xvals
end

"Parse the `kpoint_path` block of a `.win` into (start, end) fractional segment tuples."
function kpath_segments(win::WinInput)
    haskey(win.blocks, "kpoint_path") ||
        error("kpath: the .win has no kpoint_path block")
    segs = Tuple{SVector{3,Float64},SVector{3,Float64}}[]
    for ln in win.blocks["kpoint_path"]
        t = split(ln)
        length(t) >= 8 || continue
        push!(segs, (SVector{3,Float64}(parse_f64.(t[2:4])...),
                     SVector{3,Float64}(parse_f64.(t[6:8])...)))
    end
    return segs
end

"""
    kpath(m; segments, num_points=100, tasks=(:bands,), bands_colour=:none, fermi_energy=nothing,
          curv_unit=:ang2, polar=0.0, azimuth=0.0, spin=nothing,
          γ=3, α=1, β=2, smr_width=0.0, eigval_max=Inf) -> NamedTuple

Quantities along a high-symmetry path. `m` is the model matching the most demanding task:
a `BerryModel` (bands, curv), `MorbModel` (morb + curv), or `ShcModel` (shc / shc colouring);
`spin` supplies a `SpinModel` for `bands_colour = :spin`. Returns
`(; kpts, xvals, bands, colour, curv, morb, shc)` with postw90's conventions: `curv` is the
**negative** Berry curvature (Ų, or bohr² with `curv_unit = :bohr2`), `morb` the LVTS12
integrand −(G + H − 2E_F·F)/2 in eV·Ų (never unit-converted), `shc` the k-resolved
Ω^{spin γ}_{αβ} Fermi sum. SHC quantities use fixed smearing `smr_width` (as the reference
requires along a path).
"""
function kpath(m::Union{BerryModel,MorbModel,ShcModel};
               segments::Vector{<:Tuple}, num_points::Int=100,
               tasks=(:bands,), bands_colour::Symbol=:none,
               fermi_energy::Union{Nothing,Float64}=nothing,
               curv_unit::Symbol=:ang2, polar::Float64=0.0, azimuth::Float64=0.0,
               spin::Union{Nothing,SpinModel}=nothing,
               γ::Int=3, α::Int=1, β::Int=2, smr_width::Float64=0.0,
               eigval_max::Float64=Inf)
    bm = m isa BerryModel ? m : m.bm
    all(t -> t in (:bands, :curv, :morb, :shc), tasks) ||
        error("kpath: tasks must be a subset of (:bands, :curv, :morb, :shc)")
    (:curv in tasks || :morb in tasks || :shc in tasks || bands_colour != :none) &&
        fermi_energy === nothing && error("kpath: this task needs fermi_energy")
    :morb in tasks && !(m isa MorbModel) && error("kpath: morb needs a MorbModel")
    (:shc in tasks || bands_colour == :shc) && !(m isa ShcModel) &&
        error("kpath: shc needs an ShcModel")
    bands_colour == :spin && spin === nothing && error("kpath: spin colouring needs spin=SpinModel")

    kpts, xvals = kpath_points(segments, bm.lattice; num_points=num_points)
    np = length(kpts)
    nw = num_wann(bm)
    ucnv = curv_unit == :bohr2 ? 1.0 / BOHR^2 : 1.0

    bands = :bands in tasks ? Matrix{Float64}(undef, nw, np) : nothing
    colour = bands_colour != :none ? Matrix{Float64}(undef, nw, np) : nothing
    curv = :curv in tasks ? Matrix{Float64}(undef, 3, np) : nothing
    morb = :morb in tasks ? Matrix{Float64}(undef, 3, np) : nothing
    shc = :shc in tasks ? Vector{Float64}(undef, np) : nothing

    ef = fermi_energy === nothing ? 0.0 : fermi_energy
    threaded_ksum(
        (st, p) -> _kpath_kpoint!(bands, colour, curv, morb, shc, p, m, bm, kpts[p],
                                  bands_colour, ef, ucnv, polar, azimuth, spin, γ, α, β,
                                  smr_width, eigval_max, st.work),
        () -> (work=BerryKWork(nw),),
        np)
    return (; kpts, xvals, bands, colour, curv, morb, shc)
end

# Per-k body kept as a top-level function: large inlined @threads bodies with many captured
# locals are a closure-boxing hazard on Julia 1.12 (see orbital_magnetisation).
function _kpath_kpoint!(bands, colour, curv, morb, shc, p::Int, m, bm::BerryModel,
                        kf::SVector{3,Float64}, bands_colour::Symbol, ef::Float64,
                        ucnv::Float64, polar::Float64, azimuth::Float64, spin,
                        γ::Int, α::Int, β::Int, smr_width::Float64, eigval_max::Float64,
                        w::BerryKWork=BerryKWork(num_wann(bm)))
    if bands !== nothing
        E, _, U = eig_deleig_vec!(w, bm, kf; deriv=false)
        bands[:, p] = E
        if bands_colour == :spin
            spn = _spin_diag(spin, kf, U; polar=polar, azimuth=azimuth)
            colour[:, p] = clamp.(spn, -1.0 + 1e-8, 1.0 - 1e-8)
        end
    end
    if morb !== nothing
        f, g, h = _imfgh_kdata(m, kf, ef, w)
        morb[:, p] = -(g .+ h .- 2.0 * ef .* f) ./ 2.0
        curv !== nothing && (curv[:, p] = -ucnv .* f)
    elseif curv !== nothing
        curv[:, p] = -ucnv .* _imf_kdata!(w, _berry_kdata!(w, bm, kf), ef)
    end
    if shc !== nothing || bands_colour == :shc
        E, ω = _shc_k_band(m, kf, 0.0, γ, α, β, false, 0.0, 0.0, smr_width, eigval_max, w)
        shc !== nothing &&
            (shc[p] = sum(Float64(E[n] < ef) * ω[n] for n in 1:length(E)) * ucnv)
        bands_colour == :shc && (colour[:, p] = ω .* ucnv)
    end
    return nothing
end

"""
    write_kpath(seedname, res) -> seedname

Write the postw90 kpath output files for whichever quantities `res` (from [`kpath`](@ref))
contains: `-path.kpt` + `-bands.dat` (band-major, colour column when present), `-curv.dat`,
`-morb.dat`, and `-shc.dat`.
"""
function write_kpath(seedname::AbstractString, res)
    if res.bands !== nothing
        open(seedname * "-path.kpt", "w") do io
            @printf(io, "%12d\n", length(res.kpts))
            for k in res.kpts
                @printf(io, "%12.6f%12.6f%12.6f   %4.1f\n", k[1], k[2], k[3], 1.0)
            end
        end
        open(seedname * "-bands.dat", "w") do io
            nw, np = size(res.bands)
            for n in 1:nw
                for p in 1:np
                    print(io, fortran_e(res.xvals[p], 16, 8), fortran_e(res.bands[n, p], 16, 8))
                    res.colour !== nothing && print(io, fortran_e(res.colour[n, p], 16, 8))
                    println(io)
                end
                println(io, "  ")
            end
        end
    end
    for (tag, dat) in (("curv", res.curv), ("morb", res.morb))
        dat === nothing && continue
        open(seedname * "-" * tag * ".dat", "w") do io
            for p in 1:size(dat, 2)
                println(io, fortran_e(res.xvals[p], 16, 8), fortran_e(dat[1, p], 16, 8),
                        fortran_e(dat[2, p], 16, 8), fortran_e(dat[3, p], 16, 8))
            end
            println(io, "  ")
        end
    end
    if res.shc !== nothing
        open(seedname * "-shc.dat", "w") do io
            for p in 1:length(res.shc)
                println(io, fortran_e(res.xvals[p], 16, 8), fortran_e(res.shc[p], 16, 8))
            end
            println(io, "  ")
        end
    end
    return seedname
end
