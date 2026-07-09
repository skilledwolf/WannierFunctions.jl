# kslice — quantities on a 2-D slice through the BZ (postw90's kslice module): band energies
# (for Fermi-line plots), the occupied-manifold Berry curvature, the orbital-magnetisation
# integrand, and the k-resolved spin Hall conductivity at each point of an (N1+1)×(N2+1)
# inclusive grid spanned by two fractional vectors from a corner.

using Printf
using StaticArrays

"""
    kslice(m; corner, b1, b2, mesh=(50,50), fermi_energy=nothing, tasks=(:bands,),
           curv_unit=:ang2, γ=3, α=1, β=2, smr_width=0.0, eigval_max=Inf, curvature=false)
        -> (; kpts, coords, bands, curv, morb, shc)

Grid k(i,j) = corner + (i/N1)·b1 + (j/N2)·b2 (fractional; i = 0..N1 fastest). `coords` are the
2-D Cartesian coordinates in the slice plane (Å⁻¹). `m` is a `BerryModel` (bands, curv),
`MorbModel` (morb), or `ShcModel` (shc). Conventions match postw90's kslice files: `curv` is
**−**(J0+J1+J2), `morb` the LVTS12 integrand −(G + H − 2E_F·F)/2 in eV·Ų, `shc` the
fermi-summed Ω^{spin γ}_{αβ} with fixed smearing `smr_width`. `curvature=true` is a legacy
alias for adding `:curv` to `tasks`.
"""
function kslice(m::Union{BerryModel,MorbModel,ShcModel};
                corner::AbstractVector=[0.0, 0.0, 0.0],
                b1::AbstractVector, b2::AbstractVector, mesh::NTuple{2,Int}=(50, 50),
                fermi_energy::Union{Nothing,Float64}=nothing, tasks=(:bands,),
                curv_unit::Symbol=:ang2, γ::Int=3, α::Int=1, β::Int=2,
                smr_width::Float64=0.0, eigval_max::Float64=Inf, curvature::Bool=false)
    tasks = curvature && !(:curv in tasks) ? (tasks..., :curv) : tasks
    all(t -> t in (:bands, :curv, :morb, :shc), tasks) ||
        error("kslice: tasks must be a subset of (:bands, :curv, :morb, :shc)")
    (:curv in tasks || :morb in tasks || :shc in tasks) && fermi_energy === nothing &&
        error("kslice: this task needs fermi_energy")
    :morb in tasks && !(m isa MorbModel) && error("kslice: morb needs a MorbModel")
    :shc in tasks && !(m isa ShcModel) && error("kslice: shc needs an ShcModel")

    bm = m isa BerryModel ? m : m.bm
    c = SVector{3,Float64}(corner...)
    v1 = SVector{3,Float64}(b1...)
    v2 = SVector{3,Float64}(b2...)
    N1, N2 = mesh
    kpts = [c + (i / N1) * v1 + (j / N2) * v2 for j in 0:N2 for i in 0:N1]   # i fastest

    # orthonormal in-plane basis for the 2-D coordinates
    e1c = bm.lattice.B * v1
    e2c = bm.lattice.B * v2
    ê1 = e1c / norm(e1c)
    e2p = e2c - dot(e2c, ê1) * ê1
    ê2 = e2p / norm(e2p)
    coords = [(dot(bm.lattice.B * (k - c), ê1), dot(bm.lattice.B * (k - c), ê2)) for k in kpts]

    np = length(kpts)
    nw = num_wann(bm)
    ucnv = curv_unit == :bohr2 ? 1.0 / BOHR^2 : 1.0
    bands = :bands in tasks ? Matrix{Float64}(undef, nw, np) : nothing
    curv = :curv in tasks ? Matrix{Float64}(undef, 3, np) : nothing
    morb = :morb in tasks ? Matrix{Float64}(undef, 3, np) : nothing
    shc = :shc in tasks ? Vector{Float64}(undef, np) : nothing
    ef = fermi_energy === nothing ? 0.0 : fermi_energy
    threaded_ksum(
        (st, p) -> _kslice_kpoint!(bands, curv, morb, shc, p, m, bm, kpts[p], ef, ucnv,
                                   γ, α, β, smr_width, eigval_max, st.work),
        () -> (work=BerryKWork(nw),),
        np)
    return (; kpts, coords, bands, curv, morb, shc)
end

# Top-level per-k body (closure-boxing hazard in large @threads bodies on Julia 1.12).
function _kslice_kpoint!(bands, curv, morb, shc, p::Int, m, bm::BerryModel,
                         kf::SVector{3,Float64}, ef::Float64, ucnv::Float64,
                         γ::Int, α::Int, β::Int, smr_width::Float64, eigval_max::Float64,
                         w::BerryKWork=BerryKWork(num_wann(bm)))
    if bands !== nothing
        E, _, _ = eig_deleig_vec!(w, bm, kf; deriv=false)
        bands[:, p] = E
    end
    if morb !== nothing
        f, g, h = _imfgh_kdata(m, kf, ef, w)
        morb[:, p] = -(g .+ h .- 2.0 * ef .* f) ./ 2.0
        curv !== nothing && (curv[:, p] = -ucnv .* f)
    elseif curv !== nothing
        curv[:, p] = -ucnv .* _imf_kdata!(w, _berry_kdata!(w, bm, kf), ef)
    end
    if shc !== nothing
        E, ω = _shc_k_band(m, kf, 0.0, γ, α, β, false, 0.0, 0.0, smr_width, eigval_max, w)
        shc[p] = sum(Float64(E[n] < ef) * ω[n] for n in 1:length(E)) * ucnv
    end
    return nothing
end

"""
    write_kslice(seedname, coords, bands; curv=nothing, morb=nothing, shc=nothing)

Write the postw90 kslice output files: `-kslice-coord.dat`, `-kslice-bands.dat`, and (when
given) `-kslice-curv.dat`, `-kslice-morb.dat`, `-kslice-shc.dat`.
"""
function write_kslice(seedname::AbstractString, coords, bands; curv=nothing, morb=nothing,
                      shc=nothing)
    open(seedname * "-kslice-coord.dat", "w") do io
        for (x, y) in coords
            println(io, fortran_e(x, 16, 8), fortran_e(y, 16, 8))
        end
        println(io, " ")
    end
    if bands !== nothing
        open(seedname * "-kslice-bands.dat", "w") do io
            for p in 1:size(bands, 2), n in 1:size(bands, 1)
                println(io, fortran_e(bands[n, p], 16, 8))
            end
            println(io, " ")
        end
    end
    for (tag, dat) in (("curv", curv), ("morb", morb))
        dat === nothing && continue
        open(seedname * "-kslice-" * tag * ".dat", "w") do io
            for p in 1:size(dat, 2)
                println(io, fortran_e(dat[1, p], 16, 8), fortran_e(dat[2, p], 16, 8),
                        fortran_e(dat[3, p], 16, 8))
            end
            println(io, " ")
        end
    end
    if shc !== nothing
        open(seedname * "-kslice-shc.dat", "w") do io
            for p in 1:length(shc)
                println(io, fortran_e(shc[p], 16, 8))
            end
            println(io, " ")
        end
    end
    return seedname
end
