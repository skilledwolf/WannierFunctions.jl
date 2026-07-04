"""
    Wannier90

A modern Julia reimplementation of the Wannier90 core (`wannier90.x`): construction of
maximally-localised Wannier functions and Wannier interpolation, with drop-in compatibility for
the standard `.win` / `.amn` / `.mmn` / `.eig` file formats.

This is an independent, from-scratch reimplementation; it is not affiliated with the official
Wannier90 project.
"""
module Wannier90

using LinearAlgebra
using StaticArrays

include("constants.jl")
include("types.jl")
include("io.jl")
include("kmesh.jl")
include("gauge.jl")
include("spread.jl")

export Lattice, KGrid, BVectors, Model
export read_win, read_amn, read_mmn, read_eig
export read_model, initial_gauge, rotate_overlaps, compute_spread, SpreadResult

"""
    read_model(seedname) -> Model

Load a full wannierisation problem from `seedname.win`, `seedname.amn`, `seedname.mmn`, and (if
present) `seedname.eig`. Builds the lattice, k-grid, and B1 neighbour weights.
"""
function read_model(seedname::AbstractString)
    win = read_win(seedname * ".win")
    A, nb_a, nk_a, nw_a = read_amn(seedname * ".amn")
    M, kpb, gpb, nb_m, nk_m, nntot = read_mmn(seedname * ".mmn")

    nb_a == nb_m || error("num_bands mismatch between .amn ($nb_a) and .mmn ($nb_m)")
    nk_a == nk_m || error("num_kpts mismatch between .amn ($nk_a) and .mmn ($nk_m)")
    nw_a == win.num_wann || error("num_wann mismatch between .amn ($nw_a) and .win ($(win.num_wann))")
    nk_a == length(win.kpoints) || error("num_kpts mismatch between .amn ($nk_a) and .win kpoints ($(length(win.kpoints)))")

    lattice = Lattice(win.unit_cell)
    kgrid = KGrid(win.kpoints, win.mp_grid)
    bvectors = build_bvectors(kgrid, lattice, kpb, gpb; kmesh_tol=win.kmesh_tol)

    eig = isfile(seedname * ".eig") ? read_eig(seedname * ".eig") : nothing

    return Model(lattice, kgrid, bvectors, nb_m, win.num_wann, M, A, eig,
                 String(basename(seedname)))
end

"""
    initial_spread(model) -> SpreadResult

Wannier centres and spread of the initial (Löwdin-projected) gauge — the "Initial State" that
Wannier90 reports before any localisation iteration.
"""
function initial_spread(model::Model)
    U = initial_gauge(model.A)
    Mrot = rotate_overlaps(model.M, U, model.bvectors.kpb)
    return compute_spread(Mrot, model.bvectors)
end

export initial_spread

end # module
