"""
    Wannier90

A modern Julia reimplementation of the Wannier90 core (`wannier90.x`): construction of
maximally-localised Wannier functions and Wannier interpolation, with drop-in compatibility for
the standard `.win` / `.amn` / `.mmn` / `.eig` file formats.

This is an independent, from-scratch reimplementation; it is not affiliated with the official
Wannier90 project.
"""
module WannierFunctions

using LinearAlgebra
using StaticArrays

include("constants.jl")
include("types.jl")
include("known_keywords.jl")
include("io.jl")
include("kmesh.jl")
include("gauge.jl")
include("spread.jl")
include("wannierise.jl")
include("gamma.jl")
include("interpolate.jl")
include("disentangle.jl")
include("pipeline.jl")
include("operator.jl")
include("output.jl")
include("nnkp.jl")
include("chk.jl")
include("plot.jl")
include("scdm.jl")
include("berry.jl")
include("kubo.jl")
include("spin.jl")
include("dos.jl")
include("boltzwann.jl")
include("shc.jl")
include("tetrahedron_kernels.jl")
include("tetrahedron.jl")
include("morb.jl")
include("gyrotropic.jl")
include("kslice.jl")
include("geninterp.jl")
include("kdotp.jl")
include("shiftcurrent.jl")
include("kpath.jl")
include("extras.jl")
include("tbmodel.jl")
include("transport.jl")
include("postw90.jl")
include("symmetry.jl")
include("injection.jl")
include("dftk.jl")
include("slwf.jl")
include("ss.jl")
include("sitesym.jl")
include("wout.jl")
include("cli.jl")
include("show.jl")

export Lattice, KGrid, BVectors, Model
export read_win, read_amn, read_mmn, read_eig
export read_model, initial_gauge, rotate_overlaps, compute_spread, SpreadResult
export wannierise, localize, WannieriseResult
export disentangle, DisentangleResult, WindowData
export run_wannier, WannierResult, interpolate
export wigner_seitz, build_hr, interpolate_hk, interpolate_bands
export ws_translate_dist, interpolate_bands_ws
export write_hr, read_hr, write_tb, write_band_dat, write_band_kpt, write_labelinfo
export write_wout, generate_kpath, main, install_cli
export generate_nnkp, write_nnkp, parse_projections, parse_exclude_bands, Projection
export Checkpoint, read_chk, write_chk, read_chk_fmt, write_chk_fmt
export read_unk, plot_wannier_functions, write_xsf, parse_range_list, wannier_function_grid
export scdm_projections, scdm_amn, write_amn
export BerryModel, berry_curvature_k, anomalous_hall, ahc_fermiscan
export geninterp, eig_deleig, read_geninterp_kpt
export optical_conductivity, KuboResult, kubo_S, kubo_A
export density_of_states, write_dos
export SpinModel, spin_moment, spin_expectation
export kpath, kpath_points, kpath_segments, write_kpath
export gyrotropic, write_gyrotropic
export kdotp, write_kdotp
export shift_current, write_shift_current
export write_rmn, write_bxsf, write_cube, parse_atoms
export hr_diagonal, write_hr_diag, write_xyz, translate_home
export tabulate_3d, write_frmsf
export read_tb, tb_model
export transport_bulk, transport_from_tb, tran_transfer, tran_green, read_ht, write_ht,
       write_transport, run_transport, translate_centres_home
export postw90_main, fortran_g, write_kubo, write_fermiscan, write_boltzwann, write_boltzdos
export SymmetryOps, read_sym, nsym, irreducible_kmesh, cubic_point_group, anomalous_hall_sym,
       orbital_magnetisation_sym, density_of_states_sym
export SLWF, slwf_omega, slwf_gradient
export injection_current, wannier_model
export Sitesym, read_dmn, symmetrize_u!, symmetrize_gradient!
export boltzwann, BoltzWannResult
export ShcModel, ShcRyooModel, shc_fermiscan, shc_freqscan, read_spn, read_shu, write_shc
export shc_tetra, shc_imjv
export kslice, write_kslice
export MorbModel, orbital_magnetisation, read_uhu
export TBOperator, hamiltonian_operator, position_operator, bands, fourier_to_R

"""
    read_model(seedname) -> Model

Load a full wannierisation problem from `seedname.win`, `seedname.amn`, `seedname.mmn`, and (if
present) `seedname.eig`. Builds the lattice, k-grid, and B1 neighbour weights.
"""
function read_model(seedname::AbstractString)
    win = read_win(seedname * ".win")
    A, nb_a, nk_a, nw_a = read_amn(seedname * ".amn")
    M, kpb, gpb, nb_m, nk_m, nntot = read_mmn(seedname * ".mmn")

    # select_projections: keep only the named .amn columns as the WFs (num_proj > num_wann).
    if haskey(win.raw, "select_projections") && nw_a > win.num_wann
        sel = parse_range_list(replace(strip(win.raw["select_projections"]), r"\s+" => ","))
        length(sel) == win.num_wann ||
            error("select_projections picks $(length(sel)) columns but num_wann is $(win.num_wann)")
        all(c -> 1 <= c <= nw_a, sel) ||
            error("select_projections indices out of range 1:$nw_a")
        A = A[:, sel, :]
        nw_a = win.num_wann
    end

    nb_a == nb_m || error("num_bands mismatch between .amn ($nb_a) and .mmn ($nb_m)")
    nk_a == nk_m || error("num_kpts mismatch between .amn ($nk_a) and .mmn ($nk_m)")
    nw_a == win.num_wann || error("num_wann mismatch between .amn ($nw_a) and .win ($(win.num_wann))")
    nk_a == length(win.kpoints) || error("num_kpts mismatch between .amn ($nk_a) and .win kpoints ($(length(win.kpoints)))")

    # Γ-only inputs store only half the b-vectors (the −b partners are implied by
    # M(−b) = M(b)†, exactly — at Γ, u at the k+b image is e^{−ib·r}u). Expand to the closed
    # full set so the spread AND its gradient are exact with the general machinery; the B1
    # solve then recovers the standard full-set weights automatically. Output paths that must
    # match the reference file convention (.chk, .nnkp) write the first (file) half back.
    if win.gamma_only
        nk_m == 1 || error("gamma_only requires a single k-point (got $nk_m)")
        M2 = Array{ComplexF64,4}(undef, nb_m, nb_m, 2 * nntot, 1)
        kpb2 = ones(Int, 2 * nntot, 1)
        gpb2 = Array{Int,3}(undef, 3, 2 * nntot, 1)
        for b in 1:nntot
            M2[:, :, b, 1] = M[:, :, b, 1]
            M2[:, :, nntot+b, 1] = (@view M[:, :, b, 1])'
            gpb2[:, b, 1] = gpb[:, b, 1]
            gpb2[:, nntot+b, 1] = .-gpb[:, b, 1]
        end
        M, kpb, gpb, nntot = M2, kpb2, gpb2, 2 * nntot
    end

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
