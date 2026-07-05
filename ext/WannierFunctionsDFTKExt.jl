# Package extension: build a WannierFunctions `Model` from a DFTK.jl SCF result, so an
# all-Julia DFT → Wannier pipeline runs in memory with no file round-trip or external binary.
# Loaded automatically when both WannierFunctions and DFTK are present.
#
# The neighbour b-vector list comes from WannierFunctions' own kmesh search (the `-pp` code
# path, byte-compatible with wannier90), and the overlap/projection matrix elements from
# DFTK's wannier_shared primitives (`overlap_Mmn_k_kpb`, `compute_amn_kpoint`).
module WannierFunctionsDFTKExt

using WannierFunctions
const WF = WannierFunctions
using DFTK
using StaticArrays

"""
    WannierFunctions.wannier_model(scfres, projections; num_wann, num_bands, kmesh_tol,
                                   search_shells) -> Model

Assemble a wannierisation `Model` from a DFTK `scfres` and a vector of trial `projections`
(e.g. `DFTK.GaussianWannierProjection`). The SCF must have been run on the full
(symmetry-unreduced) Monkhorst–Pack grid. `num_bands` defaults to `num_wann` (isolated
manifold); energies are converted to eV.
"""
function WannierFunctions.wannier_model(scfres::NamedTuple, projections::AbstractVector;
                                        num_wann::Int, num_bands::Int=num_wann,
                                        kmesh_tol::Float64=1e-6, search_shells::Int=36)
    basis = scfres.basis
    dftkmodel = basis.model
    dftkmodel.n_spin_components == 1 ||
        error("wannier_model: only spinless DFTK models are supported")
    unit_cell = Array(dftkmodel.lattice) .* WF.BOHR       # DFTK lattice is in Bohr
    kpoints = [SVector{3,Float64}(kpt.coordinate...) for kpt in basis.kpoints]
    nk = length(kpoints)
    mp_grid = Tuple(Int.(basis.kgrid.kgrid_size))
    prod(mp_grid) == nk ||
        error("wannier_model: the SCF k-grid is symmetry-reduced ($nk of $(prod(mp_grid)) " *
              "k-points); rerun DFTK with symmetries = false")

    eig = zeros(num_bands, nk)
    for ik in 1:nk
        eig[:, ik] = [DFTK.auconvert(DFTK.Unitful.eV, ε).val
                      for ε in scfres.eigenvalues[ik][1:num_bands]]
    end

    # b-vector neighbour list from our own kmesh search (identical to the .nnkp convention)
    lattice = WF.Lattice(SMatrix{3,3,Float64}(unit_cell))
    cells = WF.supercell_cells(lattice)
    kcart = [lattice.B * k for k in kpoints]
    dnn, multi = WF.shell_distances(kcart, lattice, cells;
                                    search_shells=search_shells, tol=kmesh_tol)
    shell_list, _ = WF.select_shells(kcart, lattice, cells, dnn, multi; kmesh_tol=kmesh_tol)
    nnlist, nncell, nntot = WF.build_nnlist(kpoints, lattice, cells, dnn, shell_list, multi;
                                            tol=kmesh_tol)
    kpb = Matrix{Int}(undef, nntot, nk)
    gpb = Array{Int,3}(undef, 3, nntot, nk)
    for k in 1:nk, nn in 1:nntot
        kpb[nn, k] = nnlist[k, nn]
        gpb[:, nn, k] = nncell[:, k, nn]
    end

    # overlaps and projections from DFTK's plane-wave machinery
    M = Array{ComplexF64,4}(undef, num_bands, num_bands, nntot, nk)
    for ik in 1:nk, nn in 1:nntot
        M[:, :, nn, ik] = DFTK.overlap_Mmn_k_kpb(basis, scfres.ψ, ik, kpb[nn, ik],
                                                 gpb[:, nn, ik], num_bands)
    end
    A = Array{ComplexF64,3}(undef, num_bands, num_wann, nk)
    for ik in 1:nk
        A[:, :, ik] = DFTK.compute_amn_kpoint(basis, basis.kpoints[ik], scfres.ψ[ik],
                                              projections, num_bands)
    end

    return WannierFunctions.wannier_model(; unit_cell=unit_cell, kpoints=kpoints,
                                          mp_grid=mp_grid, num_wann=num_wann, M=M, A=A,
                                          kpb=kpb, gpb=gpb, eig=eig, kmesh_tol=kmesh_tol,
                                          seedname="dftk")
end

end # module
