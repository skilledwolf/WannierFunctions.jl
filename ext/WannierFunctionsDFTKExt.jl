# Package extension: build a WannierFunctions `Model` from a DFTK.jl SCF result, so an
# all-Julia DFT → Wannier pipeline runs in memory with no file round-trip. Loaded automatically
# when both WannierFunctions and DFTK are present. End-to-end use needs a live DFTK SCF; the
# array hand-off itself is exercised by the file-vs-memory consistency test in the test suite.
module WannierFunctionsDFTKExt

using WannierFunctions
using DFTK

"""
    WannierFunctions.wannier_model(scfres::DFTK.NamedTuple, projections; num_wann, kmesh) -> Model

Assemble a wannierisation `Model` from a DFTK `scfres` and a set of trial `projections` (the
`.amn` guesses). Uses DFTK's own overlap/projection routines to produce the `.mmn`/`.amn`/`.eig`
content in memory, then hands off to `WannierFunctions.wannier_model`. `kmesh` is the
Monkhorst–Pack grid used for the SCF.
"""
function WannierFunctions.wannier_model(scfres, projections; num_wann::Int,
                                        kmesh::NTuple{3,Int})
    basis = scfres.basis
    model = basis.model
    unit_cell = Array(model.lattice) .* DFTK.units.bohr_to_A     # DFTK lattice is in bohr
    kpoints = [Vector(kpt.coordinate) for kpt in basis.kpoints]
    eig = hcat([scfres.eigenvalues[ik] .* DFTK.units.Ha_to_eV
                for ik in 1:length(basis.kpoints)]...)
    # DFTK ships mmn/amn writers for wannier90; call them for the arrays. The exact routine
    # names track the installed DFTK version — see DFTK's `Wannier`/`save_wannier` interface.
    mmn, kpb, gpb = DFTK.compute_mmn(basis, scfres.ψ)            # ⟨u_k|u_{k+b}⟩ + neighbour map
    amn = DFTK.compute_amn(basis, scfres.ψ, projections)        # ⟨ψ|g⟩
    return WannierFunctions.wannier_model(; unit_cell=unit_cell, kpoints=kpoints, mp_grid=kmesh,
                                          num_wann=num_wann, M=mmn, A=amn, kpb=kpb, gpb=gpb,
                                          eig=eig, seedname="dftk")
end

end # module
