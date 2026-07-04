# High-level pipeline: run the full wannierisation for a model, auto-selecting the isolated-bands
# path (num_bands == num_wann) or the disentanglement path (num_bands > num_wann), then expose a
# uniform result that the interpolation routines consume.

"""
    WannierResult

Outcome of a full wannierisation. `U` is the final square gauge (num_wann × num_wann × nkpt);
`eig_interp` are the per-k energies to Fourier-interpolate (the DFT eigenvalues for the isolated
case, or the disentangled subspace eigenvalues), both num_wann × nkpt.
"""
struct WannierResult
    U::Array{ComplexF64,3}
    eig_interp::Union{Nothing,Matrix{Float64}}
    spread::SpreadResult
    disentangled::Bool
    omega_I::Float64
    niter::Int
    dis::Union{Nothing,DisentangleResult}
end

"""
    run_wannier(model, win; verbose=false) -> WannierResult

Full pipeline. Disentangles first when `num_bands > num_wann`, then runs Marzari–Vanderbilt
localisation. `win` supplies the iteration counts and (for disentanglement) the energy windows.
"""
function run_wannier(model::Model, win::WinInput; verbose::Bool=false)
    if model.num_bands > model.num_wann
        dis = disentangle(model, win; verbose=verbose)
        res = localize(dis.U0, dis.Mrot0, model.bvectors; num_iter=win.num_iter, verbose=verbose)
        return WannierResult(res.U, dis.eigval_opt, res.spread, true, dis.omega_I, res.niter, dis)
    else
        res = wannierise(model; num_iter=win.num_iter, verbose=verbose)
        return WannierResult(res.U, model.eig, res.spread, false, res.spread.ΩI, res.niter, nothing)
    end
end

"""
    run_wannier(seedname; verbose=false) -> (model, WinInput, WannierResult)

Convenience: read `seedname.{win,amn,mmn,eig}` and run the full pipeline.
"""
function run_wannier(seedname::AbstractString; verbose::Bool=false)
    model = read_model(seedname)
    win = read_win(seedname * ".win")
    return model, win, run_wannier(model, win; verbose=verbose)
end

"""
    interpolate(model, result, kpts; use_ws_distance=false) -> energies

Interpolate band energies (num_wann × npts) at fractional k-points `kpts`, using the wannierised
gauge and (subspace) eigenvalues in `result`. Requires `result.eig_interp !== nothing`. With
`use_ws_distance=true` the per-Wannier-pair minimal-image improvement (the reference default) is
applied — slightly more accurate near cell boundaries.
"""
function interpolate(model::Model, result::WannierResult, kpts::Vector{SVector{3,Float64}};
                     use_ws_distance::Bool=false)
    result.eig_interp !== nothing ||
        error("no band energies available for interpolation (isolated case without .eig)")
    irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
    Hr, _ = build_hr(result.U, result.eig_interp, model.kgrid, irvec)
    if use_ws_distance
        return interpolate_bands_ws(Hr, irvec, ndegen, result.spread.centres,
                                    model.lattice, model.kgrid.mp_grid, kpts)
    end
    return interpolate_bands(Hr, irvec, ndegen, kpts)
end
