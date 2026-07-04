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
    converged::Bool
    dis::Union{Nothing,DisentangleResult}
end

"""
    run_wannier(model; win_min=-Inf, win_max=Inf, froz_min=-Inf, froz_max=nothing,
                dis_num_iter=200, dis_mix_ratio=0.5, num_iter=100,
                algorithm=:rcg, verbose=false) -> WannierResult

Keyword-first full pipeline. Disentangles first when `num_bands > num_wann` (using the given
energy windows, eV), then localises. `algorithm` selects the spread minimiser: `:rcg` (Riemannian
conjugate gradient with a true convergence criterion — the native default) or `:w90` (the
reference-faithful Wannier90 optimiser, fixed `num_iter` sweeps).

The `run_wannier(model, win::WinInput)` / `run_wannier(seedname)` methods drive everything from a
parsed `.win` instead and default to `:w90` for drop-in parity with `wannier90.x`.
"""
function run_wannier(model::Model;
                     win_min::Float64=-Inf, win_max::Float64=Inf,
                     froz_min::Float64=-Inf, froz_max::Union{Nothing,Float64}=nothing,
                     dis_num_iter::Int=200, dis_mix_ratio::Float64=0.5,
                     num_iter::Int=100, algorithm::Symbol=:rcg, verbose::Bool=false)
    if model.num_bands > model.num_wann
        dis = disentangle(model; win_min=win_min, win_max=win_max,
                          froz_min=froz_min, froz_max=froz_max,
                          num_iter=dis_num_iter, mix_ratio=dis_mix_ratio, verbose=verbose)
        res = localize(dis.U0, dis.Mrot0, model.bvectors;
                       num_iter=num_iter, algorithm=algorithm, verbose=verbose)
        return WannierResult(res.U, dis.eigval_opt, res.spread, true, dis.omega_I, res.niter, res.converged, dis)
    else
        res = wannierise(model; num_iter=num_iter, algorithm=algorithm, verbose=verbose)
        return WannierResult(res.U, model.eig, res.spread, false, res.spread.ΩI, res.niter, res.converged, nothing)
    end
end

"Compat method: drive the pipeline from a parsed `.win` (reference-faithful `:w90` optimiser)."
function run_wannier(model::Model, win::WinInput; verbose::Bool=false)
    return run_wannier(model;
        win_min=win.dis_win_min, win_max=win.dis_win_max,
        froz_min=win.dis_froz_min,
        froz_max=(win.dis_froz_max == -Inf ? nothing : win.dis_froz_max),
        dis_num_iter=win.dis_num_iter, dis_mix_ratio=win.dis_mix_ratio,
        num_iter=win.num_iter, algorithm=:w90, verbose=verbose)
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
