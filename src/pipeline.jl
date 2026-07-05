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
    Mrot::Array{ComplexF64,4}      # final-gauge num_wann overlaps (position operator, spreads)
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
        return WannierResult(res.U, res.Mrot, dis.eigval_opt, res.spread, true, dis.omega_I, res.niter, res.converged, dis)
    else
        res = wannierise(model; num_iter=num_iter, algorithm=algorithm, verbose=verbose)
        return WannierResult(res.U, res.Mrot, model.eig, res.spread, false, res.spread.ΩI, res.niter, res.converged, nothing)
    end
end

"Compat method: drive the pipeline from a parsed `.win` (reference-faithful `:w90` optimiser)."
function run_wannier(model::Model, win::WinInput; verbose::Bool=false, sitesym=nothing)
    # Optional :w90-path localisation controls from the .win.
    lopts = Dict{Symbol,Any}()
    sitesym === nothing || (lopts[:sitesym] = sitesym)
    if _getbool(win.raw, "guiding_centres", false)
        # Branch-cut guides initialised from the projection centres (Cartesian Å), in WF order.
        projs = parse_projections(win)
        if length(projs) == model.num_wann
            lopts[:guides] = hcat([model.lattice.A * p.site for p in projs]...)
        end
    end
    if _getbool(win.raw, "precond", false)
        irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
        Rcart = [model.lattice.A * SVector{3,Float64}(r...) for r in irvec]
        lopts[:precond] = (kfrac=model.kgrid.frac, irvec=irvec, ndegen=ndegen, Rcart=Rcart)
    end
    # SLWF+C selective localisation (isolated case): route through the :rcg path which carries
    # the Ω_C objective + gradient. Active when slwf_num < num_wann.
    slwf = _slwf_from_win(win, model)
    if slwf !== nothing || sitesym !== nothing
        algo, ni = :rcg, max(win.num_iter, 2000)      # SLWF+C / site_symmetry use the :rcg path
    else
        algo, ni = :w90, win.num_iter
    end
    sitesym === nothing || model.num_bands == model.num_wann ||
        error("site_symmetry with disentanglement (num_bands > num_wann) is not yet supported; " *
              "the localisation-phase symmetrisation is implemented for the isolated case")
    if model.num_bands > model.num_wann
        # win-aware disentanglement (honours dis_spheres / dis_froz_proj / dis_proj_* / windows)
        dis = disentangle(model, win; verbose=verbose)
        res = localize(dis.U0, dis.Mrot0, model.bvectors;
                       num_iter=ni, algorithm=algo, verbose=verbose, slwf=slwf, lopts...)
        return WannierResult(res.U, res.Mrot, dis.eigval_opt, res.spread, true, dis.omega_I,
                             res.niter, res.converged, dis)
    else
        res = wannierise(model; num_iter=ni, algorithm=algo, verbose=verbose, slwf=slwf, lopts...)
        return WannierResult(res.U, res.Mrot, model.eig, res.spread, false, res.spread.ΩI,
                             res.niter, res.converged, nothing)
    end
end

"Build an `SLWF` from the .win keywords (slwf_num/slwf_constrain/slwf_lambda/slwf_centres), or nothing."
function _slwf_from_win(win::WinInput, model::Model)
    num = _getint(win.raw, "slwf_num", model.num_wann)
    num >= model.num_wann && return nothing        # slwf_num == num_wann → selective_loc off
    constrain = _getbool(win.raw, "slwf_constrain", false)
    lambda = _getfloat(win.raw, "slwf_lambda", 1.0)
    centres = zeros(3, num)
    if constrain && haskey(win.blocks, "slwf_centres")
        for ln in win.blocks["slwf_centres"]
            t = split(ln)
            length(t) >= 4 || continue
            iw = parse(Int, t[1])
            1 <= iw <= num || continue
            centres[:, iw] = model.lattice.A * SVector{3,Float64}(parse_f64.(t[2:4])...)
        end
    end
    return SLWF(num, lambda, constrain, centres)
end

"""
    run_wannier(seedname; verbose=false) -> (model, WinInput, WannierResult)

Convenience: read `seedname.{win,amn,mmn,eig}` and run the full pipeline.
"""
function run_wannier(seedname::AbstractString; verbose::Bool=false)
    model = read_model(seedname)
    win = read_win(seedname * ".win")
    # site_symmetry: load the .dmn (symmetry-adapted Wannier functions).
    ss = nothing
    if _getbool(win.raw, "site_symmetry", false) && isfile(seedname * ".dmn")
        ss = read_dmn(seedname * ".dmn", model.num_bands, model.num_wann)
    end
    return model, win, run_wannier(model, win; verbose=verbose, sitesym=ss)
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
    if use_ws_distance
        irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
        Hr, _ = build_hr(result.U, result.eig_interp, model.kgrid, irvec)
        return interpolate_bands_ws(Hr, irvec, ndegen, result.spread.centres,
                                    model.lattice, model.kgrid.mp_grid, kpts)
    end
    return bands(hamiltonian_operator(model, result), kpts)
end
