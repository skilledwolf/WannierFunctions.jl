# Drop-in command-line driver: read `seedname.win` (+ .amn/.mmn/.eig), run the full
# wannierisation + interpolation pipeline, and write the standard output files, mirroring
# `wannier90.x`. Entry point `Wannier90.main(seedname)`; see bin/wannier90.jl for the CLI wrapper.

using Printf

_winflag(win::WinInput, key, default::Bool) = _getbool(win.raw, key, default)
_winint(win::WinInput, key, default::Int) = _getint(win.raw, key, default)

"""
    main(seedname; write_files=true, verbose=true) -> (model, win, result)

Run wannierisation for `seedname` and write `.wout` plus, when requested in the `.win`,
`_hr.dat` (`write_hr`/`hr_plot`), `_tb.dat` (`write_tb`), and the band-structure files
(`bands_plot`). Returns the model, parsed input, and result.
"""
function main(seedname::AbstractString; write_files::Bool=true, verbose::Bool=true)
    seedname = replace(String(seedname), r"\.win$" => "")
    model = read_model(seedname)
    win = read_win(seedname * ".win")

    verbose && @info "Wannier90.jl" seedname num_wann=model.num_wann num_bands=model.num_bands nkpt=nkpt(model.kgrid)
    _winflag(win, "guiding_centres", false) &&
        @warn "guiding_centres is set but not yet supported; using the principal Im-ln branch " *
              "(results may differ for poorly-localised initial projections)"
    result = run_wannier(model, win; verbose=verbose)
    verbose && @info("done",
        disentangled=result.disentangled, Ω=result.spread.Ω,
        ΩI=result.spread.ΩI, ΩD=result.spread.ΩD, ΩOD=result.spread.ΩOD)

    write_files || return model, win, result

    write_wout(seedname * ".wout", model, win, result; dis=result.dis)

    # Real-space Hamiltonian outputs and band interpolation need per-k energies.
    need_hr = _winflag(win, "write_hr", false) || _winflag(win, "hr_plot", false)
    need_tb = _winflag(win, "write_tb", false)
    need_bands = _winflag(win, "bands_plot", false)
    if (need_hr || need_tb || need_bands)
        if result.eig_interp === nothing
            @warn "skipping H(R)/bands: no band energies (.eig) available"
        else
            irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
            Hr, _ = build_hr(result.U, result.eig_interp, model.kgrid, irvec)
            need_hr && write_hr(seedname * "_hr.dat", model.num_wann, irvec, ndegen, Hr)
            need_tb && write_tb(seedname * "_tb.dat", model.lattice, model.num_wann, irvec, ndegen, Hr)
            if need_bands
                npts = _winint(win, "bands_num_points", 100)
                kpts, xvals, labels, lidx = generate_kpath(win, model.lattice; bands_num_points=npts)
                if !isempty(kpts)
                    # use_ws_distance defaults to true in Wannier90; honour it for band output.
                    uws = _winflag(win, "use_ws_distance", true)
                    E = uws ?
                        interpolate_bands_ws(Hr, irvec, ndegen, result.spread.centres,
                                             model.lattice, model.kgrid.mp_grid, kpts) :
                        interpolate_bands(Hr, irvec, ndegen, kpts)
                    write_band_dat(seedname * "_band.dat", xvals, E)
                    write_band_kpt(seedname * "_band.kpt", kpts)
                    write_labelinfo(seedname * "_band.labelinfo.dat",
                                    labels, lidx, xvals[lidx], kpts[lidx])
                end
            end
        end
    end
    return model, win, result
end
