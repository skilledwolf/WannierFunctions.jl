# Drop-in command-line driver: read `seedname.win` (+ .amn/.mmn/.eig), run the full
# wannierisation + interpolation pipeline, and write the standard output files, mirroring
# `wannier90.x`. Entry point `WannierFunctions.main(seedname)`; see bin/wannier90.jl for the CLI wrapper.

using Printf

_winflag(win::WinInput, key, default::Bool) = _getbool(win.raw, key, default)
_winint(win::WinInput, key, default::Int) = _getint(win.raw, key, default)

"""
    main(seedname; pp=false, write_files=true, verbose=true) -> (model, win, result)

Run wannierisation for `seedname` and write `.wout` plus, when requested in the `.win`,
`_hr.dat` (`write_hr`/`hr_plot`), `_tb.dat` (`write_tb`), and the band-structure files
(`bands_plot`). Returns the model, parsed input, and result.

With `pp=true` (the `-pp` flag, or `postproc_setup = .true.` in the `.win`) only the
post-processing setup runs: the k-mesh is generated from the `.win` alone and `seedname.nnkp`
is written for the DFT interface. No `.amn/.mmn` needed.
"""
function main(seedname::AbstractString; pp::Bool=false, write_files::Bool=true, verbose::Bool=true)
    seedname = replace(String(seedname), r"\.win$" => "")
    if pp || _getbool(read_win(seedname * ".win").raw, "postproc_setup", false)
        out, info = generate_nnkp(seedname)
        verbose && @info "wrote $out" nntot=info.nntot shells=info.shells weights=info.weights
        return nothing
    end
    model = read_model(seedname)
    win = read_win(seedname * ".win")

    verbose && @info "WannierFunctions.jl" seedname num_wann=model.num_wann num_bands=model.num_bands nkpt=nkpt(model.kgrid)
    _winflag(win, "guiding_centres", false) &&
        @warn "guiding_centres is set but not yet supported; using the principal Im-ln branch " *
              "(results may differ for poorly-localised initial projections)"
    result = run_wannier(model, win; verbose=verbose)
    verbose && @info("done",
        disentangled=result.disentangled, Ω=result.spread.Ω,
        ΩI=result.spread.ΩI, ΩD=result.spread.ΩD, ΩOD=result.spread.ΩOD)

    write_files || return model, win, result

    write_wout(seedname * ".wout", model, win, result; dis=result.dis)
    # Checkpoint, as wannier90.x writes at the end of a run — this is what lets postw90.x
    # (or a later restart) consume our result at full precision.
    write_chk(seedname * ".chk", Checkpoint(model, win, result))

    if _winflag(win, "berry", false)
        task = lowercase(get(win.raw, "berry_task", ""))
        if occursin("ahc", task) && haskey(win.raw, "fermi_energy")
            ef = parse_f64(win.raw["fermi_energy"])
            km = haskey(win.raw, "berry_kmesh") ?
                 Tuple(parse.(Int, split(win.raw["berry_kmesh"]))) : (25, 25, 25)
            length(km) == 1 && (km = (km[1], km[1], km[1]))
            bm = BerryModel(Checkpoint(model, win, result), model.eig, model.bvectors,
                            model.kgrid, model.lattice)
            ahc = anomalous_hall(bm; fermi_energy=ef, kmesh=km)
            verbose && @info "AHC (S/cm)" x = ahc[1] y = ahc[2] z = ahc[3] kmesh = km
            open(seedname * ".wpout", "a") do io
                @printf(io, "\n AHC (S/cm)       x          y          z\n")
                @printf(io, " ==========%11.4f%11.4f%11.4f\n", ahc[1], ahc[2], ahc[3])
            end
        else
            @warn "berry=true: only berry_task=ahc with fermi_energy set is supported" task
        end
    end

    if _winflag(win, "wannier_plot", false)
        try
            paths = plot_wannier_functions(model, win, result; seedname=seedname)
            verbose && @info "wannier_plot: wrote $(length(paths)) .xsf file(s)"
        catch err
            @warn "wannier_plot failed" error = sprint(showerror, err)
        end
    end

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
            if need_tb
                rop = position_operator(model, result)
                write_tb(seedname * "_tb.dat", model.lattice, model.num_wann, irvec, ndegen, Hr;
                         pos=rop.data)
            end
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
    # Ballistic transport (transport_mode = bulk): Landauer T(E) + DOS of the 1D chain.
    if _winflag(win, "transport", false)
        if result.eig_interp === nothing
            @warn "skipping transport: no band energies (.eig) available"
        else
            run_transport(model, win, result; seedname=seedname)
        end
    end
    return model, win, result
end
