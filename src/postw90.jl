# Drop-in `postw90.x` driver: read the postw90 keywords from `seedname.win`, build the
# interpolation models from `seedname.chk(.fmt)` + companion files, dispatch every requested
# module, and write the reference-named output files in the reference formats.
#
# The physics lives in the individual modules (berry/kubo/dos/shc/…, each oracle-validated);
# this file is the keyword → kwargs mapping, the output files that had no writer yet
# (kubo_S/A, jdos, ahc/morb fermi scans, BoltzWann set), and the Fortran G-format needed for
# byte parity of the geninterp/BoltzWann outputs.

using LinearAlgebra
using Printf
using StaticArrays

"""
    fortran_g(x, w, d) -> String

Fortran `Gw.d` edit descriptor: numbers with `0.1 ≤ |x| < 10^d` (after rounding to d
significant digits) print as `F(w-4).(d-n)` followed by four blanks (n = digits before the
point, 0 for |x| < 1); everything else falls back to `Ew.d`.
"""
function fortran_g(x::Real, w::Int, d::Int)
    if x == 0
        return lpad(@sprintf("%.*f", d - 1, 0.0), w - 4) * "    "
    end
    ax = abs(x)
    r = round(ax, sigdigits=d)
    if 0.1 <= r < 10.0^d
        n = r < 1.0 ? 0 : floor(Int, log10(r)) + 1
        return lpad(@sprintf("%.*f", d - n, Float64(x)), w - 4) * "    "
    end
    return fortran_e(Float64(x), w, d)
end

# ---------------------------------------------------------------------------
# keyword parsing helpers
# ---------------------------------------------------------------------------

_pw_tokens(s) = [lowercase(t) for t in split(s, r"[,;+\s]+") if !isempty(t)]

"Parse `key = n` or `key = n1 n2 n3` into an NTuple{3,Int} k-mesh."
function _pw_kmesh(raw, key, default::NTuple{3,Int})
    haskey(raw, key) || return default
    v = parse.(Int, split(raw[key]))
    return length(v) == 1 ? (v[1], v[1], v[1]) : (v[1], v[2], v[3])
end

"The postw90 Fermi-level list: fermi_energy, or fermi_energy_min/max(/step)."
function _pw_fermi_list(raw)
    if haskey(raw, "fermi_energy_min")
        fmin = parse_f64(split(raw["fermi_energy_min"])[1])
        fmax = haskey(raw, "fermi_energy_max") ? parse_f64(split(raw["fermi_energy_max"])[1]) : fmin
        step = haskey(raw, "fermi_energy_step") ? parse_f64(split(raw["fermi_energy_step"])[1]) : 0.01
        n = floor(Int, (fmax - fmin) / step + 1e-8) + 1
        return [fmin + (i - 1) * step for i in 1:n]
    end
    haskey(raw, "fermi_energy") && return [parse_f64(split(raw["fermi_energy"])[1])]
    return Float64[]
end

"kubo/gyrotropic-style frequency list from `<prefix>_freq_min/max/step`."
function _pw_freq_list(raw, prefix, fmax_default)
    fmin = _getfloat(raw, prefix * "_freq_min", 0.0)
    fmax = _getfloat(raw, prefix * "_freq_max", fmax_default)
    step = _getfloat(raw, prefix * "_freq_step", 0.01)
    n = floor(Int, (fmax - fmin) / step + 1e-8) + 1
    return [fmin + (i - 1) * step for i in 1:n]
end

"postw90's kubo_eigval_max default: the frozen-window top + 2/3, else max eigenvalue + 2/3."
function _pw_eigval_max(raw, win::WinInput, eigmax::Float64)
    haskey(raw, "kubo_eigval_max") && return parse_f64(split(raw["kubo_eigval_max"])[1])
    win.dis_froz_max != -Inf && return win.dis_froz_max + 2.0 / 3.0
    return eigmax + 2.0 / 3.0
end

# ---------------------------------------------------------------------------
# writers that had no file counterpart yet (reference formats from berry.F90 /
# boltzwann.F90 / geninterp.F90)
# ---------------------------------------------------------------------------

"Write `-kubo_S_ab.dat` / `-kubo_A_ab.dat` / `-jdos.dat` (3E16.8 / 2E16.8, berry.F90)."
function write_kubo(seedname::AbstractString, res::KuboResult)
    dirs = ("x", "y", "z")
    pairs_S = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))
    pairs_A = ((1, 2), (2, 3), (3, 1))
    for (i, j) in pairs_S
        open(seedname * "-kubo_S_" * dirs[i] * dirs[j] * ".dat", "w") do io
            for f in 1:length(res.freqs)
                z = kubo_S(res, i, j, f)
                println(io, fortran_e(res.freqs[f], 16, 8), fortran_e(real(z), 16, 8),
                        fortran_e(imag(z), 16, 8))
            end
        end
    end
    for (i, j) in pairs_A
        open(seedname * "-kubo_A_" * dirs[i] * dirs[j] * ".dat", "w") do io
            for f in 1:length(res.freqs)
                z = kubo_A(res, i, j, f)
                println(io, fortran_e(res.freqs[f], 16, 8), fortran_e(real(z), 16, 8),
                        fortran_e(imag(z), 16, 8))
            end
        end
    end
    open(seedname * "-jdos.dat", "w") do io
        for f in 1:length(res.freqs)
            println(io, fortran_e(res.freqs[f], 16, 8), fortran_e(res.jdos[f], 16, 8))
        end
    end
    return seedname
end

"Write `-ahc-fermiscan.dat` or `-morb-fermiscan.dat` rows `(4(F12.6,1x))` (berry.F90)."
function write_fermiscan(path::AbstractString, fermis::Vector{Float64}, vals)
    open(path, "w") do io
        for (i, ef) in enumerate(fermis)
            @printf(io, "%12.6f %12.6f %12.6f %12.6f \n", ef, vals[i][1], vals[i][2], vals[i][3])
        end
    end
    return path
end

"""
Write the BoltzWann output set (`_tdf/_elcond/_sigmas/_seebeck/_kappa.dat`, boltzwann.F90
headers and G18.10 rows). Component order xx xy yy xz yz zz; Seebeck is the full 3×3.
"""
function write_boltzwann(seedname::AbstractString, r::BoltzWannResult)
    voigt = ((1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3))
    g(x) = fortran_g(x, 18, 10)
    open(seedname * "_tdf.dat", "w") do io
        println(io, "# Written by the BoltzWann module of the Wannier90 code.")
        println(io, "# Transport distribution function (in units of 1/hbar^2 * eV * fs / angstrom)" *
                    " vs energy in eV")
        println(io, "# Content of the columns:")
        println(io, "# Energy TDF_xx TDF_xy TDF_yy TDF_xz TDF_yz TDF_zz")
        println(io, "#   (if spin decomposition is required, 12 further columns are provided, with the 6")
        println(io, "#    components of the TDF for the spin up, followed by those for the spin down)")
        for i in 1:length(r.energies)
            println(io, g(r.energies[i]), (g(r.tdf[c, i]) for c in 1:6)...)
        end
    end
    blocks = (("_elcond.dat", "# [Electrical conductivity in SI units, i.e. in 1/Ohm/m]",
               "# Mu(eV) Temp(K) ElCond_xx ElCond_xy ElCond_yy ElCond_xz ElCond_yz ElCond_zz",
               r.elcond, false),
              ("_sigmas.dat",
               "# [(Electrical conductivity * Seebeck coefficient) in SI units, i.e. in Ampere/m/K]",
               "# Mu(eV) Temp(K) (Sigma*S)_xx (Sigma*S)_xy (Sigma*S)_yy (Sigma*S)_xz (Sigma*S)_yz (Sigma*S)_zz",
               r.sigmas, false),
              ("_seebeck.dat", "# [Seebeck coefficient in SI units, i.e. in V/K]",
               "# Mu(eV) Temp(K) Seebeck_xx Seebeck_xy Seebeck_xz Seebeck_yx Seebeck_yy Seebeck_yz Seebeck_zx Seebeck_zy Seebeck_zz",
               r.seebeck, true),
              ("_kappa.dat",
               "# [K coefficient in SI units, i.e. in W/m/K]\n" *
               "# [the K coefficient is defined in the documentation, and is an ingredient of\n" *
               "#  the thermal conductivity. See the docs for further information.]",
               "# Mu(eV) Temp(K) Kappa_xx Kappa_xy Kappa_yy Kappa_xz Kappa_yz Kappa_zz",
               r.kappa, false))
    for (suffix, unitline, colline, arr, full) in blocks
        open(seedname * suffix, "w") do io
            println(io, "# Written by the BoltzWann module of the Wannier90 code.")
            println(io, unitline)
            println(io, colline)
            for (iμ, μ) in enumerate(r.mus), (iT, T) in enumerate(r.temps)
                comps = full ? vec(permutedims(arr[:, :, iμ, iT])) :   # row-major: xx xy xz yx …
                        [arr[i, j, iμ, iT] for (i, j) in voigt]
                println(io, g(μ), g(T), (g(c) for c in comps)...)
            end
        end
    end
    return seedname
end

"Write `_boltzdos.dat` (fixed-width / unsmeared variant; boltzwann.F90:1180)."
function write_boltzdos(seedname::AbstractString, es::Vector{Float64}, dos::Vector{Float64};
                        smr_width::Float64=0.0, volume::Float64)
    g(x) = fortran_g(x, 18, 10)
    binwidth = length(es) > 1 ? es[2] - es[1] : 1.0
    open(seedname * "_boltzdos.dat", "w") do io
        println(io, "# Written by the BoltzWann module of the Wannier90 code.")
        println(io, "# The first column.")
        if smr_width / binwidth < 2.0
            println(io, "# The second column is the unsmeared DOS.")
        else
            println(io, "# The second column is the DOS for a fixed smearing of ",
                    strip(fortran_g(smr_width, 14, 6)), " eV.")
        end
        println(io, "# Cell volume (ang^3):  ", fortran_g(volume, 14, 6))   # '(A,1X,G14.6)'
        println(io, "# Energy(eV) DOS [DOS DOS ...]")
        for i in 1:length(es)                       # '(1X,2G18.10)'
            println(io, " ", g(es[i]), g(dos[i]))
        end
    end
    return seedname
end

# ---------------------------------------------------------------------------
# the driver
# ---------------------------------------------------------------------------

"""
    postw90_main(seedname; verbose=true)

Drop-in `postw90.x`: dispatch every module requested in `seedname.win` (berry tasks
ahc/morb/kubo/sc/shc/kdotp, gyrotropic, dos, kpath, kslice, geninterp, boltzwann,
spin_moment) and write the reference-named output files.
"""
function postw90_main(seedname::AbstractString; verbose::Bool=true)
    win = read_win(seedname * ".win")
    raw = win.raw
    say(x) = verbose && println(x)

    eig = read_eig(seedname * ".eig")
    eigmax = maximum(eig)
    fermis = _pw_fermi_list(raw)
    # postw90's global interpolation mesh keyword is `kmesh`; module-specific ones override.
    kmesh_glob = _pw_kmesh(raw, "kmesh", win.mp_grid)
    kmesh_def = _pw_kmesh(raw, "berry_kmesh", kmesh_glob)

    # lazy, shared model instances (typed small-Union Refs: `Any` here would make every
    # downstream `bm()`/`mm()` call dynamically dispatched)
    _bm = Ref{Union{Nothing,BerryModel}}(nothing)
    bm() = (_bm[] === nothing && (_bm[] = BerryModel(seedname)); _bm[]::BerryModel)
    _mm = Ref{Union{Nothing,MorbModel}}(nothing)
    mm() = (_mm[] === nothing &&
            (_mm[] = MorbModel(seedname; transl_inv_full=_getbool(raw, "transl_inv_full", false)));
            _mm[]::MorbModel)
    function shcmodel()
        method = lowercase(get(raw, "shc_method", "qiao"))
        kw = Dict{Symbol,Any}()
        if _getbool(raw, "shc_bandshift", false)
            kw[:scissors_shift] = _getfloat(raw, "shc_bandshift_energyshift", 0.0)
            kw[:num_valence_bands] = _getint(raw, "shc_bandshift_firstband", 1) - 1
        end
        occursin("ryoo", method) &&
            return ShcRyooModel(seedname; transl_inv_full=_getbool(raw, "transl_inv_full", false))
        return ShcModel(seedname; kw...)
    end
    smr_kubo = (adaptive=_getbool(raw, "kubo_adpt_smr", true),
                adpt_fac=_getfloat(raw, "kubo_adpt_smr_fac", sqrt(2.0)),
                adpt_max=_getfloat(raw, "kubo_adpt_smr_max", 1.0),
                smr_width=_getfloat(raw, "kubo_smr_fixed_en_width", 0.0))

    # ---------------- berry ----------------
    if _getbool(raw, "berry", false)
        # the reference matches berry_task by substring (index()), e.g. "eval_shc" → shc
        taskstr = lowercase(get(raw, "berry_task", ""))
        tasks = [t for t in ("ahc", "morb", "kubo", "sc", "shc", "kdotp") if occursin(t, taskstr)]
        if "ahc" in tasks
            isempty(fermis) && error("berry_task=ahc needs fermi_energy(_min/max)")
            ahc = ahc_fermiscan(bm(); fermi_energies=fermis, kmesh=kmesh_def,
                                adpt_kmesh=_getint(raw, "berry_curv_adpt_kmesh", 1),
                                adpt_thresh=_getfloat(raw, "berry_curv_adpt_kmesh_thresh", 100.0))
            vals = collect(eachcol(ahc))                  # 3 × nf matrix → per-Fermi vectors
            if length(fermis) > 1
                write_fermiscan(seedname * "-ahc-fermiscan.dat", fermis, vals)
                say("  * $(seedname)-ahc-fermiscan.dat")
            end
            for (ef, v) in zip(fermis, vals)
                say(@sprintf("  AHC (S/cm)  E_f=%10.4f   x: %10.4f  y: %10.4f  z: %10.4f",
                             ef, v[1], v[2], v[3]))
            end
        end
        if "morb" in tasks
            isempty(fermis) && error("berry_task=morb needs fermi_energy(_min/max)")
            vals = [orbital_magnetisation(mm(); fermi_energy=ef, kmesh=kmesh_def) for ef in fermis]
            if length(fermis) > 1
                write_fermiscan(seedname * "-morb-fermiscan.dat", fermis, vals)
                say("  * $(seedname)-morb-fermiscan.dat")
            end
            for (ef, v) in zip(fermis, vals)
                say(@sprintf("  M_orb (μ_B/cell)  E_f=%10.4f   x: %10.4f  y: %10.4f  z: %10.4f",
                             ef, v[1], v[2], v[3]))
            end
        end
        if "kubo" in tasks
            isempty(fermis) && error("berry_task=kubo needs fermi_energy")
            freqs = _pw_freq_list(raw, "kubo", eigmax - minimum(eig))
            res = optical_conductivity(bm(); fermi_energy=fermis[1], kmesh=kmesh_def,
                                       freqs=freqs,
                                       eigval_max=_pw_eigval_max(raw, win, eigmax), smr_kubo...)
            write_kubo(seedname, res)
            say("  * $(seedname)-kubo_S_*.dat / -kubo_A_*.dat / -jdos.dat")
        end
        if "sc" in tasks
            isempty(fermis) && error("berry_task=sc needs fermi_energy")
            freqs = _pw_freq_list(raw, "kubo", eigmax - minimum(eig))
            sc = shift_current(seedname; fermi_energy=fermis[1], freqs=freqs, kmesh=kmesh_def,
                               phase_conv=_getint(raw, "sc_phase_conv", 1),
                               sc_eta=_getfloat(raw, "sc_eta", 0.04),
                               w_thr=_getfloat(raw, "sc_w_thr", 5.0),
                               eigval_max=_pw_eigval_max(raw, win, eigmax), smr_kubo...)
            write_shift_current(seedname, freqs, sc)
            say("  * $(seedname)-sc_*.dat")
        end
        if "shc" in tasks
            isempty(fermis) && error("berry_task=shc needs fermi_energy(_min/max)")
            sm = shcmodel()
            γ = _getint(raw, "shc_gamma", 3)
            α = _getint(raw, "shc_alpha", 1)
            β = _getint(raw, "shc_beta", 2)
            emax = _pw_eigval_max(raw, win, eigmax)
            if _getbool(raw, "tetrahedron_method", false)
                if _getbool(raw, "shc_freq_scan", false)
                    freqs = _pw_freq_list(raw, "shc", 8.0)
                    out = shc_tetra(sm; kmesh=kmesh_def, freqs=freqs, fermi_energy=fermis[1],
                                    γ=γ, α=α, β=β, eigval_max=emax)
                    write_shc(seedname * "-shc-freqscan.dat", freqs, out; freq_scan=true)
                else
                    out = shc_tetra(sm; kmesh=kmesh_def, fermi_energies=fermis,
                                    γ=γ, α=α, β=β, eigval_max=emax)
                    write_shc(seedname * "-shc-fermiscan.dat", fermis, out)
                end
            elseif _getbool(raw, "shc_freq_scan", false)
                freqs = _pw_freq_list(raw, "shc", 8.0)
                out = shc_freqscan(sm; freqs=freqs, fermi_energy=fermis[1], kmesh=kmesh_def,
                                   γ=γ, α=α, β=β, eigval_max=emax, smr_kubo...)
                write_shc(seedname * "-shc-freqscan.dat", freqs, out; freq_scan=true)
            else
                out = shc_fermiscan(sm; fermi_energies=fermis, kmesh=kmesh_def,
                                    γ=γ, α=α, β=β, eigval_max=emax, smr_kubo...)
                write_shc(seedname * "-shc-fermiscan.dat", fermis, out)
            end
            say("  * $(seedname)-shc-*.dat")
        end
        if "kdotp" in tasks
            kp = haskey(raw, "kdotp_kpoint") ? parse_f64.(split(raw["kdotp_kpoint"])) : [0.0, 0.0, 0.0]
            bands = haskey(raw, "kdotp_bands") ?
                    parse_range_list(replace(strip(raw["kdotp_bands"]), r"\s+" => ",")) :
                    collect(1:_getint(raw, "kdotp_num_bands", 0))
            res = kdotp(bm(); kpoint=kp, bands=bands)
            write_kdotp(seedname, res)
            say("  * $(seedname)-kdotp_*.dat")
        end
    end

    # ---------------- gyrotropic ----------------
    if _getbool(raw, "gyrotropic", false)
        # task string may be glued ("-C-dos-D0"); the reference matches by substring
        ts = lowercase(get(raw, "gyrotropic_task", "-d0-dw-c-k-noa-dos"))
        gtasks = Tuple(sym for (key, sym) in (("-d0", :D0), ("-dw", :Dw), ("-c", :C),
                                              ("-k", :K), ("-noa", :NOA), ("-dos", :dos))
                       if occursin(key, ts))
        m = isfile(seedname * ".uHu") ? mm() : bm()
        # integration box: b1/b2/b3 rows (fractional) around gyrotropic_box_center
        box = Matrix{Float64}(I, 3, 3)
        for (r, key) in enumerate(("gyrotropic_box_b1", "gyrotropic_box_b2", "gyrotropic_box_b3"))
            haskey(raw, key) && (box[r, :] = parse_f64.(split(raw[key])))
        end
        center = haskey(raw, "gyrotropic_box_center") ?
                 parse_f64.(split(raw["gyrotropic_box_center"])) : [0.5, 0.5, 0.5]
        corner = center - 0.5 .* vec(sum(box; dims=1))
        res = gyrotropic(m;
                         tasks=gtasks, fermi_energies=fermis,
                         freqs=_pw_freq_list(raw, "gyrotropic", 0.0),
                         kmesh=_pw_kmesh(raw, "gyrotropic_kmesh", kmesh_def),
                         smr_width=_getfloat(raw, "gyrotropic_smr_fixed_en_width", 0.1),
                         smr_max_arg=_getfloat(raw, "gyrotropic_smr_max_arg", 5.0),
                         box=box, box_corner=corner,
                         degen_thresh=_getfloat(raw, "gyrotropic_degen_thresh", 0.0),
                         eigval_max=_getfloat(raw, "gyrotropic_eigval_max", Inf),
                         band_list=(haskey(raw, "gyrotropic_band_list") ?
                                    parse_range_list(replace(strip(raw["gyrotropic_band_list"]),
                                                             r"\s+" => ",")) : nothing))
        write_gyrotropic(seedname, res; tasks=gtasks, spin=isfile(seedname * ".spn"))
        say("  * $(seedname)-gyrotropic-*.dat")
    end

    # ---------------- dos ----------------
    if _getbool(raw, "dos", false)
        emin = _getfloat(raw, "dos_energy_min", minimum(eig) - 0.6667)
        emax_ = _getfloat(raw, "dos_energy_max", eigmax + 0.6667)
        step = _getfloat(raw, "dos_energy_step", 0.01)
        n = floor(Int, (emax_ - emin) / step + 1e-8) + 1
        es = [emin + (i - 1) * step for i in 1:n]
        spinm = _getbool(raw, "spin_decomp", false) ? SpinModel(seedname) : nothing
        project = haskey(raw, "dos_project") ?
                  parse_range_list(replace(strip(raw["dos_project"]), r"\s+" => ",")) : nothing
        out = density_of_states(bm(); energies=es,
                                kmesh=_pw_kmesh(raw, "dos_kmesh", kmesh_def),
                                adaptive=_getbool(raw, "dos_adpt_smr", true),
                                adpt_fac=_getfloat(raw, "dos_adpt_smr_fac", sqrt(2.0)),
                                adpt_max=_getfloat(raw, "dos_adpt_smr_max", 1.0),
                                smr_width=_getfloat(raw, "dos_smr_fixed_en_width", 0.0),
                                spin=spinm, project=project,
                                elec_per_state=(_getbool(raw, "spinors", false) ? 1 : 2))
        write_dos(seedname * "-dos.dat", out[1], out[2:end]...)
        say("  * $(seedname)-dos.dat")
    end

    # ---------------- kpath ----------------
    if _getbool(raw, "kpath", false)
        tokens = _pw_tokens(get(raw, "kpath_task", "bands"))
        ktasks0 = [Symbol(t) for t in tokens if t in ("bands", "curv", "morb", "shc")]
        if :morb in ktasks0 && !isfile(seedname * ".uHu")
            @warn "kpath morb task skipped: no $(seedname).uHu"
            filter!(!=(:morb), ktasks0)
        end
        ktasks = Tuple(ktasks0)
        m = :shc in ktasks ? shcmodel() : (:morb in ktasks ? mm() : bm())
        res = kpath(m;
                    segments=kpath_segments(win),
                    num_points=_getint(raw, "kpath_num_points", 100),
                    tasks=ktasks,
                    bands_colour=Symbol(lowercase(get(raw, "kpath_bands_colour", "none"))),
                    fermi_energy=(isempty(fermis) ? nothing : fermis[1]),
                    curv_unit=Symbol(lowercase(get(raw, "berry_curv_unit", "ang2"))))
        write_kpath(seedname, res)
        say("  * $(seedname)-path.kpt / -bands.dat / task files")
    end

    # ---------------- kslice ----------------
    if _getbool(raw, "kslice", false)
        tokens = _pw_tokens(get(raw, "kslice_task", "bands"))
        # fermi_lines is a plot-script task; its data is the bands grid we always write
        "fermi_lines" in tokens && !("bands" in tokens) && push!(tokens, "bands")
        stasks0 = [Symbol(t) for t in tokens if t in ("bands", "curv", "morb", "shc")]
        if :morb in stasks0 && !isfile(seedname * ".uHu")
            @warn "kslice morb task skipped: no $(seedname).uHu"
            filter!(!=(:morb), stasks0)
        end
        stasks = Tuple(stasks0)
        m = :shc in stasks ? shcmodel() : (:morb in stasks ? mm() : bm())
        n2d = haskey(raw, "kslice_2dkmesh") ? parse.(Int, split(raw["kslice_2dkmesh"])) : [50]
        mesh2 = length(n2d) == 1 ? (n2d[1], n2d[1]) : (n2d[1], n2d[2])
        res = kslice(m;
                     corner=(haskey(raw, "kslice_corner") ?
                             parse_f64.(split(raw["kslice_corner"])) : [0.0, 0.0, 0.0]),
                     b1=parse_f64.(split(raw["kslice_b1"])),
                     b2=parse_f64.(split(raw["kslice_b2"])),
                     mesh=mesh2, tasks=stasks,
                     fermi_energy=(isempty(fermis) ? nothing : fermis[1]),
                     curv_unit=Symbol(lowercase(get(raw, "berry_curv_unit", "ang2"))))
        write_kslice(seedname, res.coords, res.bands;
                     curv=(:curv in stasks ? res.curv : nothing),
                     morb=(:morb in stasks ? res.morb : nothing),
                     shc=(:shc in stasks ? res.shc : nothing))
        say("  * $(seedname)-kslice-*.dat")
    end

    # ---------------- geninterp ----------------
    if _getbool(raw, "geninterp", false)
        geninterp(bm(), seedname; alsofirstder=_getbool(raw, "geninterp_alsofirstder", false))
        say("  * $(seedname)_geninterp.dat")
    end

    # ---------------- boltzwann ----------------
    if _getbool(raw, "boltzwann", false)
        μmin = _getfloat(raw, "boltz_mu_min", NaN)
        μmax = _getfloat(raw, "boltz_mu_max", μmin)
        μstep = _getfloat(raw, "boltz_mu_step", 1.0)
        isnan(μmin) && error("boltzwann needs boltz_mu_min")
        nμ = max(1, floor(Int, (μmax - μmin) / μstep + 1e-8) + 1)
        mus = [μmin + (i - 1) * μstep for i in 1:nμ]
        Tmin = _getfloat(raw, "boltz_temp_min", 300.0)
        Tmax = _getfloat(raw, "boltz_temp_max", Tmin)
        Tstep = _getfloat(raw, "boltz_temp_step", 1.0)
        nT = max(1, floor(Int, (Tmax - Tmin) / Tstep + 1e-8) + 1)
        temps = [Tmin + (i - 1) * Tstep for i in 1:nT]
        # the TDF energy window is the disentanglement window (eig bounds when unset)
        wmin = win.dis_win_min == -Inf ? minimum(eig) : win.dis_win_min
        wmax = win.dis_win_max == Inf ? eigmax : win.dis_win_max
        r = boltzwann(bm();
                      kmesh=_pw_kmesh(raw, "boltz_kmesh", kmesh_glob),
                      relax_time=_getfloat(raw, "boltz_relax_time", 10.0),
                      mus=mus, temps=temps,
                      tdf_energy_step=_getfloat(raw, "boltz_tdf_energy_step", 0.001),
                      tdf_smr_width=_getfloat(raw, "boltz_tdf_smr_fixed_en_width", 0.0),
                      win=(wmin, wmax),
                      elec_per_state=(_getbool(raw, "spinors", false) ? 1 : 2))
        write_boltzwann(seedname, r)
        say("  * $(seedname)_tdf/_elcond/_sigmas/_seebeck/_kappa.dat")
        if _getbool(raw, "boltz_calc_also_dos", false)
            emin = _getfloat(raw, "boltz_dos_energy_min", minimum(eig) - 0.6667)
            emax_ = _getfloat(raw, "boltz_dos_energy_max", eigmax + 0.6667)
            step = _getfloat(raw, "boltz_dos_energy_step", 0.001)
            n = floor(Int, (emax_ - emin) / step + 1e-8) + 1
            es = [emin + (i - 1) * step for i in 1:n]
            smr = _getfloat(raw, "boltz_dos_smr_fixed_en_width", 0.0)
            _, d = density_of_states(bm(); energies=es,
                                     kmesh=_pw_kmesh(raw, "boltz_kmesh", kmesh_def),
                                     adaptive=_getbool(raw, "boltz_dos_adpt_smr", true),
                                     adpt_fac=_getfloat(raw, "boltz_dos_adpt_smr_fac", sqrt(2.0)),
                                     smr_width=smr,
                                     elec_per_state=(_getbool(raw, "spinors", false) ? 1 : 2))
            write_boltzdos(seedname, es, d; smr_width=smr,
                           volume=cell_volume(bm().lattice))
            say("  * $(seedname)_boltzdos.dat")
        end
    end

    # ---------------- spin moment ----------------
    if _getbool(raw, "spin_moment", false)
        isempty(fermis) && error("spin_moment needs fermi_energy")
        r = spin_moment(SpinModel(seedname); fermi_energy=fermis[1], kmesh=kmesh_def)
        say(@sprintf("  Spin moment (μ_B/cell): (%.6f %.6f %.6f)  |S| = %.6f",
                     r.moment..., norm(r.moment)))
    end

    return nothing
end
