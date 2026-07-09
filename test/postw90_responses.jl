# Split from the original monolithic runtests.jl — included by runtests.jl (common.jl
# provides REFROOT/DATAROOT, the shared GaAs/diamond models, and helpers).
@testset "spin + kpath + gyrotropic + Ryoo SHC + pdos" begin
    # Stage a reference test dir into a tmpdir: plain files copied, .bz2 decompressed.
    function stage(dir, plain, zipped)
        tmp = mktempdir()
        for f in plain
            cp(joinpath(dir, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in zipped
            run(pipeline(`bunzip2 -kc $(joinpath(dir, f * ".bz2"))`, stdout = joinpath(tmp, f)))
        end
        tmp
    end

    # Spin moment + spin-decomposed DOS: Fe (spn is plain text here; no .mmn needed)
    sd = joinpath(REFROOT, "testpostw90_fe_spin")
    if isfile(joinpath(sd, "Fe.spn")) && Sys.which("bunzip2") !== nothing
        tmp = stage(sd, ("Fe.win", "Fe.eig", "Fe.spn"), ("Fe.chk.fmt",))
        sm = SpinModel(joinpath(tmp, "Fe"))
        r = spin_moment(sm; fermi_energy = 12.6279, kmesh = (4, 4, 4))
        @test r.moment[3] ≈ 3.090787 atol = 1e-4      # oracle wpout (tol 1e-3)
        @test abs(r.moment[1]) < 1e-4
        @test r.phi ≈ 5.831720 atol = 1e-2
        es, dtot, dup, ddn = density_of_states(sm.bm; energies = range(10.0, 13.0; length = 16),
                                               kmesh = (4, 4, 4), adaptive = false,
                                               smr_width = 0.5, spin = sm, elec_per_state = 1)
        @test dtot[1] ≈ 0.83154571 atol = 1e-4        # oracle Fe-dos.dat row 1
        @test dup[1] ≈ 0.11547299 atol = 1e-4
        @test ddn[1] ≈ 0.71607272 atol = 1e-4
        @test maximum(abs.(dtot .- dup .- ddn)) < 1e-10
    else
        @test_skip false
    end

    # Projected DOS: copper d-manifold (WF 1:5)
    pd = joinpath(REFROOT, "testpostw90_example04_pdos")
    if isfile(joinpath(pd, "copper.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(pd, ("copper.win", "copper.eig"), ("copper.chk.fmt",))
        bm = BerryModel(joinpath(tmp, "copper"))
        es, d = density_of_states(bm; energies = range(8.0, 10.0; length = 9),
                                  kmesh = (10, 10, 10), project = collect(1:5))
        @test d[1] ≈ 1.6146066 atol = 1e-4            # oracle copper-dos.dat rows 1, 9
        @test d[9] ≈ 2.8830387 atol = 1e-4
    else
        @test_skip false
    end

    # kpath: Fe bands+morb+curv (byte-level parity checked against staged oracle files in
    # scratch/ during development; here spot values, incl. the xval construction)
    kd = joinpath(REFROOT, "testpostw90_fe_kpathcurv")
    if isfile(joinpath(kd, "Fe.uHu.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(kd, ("Fe.win", "Fe.eig"), ("Fe.chk.fmt", "Fe.mmn", "Fe.uHu"))
        mm = MorbModel(joinpath(tmp, "Fe"))
        win = read_win(joinpath(tmp, "Fe.win"))
        res = kpath(mm; segments = kpath_segments(win), num_points = 10,
                    tasks = (:bands, :morb, :curv), fermi_energy = 12.6279)
        @test length(res.kpts) == 20
        @test res.xvals[2] ≈ 0.21892688 atol = 1e-7   # current-segment step
        @test res.xvals[11] ≈ 2.1810044 atol = 1e-6   # vertex xval trap (NOT cumulative len)
        @test res.bands[1, 1] ≈ 4.4341400 atol = 1e-4
        @test res.curv[3, 3] ≈ -5.6272413 atol = 1e-3  # −Ω convention
        @test res.morb[3, 1] ≈ 0.32522033 atol = 1e-3  # −(G+H−2E_F F)/2, eV·Å²
    else
        @test_skip false
    end

    # kpath SHC colouring: Pt (fixed smearing 1 eV, Å² band-resolved term)
    pk = joinpath(REFROOT, "testpostw90_pt_kpathshc")
    if isfile(joinpath(pk, "Pt.spn.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(pk, ("Pt.win", "Pt.eig"), ("Pt.chk.fmt", "Pt.mmn", "Pt.spn"))
        sm = ShcModel(joinpath(tmp, "Pt"))
        win = read_win(joinpath(tmp, "Pt.win"))
        res = kpath(sm; segments = kpath_segments(win), num_points = 10, tasks = (:shc,),
                    fermi_energy = 17.9919, smr_width = 1.0, eigval_max = 30.6667)
        @test length(res.shc) == 60
        @test res.shc[1] ≈ 0.073385059 atol = 1e-4    # oracle Pt-shc.dat row 1 (tol 0.1)
    else
        @test_skip false
    end

    # Gyrotropic: Te, all tasks (oracle: five of six files byte-identical)
    gd = joinpath(REFROOT, "testpostw90_te_gyrotropic")
    if isfile(joinpath(gd, "Te.uHu.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(gd, ("Te.win", "Te.eig", "Te.mmn"), ("Te.chk.fmt", "Te.uHu"))
        mm = MorbModel(joinpath(tmp, "Te"))
        res = gyrotropic(mm; tasks = (:D0, :Dw, :C, :K, :NOA, :dos),
                         fermi_energies = [2.0, 4.0, 6.0, 8.0, 10.0], freqs = [0.0, 0.05, 0.1],
                         kmesh = (5, 5, 5), smr_width = 0.1, degen_thresh = 0.001,
                         box = [0.2 0.0 0.0; 0.0 0.2 0.0; 0.0 0.0 0.2],
                         box_corner = [0.23333, 0.23333, 0.4], eigval_max = 8.6667)
        @test res.C[1, 1, 1] ≈ 0.361959e1 rtol = 1e-4      # oracle C.dat, E_f=2
        @test res.D[1, 1, 1] ≈ 0.472879e-2 rtol = 1e-4     # oracle D.dat
        @test res.Dw[1, 1, 1, 1] ≈ 0.349534e-3 rtol = 1e-4 # oracle tildeD.dat, ω=0
        @test res.K_orb[1, 1, 1] ≈ -0.825797e-7 rtol = 1e-4
        @test res.NOA_orb[1, 1, 1, 1] ≈ -0.443840e2 rtol = 1e-4
        @test res.dos[1] ≈ 0.277695e-3 rtol = 1e-4
    else
        @test_skip false
    end

    # Ryoo SHC (.sHu/.sIu) + transl_inv_full: Pt frequency scan
    rd = joinpath(REFROOT, "testpostw90_pt_shc_ryoo")
    if isfile(joinpath(rd, "Pt.sHu.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(rd, ("Pt.win", "Pt.eig"),
                    ("Pt.chk.fmt", "Pt.mmn", "Pt.spn", "Pt.sHu", "Pt.sIu"))
        sm = ShcRyooModel(joinpath(tmp, "Pt"))
        shc = shc_freqscan(sm; freqs = [0.0, 7.0], fermi_energy = 18.3823, kmesh = (9, 9, 9),
                           γ = 1, α = 3, β = 2, adaptive = false, smr_width = 0.1,
                           eigval_max = 1000.0)
        @test real(shc[1]) ≈ -2576.0024 rtol = 1e-5        # oracle freqscan rows 1, 71
        @test imag(shc[1]) ≈ 0.0 atol = 1e-9
        @test real(shc[2]) ≈ -67.676662 rtol = 1e-4
        @test imag(shc[2]) ≈ -377.14785 rtol = 1e-4
    else
        @test_skip false
    end
    rt = joinpath(REFROOT, "testpostw90_pt_shc_ryoo_transl_inv")
    if isfile(joinpath(rt, "Pt.sHu.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(rt, ("Pt.win", "Pt.eig"),
                    ("Pt.chk.fmt", "Pt.mmn", "Pt.spn", "Pt.sHu", "Pt.sIu"))
        sm = ShcRyooModel(joinpath(tmp, "Pt"); transl_inv_full = true)
        shc = shc_freqscan(sm; freqs = [0.0], fermi_energy = 18.3823, kmesh = (9, 9, 9),
                           γ = 1, α = 3, β = 2, adaptive = false, smr_width = 0.1,
                           eigval_max = 1000.0)
        @test real(shc[1]) ≈ -2130.4230 rtol = 1e-5        # transl_inv_full + ws_distance
    else
        @test_skip false
    end

    # GaAs ac-SHC: Qiao + frequency scan + scissors shift
    ga = joinpath(REFROOT, "testpostw90_gaas_shc")
    if isfile(joinpath(ga, "GaAs.spn.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(ga, ("GaAs.win", "GaAs.eig"), ("GaAs.chk.fmt", "GaAs.mmn", "GaAs.spn"))
        emax = maximum(read_eig(joinpath(tmp, "GaAs.eig"))) + 0.6667
        sm = ShcModel(joinpath(tmp, "GaAs"); scissors_shift = 1.117, num_valence_bands = 8)
        shc = shc_freqscan(sm; freqs = [0.0, 8.0], fermi_energy = 7.9366, kmesh = (10, 10, 10),
                           adaptive = false, smr_width = 0.05, eigval_max = emax)
        @test real(shc[1]) ≈ -428.20457 rtol = 1e-5        # oracle freqscan rows 1, 801
        @test real(shc[2]) ≈ 404.66629 rtol = 1e-4
    else
        @test_skip false
    end

    # morb with transl_inv_full: Fe (oracle wpout: M_z = 0.0415 vs 0.0431 plain)
    mt = joinpath(REFROOT, "testpostw90_fe_morb_transl_inv")
    if isfile(joinpath(mt, "Fe.uHu.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(mt, ("Fe.win", "Fe.eig"), ("Fe.chk.fmt", "Fe.mmn", "Fe.uHu"))
        mm = MorbModel(joinpath(tmp, "Fe"); transl_inv_full = true)
        M = orbital_magnetisation(mm; fermi_energy = 12.6279, kmesh = (10, 10, 10))
        @test M[3] ≈ 0.0415 atol = 1e-3
        @test abs(M[1]) < 1e-3
    else
        @test_skip false
    end
end

@testset "shift current + kdotp + w90 output extras" begin
    # kdotp: GaAs at L, bands 4-5 (oracle checks only order 0 — orders 1/2 are gauge-sensitive)
    kd = joinpath(REFROOT, "testpostw90_gaas_kdotp")
    if isfile(joinpath(kd, "gaas.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        cp(joinpath(kd, "gaas.win"), joinpath(tmp, "gaas.win"); follow_symlinks = true)
        cp(joinpath(kd, "gaas.eig"), joinpath(tmp, "gaas.eig"); follow_symlinks = true)
        run(pipeline(`bunzip2 -kc $(joinpath(kd, "gaas.chk.fmt.bz2"))`, stdout = joinpath(tmp, "gaas.chk.fmt")))
        bm = BerryModel(joinpath(tmp, "gaas"))
        res = kdotp(bm; kpoint = [0.0, 0.5, 0.0], bands = [4, 5])
        @test real(res.T0[1, 1]) ≈ 6.4827435 atol = 1e-4      # oracle gaas-kdotp_0.dat
        @test real(res.T0[2, 2]) ≈ 8.6209080 atol = 1e-4
        @test abs(res.T0[1, 2]) < 1e-10
    else
        @test_skip false
    end

    # shift current: GaAs σ_xyz, TB phase convention, no eta correction (25³ in the oracle —
    # here the same physics on the benchmark mesh but only spot frequencies via coarse 10³ run)
    sc = joinpath(REFROOT, "testpostw90_gaas_sc_xyz_ws")
    if isfile(joinpath(sc, "gaas.mmn.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("gaas.win", "gaas.eig")
            cp(joinpath(sc, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in ("gaas.chk.fmt", "gaas.mmn")
            if isfile(joinpath(sc, f))               # already decompressed in place
                cp(joinpath(sc, f), joinpath(tmp, f); follow_symlinks = true)
            else
                run(pipeline(`bunzip2 -kc $(joinpath(sc, f * ".bz2"))`, stdout = joinpath(tmp, f)))
            end
        end
        emax = maximum(read_eig(joinpath(tmp, "gaas.eig"))) + 0.6667
        freqs = collect(range(0.0, 10.0; length = 201))
        res = shift_current(joinpath(tmp, "gaas"); fermi_energy = 7.7414, freqs = freqs,
                            kmesh = (10, 10, 10), phase_conv = 1, sc_eta = 0.04,
                            eta_corr = true, eigval_max = emax)
        # oracle gaas-sc_xyz.dat (ws variant, 10³ mesh): σ_xyz at ω = 1.65 eV (row 34)
        ixyz = 6                                    # bc packing: xx yy zz xy xz yz
        @test res.sc[1, ixyz, 34] ≈ 0.72455917e-5 atol = 1e-6
    else
        @test_skip false
    end

    # _r.dat: diamond (reference-exact gauge trajectory via the :w90 optimiser)
    rd = joinpath(REFROOT, "testw90_rmn")
    if isfile(joinpath(rd, "diamond.mmn")) || isfile(joinpath(rd, "diamond.mmn.bz2"))
        tmp = mktempdir()
        for f in ("diamond.win", "diamond.eig", "diamond.amn", "diamond.mmn")
            src = joinpath(rd, f)
            if isfile(src)
                cp(src, joinpath(tmp, f); follow_symlinks = true)
            elseif isfile(src * ".bz2") && Sys.which("bunzip2") !== nothing
                run(pipeline(`bunzip2 -kc $(src * ".bz2")`, stdout = joinpath(tmp, f)))
            end
        end
        model = read_model(joinpath(tmp, "diamond"))
        win = read_win(joinpath(tmp, "diamond.win"))
        result = run_wannier(model, win)
        out = write_rmn(joinpath(tmp, "diamond"), model, result.Mrot)
        lines = readlines(out)
        @test strip(lines[2]) == "4"
        @test strip(lines[3]) == "93"
        @test length(lines) == 3 + 93 * 16
    else
        @test_skip false
    end
end

@testset "TB input + disentangle/wannierise options + FermiSurfer" begin
    # dis_spheres: LaVO3 (Ω_total 7.508128029, Ω_I 7.457463597) — plain inputs, no bz2
    dd = joinpath(REFROOT, "testw90_lavo3_dissphere")
    if isfile(joinpath(dd, "LaVO3.win"))
        m = read_model(joinpath(dd, "LaVO3")); w = read_win(joinpath(dd, "LaVO3.win"))
        r = run_wannier(m, w)
        @test r.spread.Ω ≈ 7.508128029 atol = 1e-5
        @test r.omega_I ≈ 7.457463597 atol = 1e-5
    else
        @test_skip false
    end

    # PDWF (dis_froz_proj + dis_proj_min/max): graphene (Ω_total 15.803349910)
    pd = joinpath(REFROOT, "testw90_graphene_pdwf")
    if isfile(joinpath(pd, "graphene.win"))
        m = read_model(joinpath(pd, "graphene")); w = read_win(joinpath(pd, "graphene.win"))
        r = run_wannier(m, w)
        @test r.spread.Ω ≈ 15.803349910 atol = 1e-5
        @test r.omega_I ≈ 7.962090079 atol = 1e-5
    else
        @test_skip false
    end

    # guiding_centres + select_projections: silicon (22.738496505 Bohr² = Ω/bohr²)
    gs = joinpath(REFROOT, "testw90_guidingcentre_selectproj")
    if isfile(joinpath(gs, "silicon.win"))
        m = read_model(joinpath(gs, "silicon")); w = read_win(joinpath(gs, "silicon.win"))
        r = run_wannier(m, w)
        bohr = 0.52917720859
        @test r.spread.Ω / bohr^2 ≈ 22.738496505 atol = 1e-4
        @test all(s -> isapprox(s / bohr^2, 5.68462413; atol = 1e-4), r.spread.spreads)
    else
        @test_skip false
    end

    # preconditioned CG: GaAs reaches the same minimum (4.466880976)
    pc = joinpath(REFROOT, "testw90_precond_1")
    if isfile(joinpath(pc, "gaas1.win"))
        m = read_model(joinpath(pc, "gaas1")); w = read_win(joinpath(pc, "gaas1.win"))
        r = run_wannier(m, w)
        @test r.spread.Ω ≈ 4.466880976 atol = 1e-6
    else
        @test_skip false
    end

    # TB-model input round-trip: build a BerryModel, write/read a _tb.dat, reproduce it
    fe = joinpath(REFROOT, "testpostw90_fe_ahc")
    if isfile(joinpath(fe, "Fe.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("Fe.win", "Fe.eig")
            cp(joinpath(fe, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in ("Fe.chk.fmt", "Fe.mmn")
            run(pipeline(`bunzip2 -kc $(joinpath(fe, f * ".bz2"))`, stdout = joinpath(tmp, f)))
        end
        bm = BerryModel(joinpath(tmp, "Fe"))
        pos = reshape(bm.Ar, size(bm.Ar, 1), size(bm.Ar, 2), size(bm.Ar, 3), 3)
        tbf = joinpath(tmp, "RT_tb.dat")
        write_tb(tbf, bm.lattice, size(bm.Hr, 1), bm.irvec, bm.ndegen, bm.Hr; pos = pos)
        tb = tb_model(tbf)
        @test tb.irvec == bm.irvec
        # bands identical to the source model
        k = SVector(0.1, 0.2, 0.3)
        E1, _ = eig_deleig(bm, k; deriv = false)
        E2, _ = eig_deleig(tb, k; deriv = false)
        @test maximum(abs.(E1 .- E2)) < 1e-5
        # AHC reproduced from the TB model alone
        a1 = anomalous_hall(bm; fermi_energy = 12.6279, kmesh = (8, 8, 8))
        a2 = anomalous_hall(tb; fermi_energy = 12.6279, kmesh = (8, 8, 8))
        @test a2[3] ≈ a1[3] rtol = 1e-4

        # FermiSurfer .frmsf: layout + energies match interpolation
        E, _ = tabulate_3d(tb; mesh = (4, 4, 4))
        @test size(E) == (size(bm.Hr, 1), 4, 4, 4)
        ffile = joinpath(tmp, "t.frmsf")
        write_frmsf(ffile, tb.lattice, E; fermi_energy = 12.0)
        lines = readlines(ffile)
        @test length(lines) == 6 + size(E, 1) * 64      # header(6) + nband*4^3 energies
        Ek, _ = eig_deleig(tb, [0.0, 0.0, 0.0]; deriv = false)
        @test E[1, 1, 1, 1] ≈ Ek[1] atol = 1e-8
    else
        @test_skip false
    end
end

@testset "tetrahedron SHC" begin
    # Pt tetrahedron Fermi scan (Ghim-Park + Kawamura), qiao operators — 10^3 mesh (~7s)
    td = joinpath(REFROOT, "testpostw90_pt_tetra_shcfermi")
    if isfile(joinpath(td, "Pt.spn.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("Pt.win", "Pt.eig")
            cp(joinpath(td, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in ("Pt.chk.fmt", "Pt.mmn", "Pt.spn")
            run(pipeline(`bunzip2 -kc $(joinpath(td, f * ".bz2"))`, stdout = joinpath(tmp, f)))
        end
        sm = ShcModel(joinpath(tmp, "Pt"))
        efs = collect(13.0:0.5:23.0)
        shc = shc_tetra(sm; kmesh = (10, 10, 10), fermi_energies = efs, γ = 3, α = 1, β = 2,
                        cutoff = 1e-1, avoid_deg = 3e-4)
        @test shc[1] ≈ -1047.7881 atol = 0.1        # oracle fermiscan rows 1, 11, 21
        @test shc[11] ≈ 1863.5104 atol = 0.1
        @test shc[21] ≈ 522.79635 atol = 0.1
    else
        @test_skip false
    end
end

