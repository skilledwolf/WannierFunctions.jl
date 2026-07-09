# Split from the original monolithic runtests.jl — included by runtests.jl (common.jl
# provides REFROOT/DATAROOT, the shared GaAs/diamond models, and helpers).
@testset "symmetry: BZ reduction + irreducible AHC" begin
    # k-mesh reduction: eigenvalue multiset on the irreducible wedge must reproduce the full mesh
    gs = joinpath(REFROOT, "testw90_example21_As_sp")
    if isfile(joinpath(gs, "GaAs.sym"))
        m = read_model(joinpath(gs, "GaAs")); w = read_win(joinpath(gs, "GaAs.win"))
        r = run_wannier(m, w)
        irvec, ndeg = wigner_seitz(m.lattice, m.kgrid.mp_grid)
        Hr, _ = build_hr(r.U, r.eig_interp, m.kgrid, irvec)
        bm = BerryModel(m.lattice, irvec, ndeg, Hr)
        sym = read_sym(joinpath(gs, "GaAs.sym"))
        @test nsym(sym) == 24
        reps, wts, _ = irreducible_kmesh((6, 6, 6), sym; kaction = :cart, lattice = m.lattice)
        @test sum(wts) == 216
        full = Float64[]
        for i in 0:5, j in 0:5, k in 0:5
            E, _ = eig_deleig(bm, [i / 6, j / 6, k / 6]; deriv = false); append!(full, E)
        end
        irr = Float64[]
        for (ki, kf) in enumerate(reps)
            E, _ = eig_deleig(bm, kf; deriv = false)
            for _ in 1:wts[ki], e in E
                push!(irr, e)
            end
        end
        @test maximum(abs.(sort(irr) .- sort(full))) < 1e-6   # exact multiset match
    else
        @test_skip false
    end

    # irreducible-wedge orbital magnetisation + DOS == full BZ (Fe; magnetic subgroup)
    fm = joinpath(REFROOT, "testpostw90_fe_morb")
    if isfile(joinpath(fm, "Fe.uHu.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("Fe.win", "Fe.eig")
            cp(joinpath(fm, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in ("Fe.chk.fmt", "Fe.mmn", "Fe.uHu")
            if isfile(joinpath(fm, f))
                cp(joinpath(fm, f), joinpath(tmp, f); follow_symlinks = true)
            else
                run(pipeline(`bunzip2 -kc $(joinpath(fm, f * ".bz2"))`, stdout = joinpath(tmp, f)))
            end
        end
        mm = MorbModel(joinpath(tmp, "Fe"))
        ef = 12.6279
        symg = WannierFunctions._pseudovector_subgroup(mm.bm, cubic_point_group(), (4, 4, 4), ef)
        @test nsym(symg) == 8
        M_full = orbital_magnetisation(mm; fermi_energy = ef, kmesh = (8, 8, 8))
        M_irr, ninfo = orbital_magnetisation_sym(mm, symg; fermi_energy = ef, kmesh = (8, 8, 8))
        @test ninfo[1] < ninfo[2] ÷ 6                 # genuine wedge reduction
        @test maximum(abs.(M_full .- M_irr)) < 1e-5
        es = collect(range(10.0, 15.0; length = 21))
        _, d_full = density_of_states(mm.bm; energies = es, kmesh = (8, 8, 8),
                                      adaptive = false, smr_width = 0.3)
        _, d_irr, _ = density_of_states_sym(mm.bm, symg; energies = es, kmesh = (8, 8, 8),
                                            adaptive = false, smr_width = 0.3)
        @test maximum(abs.(d_full .- d_irr)) < 1e-10  # scalar wedge sum is exact
    else
        @test_skip false
    end

    # irreducible-wedge symmetrised AHC == full-BZ AHC (Fe; magnetic subgroup filtered from Oₕ)
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
        Ef = 12.6279
        full = anomalous_hall(bm; fermi_energy = Ef, kmesh = (10, 10, 10))
        mg = WannierFunctions._pseudovector_subgroup(bm, cubic_point_group(), (10, 10, 10), Ef; tol = 1e-4)
        @test nsym(mg) >= 1
        sym_ahc, (nirr, nfull) = anomalous_hall_sym(bm, mg; fermi_energy = Ef, kmesh = (10, 10, 10))
        @test nirr < nfull                                      # genuine reduction
        @test maximum(abs.(collect(sym_ahc) .- collect(full))) < 1e-3   # matches full-BZ
    else
        @test_skip false
    end
end

@testset "SLWF+C selective localisation" begin
    # example26: GaAs, slwf_num=1, constrain WF1 to (0.25,0.25,0.25); objective Ω_C=1.634087565
    ed = joinpath(REFROOT, "testw90_example26")
    if isfile(joinpath(ed, "gaas.win")) && isfile(joinpath(ed, "gaas.amn"))
        m = read_model(joinpath(ed, "gaas")); w = read_win(joinpath(ed, "gaas.win"))
        # objective on the reference gauge (if a chk is present), else via the optimiser
        r = run_wannier(m, w)
        @test r.spread.Ω ≈ 1.634087565 atol = 1e-6     # Ω_C (compute_spread returns it for SLWF)
        @test r.converged
    else
        @test_skip false
    end
end

@testset "symmetry-adapted WFs (site_symmetry, localisation)" begin
    # GaAs example21: isolated, site_symmetry; symmetric Ω=10.136492662, U symmetry-adapted
    gd = joinpath(REFROOT, "testw90_example21_As_sp")
    if isfile(joinpath(gd, "GaAs.dmn")) && isfile(joinpath(gd, "GaAs.amn"))
        m = read_model(joinpath(gd, "GaAs"))
        ss = read_dmn(joinpath(gd, "GaAs.dmn"), m.num_bands, m.num_wann)
        @test ss.nsym == 24 && ss.nkptirr == 10
        res = wannierise(m; num_iter = 2000, algorithm = :rcg, sitesym = ss)
        @test res.spread.Ω ≈ 10.136492662 atol = 1e-5
        @test res.converged
        # U must be symmetry-adapted: U(Rk) = d_band·U(k)·d_wann†
        dev = 0.0
        for ir in 1:ss.nkptirr, is in 1:ss.nsym
            ik = ss.ir2ik[ir]; irk = ss.kptsym[is, ir]
            pred = ss.d_band[:, :, is, ir] * res.U[:, :, ik] * ss.d_wann[:, :, is, ir]'
            dev = max(dev, maximum(abs.(res.U[:, :, irk] - pred)))
        end
        @test dev < 1e-6
    else
        @test_skip false
    end
end

@testset "symmetry-adapted disentanglement (dis_extract_symmetry, H3S)" begin
    # testw90_disentanglement_sawfs: 20 bands → 12 WFs with site_symmetry. The benchmark runs
    # exactly 10 constrained Ω_I iterations (mix_ratio 0.2), so the trajectory itself is the
    # oracle; the localisation Ω is checked against a converged wannier90.x run (num_iter=5000,
    # Ω_total = 6.301957261).
    gd = joinpath(REFROOT, "testw90_disentanglement_sawfs")
    if isfile(joinpath(gd, "H3S.dmn")) && isfile(joinpath(gd, "H3S.mmn"))
        m = read_model(joinpath(gd, "H3S"))
        win = read_win(joinpath(gd, "H3S.win"))
        ss = read_dmn(joinpath(gd, "H3S.dmn"), m.num_bands, m.num_wann; eps = 1e-8)
        dis = disentangle(m, win; sitesym = ss)
        # Ω_I(i-1)/Ω_I(i) trajectory, all 10 iterations of the benchmark
        bench = [(3.61187069, 3.46062472), (3.58098446, 3.44855662), (3.55399227, 3.43895830),
                 (3.53058214, 3.43131533), (3.51040746, 3.42522424), (3.49311461, 3.42036598),
                 (3.47836056, 3.41648814), (3.46582260, 3.41338795), (3.45520495, 3.41090822),
                 (3.44624098, 3.40892357)]
        for (i, (b1, b2)) in enumerate(bench)
            @test dis.omega_I_trace[i][2] ≈ b1 atol = 5e-7
            @test dis.omega_I_trace[i][3] ≈ b2 atol = 5e-7
        end
        @test dis.omega_I ≈ 3.408923571 atol = 1e-8
        # Full pipeline: Ω_I preserved, converged symmetric Ω_total, U symmetry-adapted
        _, _, res = run_wannier(joinpath(gd, "H3S"))
        @test res.omega_I ≈ 3.408923571 atol = 1e-8
        @test res.spread.Ω ≈ 6.301957261 atol = 1e-7
        @test res.converged
        dev = 0.0
        for ir in 1:ss.nkptirr, is in 1:ss.nsym
            ik = ss.ir2ik[ir]; irk = ss.kptsym[is, ir]
            D = ss.d_wann[:, :, is, ir]
            dev = max(dev, maximum(abs.(res.U[:, :, irk] - D * res.U[:, :, ik] * D')))
        end
        @test dev < 1e-9
    else
        @test_skip false
    end
end

