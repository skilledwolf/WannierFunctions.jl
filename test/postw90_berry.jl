# Split from the original monolithic runtests.jl — included by runtests.jl (common.jl
# provides REFROOT/DATAROOT, the shared GaAs/diamond models, and helpers).
@testset "Berry curvature / AHC" begin
    # Time-reversal-symmetric insulator: AHC must vanish (CI-safe physics invariant).
    if DIAMOND_MODEL !== nothing
        win = read_win(DIAMOND_SEED * ".win")
        res = run_wannier(DIAMOND_MODEL, win)
        chk = Checkpoint(DIAMOND_MODEL, win, res)
        bm = BerryModel(chk, DIAMOND_MODEL.eig, DIAMOND_MODEL.bvectors,
                        DIAMOND_MODEL.kgrid, DIAMOND_MODEL.lattice)
        ahc = anomalous_hall(bm; fermi_energy = 10.0, kmesh = (4, 4, 4))  # gap: 4 bands filled
        # Time reversal ⇒ AHC = 0 exactly only for the full BZ integral; a coarse 4³ sum
        # leaves cancellation residuals of a few µS/cm (vs Fe's ~1222 S/cm signal).
        @test maximum(abs.(ahc)) < 1e-4
    else
        @test_skip false
    end

    # Fe (spinor, disentangled, magnetic): match the postw90 benchmark on its 10³ mesh.
    fed = joinpath(REFROOT, "testpostw90_fe_ahc")
    if isfile(joinpath(fed, "Fe.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("Fe.win", "Fe.eig")
            cp(joinpath(fed, f), joinpath(tmp, f); follow_symlinks = true)
        end
        run(pipeline(`bunzip2 -kc $(joinpath(fed, "Fe.chk.fmt.bz2"))`,
                     stdout = joinpath(tmp, "Fe.chk.fmt")))
        run(pipeline(`bunzip2 -kc $(joinpath(fed, "Fe.mmn.bz2"))`,
                     stdout = joinpath(tmp, "Fe.mmn")))
        bm = BerryModel(joinpath(tmp, "Fe"))
        ahc = anomalous_hall(bm; fermi_energy = 12.6279, kmesh = (10, 10, 10))
        # harness tolerances for AHC: abs 1e-3, rel 2e-3
        @test ahc[1] ≈ 0.0334 atol = 1e-3
        @test ahc[2] ≈ 0.0572 atol = 1e-3
        @test ahc[3] ≈ 1222.1510 rtol = 2e-3
    else
        @test_skip false
    end
end

# =========================================================================
# (14) ADAPTIVE AHC FERMISCAN + GENINTERP
# =========================================================================
@testset "Adaptive fermiscan + geninterp" begin
    fed = joinpath(REFROOT, "testpostw90_fe_ahc_adaptandfermi")
    if isfile(joinpath(fed, "Fe.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("Fe.win", "Fe.eig")
            cp(joinpath(fed, f), joinpath(tmp, f); follow_symlinks = true)
        end
        run(pipeline(`bunzip2 -kc $(joinpath(fed, "Fe.chk.fmt.bz2"))`, stdout = joinpath(tmp, "Fe.chk.fmt")))
        run(pipeline(`bunzip2 -kc $(joinpath(fed, "Fe.mmn.bz2"))`, stdout = joinpath(tmp, "Fe.mmn")))
        bm = BerryModel(joinpath(tmp, "Fe"))
        # two Fermi levels, small mesh, adaptive on: consistency + reproducibility
        out = ahc_fermiscan(bm; fermi_energies = [12.4279, 12.6279], kmesh = (4, 4, 4),
                            adpt_kmesh = 3, adpt_thresh = 10.0)
        @test size(out) == (3, 2)
        @test all(isfinite, out)
        # single-level entry point must agree with the scan's column
        one = anomalous_hall(bm; fermi_energy = 12.6279, kmesh = (4, 4, 4))
        scan = ahc_fermiscan(bm; fermi_energies = [12.6279], kmesh = (4, 4, 4))
        @test maximum(abs.(one .- scan[:, 1])) < 1e-10
    else
        @test_skip false
    end

    gid = joinpath(REFROOT, "testpostw90_si_geninterp")
    if isfile(joinpath(gid, "silicon.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("silicon.win", "silicon.eig", "silicon_geninterp.kpt")
            cp(joinpath(gid, f), joinpath(tmp, f); follow_symlinks = true)
        end
        run(pipeline(`bunzip2 -kc $(joinpath(gid, "silicon.chk.fmt.bz2"))`, stdout = joinpath(tmp, "silicon.chk.fmt")))
        run(pipeline(`bunzip2 -kc $(joinpath(gid, "silicon.mmn.bz2"))`, stdout = joinpath(tmp, "silicon.mmn")))
        bm = BerryModel(joinpath(tmp, "silicon"))
        out = geninterp(bm, joinpath(tmp, "silicon"))
        @test isfile(out)
        rows = [parse.(Float64, split(l)) for l in eachline(out) if !startswith(l, "#")]
        @test length(rows) == 24                              # 3 k-points × 8 bands
        # Time reversal relates E(k) and E(−k); for a complex (not reality-constrained) Wannier
        # gauge the *interpolated* values obey it only to interpolation accuracy off-grid
        # (~1e-5 eV here), exactly on the original grid.
        E1, dE1 = eig_deleig(bm, [0.1, 0.2, 0.3])
        E2, dE2 = eig_deleig(bm, [-0.1, -0.2, -0.3])
        @test maximum(abs.(E1 .- E2)) < 1e-3
        @test maximum(abs.(dE1 .+ dE2)) < 5e-2
    else
        @test_skip false
    end
end

# =========================================================================
# (15) KUBO OPTICAL CONDUCTIVITY + ORBITAL MAGNETISATION (Fe oracle values)
# =========================================================================
@testset "Kubo + orbital magnetisation" begin
    kd = joinpath(REFROOT, "testpostw90_fe_kubo_Axy")
    if isfile(joinpath(kd, "Fe.uHu.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("Fe.win", "Fe.eig")
            cp(joinpath(kd, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in ("Fe.chk.fmt", "Fe.mmn", "Fe.uHu")
            run(pipeline(`bunzip2 -kc $(joinpath(kd, f * ".bz2"))`, stdout = joinpath(tmp, f)))
        end
        bm = BerryModel(joinpath(tmp, "Fe"))
        res = optical_conductivity(bm; fermi_energy = 12.6279, kmesh = (10, 10, 10),
                                   freqs = 0.0:0.5:7.0, eigval_max = 30.0 + 2.0 / 3.0)
        # postw90 oracle values (Fe-kubo_A_xy.dat / Fe-jdos.dat, 10³ mesh)
        @test real(kubo_A(res, 1, 2, 1)) ≈ 304.66638 atol = 1e-2      # ω=0, Re σ_A,xy
        @test real(kubo_A(res, 1, 2, 2)) ≈ 119.20666 atol = 1e-2      # ω=0.5
        @test imag(kubo_A(res, 1, 2, 2)) ≈ 156.15522 atol = 1e-2
        @test abs(imag(kubo_A(res, 1, 2, 1))) < 1e-8                   # Im σ_A(ω=0) = 0
        # morb: same staged data
        mm = MorbModel(joinpath(tmp, "Fe"))
        M = orbital_magnetisation(mm; fermi_energy = 12.6279, kmesh = (10, 10, 10))
        @test abs(M[1]) < 1e-3
        @test abs(M[2]) < 1e-3
        @test M[3] ≈ 0.0431 atol = 1e-3                                # benchmark, μ_B/cell
    else
        @test_skip false
    end
end

# =========================================================================
# (16) DOS, BOLTZWANN, SHC, KSLICE (oracle benchmark values)
# =========================================================================
@testset "DOS + BoltzWann + SHC + kslice" begin
    # DOS: copper (chk+eig only — exercises the H(R)-only BerryModel)
    dd = joinpath(REFROOT, "testpostw90_example04_dos")
    if isfile(joinpath(dd, "copper.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("copper.win", "copper.eig")
            cp(joinpath(dd, f), joinpath(tmp, f); follow_symlinks = true)
        end
        run(pipeline(`bunzip2 -kc $(joinpath(dd, "copper.chk.fmt.bz2"))`, stdout = joinpath(tmp, "copper.chk.fmt")))
        bm = BerryModel(joinpath(tmp, "copper"))
        es, d = density_of_states(bm; energies = 8.0:0.25:10.0, kmesh = (10, 10, 10))
        @test d[6] ≈ 4.7047 atol = 1e-3          # oracle copper-dos.dat at E=9.25 (peak)
        @test d[1] ≈ 1.9268 atol = 1e-3
    else
        @test_skip false
    end

    # BoltzWann: silicon (ws_distance=true path)
    bd = joinpath(REFROOT, "testpostw90_boltzwann")
    if isfile(joinpath(bd, "silicon.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("silicon.win", "silicon.eig")
            cp(joinpath(bd, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in ("silicon.chk.fmt", "silicon.mmn")
            run(pipeline(`bunzip2 -kc $(joinpath(bd, f * ".bz2"))`, stdout = joinpath(tmp, f)))
        end
        bm = BerryModel(joinpath(tmp, "silicon"))
        @test bm.wsdist !== nothing               # use_ws_distance honoured
        eig = read_eig(joinpath(tmp, "silicon.eig"))
        res = boltzwann(bm; kmesh = (20, 20, 20), relax_time = 10.0, mus = [5.0],
                        temps = [300.0], win = (minimum(eig), maximum(eig)))
        # oracle silicon_elcond.dat (harness tol: abs 10, rel 1e-4)
        @test res.elcond[1, 1, 1, 1] ≈ 6.504319e6 rtol = 1e-4
        @test res.elcond[3, 3, 1, 1] ≈ res.elcond[1, 1, 1, 1] rtol = 1e-2   # cubic symmetry
    else
        @test_skip false
    end

    # SHC: Pt (spinor, Qiao method, ws_distance default-true)
    sd = joinpath(REFROOT, "testpostw90_pt_shc")
    if isfile(joinpath(sd, "Pt.spn.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("Pt.win", "Pt.eig")
            cp(joinpath(sd, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in ("Pt.chk.fmt", "Pt.mmn", "Pt.spn")
            run(pipeline(`bunzip2 -kc $(joinpath(sd, f * ".bz2"))`, stdout = joinpath(tmp, f)))
        end
        sm = ShcModel(joinpath(tmp, "Pt"))
        # 3 spot Fermi energies on a smaller mesh would differ from the benchmark; use the
        # benchmark mesh (15³) at 3 levels only — cheap (per-level reuse of k-data).
        out = shc_fermiscan(sm; fermi_energies = [6.0, 17.9, 26.0], kmesh = (15, 15, 15),
                            eigval_max = 30.0 + 2.0 / 3.0)
        @test abs(out[1]) < 1e-6                   # empty bands → zero
        @test out[2] ≈ 1812.32 rtol = 1e-4         # oracle fermiscan at E_F=17.9
        @test out[3] ≈ 244.09181 rtol = 1e-4       # at E_F=26.0
    else
        @test_skip false
    end

    # kslice: Fe curv+bands (5×5 slice)
    kd2 = joinpath(REFROOT, "testpostw90_fe_kslicecurv")
    if isfile(joinpath(kd2, "Fe.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = mktempdir()
        for f in ("Fe.win", "Fe.eig")
            cp(joinpath(kd2, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in ("Fe.chk.fmt", "Fe.mmn")
            run(pipeline(`bunzip2 -kc $(joinpath(kd2, f * ".bz2"))`, stdout = joinpath(tmp, f)))
        end
        bm = BerryModel(joinpath(tmp, "Fe"))
        ks = kslice(bm; corner = [0.0, 0.0, 0.0], b1 = [0.5, -0.5, -0.5],
                    b2 = [0.5, 0.5, 0.5], mesh = (5, 5),
                    fermi_energy = 12.6279, tasks = (:bands, :curv))
        @test length(ks.kpts) == 36
        @test ks.bands[1, 1] ≈ 4.43414 atol = 1e-4  # oracle Fe-kslice-bands.dat first entries
        @test ks.bands[2, 1] ≈ 4.5557612 atol = 1e-4
        @test ks.curv[3, 2] ≈ -5.6272413 atol = 1e-4    # oracle curv z at point 2
    else
        @test_skip false
    end
end

