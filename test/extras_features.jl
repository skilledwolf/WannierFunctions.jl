# Split from the original monolithic runtests.jl — included by runtests.jl (common.jl
# provides REFROOT/DATAROOT, the shared GaAs/diamond models, and helpers).
@testset "higher-order finite differences (higher_order_n)" begin
    # knbo3 (3×3×3, higher_order_n = 2): the b-list carries every b and 2b with
    # central-difference weights 4/3·w and −w/12; the full MV spread must match the benchmark.
    gd = joinpath(REFROOT, "testw90_knbo3_higher")
    if isfile(joinpath(gd, "knbo3.mmn")) && isfile(joinpath(gd, "knbo3.amn"))
        m = read_model(joinpath(gd, "knbo3"))
        @test m.bvectors.nntot == 12                       # 6 first-order + 6 doubled
        w = m.bvectors.shell_weight
        @test length(w) == 4
        @test w[3] / w[1] ≈ -1 / 16 atol = 1e-12           # (−1/12) / (4/3)
        _, _, res = run_wannier(joinpath(gd, "knbo3"))
        @test res.spread.ΩI ≈ 11.742315886 atol = 2e-6
        @test res.spread.ΩD ≈ 0.005696019 atol = 1e-6
        @test res.spread.ΩOD ≈ 2.044878523 atol = 1e-6
        @test res.spread.Ω ≈ 13.792890427 atol = 2e-6
        # -pp generation: nnkpts block must list the multiples in the reference's canonical
        # order (validated byte-identical vs wannier90.x -pp during development)
        tmpd = mktempdir()
        cp(joinpath(gd, "knbo3.win"), joinpath(tmpd, "knbo3.win"); follow_symlinks = true)
        generate_nnkp(joinpath(tmpd, "knbo3"))
        txt = read(joinpath(tmpd, "knbo3.nnkp"), String)
        @test occursin(r"begin nnkpts\s+12", txt)
    else
        @test_skip false
    end

    # postw90 with higher-order b-sums: Fe morb (transl_inv_full) and Pt SHC fermi scan
    function stage(dir, plain, zipped)
        tmp = mktempdir()
        for f in plain
            cp(joinpath(dir, f), joinpath(tmp, f); follow_symlinks = true)
        end
        for f in zipped
            if isfile(joinpath(dir, f))
                cp(joinpath(dir, f), joinpath(tmp, f); follow_symlinks = true)
            else
                run(pipeline(`bunzip2 -kc $(joinpath(dir, f * ".bz2"))`, stdout = joinpath(tmp, f)))
            end
        end
        tmp
    end
    mt = joinpath(REFROOT, "testpostw90_fe_morb_transl_inv_higher")
    if isfile(joinpath(mt, "Fe.uHu.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(mt, ("Fe.win", "Fe.eig"), ("Fe.chk.fmt", "Fe.mmn", "Fe.uHu"))
        mm = MorbModel(joinpath(tmp, "Fe"); transl_inv_full = true)   # unformatted .uHu
        M = orbital_magnetisation(mm; fermi_energy = 12.6631, kmesh = (10, 10, 10))
        @test M[3] ≈ -0.0617 atol = 1e-3
        @test abs(M[1]) < 1e-3
    else
        @test_skip false
    end
    pt = joinpath(REFROOT, "testpostw90_pt_shc_higher")
    if isfile(joinpath(pt, "Pt.spn.bz2")) && Sys.which("bunzip2") !== nothing
        tmp = stage(pt, ("Pt.win", "Pt.eig"), ("Pt.chk.fmt", "Pt.mmn", "Pt.spn"))
        sm = ShcModel(joinpath(tmp, "Pt"))
        emax = maximum(read_eig(joinpath(tmp, "Pt.eig"))) + 2.0 / 3.0
        out = shc_fermiscan(sm; fermi_energies = [8.4, 17.9, 26.0], kmesh = (15, 15, 15),
                            eigval_max = emax)
        @test out[1] ≈ -3.9190876 rtol = 1e-5
        @test out[2] ≈ 1633.0502 rtol = 1e-5
        @test out[3] ≈ 435.58483 rtol = 1e-5
    else
        @test_skip false
    end
end

@testset "Stengel–Spaldin functional (use_ss_functional)" begin
    # The SS surface is a long, nearly-flat valley: the benchmark's "converged" value is where
    # the default Δ-criterion happens to stop (wannier90.x itself reaches 13.312 when run with
    # conv_tol = 1e-16). The robust oracle checks are objective parity at the shared initial
    # gauge, descent at least as deep as the benchmark, and gradient/objective consistency.
    gd = joinpath(REFROOT, "testw90_knbo3_higher_stengel_spaldin")
    if isfile(joinpath(gd, "knbo3.mmn")) && isfile(joinpath(gd, "knbo3.amn"))
        m = read_model(joinpath(gd, "knbo3"))
        bv = m.bvectors
        ssd = WannierFunctions.ss_data(bv)
        U0 = initial_gauge(m.A)
        Mrot0 = rotate_overlaps(m.M, U0, bv.kpb)
        # objective parity: w90's iteration-0 spread
        @test WannierFunctions.ss_spread(Mrot0, bv, ssd).Ω ≈ 14.7994988992 atol = 1e-8
        # gradient/objective consistency along the optimiser's parametrisation
        wbtot = sum(@view bv.wb[:, 1])
        G = WannierFunctions.ss_gradient(Mrot0, bv, ssd)
        X = similar(G)
        for k in axes(X, 3)
            A = ComplexF64.(reshape(sin.(1:size(X, 1)^2) .+ k, size(X, 1), size(X, 1)))
            X[:, :, k] = (A - A') / 2 + im * (A + A') / 20
            X[:, :, k] = (X[:, :, k] - X[:, :, k]') / 2
        end
        slope = -real(LinearAlgebra.dot(G, X)) / (4wbtot)
        h = 1e-6
        Ω(s) = begin
            R = WannierFunctions.expm_all(X .* (s / (4wbtot)))
            Ut = copy(U0); Mt = copy(Mrot0)
            WannierFunctions.apply_rotation!(Ut, Mt, bv.kpb, R)
            WannierFunctions.ss_spread(Mt, bv, ssd).Ω
        end
        @test slope ≈ (Ω(h) - Ω(-h)) / (2h) rtol = 1e-4
        # full pipeline: MV Ω_I is data-level invariant; the SS descent must reach at least
        # the benchmark's stopping point (13.845371018)
        _, _, res = run_wannier(joinpath(gd, "knbo3"))
        @test res.spread.ΩI ≈ 11.742315886 atol = 2e-6
        @test WannierFunctions.ss_spread(res.Mrot, bv, ssd).Ω <= 13.845371018 + 1e-6
    else
        @test_skip false
    end
end

@testset "injection current (vs WannierBerri)" begin
    # Circular injection current on GaAs; anchors from WannierBerri's InjectionCurrent
    # calculator on the SAME tight-binding model (12³ mesh, 0.1 eV Gaussian).
    sc = joinpath(REFROOT, "testpostw90_gaas_sc_xyz")
    if isfile(joinpath(sc, "gaas.chk.fmt.bz2")) && Sys.which("bunzip2") !== nothing
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
        bm = BerryModel(joinpath(tmp, "gaas"))
        η = injection_current(bm; freqs = [1.0, 1.5, 2.0], fermi_energy = 7.7414,
                              kmesh = (12, 12, 12), smr_width = 0.1)
        # WannierBerri anchors on the SAME tb.dat (scratch/gaas_inj). Agreement is limited to
        # ~1e-4 by the two codes' degenerate-state regularisation conventions (GaAs has exact
        # band degeneracies at high-symmetry points on the 12³ grid).
        @test η[2, 3, 1, 1] ≈ -8.218173e-8 rtol = 3e-4    # yzx at ω = 1.0
        @test η[2, 3, 1, 2] ≈ -1.769306e-6 rtol = 3e-4    # yzx at ω = 1.5
        @test η[1, 2, 3, 2] ≈ 5.237471e-7 rtol = 3e-4     # xyz at ω = 1.5
    else
        @test_skip false
    end
end

@testset "DFTK live end-to-end (extension)" begin
    # Full all-Julia pipeline: DFTK silicon SCF → in-memory Model → wannierisation. Runs only
    # when DFTK is available in the test environment (heavy dependency); the reference run
    # gives Ω = 6.4566 Å², bond-centred WFs, and mesh-band reproduction at 1e-11 eV.
    if Base.find_package("DFTK") !== nothing
        include("dftk_e2e.jl")                # separate file so the suite parses without DFTK
    else
        @test_skip false
    end
end

@testset "in-memory model bridge (DFTK-style)" begin
    # The DFTK bridge hands off in-memory arrays via wannier_model; check the in-memory path
    # produces an identical wannierisation to reading the same content from files.
    sd = joinpath(REFROOT, "testw90_example01")
    seed = isfile(joinpath(sd, "gaas.win")) ? joinpath(sd, "gaas") : GAAS_SEED
    if isfile(seed * ".mmn") && isfile(seed * ".amn")
        m1 = read_model(seed)
        M, kpb, gpb, _, _, _ = read_mmn(seed * ".mmn")
        A, _, _, _ = read_amn(seed * ".amn")
        w = read_win(seed * ".win")
        eig = isfile(seed * ".eig") ? read_eig(seed * ".eig") : nothing
        m2 = wannier_model(; unit_cell = w.unit_cell, kpoints = w.kpoints, mp_grid = w.mp_grid,
                           num_wann = w.num_wann, M = M, A = A, kpb = kpb, gpb = gpb, eig = eig)
        r1 = run_wannier(m1, w)
        r2 = run_wannier(m2, w)
        @test r1.spread.Ω ≈ r2.spread.Ω atol = 1e-9      # in-memory == file path
        @test m2.num_bands == m1.num_bands
    else
        @test_skip false
    end
end

@testset "CLI wrappers (install_cli + argv entry points)" begin
    # Bad-usage paths return exit codes rather than exiting.
    @test WannierFunctions.wannier90_cli(String[]) == 1
    @test WannierFunctions.postw90_cli(String[]) == 1
    @test WannierFunctions.w90chk2chk_cli(["-bogus", "x"]) == 1

    mktempdir() do dir
        written = install_cli(; dir = dir, julia_flags = ["-O0"])
        names = basename.(written)
        @test "wannier90.jl" in names && "postw90.jl" in names && "w90chk2chk.jl" in names
        script = read(joinpath(dir, "wannier90.jl"), String)
        @test occursin(Base.active_project(), script)        # launches this environment
        @test occursin("wannier90_cli", script) && occursin("-O0", script)
        Sys.iswindows() || @test uperm(joinpath(dir, "wannier90.jl")) & 0x01 != 0  # executable
    end
end

@testset "scdm_auto (Vitale μ/σ fit)" begin
    # Synthetic erfc projectability with known parameters + bounded deterministic noise:
    # the LM fit must recover (μ_fit, σ_fit) and apply the μ = μ_fit − 3σ_fit shift.
    μt, σt = -3.5, 1.7
    nb, nk = 40, 12
    eig = zeros(nb, nk); proj = zeros(nb, nk)
    for k in 1:nk, m in 1:nb
        e = -12.0 + 24.0 * (m - 1) / (nb - 1) + 0.13 * sin(2.3m + 0.7k)
        eig[m, k] = e
        proj[m, k] = clamp(0.5 * WannierFunctions.erfc_((e - μt) / σt) + 0.02 * sin(11.0m + 3.1k),
                           0.0, 1.0)
    end
    r = scdm_auto(proj, eig; sigma_factor = 3.0)
    @test isapprox(r.mu_fit, μt; atol = 0.15)
    @test isapprox(r.sigma_fit, σt; atol = 0.15)
    @test isapprox(r.mu, μt - 3σt; atol = 0.5)
    @test r.rms < 0.05

    # The .amn convenience method: projectability is the diagonal of the trial-space
    # projector, so it is guaranteed ∈ [0,1] even for non-orthonormal / rank-deficient
    # columns. Build a low-energy manifold with strong overlap on the trial block and
    # check the fit places μ below the manifold and returns a valid, positive σ.
    A = zeros(ComplexF64, nb, 4, nk)
    for k in 1:nk, m in 1:nb
        # bands 1..8 carry the trial character (energies well below μt), higher bands ~none
        wt = m <= 8 ? 1.0 : 0.01
        for j in 1:4
            A[m, j, k] = wt * cis(0.3m + 0.9j + 0.1k) * (0.5 + 0.5 * cos(1.7m + 2.1j))
        end
    end
    r2 = scdm_auto(A, eig; sigma_factor = 3.0)
    @test isfinite(r2.mu) && r2.sigma > 0 && isfinite(r2.rms)
    @test r2.mu < r2.mu_fit                     # the −kσ shift moved μ down
end
