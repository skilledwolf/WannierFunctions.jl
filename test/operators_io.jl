# Split from the original monolithic runtests.jl — included by runtests.jl (common.jl
# provides REFROOT/DATAROOT, the shared GaAs/diamond models, and helpers).
@testset "TBOperator + position operator" begin
    if DIAMOND_MODEL !== nothing
        res = run_wannier(DIAMOND_MODEL)                      # native path, :rcg
        H = hamiltonian_operator(DIAMOND_MODEL, res)
        # operator evaluation must agree with the array-based interpolation
        kf = DIAMOND_MODEL.kgrid.frac
        E_op = bands(H, kf)
        irvec, ndegen = wigner_seitz(DIAMOND_MODEL.lattice, DIAMOND_MODEL.kgrid.mp_grid)
        Hr, _ = build_hr(res.U, res.eig_interp, DIAMOND_MODEL.kgrid, irvec)
        E_arr = interpolate_bands(Hr, irvec, ndegen, kf)
        @test maximum(abs.(E_op .- E_arr)) < 1e-12

        # position operator: diag(r(R=0)) must equal the Wannier centres, imag ≈ 0
        rop = position_operator(DIAMOND_MODEL, res)
        ir0 = findfirst(==((0, 0, 0)), rop.irvec)
        @test ir0 !== nothing
        for n in 1:DIAMOND_MODEL.num_wann, c in 1:3
            @test real(rop.data[n, n, ir0, c]) ≈ res.spread.centres[c, n] atol = 1e-8
            @test abs(imag(rop.data[n, n, ir0, c])) < 1e-8
        end
        # NB: r(-R) = r(R)† does NOT hold exactly for this operator — the finite-difference
        # Berry connection i·Σ_b w_b b M is a first-order stencil whose Hermitian defect is a
        # k-mesh artifact (~0.07 Å on this 4×4×4 mesh). The reference writes the same
        # non-symmetrised object to _tb.dat (verified: max element diff vs wannier90.x is
        # 1.3e-8 Å, the E15.8 file precision), so we assert only the guaranteed invariants
        # above and finiteness here.
        @test all(isfinite, rop.data)
    else
        @test_skip false
    end
end

# =========================================================================
# (9) CHECKPOINT (.chk) — binary Fortran-record interchange
# =========================================================================
@testset "Checkpoint read/write" begin
    if DIAMOND_MODEL !== nothing
        win = read_win(DIAMOND_SEED * ".win")
        res = run_wannier(DIAMOND_MODEL, win)
        c = Checkpoint(DIAMOND_MODEL, win, res)
        p = joinpath(mktempdir(), "diamond.chk")
        write_chk(p, c)
        c2 = read_chk(p)
        @test c2.checkpoint == "postwann"
        @test !c2.have_disentangled
        @test c2.u_matrix ≈ c.u_matrix
        @test c2.m_matrix ≈ c.m_matrix
        @test c2.centres ≈ c.centres
        @test c2.spreads ≈ c.spreads
        @test c2.mp_grid == DIAMOND_MODEL.kgrid.mp_grid
        @test c2.real_lattice ≈ Matrix(DIAMOND_MODEL.lattice.A)
        # physics consistency: spread recomputed from the checkpoint's M equals the stored one
        sr = compute_spread(c2.m_matrix, DIAMOND_MODEL.bvectors)
        @test sr.Ω ≈ sum(c2.spreads) atol = 1e-9
    else
        @test_skip false
    end
end

# =========================================================================
# (10) FROZEN-LOCKING ORTHO-FIX + exclude_bands
# =========================================================================
@testset "dis_proj_froz ortho-fix" begin
    # Synthetic near-degenerate case: 4 window states, bands 1–2 frozen, 3 WFs. The third trial
    # column leaks only ε into the non-frozen space, so the required third QPQ eigenvalue is
    # ε² < eps8 and its eigenvector is degenerate with the frozen null space — the case the
    # reference ortho-fix exists for. The selected non-frozen columns must come out orthogonal
    # to the frozen states and orthonormal.
    nd, nf, nwann = 4, 2, 3
    ε = 1e-6
    wd = WannierFunctions.WindowData([1], [nd], [nf], [[1, 2]], [[3, 4]],
                                     [[true, true, false, false]], true)
    A = zeros(ComplexF64, nd, nwann)
    A[1, 1] = 1                       # trial 1 = frozen band 1
    A[2, 2] = 1                       # trial 2 = frozen band 2
    A[1, 3] = 1; A[3, 3] = ε          # trial 3 ≈ frozen band 1 + tiny non-frozen leak
    Uopt = [WannierFunctions.svd_orthonormalize(A)]
    WannierFunctions.dis_proj_froz!(Uopt, wd, nwann)
    U = Uopt[1]
    @test maximum(abs.(U' * U - I(nwann))) < 1e-10            # orthonormal columns
    @test U[1, 1] == 1 && U[2, 2] == 1                        # frozen unit vectors
    for l in nf+1:nwann, ifz in (1, 2)
        @test abs(U[ifz, l]) ≤ 1e-8                           # non-frozen ⊥ frozen states
    end
end

@testset "exclude_bands parsing" begin
    winpath = joinpath(mktempdir(), "x.win")
    write(winpath, """
    num_wann = 4
    exclude_bands = 1-3, 7, 10-11
    mp_grid : 1 1 1
    begin unit_cell_cart
    1 0 0
    0 1 0
    0 0 1
    end unit_cell_cart
    begin kpoints
    0 0 0
    end kpoints
    """)
    win = read_win(winpath)
    @test parse_exclude_bands(win) == [1, 2, 3, 7, 10, 11]
end

# =========================================================================
# (11) Γ-ONLY, SPINOR PROJECTIONS, WF PLOTTING
# =========================================================================
@testset "Gamma-only (silane/benzene)" begin
    si = joinpath(REFROOT, "testw90_example07", "silane")
    bz = joinpath(REFROOT, "testw90_benzene_gamma_val", "benzene")
    if isfile(si * ".mmn")
        model = read_model(si)
        @test model.bvectors.nntot == 6            # half set expanded to the closed full set
        res = run_wannier(model; num_iter = 500)
        @test res.converged
        @test res.spread.Ω ≈ 4.04498078 atol = 2e-6
        @test res.spread.ΩI ≈ 3.920640338 atol = 1e-6
    else
        @test_skip false
    end
    if isfile(bz * ".mmn")
        model = read_model(bz)
        res = run_wannier(model; num_iter = 2000)
        @test res.converged
        @test res.spread.ΩI ≈ 10.455472666 atol = 1e-6
        # We minimise over the full unitary group; the reference Γ algorithm restricts to real
        # orthogonal gauges, so our converged Ω is allowed to be marginally LOWER (verified:
        # the reference is stationary at 12.958338012 even after 20000 sweeps).
        @test res.spread.Ω ≤ 12.958338012 + 1e-6
        @test res.spread.Ω ≈ 12.958338012 atol = 5e-6
    else
        @test_skip false
    end
end

@testset "Spinor projections + orbital ordering" begin
    winpath = joinpath(mktempdir(), "s.win")
    write(winpath, """
    num_wann = 18
    spinors = true
    mp_grid : 1 1 1
    begin unit_cell_cart
    2.0 0 0
    0 2.0 0
    0 0 2.0
    end unit_cell_cart
    begin atoms_frac
    Pt 0.0 0.0 0.0
    end atoms_frac
    begin projections
    Pt: d;s;p
    end projections
    begin kpoints
    0 0 0
    end kpoints
    """)
    win = read_win(winpath)
    projs = parse_projections(win)
    @test length(projs) == 18                                  # 9 spatial × 2 spins
    @test [p.s for p in projs[1:4]] == [1, -1, 1, -1]          # up/down interleaved
    # ascending-l emission regardless of "d;s;p" spec order: s, p×3, d×5 (each doubled)
    @test [p.l for p in projs[1:2:end]] == [0, 1, 1, 1, 2, 2, 2, 2, 2]
end

@testset "Wannier-function plotting (xsf)" begin
    gd = joinpath(REFROOT, "testw90_example01")
    if isfile(joinpath(gd, "UNK00001.1")) && GAAS_MODEL !== nothing
        win = read_win(joinpath(gd, "gaas.win"))
        res = run_wannier(GAAS_MODEL, win)
        tmp = mktempdir()
        paths = plot_wannier_functions(GAAS_MODEL, win, res;
                                       seedname = joinpath(tmp, "gaas"), dir = gd)
        @test length(paths) == 4
        @test all(isfile, paths)
        txt = read(paths[1], String)
        @test occursin("PRIMVEC", txt) && occursin("BEGIN_DATAGRID_3D", txt)
        # the WF grid is normalised like the UNK inputs; peak amplitude is O(1)
        w, ng, los = WannierFunctions.wannier_function_grid(GAAS_MODEL, win, res; list = [1], dir = gd)
        @test ng == (20, 20, 20)
        @test 0.5 < maximum(abs.(w)) < 20
        # phase fixing: the max-|w| point is real positive
        v = w[:, :, :, 1][argmax(abs2.(w[:, :, :, 1]))]
        @test abs(imag(v)) < 1e-10 && real(v) > 0
    else
        @test_skip false
    end
end

# =========================================================================
# (12) SCDM PROJECTIONS + FORMATTED CHECKPOINT
# =========================================================================
@testset "SCDM projections (GaAs)" begin
    gd = joinpath(REFROOT, "testw90_example01")
    if isfile(joinpath(gd, "UNK00001.1")) && GAAS_MODEL !== nothing
        model = read_model(joinpath(gd, "gaas"))
        A = scdm_projections(model; dir = gd)
        @test size(A) == (4, 4, 8)
        model.A = A
        res = run_wannier(model)
        @test res.converged
        # SCDM start must reach the same gauge-invariant minimum as the shipped projections
        @test res.spread.Ω ≈ 4.466880976 atol = 2e-6
        # .amn writer round-trip
        p = joinpath(mktempdir(), "scdm.amn")
        write_amn(p, A)
        A2, nb, nk, nw = read_amn(p)
        @test A2 ≈ A atol = 1e-11
    else
        @test_skip false
    end
end

@testset "Formatted checkpoint (.chk.fmt)" begin
    if DIAMOND_MODEL !== nothing
        win = read_win(DIAMOND_SEED * ".win")
        res = run_wannier(DIAMOND_MODEL, win)
        c = Checkpoint(DIAMOND_MODEL, win, res)
        p = joinpath(mktempdir(), "d.chk.fmt")
        write_chk_fmt(p, c)
        c2 = read_chk_fmt(p)
        @test c2.u_matrix ≈ c.u_matrix
        @test c2.m_matrix ≈ c.m_matrix
        @test c2.centres ≈ c.centres
        @test c2.mp_grid == c.mp_grid
    else
        @test_skip false
    end
end

# =========================================================================
# (13) BERRY CURVATURE / AHC
# =========================================================================
