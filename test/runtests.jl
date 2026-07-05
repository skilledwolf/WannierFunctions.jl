using Test
using LinearAlgebra
using StaticArrays
using WannierFunctions

# ---------------------------------------------------------------------------
# Reference-tree location. All validation (and most unit) tests read the GaAs
# and diamond inputs from the vendored Wannier90 test-suite. If the reference
# tree is absent (e.g. lean CI), those tests are skipped rather than errored.
# ---------------------------------------------------------------------------
const REFROOT = joinpath(@__DIR__, "..", "reference", "wannier90",
                         "test-suite", "tests")
const DATAROOT = joinpath(@__DIR__, "..", "examples", "data")

# Prefer the vendored reference test-suite inputs; fall back to the identical files shipped
# under examples/data so CI without the reference clone still validates the physics.
_seed(refdir, name) = isfile(joinpath(REFROOT, refdir, name * ".win")) ?
                      joinpath(REFROOT, refdir, name) : joinpath(DATAROOT, name)
const GAAS_SEED    = _seed("testw90_example01", "gaas")
const DIAMOND_SEED = _seed("testw90_example05", "diamond")

has_gaas()    = isfile(GAAS_SEED * ".win") && isfile(GAAS_SEED * ".amn") &&
                isfile(GAAS_SEED * ".mmn")
has_diamond() = isfile(DIAMOND_SEED * ".win") && isfile(DIAMOND_SEED * ".amn") &&
                isfile(DIAMOND_SEED * ".mmn") && isfile(DIAMOND_SEED * ".eig")

# Build models once and reuse (read_model is the expensive I/O step).
const GAAS_MODEL    = has_gaas()    ? read_model(GAAS_SEED)    : nothing
const DIAMOND_MODEL = has_diamond() ? read_model(DIAMOND_SEED) : nothing

"Σ_b w_b b_α b_β at a single k-point (the per-k B1 completeness matrix)."
function b1_matrix(bv::WannierFunctions.BVectors, k::Int)
    S = zeros(3, 3)
    for b in 1:bv.nntot
        w = bv.wb[b, k]
        bb = SVector{3,Float64}(bv.bvec[1, b, k], bv.bvec[2, b, k], bv.bvec[3, b, k])
        S .+= w .* (bb * bb')
    end
    return S
end

@testset "Wannier90 test suite" begin

# =========================================================================
# (1) UNIT TESTS
# =========================================================================
@testset "Unit tests" begin

    # ---- .amn / .mmn / .eig parse dimensions ----
    @testset "I/O parse dimensions" begin
        if has_gaas()
            A, nb_a, nk_a, nw_a = read_amn(GAAS_SEED * ".amn")
            @test size(A) == (nb_a, nw_a, nk_a)
            @test (nb_a, nw_a, nk_a) == (4, 4, 8)      # GaAs 4 bands / 4 WF / 2×2×2

            M, kpb, gpb, nb_m, nk_m, nntot = read_mmn(GAAS_SEED * ".mmn")
            @test size(M) == (nb_m, nb_m, nntot, nk_m)
            @test size(kpb) == (nntot, nk_m)
            @test size(gpb) == (3, nntot, nk_m)
            @test nb_m == nb_a
            @test nk_m == nk_a

            # GaAs ships no .eig
            @test !isfile(GAAS_SEED * ".eig")
        else
            @info "GaAs reference inputs missing; skipping GaAs I/O dimension tests" GAAS_SEED
            @test_skip false
        end

        if has_diamond()
            A, nb_a, nk_a, nw_a = read_amn(DIAMOND_SEED * ".amn")
            @test size(A) == (nb_a, nw_a, nk_a)
            @test (nb_a, nw_a, nk_a) == (4, 4, 64)     # diamond 4/4, 4×4×4

            M, kpb, gpb, nb_m, nk_m, nntot = read_mmn(DIAMOND_SEED * ".mmn")
            @test size(M) == (nb_m, nb_m, nntot, nk_m)
            @test nk_m == nk_a
            @test nb_m == nb_a

            eig = read_eig(DIAMOND_SEED * ".eig")
            @test size(eig) == (nb_a, nk_a)            # (num_bands × num_kpts)
        else
            @info "diamond reference inputs missing; skipping diamond I/O dimension tests" DIAMOND_SEED
            @test_skip false
        end
    end

    # ---- kmesh B1 completeness: Σ_b w_b b⊗b ≈ I (per k-point) ----
    @testset "kmesh B1 completeness" begin
        for (name, model) in (("GaAs", GAAS_MODEL), ("diamond", DIAMOND_MODEL))
            if model === nothing
                @info "$name reference missing; skipping B1 completeness" name
                @test_skip false
                continue
            end
            bv = model.bvectors
            nk = WannierFunctions.nkpt(model.kgrid)
            # Completeness is a PER-K relation; check every k-point.
            for k in 1:nk
                @test b1_matrix(bv, k) ≈ Matrix(I, 3, 3) atol = 1e-8
            end
        end
    end

    # ---- Löwdin gauge unitarity: U_k† U_k ≈ I ----
    @testset "Löwdin gauge unitarity" begin
        for (name, model) in (("GaAs", GAAS_MODEL), ("diamond", DIAMOND_MODEL))
            if model === nothing
                @info "$name reference missing; skipping Löwdin unitarity" name
                @test_skip false
                continue
            end
            U = initial_gauge(model.A)
            nw = model.num_wann
            for k in 1:size(U, 3)
                Uk = U[:, :, k]
                @test Uk' * Uk ≈ Matrix(I, nw, nw) atol = 1e-10
            end
        end
    end

    # ---- Wigner–Seitz sum rule: Σ 1/ndegen == ∏ mp_grid ----
    @testset "Wigner–Seitz sum rule" begin
        for (name, model) in (("GaAs", GAAS_MODEL), ("diamond", DIAMOND_MODEL))
            if model === nothing
                @info "$name reference missing; skipping WS sum rule" name
                @test_skip false
                continue
            end
            irvec, ndegen = wigner_seitz(model.lattice, model.kgrid.mp_grid)
            @test length(irvec) == length(ndegen)
            @test all(>(0), ndegen)
            @test sum(1.0 / d for d in ndegen) ≈ prod(model.kgrid.mp_grid) atol = 1e-8
        end
    end

    # ---- spread self-consistency: Ω ≈ ΩI+ΩOD+ΩD ≈ Σ spreads ----
    @testset "spread self-consistency" begin
        for (name, model) in (("GaAs", GAAS_MODEL), ("diamond", DIAMOND_MODEL))
            if model === nothing
                @info "$name reference missing; skipping spread self-consistency" name
                @test_skip false
                continue
            end
            sr = initial_spread(model)
            # The load-bearing check: decomposition sums to the total.
            @test sr.Ω ≈ sr.ΩI + sr.ΩOD + sr.ΩD atol = 1e-10
            # Ω is assembled as sum(spreads) in compute_spread — machine-exact,
            # but assert it anyway as a guard against future refactors.
            @test sr.Ω ≈ sum(sr.spreads) atol = 1e-12
        end
    end
end

# =========================================================================
# (2) VALIDATION TESTS — exact reference numbers, test-suite tolerances.
#
# Tolerances (docs/reference-notes/test-suite-and-targets.md §1.4):
#   omegaI / omegaOD / omegaTotal : abs 1e-6
#   omegaD                        : rel 5e-6 (abs 1e-6)
#   centres                       : abs 1e-5
#   spreads                       : abs 3e-6
# =========================================================================
@testset "Validation vs reference benchmarks" begin

    # ---- GaAs testw90_example01 ----
    @testset "GaAs example01" begin
        if GAAS_MODEL === nothing
            @info "GaAs reference missing; skipping validation" GAAS_SEED
            @test_skip false
        else
            # The GaAs .win cell is in BOHR; converting to Å introduces a
            # ~4e-8 constant offset in the Ω components relative to the golden
            # .wout. That is well inside the 1e-6 relative tolerance, but to
            # stay on the safe side of the strict "< abs" comparison we give a
            # tiny 2e-6 absolute margin on the Ω components (still tighter than
            # the reference rel tol at these magnitudes).
            ω_abs = 2e-6

            # -- Initial state --
            sr0 = initial_spread(GAAS_MODEL)
            @test sr0.Ω  ≈ 4.4688121 atol = ω_abs   # initial Omega_Total
            @test sr0.ΩD ≈ 0.0083198 atol = 1e-6    # Iter 0 O_D (7-digit ref)
            @test sr0.ΩOD ≈ 0.5036294 atol = 1e-6   # Iter 0 O_OD

            # -- Final state after 20 iterations --
            res = wannierise(GAAS_MODEL; num_iter = 20)
            sf = res.spread
            @test sf.Ω   ≈ 4.466880976 atol = ω_abs
            @test sf.ΩI  ≈ 3.956862958 atol = ω_abs
            @test sf.ΩD  ≈ 0.008030049 atol = ω_abs
            @test sf.ΩD  ≈ 0.008030049 rtol = 5e-6   # omegaD relative tol
            @test sf.ΩOD ≈ 0.501987969 atol = ω_abs

            # Per-WF spreads: all four equal 1.11672024 (abs 3e-6).
            for n in 1:GAAS_MODEL.num_wann
                @test sf.spreads[n] ≈ 1.11672024 atol = 3e-6
            end

            # Centres (Å, abs 1e-5). Reference final-state WF centres.
            ref_centres = [
                -0.866253  1.973841  1.973841;
                -0.866253  0.866253  0.866253;
                -1.973841  1.973841  0.866253;
                -1.973841  0.866253  1.973841;
            ]  # rows = WF, cols = x,y,z
            # centres in SpreadResult are (3 × nw): compare as a set (WF order
            # from Löwdin need not match the .wout print order), so match each
            # reference centre to the nearest computed centre.
            computed = [SVector{3,Float64}(sf.centres[:, n]) for n in 1:GAAS_MODEL.num_wann]
            for r in 1:size(ref_centres, 1)
                rc = SVector{3,Float64}(ref_centres[r, :])
                dmin = minimum(norm(c - rc) for c in computed)
                @test dmin ≤ 1e-5
            end
        end
    end

    # ---- diamond testw90_example05 ----
    @testset "diamond example05" begin
        if DIAMOND_MODEL === nothing
            @info "diamond reference missing; skipping validation" DIAMOND_SEED
            @test_skip false
        else
            res = wannierise(DIAMOND_MODEL; num_iter = 20)
            sf = res.spread
            @test sf.ΩI  ≈ 1.954619860 atol = 1e-6
            @test sf.ΩD  ≈ 0.0          atol = 1e-9   # exactly zero by symmetry
            @test sf.ΩOD ≈ 0.366285055 atol = 1e-6
            @test sf.Ω   ≈ 2.320904915 atol = 1e-6

            # -- M2 exact-at-grid-points invariant --
            # Interpolated eigenvalues at the mp_grid k-points must equal the
            # sorted input .eig to 1e-8 (H(k) = U†diag(ε)U is unitarily similar
            # to diag(ε); the WS/ndegen set is dual to the k-grid).
            @test DIAMOND_MODEL.eig !== nothing
            irvec, ndegen = wigner_seitz(DIAMOND_MODEL.lattice, DIAMOND_MODEL.kgrid.mp_grid)
            Hr, _ = build_hr(res.U, DIAMOND_MODEL.eig, DIAMOND_MODEL.kgrid, irvec)
            kfrac = DIAMOND_MODEL.kgrid.frac
            interp = interpolate_bands(Hr, irvec, ndegen, kfrac)   # (nw × nk), ascending
            nk = length(kfrac)
            maxerr = 0.0
            for k in 1:nk
                ref = sort(DIAMOND_MODEL.eig[:, k])
                maxerr = max(maxerr, maximum(abs.(interp[:, k] .- ref)))
            end
            @test maxerr ≤ 1e-8
        end
    end
end

end  # top-level testset

# =========================================================================
# (3) DISENTANGLEMENT VALIDATION (M3) — silicon & copper
#     Silicon's overlap file ships bz2-compressed; both cases are staged into
#     a tempdir. Skipped cleanly if the reference tree or `bunzip2` is absent.
# =========================================================================
@testset "Disentanglement validation" begin
    tests_dir = joinpath(@__DIR__, "..", "reference", "wannier90", "test-suite", "tests")

    # -- silicon example03: 12 → 8 WF, outer + frozen window --
    @testset "silicon example03" begin
        si = joinpath(tests_dir, "testw90_example03")
        mmnbz = joinpath(@__DIR__, "..", "reference", "wannier90", "test-suite",
                         "checkpoints", "si_geninterp", "silicon.mmn.bz2")
        staged = false
        seed = ""
        if isfile(joinpath(si, "silicon.win")) && isfile(mmnbz) &&
           Sys.which("bunzip2") !== nothing
            tmp = mktempdir()
            for f in ("silicon.win", "silicon.amn", "silicon.eig")
                cp(joinpath(si, f), joinpath(tmp, f); follow_symlinks = true)
            end
            try
                run(pipeline(`bunzip2 -kc $mmnbz`, stdout = joinpath(tmp, "silicon.mmn")))
                seed = joinpath(tmp, "silicon"); staged = true
            catch
                staged = false
            end
        end
        if staged
            model, win, res = run_wannier(seed)
            @test res.disentangled
            s = res.spread
            @test s.ΩI  ≈ 11.849193709 atol = 1e-6
            @test s.ΩD  ≈ 0.105470244  atol = 1e-6
            @test s.ΩOD ≈ 2.544910550  atol = 1e-6
            @test s.Ω   ≈ 14.499574503 atol = 1e-6
            # Ω_I convergence trace first row matches the reference exactly.
            @test res.dis.omega_I_trace[1][2] ≈ 12.70775084 atol = 1e-5
        else
            @info "silicon inputs unavailable (need reference tree + bunzip2); skipping"
            @test_skip false
        end
    end

    # -- copper example04: 12 → 7 WF (metal), outer + frozen window --
    @testset "copper example04" begin
        cu = joinpath(tests_dir, "testw90_example04")
        if all(isfile(joinpath(cu, "copper." * e)) for e in ("win", "amn", "mmn", "eig"))
            model, win, res = run_wannier(joinpath(cu, "copper"))
            @test res.disentangled
            s = res.spread
            @test s.ΩI  ≈ 3.662691490 atol = 1e-6
            @test s.ΩOD ≈ 0.363454087 atol = 1e-6
            @test s.Ω   ≈ 4.028040058 atol = 1e-6
        else
            @info "copper inputs missing; skipping" cu
            @test_skip false
        end
    end
end

# =========================================================================
# (4) INTERPOLATION REFINEMENTS — use_ws_distance and k-path generation
# =========================================================================
@testset "Interpolation refinements" begin
    # use_ws_distance interpolation must still satisfy the exact-at-grid-points invariant
    # (the per-pair minimal-image shifts are multiples of mp_grid, so they add no phase at
    # grid k-points), and must run without error.
    @testset "ws_distance grid invariant" begin
        if DIAMOND_MODEL !== nothing
            res = wannierise(DIAMOND_MODEL; num_iter = 20)
            irvec, ndegen = wigner_seitz(DIAMOND_MODEL.lattice, DIAMOND_MODEL.kgrid.mp_grid)
            Hr, _ = build_hr(res.U, DIAMOND_MODEL.eig, DIAMOND_MODEL.kgrid, irvec)
            kfrac = DIAMOND_MODEL.kgrid.frac
            Ews = interpolate_bands_ws(Hr, irvec, ndegen, res.spread.centres,
                                       DIAMOND_MODEL.lattice, DIAMOND_MODEL.kgrid.mp_grid, kfrac)
            maxerr = maximum(maximum(abs.(Ews[:, k] .- sort(DIAMOND_MODEL.eig[:, k])))
                             for k in 1:length(kfrac))
            @test maxerr ≤ 1e-8
        else
            @test_skip false
        end
    end

    # generate_kpath must reprint the endpoint of a discontinuous segment (and its label).
    @testset "k-path discontinuity" begin
        winpath = joinpath(mktempdir(), "disc.win")
        write(winpath, """
        num_wann = 1
        mp_grid : 1 1 1
        begin unit_cell_cart
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
        end unit_cell_cart
        begin kpoints
        0.0 0.0 0.0
        end kpoints
        begin kpoint_path
        A 0.0 0.0 0.0  B 0.5 0.0 0.0
        B 0.5 0.0 0.0  C 0.5 0.5 0.0
        D 0.0 0.5 0.0  A 0.0 0.0 0.0
        end kpoint_path
        """)
        win = read_win(winpath)
        lat = WannierFunctions.Lattice(win.unit_cell)
        kpts, xvals, labels, lidx = generate_kpath(win, lat; bands_num_points = 10)
        # A, B, C (end of the discontinuous 2nd segment), D (start of 3rd), A → 5 labels
        @test "C" in labels
        @test "D" in labels
        @test length(labels) == 5
        @test issorted(xvals)                       # monotone non-decreasing path coordinate
        @test length(lidx) == length(labels)
    end
end

# =========================================================================
# (5) MULTI-SHELL kmesh — the FCC validation cases are all single-shell, so
#     exercise the 2-shell B1 weight solve on a synthetic tetragonal mesh.
# =========================================================================
@testset "Multi-shell kmesh (tetragonal)" begin
    mp = (2, 2, 2)
    A = SMatrix{3,3,Float64}([1.0 0 0; 0 1.0 0; 0 0 2.0])   # a=a=1, c=2 → two shells needed
    lat = WannierFunctions.Lattice(A)
    kfrac = vec([SVector(x, y, z) for x in (0.0, 0.5), y in (0.0, 0.5), z in (0.0, 0.5)])
    kg = WannierFunctions.KGrid(kfrac, mp)
    dirs = [SVector(0.5,0,0), SVector(-0.5,0,0), SVector(0,0.5,0),
            SVector(0,-0.5,0), SVector(0,0,0.5), SVector(0,0,-0.5)]
    nk = length(kfrac); nntot = length(dirs)
    kpb = Matrix{Int}(undef, nntot, nk); gpb = Array{Int,3}(undef, 3, nntot, nk)
    for k in 1:nk, (b, Δ) in enumerate(dirs)
        target = kfrac[k] + Δ
        kp = findfirst(kp -> norm(round.(target - kfrac[kp]) - (target - kfrac[kp])) < 1e-9, 1:nk)
        gpb[:, b, k] = round.(Int, target - kfrac[kp]); kpb[b, k] = kp
    end
    bv = WannierFunctions.build_bvectors(kg, lat, kpb, gpb)        # errors if B1 not satisfied
    @test length(bv.shells) == 2                            # two distinct shell radii
    @test bv.shell_weight[1] != bv.shell_weight[2]          # genuinely different weights
    # Explicit B1 completeness Σ_b w_b b⊗b = I.
    S = zeros(3, 3)
    for b in 1:bv.nntot
        w = bv.wb[b, 1]
        bb = SVector{3,Float64}(bv.bvec[1,b,1], bv.bvec[2,b,1], bv.bvec[3,b,1])
        S .+= w .* (bb * bb')
    end
    @test norm(S - I) < 1e-10
end

# =========================================================================
# (6) OPTIMIZER PARITY — :rcg (native default) and :w90 (reference-faithful)
#     must find the same spread minimum; :rcg must actually converge.
# =========================================================================
@testset "Optimizer parity (:rcg vs :w90)" begin
    if GAAS_MODEL !== nothing
        r_rcg = wannierise(GAAS_MODEL)                              # :rcg default
        r_w90 = wannierise(GAAS_MODEL; algorithm = :w90, num_iter = 20)
        @test r_rcg.converged
        @test r_rcg.spread.Ω ≈ r_w90.spread.Ω atol = 1e-8           # same minimum
        @test r_rcg.spread.Ω ≈ 4.466880976 atol = 2e-6              # reference value
    else
        @test_skip false
    end
    if DIAMOND_MODEL !== nothing
        r = wannierise(DIAMOND_MODEL)
        @test r.converged
        @test r.spread.Ω ≈ 2.320904915 atol = 1e-6
        @test r.spread.ΩD ≈ 0.0 atol = 1e-9
    else
        @test_skip false
    end
end

# =========================================================================
# (7) .nnkp GENERATION (-pp mode) — the k-mesh built from the .win alone must
#     reproduce the shells/weights derived from the .mmn connectivity, and the
#     projections block must parse.
# =========================================================================
@testset ".nnkp generation (-pp)" begin
    if GAAS_MODEL !== nothing
        win = read_win(GAAS_SEED * ".win")
        out, info = generate_nnkp(GAAS_SEED; out = joinpath(mktempdir(), "gaas.nnkp"))
        @test isfile(out)
        @test info.nntot == GAAS_MODEL.bvectors.nntot                     # 8 neighbours
        @test length(info.weights) == length(GAAS_MODEL.bvectors.shell_weight)
        @test info.weights[1] ≈ GAAS_MODEL.bvectors.shell_weight[1] atol = 1e-9
        projs = parse_projections(win)
        @test length(projs) == 4                                          # As:sp3 → 4 orbitals
        @test all(p -> p.l == -3, projs)                                  # sp3 code
        @test sort([p.mr for p in projs]) == [1, 2, 3, 4]
        # nnkpts block consistency vs the .mmn-derived connectivity: the b-vector SETS at
        # k=1 must agree (ordering may legitimately differ between the two constructions).
        lines = readlines(out)
        i0 = findfirst(==("begin nnkpts"), lines)
        nn = parse(Int, lines[i0+1])
        mine = Set{NTuple{4,Int}}()
        for j in 1:nn
            t = parse.(Int, split(lines[i0+1+j]))
            push!(mine, (t[2], t[3], t[4], t[5]))
        end
        frommmn = Set{NTuple{4,Int}}()
        for b in 1:GAAS_MODEL.bvectors.nntot
            push!(frommmn, (GAAS_MODEL.bvectors.kpb[b, 1], GAAS_MODEL.bvectors.gpb[1, b, 1],
                            GAAS_MODEL.bvectors.gpb[2, b, 1], GAAS_MODEL.bvectors.gpb[3, b, 1]))
        end
        @test mine == frommmn
    else
        @test_skip false
    end
end

# =========================================================================
# (8) OPERATOR API — TBOperator invariants
# =========================================================================
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
            run(pipeline(`bunzip2 -kc $(joinpath(sc, f * ".bz2"))`, stdout = joinpath(tmp, f)))
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
