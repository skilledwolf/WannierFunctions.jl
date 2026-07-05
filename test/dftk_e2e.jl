# Live DFTK → WannierFunctions end-to-end (included by runtests.jl only when DFTK is
# available). Mirrors examples/06_dftk_end_to_end.jl with test assertions.
using DFTK

let
    a = 10.26
    lattice = a / 2 * [[0 1 1.0]; [1 0 1.0]; [1 1 0.0]]
    Si = ElementPsp(:Si; psp = load_psp("hgh/lda/si-q4"))
    model = model_DFT(lattice, [Si, Si], [ones(3) / 8, -ones(3) / 8];
                      functionals = LDA(), symmetries = false)
    basis = PlaneWaveBasis(model; Ecut = 14, kgrid = (4, 4, 4))
    scfres = self_consistent_field(basis; tol = 1e-10)

    centers = [[1, 1, 1] / 8, [-3, 1, 1] / 8, [1, -3, 1] / 8, [1, 1, -3] / 8]
    projs = [DFTK.GaussianWannierProjection(c) for c in centers]
    wmodel = wannier_model(scfres, projs; num_wann = 4)
    res = wannierise(wmodel; num_iter = 500, algorithm = :w90, conv_tol = 1e-10,
                     conv_window = 5)
    @test res.converged
    @test res.spread.Ω ≈ 6.4566 atol = 0.01
    # WFs sit on the four Si–Si bond centres: (0,0,0), (½,0,0), (0,½,0), (0,0,½) mod 1
    Binv = inv(Matrix(wmodel.lattice.A))
    for e in ([0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5])
        @test any(1:4) do n
            d = Binv * res.spread.centres[:, n] .- e
            all(abs.(d .- round.(d)) .< 0.02)             # periodic distance
        end
    end
    # interpolation reproduces the SCF eigenvalues on the mesh
    irvec, ndegen = wigner_seitz(wmodel.lattice, wmodel.kgrid.mp_grid)
    Hr, _ = build_hr(res.U, wmodel.eig, wmodel.kgrid, irvec)
    maxdev = 0.0
    for (ik, kf) in enumerate(wmodel.kgrid.frac)
        E = interpolate_bands(Hr, irvec, ndegen, [kf])[:, 1]
        maxdev = max(maxdev, maximum(abs.(E .- wmodel.eig[:, ik])))
    end
    @test maxdev < 1e-8

    # SCDM (projection-free) path: num_wann is the only Wannier-specific input, and it must
    # reach the identical maximally-localised gauge.
    wmodel2 = wannier_model(scfres; num_wann = 4)
    res2 = wannierise(wmodel2; num_iter = 500, algorithm = :w90, conv_tol = 1e-10,
                      conv_window = 5)
    @test res2.converged
    @test res2.spread.Ω ≈ res.spread.Ω atol = 1e-6
end
