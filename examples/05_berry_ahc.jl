# Example 5 — Berry curvature and the intrinsic anomalous Hall conductivity of bcc iron.
#
# Consumes a completed Wannier90 run (checkpoint + eig + mmn) — here the Fe test case from the
# reference test suite (spinor, 28 bands disentangled to 18 WFs) — and Wannier-interpolates the
# occupied-manifold Berry curvature over a dense k-mesh. Run from the repository root:
#
#   julia -t auto --project=. examples/05_berry_ahc.jl
#
# Requires the reference tree (reference/wannier90) for the Fe input data; the script stages
# and decompresses what it needs.
using WannierFunctions
using Printf

fed = joinpath(@__DIR__, "..", "reference", "wannier90", "test-suite", "tests", "testpostw90_fe_ahc")
isfile(joinpath(fed, "Fe.chk.fmt.bz2")) || error("Fe test data not found — clone the reference " *
                                                 "wannier90 tree under reference/wannier90 first")
tmp = mktempdir()
for f in ("Fe.win", "Fe.eig")
    cp(joinpath(fed, f), joinpath(tmp, f))
end
run(pipeline(`bunzip2 -kc $(joinpath(fed, "Fe.chk.fmt.bz2"))`, stdout = joinpath(tmp, "Fe.chk.fmt")))
run(pipeline(`bunzip2 -kc $(joinpath(fed, "Fe.mmn.bz2"))`, stdout = joinpath(tmp, "Fe.mmn")))

bm = BerryModel(joinpath(tmp, "Fe"))     # H(R) + Berry connection A(R) from the checkpoint
display(bm); println()

# Berry curvature along a path through the BZ (the sharp spikes near band crossings at E_F
# are what make the AHC integral demanding — and what adaptive meshes chase).
for k in ([0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [0.5, -0.5, -0.5])
    Ω = berry_curvature_k(bm, k, 12.6279)
    @printf("  -2Im f(k=%s)  = (%10.2f, %10.2f, %10.2f) Å²\n", string(k), Ω...)
end

t = @elapsed ahc = anomalous_hall(bm; fermi_energy = 12.6279, kmesh = (10, 10, 10))
@printf("\nAHC on a 10×10×10 mesh (%.1fs, %d threads):\n", t, Threads.nthreads())
@printf("  σ = (%.4f, %.4f, %.4f) S/cm\n", ahc...)
println("  postw90.x reference: (0.0334, 0.0572, 1222.1510) S/cm")
println("  (the large z-component is the magnetisation direction of bcc Fe)")
