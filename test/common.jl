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

