# Static package-quality checks: Aqua (method ambiguities, unbound type parameters, stale
# deps, compat bounds, piracy) and, when the resolver could install it for this Julia
# version, a JET type-level error scan of the package.
using Aqua

@testset "Aqua.jl package quality" begin
    # Ambiguities are tested against our own methods only (Base/StaticArrays cross-package
    # ambiguities are not actionable here).
    Aqua.test_all(WannierFunctions; ambiguities=(recursive = false,))
end
