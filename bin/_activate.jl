# Shared project activation for the bin/ entry points. Must run BEFORE `using
# WannierFunctions`, so it lives here rather than in the package.
let proj = normpath(joinpath(@__DIR__, ".."))
    if Base.active_project() === nothing ||
       !isfile(joinpath(dirname(Base.active_project()), "src", "WannierFunctions.jl"))
        # Load Pkg lazily: with a correct --project (the common case) this costs nothing.
        # invokelatest: Pkg enters the world inside this same top-level expression.
        Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
        Base.invokelatest(Pkg.activate, proj; io=devnull)
    end
end
