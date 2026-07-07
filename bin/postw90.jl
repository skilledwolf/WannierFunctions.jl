#!/usr/bin/env julia
#
# Drop-in command-line front end: `postw90.jl <seedname>` post-processes a wannierised
# calculation from seedname.win + seedname.chk(.fmt) (+ .eig/.mmn/.spn/.uHu/... as the tasks
# require) and writes the reference-named output files, analogous to `postw90.x <seedname>`.
#
# Usage:
#   julia --project=/path/to/WannierFunctions.jl bin/postw90.jl <seedname>

let proj = normpath(joinpath(@__DIR__, ".."))
    if Base.active_project() === nothing ||
       !isfile(joinpath(dirname(Base.active_project()), "src", "WannierFunctions.jl"))
        # Load Pkg lazily: with a correct --project (the common case) this costs nothing.
        Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
        Pkg.activate(proj; io=devnull)
    end
end

using WannierFunctions

exit(WannierFunctions.postw90_cli(ARGS))
