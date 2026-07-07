#!/usr/bin/env julia
#
# Drop-in for `w90chk2chk.x`: convert between the binary `.chk` and formatted `.chk.fmt`
# checkpoint formats (either direction, auto-detected from the flag).
#
# Usage:
#   w90chk2chk.jl -export <seedname>   # seedname.chk     -> seedname.chk.fmt
#   w90chk2chk.jl -import <seedname>   # seedname.chk.fmt -> seedname.chk

let proj = normpath(joinpath(@__DIR__, ".."))
    if Base.active_project() === nothing ||
       !isfile(joinpath(dirname(Base.active_project()), "src", "WannierFunctions.jl"))
        # Load Pkg lazily: with a correct --project (the common case) this costs nothing.
        Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
        Pkg.activate(proj; io=devnull)
    end
end

using WannierFunctions

exit(WannierFunctions.w90chk2chk_cli(ARGS))
