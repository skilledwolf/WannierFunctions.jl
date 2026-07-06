#!/usr/bin/env julia
#
# Drop-in command-line front end: `wannier90.jl <seedname>` runs the full wannierisation +
# interpolation pipeline on seedname.win (+ .amn/.mmn/.eig) and writes seedname.wout and the
# requested output files, analogous to the reference `wannier90.x <seedname>`.
#
# Usage:
#   julia --project=/path/to/WannierFunctions.jl bin/wannier90.jl <seedname>

import Pkg
let proj = normpath(joinpath(@__DIR__, ".."))
    if Base.active_project() === nothing ||
       !isfile(joinpath(dirname(Base.active_project()), "src", "WannierFunctions.jl"))
        Pkg.activate(proj; io=devnull)
    end
end

using WannierFunctions

exit(WannierFunctions.wannier90_cli(ARGS))
