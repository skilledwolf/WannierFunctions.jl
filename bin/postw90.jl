#!/usr/bin/env julia
#
# Drop-in command-line front end: `postw90.jl <seedname>` post-processes a wannierised
# calculation from seedname.win + seedname.chk(.fmt) (+ .eig/.mmn/.spn/.uHu/... as the tasks
# require) and writes the reference-named output files, analogous to `postw90.x <seedname>`.
#
# Usage:
#   julia --project=/path/to/WannierFunctions.jl bin/postw90.jl <seedname>

import Pkg
let proj = normpath(joinpath(@__DIR__, ".."))
    if Base.active_project() === nothing ||
       !isfile(joinpath(dirname(Base.active_project()), "src", "WannierFunctions.jl"))
        Pkg.activate(proj; io=devnull)
    end
end

using WannierFunctions

if isempty(ARGS)
    println(stderr, "usage: postw90.jl <seedname>")
    exit(1)
end
postw90_main(replace(ARGS[1], r"\.win$" => ""))
