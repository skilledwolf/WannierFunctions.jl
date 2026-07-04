#!/usr/bin/env julia
#
# Drop-in command-line front end: `wannier90.jl <seedname>` runs the full wannierisation +
# interpolation pipeline on seedname.win (+ .amn/.mmn/.eig) and writes seedname.wout and the
# requested output files, analogous to the reference `wannier90.x <seedname>`.
#
# Usage:
#   julia --project=/path/to/Wannier90.jl bin/wannier90.jl <seedname>

import Pkg
let proj = normpath(joinpath(@__DIR__, ".."))
    if Base.active_project() === nothing ||
       !isfile(joinpath(dirname(Base.active_project()), "src", "Wannier90.jl"))
        Pkg.activate(proj; io=devnull)
    end
end

using Wannier90

args = copy(ARGS)
pp = false
filter!(args) do a
    a in ("-pp", "--pp", "-postproc") && (pp = true; return false)
    return true
end
if isempty(args)
    println(stderr, "usage: wannier90.jl [-pp] <seedname>")
    exit(1)
end
Wannier90.main(replace(args[1], r"\.win$" => ""); pp=pp)
