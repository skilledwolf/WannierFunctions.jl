#!/usr/bin/env julia
#
# Drop-in for `w90chk2chk.x`: convert between the binary `.chk` and formatted `.chk.fmt`
# checkpoint formats (either direction, auto-detected from the flag).
#
# Usage:
#   w90chk2chk.jl -export <seedname>   # seedname.chk     -> seedname.chk.fmt
#   w90chk2chk.jl -import <seedname>   # seedname.chk.fmt -> seedname.chk

import Pkg
let proj = normpath(joinpath(@__DIR__, ".."))
    if Base.active_project() === nothing ||
       !isfile(joinpath(dirname(Base.active_project()), "src", "WannierFunctions.jl"))
        Pkg.activate(proj; io=devnull)
    end
end

using WannierFunctions

function usage()
    println(stderr, "usage: w90chk2chk.jl -export|-import <seedname>")
    exit(1)
end

length(ARGS) == 2 || usage()
mode, seed = ARGS[1], replace(ARGS[2], r"\.chk(\.fmt)?$" => "")
if mode in ("-export", "--export", "-u2f")          # binary -> formatted
    chk = read_chk(seed * ".chk")
    write_chk_fmt(seed * ".chk.fmt", chk)
    println("$(seed).chk -> $(seed).chk.fmt")
elseif mode in ("-import", "--import", "-f2u")      # formatted -> binary
    chk = read_chk_fmt(seed * ".chk.fmt")
    write_chk(seed * ".chk", chk)
    println("$(seed).chk.fmt -> $(seed).chk")
else
    usage()
end
