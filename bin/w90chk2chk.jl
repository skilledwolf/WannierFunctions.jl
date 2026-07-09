#!/usr/bin/env julia
#
# Drop-in for `w90chk2chk.x`: convert between the binary `.chk` and formatted `.chk.fmt`
# checkpoint formats (either direction, auto-detected from the flag).
#
# Usage:
#   w90chk2chk.jl -export <seedname>   # seedname.chk     -> seedname.chk.fmt
#   w90chk2chk.jl -import <seedname>   # seedname.chk.fmt -> seedname.chk

include(joinpath(@__DIR__, "_activate.jl"))

using WannierFunctions

exit(WannierFunctions.w90chk2chk_cli(ARGS))
