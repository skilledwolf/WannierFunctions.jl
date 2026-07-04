# Human-readable `.wout` log writer. Reproduces the reference Wannier90 layout for the sections
# that carry the physics (nearest-neighbour shells + b-vector weights, initial/final Wannier
# centres and spreads, and the Ω_I/Ω_D/Ω_OD/Ω_Total decomposition) in the same fixed-decimal
# `F`-format the reference uses, so the numbers are directly comparable to a reference `.wout`.

using Printf

"Write the b-vector shell / weight report block."
function _wout_shells(io, model::Model)
    bv = model.bvectors
    println(io, " *---------------------------------- K-MESH ----------------------------------*")
    @printf(io, " | Number of Nearest Neighbour shells used : %3d                              |\n",
            length(bv.shells))
    println(io, " +----------------------------------------------------------------------------+")
    println(io, " |                  b_k Vectors (Ang^-1) and Weights (Ang^2)                  |")
    println(io, " |            No.         b_k(x)      b_k(y)      b_k(z)        w_b            |")
    println(io, " |            ---        --------------------------------     --------        |")
    for b in 1:bv.nntot
        @printf(io, " | %13d  %11.6f %11.6f %11.6f %13.6f     |\n",
                b, bv.bvec[1, b, 1], bv.bvec[2, b, 1], bv.bvec[3, b, 1], bv.wb[b, 1])
    end
    println(io, " +----------------------------------------------------------------------------+")
    println(io)
end

"Write a Wannier-centres-and-spreads block (used for both Initial and Final State)."
function _wout_centres(io, sr::SpreadResult)
    nw = length(sr.spreads)
    sx = sy = sz = 0.0
    for n in 1:nw
        x, y, z = sr.centres[1, n], sr.centres[2, n], sr.centres[3, n]
        @printf(io, "  WF centre and spread %4d  ( %10.6f, %10.6f, %10.6f )  %15.8f\n",
                n, x, y, z, sr.spreads[n])
        sx += x; sy += y; sz += z
    end
    @printf(io, "  Sum of centres and spreads ( %10.6f, %10.6f, %10.6f )  %15.8f\n",
            sx, sy, sz, sr.Ω)
    println(io)
end

"""
    write_wout(path, model, win, result; dis=nothing, omega_trace=nothing)

Write the `.wout` log. `result::WannierResult` supplies the final spread; `dis::DisentangleResult`
(optional) adds the disentanglement Ω_I convergence table.
"""
function write_wout(path::AbstractString, model::Model, win::WinInput, result::WannierResult;
                    dis=nothing)
    open(path, "w") do io
        println(io, " +---------------------------------------------------------------------------+")
        println(io, " |                       WannierFunctions.jl  (Julia)                         |")
        println(io, " |     A modern from-scratch reimplementation of the Wannier90 core.          |")
        println(io, " |     Independent project; not affiliated with wannier-developers.           |")
        println(io, " +---------------------------------------------------------------------------+")
        println(io)
        @printf(io, "  Number of Wannier Functions               :  %6d\n", model.num_wann)
        @printf(io, "  Number of input Bloch states              :  %6d\n", model.num_bands)
        @printf(io, "  Number of k-points                        :  %6d\n", nkpt(model.kgrid))
        @printf(io, "  Grid size                                 =  %2d x %2d x %2d\n",
                model.kgrid.mp_grid...)
        @printf(io, "  Unit cell volume                          :  %14.6f  Ang^3\n",
                cell_volume(model.lattice))
        println(io)
        println(io, "                             Lattice Vectors (Ang)")
        for i in 1:3
            @printf(io, "                 a_%d  %11.6f %11.6f %11.6f\n",
                    i, model.lattice.A[1, i], model.lattice.A[2, i], model.lattice.A[3, i])
        end
        println(io)
        _wout_shells(io, model)

        if result.disentangled && dis !== nothing
            println(io, " *------------------------- DISENTANGLE -------------------------*")
            println(io, " +---------------------------------------------------------------+<-- DIS")
            println(io, " |  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)      |<-- DIS")
            println(io, " +---------------------------------------------------------------+<-- DIS")
            for (it, wi1, wi, d) in dis.omega_I_trace
                @printf(io, " %8d   %14.8f  %14.8f   %13.3E    <-- DIS\n", it, wi1, wi, d)
            end
            println(io)
            @printf(io, "  Final Omega_I (disentanglement) = %15.9f Ang^2\n", dis.omega_I)
            println(io)
        end

        # Final localised state.
        println(io, " Final State")
        _wout_centres(io, result.spread)
        s = result.spread
        @printf(io, "         Spreads (Ang^2)       Omega I      = %15.9f\n", s.ΩI)
        @printf(io, "        ================       Omega D      = %15.9f\n", s.ΩD)
        @printf(io, "                               Omega OD     = %15.9f\n", s.ΩOD)
        @printf(io, "    Final Spread (Ang^2)       Omega Total  = %15.9f\n", s.Ω)
        println(io, " ------------------------------------------------------------------------------")
        println(io)
        println(io, " All done: WannierFunctions.jl exiting cleanly.")
    end
    return path
end
