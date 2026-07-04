# Human-friendly REPL/notebook display for the core types.

using Printf

function Base.show(io::IO, ::MIME"text/plain", m::Model)
    println(io, "WannierFunctions Model \"", m.seedname, "\"")
    @printf(io, "  bands → Wannier functions : %d → %d%s\n", m.num_bands, m.num_wann,
            m.num_bands > m.num_wann ? "   (needs disentanglement)" : "   (isolated)")
    @printf(io, "  k-mesh                    : %d × %d × %d   (%d k-points)\n",
            m.kgrid.mp_grid..., nkpt(m.kgrid))
    @printf(io, "  b-vectors                 : %d per k, %d shell%s, Σw_b b⊗b = 1 ✓\n",
            m.bvectors.nntot, length(m.bvectors.shells),
            length(m.bvectors.shells) == 1 ? "" : "s")
    @printf(io, "  cell volume               : %.4f Å³\n", cell_volume(m.lattice))
    print(io,   "  band energies (.eig)      : ", m.eig === nothing ? "absent" : "present")
end

function Base.show(io::IO, ::MIME"text/plain", s::SpreadResult)
    nw = length(s.spreads)
    println(io, "Wannier spread  Ω = ", @sprintf("%.9f", s.Ω), " Å²")
    @printf(io, "  Ω_I  = %14.9f   (gauge-invariant)\n", s.ΩI)
    @printf(io, "  Ω_OD = %14.9f   Ω_D = %.9f\n", s.ΩOD, s.ΩD)
    show_n = min(nw, 8)
    for n in 1:show_n
        @printf(io, "  WF %-2d  centre (% .5f, % .5f, % .5f) Å   spread %.6f Å²\n",
                n, s.centres[1, n], s.centres[2, n], s.centres[3, n], s.spreads[n])
    end
    nw > show_n && print(io, "  ⋮ (", nw - show_n, " more)")
end

function Base.show(io::IO, ::MIME"text/plain", r::WannieriseResult)
    println(io, "WannieriseResult: Ω = ", @sprintf("%.9f", r.spread.Ω), " Å² after ", r.niter,
            " iterations", r.converged ? " (converged)" : " (iteration limit)")
    show(io, MIME"text/plain"(), r.spread)
end

function Base.show(io::IO, ::MIME"text/plain", r::WannierResult)
    println(io, "WannierResult", r.disentangled ? " (disentangled)" : " (isolated bands)")
    @printf(io, "  Ω_I locked at %.9f Å²\n", r.omega_I)
    println(io, "  localiser: ", r.niter, " iterations",
            r.converged ? ", converged" : ", iteration limit")
    show(io, MIME"text/plain"(), r.spread)
end

function Base.show(io::IO, ::MIME"text/plain", d::DisentangleResult)
    nk = length(d.Uopt)
    println(io, "DisentangleResult: Ω_I = ", @sprintf("%.9f", d.omega_I), " Å² after ", d.niter,
            " iterations over ", nk, " k-points")
    n = min(3, length(d.omega_I_trace))
    for (it, wi1, wi, δ) in d.omega_I_trace[1:n]
        @printf(io, "  iter %-3d  Ω_I = %.8f   Δ = %+.2e\n", it, wi, δ)
    end
    length(d.omega_I_trace) > n && println(io, "  ⋮")
end
