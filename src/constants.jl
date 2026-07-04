# Physical constants and tolerances, matched to the reference Wannier90 (src/constants.F90).

"""
Bohr radius in Ångström. The reference Wannier90 defaults to **CODATA2006** (compile flags select
2010/2018/2022); we match the default so bohr-specified cells reproduce reference output exactly.
"""
const BOHR = 0.52917720859

"2π."
const TWOPI = 2 * π

# Default numerical tolerances (mirror Wannier90 defaults).
const KMESH_TOL_DEFAULT = 1.0e-6      # shell-distance / B1 degeneracy tolerance
const CONV_TOL_DEFAULT  = 1.0e-10     # spread-convergence tolerance (Ų)

"""
    @maybe_threads cond for ... end

Thread the loop only when `cond` is true (and more than one thread is available). Threading a
64-iteration loop of microsecond bodies costs more in scheduling than it saves — every threaded
k/R loop in this package is gated on a problem-size condition.
"""
macro maybe_threads(cond, ex)
    esc(quote
        if $cond && Threads.nthreads() > 1
            Threads.@threads $ex
        else
            $ex
        end
    end)
end

"Default minimum loop length before a k/R loop is worth threading."
const THREAD_MIN = 128
