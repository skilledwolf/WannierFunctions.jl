# Physical constants and tolerances, matched to the reference Wannier90 (src/constants.F90).

"""
Bohr radius in Ångström. The reference Wannier90 defaults to **CODATA2006** (compile flags select
2010/2018/2022); we match the default so bohr-specified cells reproduce reference output exactly.
"""
const BOHR = 0.52917720859

"2π."
const TWOPI = 2 * π

# SI / atomic-unit constants (CODATA2006, matching the reference default build's
# constants.F90). Every response module converts with these — keep them in one place.
const ELEM_CHARGE_SI = 1.602176487e-19    # C
const HBAR_SI        = 1.054571628e-34    # J·s
const ELEC_MASS_SI   = 9.10938215e-31     # kg
const EPS0_SI        = 8.854187817e-12    # F/m
const KB_SI          = 1.3806504e-23      # J/K
const EV_AU          = 3.674932540e-2     # eV → Hartree (constants.F90:178)
const EV_SECONDS     = 6.582119e-16       # ħ/e in eV·s — CODATA2006 set (7 digits!)

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

"""
    threaded_ksum(body!, init, n) -> Vector{typeof(init())}

Chunked parallel reduction over `1:n`: the range is split into one chunk per thread, each
chunk task gets its own state from `init()` (accumulators + scratch buffers), and
`body!(state, i)` runs sequentially within a chunk. Returns the per-chunk states for the
caller to reduce. Falls back to a single sequential chunk for small `n` or one thread —
the same size gating as [`@maybe_threads`](@ref). Compared to a `Threads.@threads` loop
with per-`i` allocations, this bounds memory at `nthreads` states (not `n`) and lets hot
per-k scratch be reused across the whole chunk.
"""
function threaded_ksum(body!::F, init::G, n::Int) where {F,G}
    nchunks = (Threads.nthreads() > 1 && n >= THREAD_MIN) ? min(Threads.nthreads(), n) : 1
    if nchunks == 1
        return [_ksum_chunk!(body!, init(), 1:n)]
    end
    chunks = collect(Iterators.partition(1:n, cld(n, nchunks)))
    states = [init() for _ in 1:length(chunks)]
    @sync for (ci, chunk) in enumerate(chunks)
        Threads.@spawn _ksum_chunk!(body!, states[ci], chunk)
    end
    return states
end

# Function barrier: keeps the @spawn closure tiny (three read-only captures) — large
# inlined threaded bodies are a closure-boxing hazard on Julia 1.12.
function _ksum_chunk!(body!::F, state, chunk) where {F}
    for i in chunk
        body!(state, i)
    end
    return state
end
