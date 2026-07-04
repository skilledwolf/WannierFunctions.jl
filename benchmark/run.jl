# Benchmark harness — plain wall-clock timings of the pipeline stages on the validation systems
# plus a dense interpolation workload. No external dependencies; run twice to compare threading:
#
#   julia -t 1 --project=. benchmark/run.jl
#   julia -t 8 --project=. benchmark/run.jl
#
# Physics invariants are printed alongside timings so a speed change can never silently be a
# correctness change.

import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")); io=devnull)
using Wannier90, StaticArrays, Printf

root = normpath(joinpath(@__DIR__, ".."))
silicon = joinpath(root, "scratch", "silicon", "silicon")
copper = joinpath(root, "scratch", "copper", "copper")
diamond = joinpath(root, "examples", "data", "diamond")

function timeit(f; warmup=true)
    warmup && f()
    t0 = time_ns()
    out = f()
    return (time_ns() - t0) / 1e9, out
end

@printf("Wannier90.jl benchmarks — %d thread(s)\n", Threads.nthreads())
@printf("%-46s %10s   %s\n", "workload", "seconds", "invariant")
println("-"^88)

# 1. Localisation, reference-faithful optimiser (copper: 12→7, 200 fixed sweeps).
if isfile(copper * ".mmn")
    mc = read_model(copper)
    wc = read_win(copper * ".win")
    t, r = timeit(() -> run_wannier(mc, wc))
    @printf("%-46s %10.3f   Ω = %.9f\n", "copper  disentangle + :w90 (200 sweeps)", t, r.spread.Ω)
    t, r = timeit(() -> run_wannier(mc; win_min=wc.dis_win_min, win_max=wc.dis_win_max,
                                    froz_min=wc.dis_froz_min, froz_max=wc.dis_froz_max,
                                    dis_num_iter=wc.dis_num_iter, dis_mix_ratio=wc.dis_mix_ratio))
    @printf("%-46s %10.3f   Ω = %.9f (converged=%s, %d its)\n",
            "copper  disentangle + :rcg (converged)", t, r.spread.Ω, r.converged, r.niter)
end

# 2. Disentanglement + localisation, silicon (12→8, frozen window).
if isfile(silicon * ".mmn")
    ms = read_model(silicon)
    ws = read_win(silicon * ".win")
    t, r = timeit(() -> run_wannier(ms, ws))
    @printf("%-46s %10.3f   Ω = %.9f\n", "silicon disentangle + :w90 (50 sweeps)", t, r.spread.Ω)
end

# 3. Dense band interpolation (diamond H on 20 000 k-points).
md = read_model(diamond)
wd = read_win(diamond * ".win")
res = run_wannier(md, wd)
H = hamiltonian_operator(md, res)
kdense = [SVector{3,Float64}(i / 29, j / 31, k / 23)
          for i in 0:28 for j in 0:30 for k in 0:22]        # 20 677 deterministic points
t, E = timeit(() -> bands(H, kdense))
@printf("%-46s %10.3f   ⟨E⟩ = %.6f eV\n", "diamond bands() on 20 677 k-points", t,
        sum(E) / length(E))

# 4. use_ws_distance interpolation (per-element phases, 2 000 k-points).
t, E2 = timeit(() -> interpolate(md, res, kdense[1:2_000]; use_ws_distance=true))
@printf("%-46s %10.3f   ⟨E⟩ = %.6f eV\n", "diamond ws-distance interp, 2 000 k-points", t,
        sum(E2) / length(E2))
