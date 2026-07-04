# Example 4 — the Julia-native API: keyword-first pipeline, operators, and -pp setup.
#
# Where examples 1–3 mirror the classic wannier90.x workflow, this one shows the surface a new
# Julia project would use directly. Run from the repository root:
#
#   julia --project=. examples/04_native_api.jl
#
using WannierFunctions
using Printf

seed = joinpath(@__DIR__, "data", "diamond")

# --- 1. The model, with a readable REPL display -------------------------------------------
model = read_model(seed)
display(model); println()

# --- 2. One keyword-first call; :rcg converges instead of running fixed sweeps ------------
res = run_wannier(model)                       # isolated case; windows would be keywords here
@printf("\nΩ = %.9f Å²  after %d :rcg iterations (converged = %s)\n",
        res.spread.Ω, res.niter, res.converged)

# --- 3. Operators: the Hamiltonian and the position operator are the same abstraction -----
H = hamiltonian_operator(model, res)
r = position_operator(model, res)
display(H); println(); display(r); println()

Hk = H([0.25, 0.0, 0.0])                       # H(k) at an arbitrary fractional k
E = bands(H, [[0.0, 0.0, 0.0], [0.5, 0.0, 0.5]])
@printf("\nvalence bands at Γ: %s eV\n", join((@sprintf("%.4f", e) for e in E[:, 1]), "  "))

# diag r(R=0) is the Wannier centres — the operator view and the spread view agree:
ir0 = findfirst(==((0, 0, 0)), r.irvec)
@printf("WF1 centre from r(0):   (% .6f, % .6f, % .6f) Å\n",
        real(r.data[1, 1, ir0, 1]), real(r.data[1, 1, ir0, 2]), real(r.data[1, 1, ir0, 3]))
@printf("WF1 centre from spread: (% .6f, % .6f, % .6f) Å\n",
        res.spread.centres[1, 1], res.spread.centres[2, 1], res.spread.centres[3, 1])

# --- 4. Starting a NEW calculation: -pp / .nnkp generation from the .win alone ------------
out, info = generate_nnkp(seed; out = joinpath(mktempdir(), "diamond.nnkp"))
@printf("\n-pp: wrote %s  (%d neighbours/k, shell weights %s Ų)\n",
        basename(out), info.nntot, string(round.(info.weights, digits = 6)))

# --- 5. Strict input validation: typos are errors, not silence ----------------------------
winpath = joinpath(mktempdir(), "typo.win")
write(winpath, replace(read(seed * ".win", String), "num_iter" => "num_itre"))
try
    read_win(winpath)
catch err
    println("\nstrict parsing: ", split(sprint(showerror, err), "\n")[1])
end
