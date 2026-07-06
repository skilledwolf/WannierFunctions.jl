# Migrating from Wannier90

If you already run `wannier90.x`, this package slots into the same place in your workflow. It
reads the **same input files** your DFT interface produces — `.win`, `.amn`, `.mmn`, `.eig` —
so nothing upstream changes. You can drive it two ways: as drop-in binaries writing the same
output files, or from Julia with results as native data structures instead of a `.wout` to
parse.

## The pipeline mapping

A typical Wannier90 run:

```
1. SCF + NSCF (DFT)                    → wavefunctions
2. wannier90.x -pp seedname            → seedname.nnkp
3. pw2wannier90 (or your interface)    → seedname.amn/.mmn/.eig (+ .spn/.uHu/.dmn as needed)
4. wannier90.x seedname                → .wout, .chk, _hr.dat, band files
5. postw90.x seedname                  → AHC, DOS, kpath, BoltzWann, … .dat files
```

Steps 1–3 are unchanged. Steps 2, 4, and 5 map one-to-one:

```bash
wannier90.jl -pp seedname     # step 2 (byte-identical .nnkp)
wannier90.jl seedname         # step 4 (same outputs, same formats)
postw90.jl seedname           # step 5 (same .dat files)
```

The commands come either from a repository clone (`julia --project=. bin/wannier90.jl …`) or,
for a `pkg> add`-installed package, from `using WannierFunctions; install_cli()`, which writes
these launchers to `~/.julia/bin`.

Or as a library, replacing steps 4–5 with data you can compute on:

```julia
using WannierFunctions

model = read_model("seedname")               # .win/.amn/.mmn/.eig
res   = run_wannier(model)                   # windows via keywords for entangled bands

H = hamiltonian_operator(model, res)         # the H(R) behind _hr.dat
E = bands(H, kpath)                          # replaces bands_plot

bm = BerryModel("seedname")                  # from the .chk you (or we) wrote
anomalous_hall(bm; fermi_energy=…, kmesh=…)  # replaces postw90's berry_task = ahc
```

Checkpoints interchange in both directions: `wannier90.x restart=plot` consumes our `.chk`,
and `BerryModel` consumes theirs. `bin/w90chk2chk.jl` converts `.chk ↔ .chk.fmt`.

## Where your `.win` keywords go

The short version: **they work**. The parser is strict (unknown keywords error with a
did-you-mean; recognised-but-unsupported ones warn once), and the supported set covers the
wannierise, disentanglement (energy windows, `dis_spheres`, PDWF projectability, symmetry-
adapted), plotting/output, transport, and all postw90 module keywords. The precise policy,
keyword semantics matched to the reference source, and the (short) list of behavioural
differences are in [Wannier90 compatibility](wannier90-compat.md).

| You used | Here |
|----------|------|
| `num_bands > num_wann` + `dis_win_*`/`dis_froz_*` | same keywords; `run_wannier` auto-selects disentanglement |
| `guiding_centres`, `precond`, `slwf_*`, `site_symmetry`, `use_ss_functional`, `gamma_only`, `higher_order_n` | same keywords, implemented and validated |
| `berry_task = ahc/morb/kubo/sc/shc/kdotp`, `gyrotropic`, `dos`, `kpath`, `kslice`, `geninterp`, `boltzwann`, `spin_moment` | `bin/postw90.jl`, or the per-module Julia API ([How-to](howto.md)) |
| `transport = true`, `transport_mode = bulk` | supported (`tran_lcr` is the one exclusion) |
| `wannier_plot`, cube/xsf, `write_rmn/_tb/_hr/_hr_diag/xyz`, `.bxsf` | supported, reference formats |

## What you get back (library route)

Instead of scraping a `.wout`:

- `res.spread.Ω, .ΩI, .ΩOD, .ΩD` — the `Omega I/D/OD/Total` lines, in Å².
- `res.spread.centres` (3 × num_wann, Cartesian Å) and `.spreads` — the per-WF lines.
- `res.U`, `res.omega_trace`, `res.niter`, `res.converged` — the gauge and convergence trace.
- Disentangled runs: `res.omega_I`, `res.dis.omega_I_trace`.

These reproduce the reference `.wout` numbers to the test-suite tolerances or better — see
[Validation](validation.md).

## Behavioural notes to expect

- **Convergence semantics differ by optimiser.** The `.win`/CLI path uses `:w90`, where —
  matching Wannier90 — the convergence check is off unless `conv_window > 1` and the loop runs
  the full `num_iter`. The Julia-native default `:rcg` has a real convergence criterion and
  `converged` means what it says. Same minima either way (asserted by tests).
- **Units.** Centres/spreads come back in Å / Å² even for a `bohr` cell, matching the `.wout`
  convention; energies are eV. The bohr constant is the reference default (CODATA2006).
- **`use_ws_distance` defaults to true** (as in postw90) — mind it when comparing against old
  runs that disabled it.
- **Γ-only** runs the reference's real-orthogonal algorithm and returns exactly real gauges.
- One step further than migration: with DFTK you can skip the file round-trip entirely —
  see [Getting started](getting-started.md), Workflow B.
