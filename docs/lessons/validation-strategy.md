# Validate on gauge-invariant quantities only

Wannier functions are gauge-dependent: two *correct* implementations can produce different U(k)
matrices and different `_hr.dat` H(R) matrix elements. Diffing those against the reference is a
phantom-bug trap. Compare only gauge-invariant quantities:

- **Interpolated band energies on a k-path** — eigenvalues of H(k). Fully gauge-invariant. The
  gold-standard oracle for the interpolation stage.
- **Total spread Ω and its Ω_I / Ω_OD / Ω_D decomposition** at convergence — deterministic when the
  iteration starts from the *reference* `.amn` (same starting gauge).
- **Wannier centres** — well-defined at the spread minimum, modulo lattice vectors.

Not U(k), not raw H(R) elements. Their *eigenvalues*, yes; the elements, no.

## Milestone ladder (each validates the previous subsystem)

- **M0** — parse `.win/.amn/.mmn/.eig`; compute the *initial* spread from projected M-matrices;
  match the reference `.wout` "Initial State" spread. Validates I/O + b-vector shells/weights +
  spread evaluation with zero iteration.
- **M1** — localization iteration → match final Ω, the Ω_I/Ω_OD/Ω_D decomposition, and centres. The
  `.wout` reports Ω per iteration: use that trace as a debugging ladder (iter-1 off → gradient/M-update
  bug; iter-1 matches, iter-50 diverges → step-size/convergence).
- **M2** — Wannier interpolation → match reference bands on the example's k-path.
- **M3** — add disentanglement (Souza–Marzari–Vanderbilt) against an entangled/metal example.
- **M4** — breadth, docs, examples.

## First target = isolated-bands example (num_bands == num_wann)

No disentanglement, no frozen window. Canonically a valence-bands case (verify which shipped example
is simplest and has benchmark outputs — don't assume the example number). Disentanglement is the
second-hardest algorithm; adding it in step one triples the debugging surface.

## Tooling notes

- **No DFT run needed**: the test-suite ships `.mmn/.amn/.eig` inputs AND benchmark outputs. Extract
  the harness's exact tolerances — that's the pass/fail threshold; don't guess it.
- **gfortran = live oracle**: build the reference binary to dump *intermediate* quantities
  (per-iteration spreads, chosen shells/weights, centres) when M1/M2 diverge. Worth it at a
  divergence, not before.
- **Remove the kmesh variable first**: b-vector shell finding + B1 finite-difference weights
  (Σ_b w_b b_α b_β = δ_αβ) is a known rabbit hole. The `.nnkp` gives k+b connections; the `.wout`
  reports chosen shells/weights. Validate kmesh against those before trusting anything downstream.

## Scope

Core milestone = `wannier90.x` (wannierization + Wannier interpolation). `postw90.x` (BoltzWann,
Berry/AHC, gyrotropic, spin Hall, geninterp) is a large separate surface — explicitly out of scope
for the core milestone.
