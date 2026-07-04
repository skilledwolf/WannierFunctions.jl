# Theory

A derivation-level overview of the physics implemented by this package: the maximally-localised
Wannier function (MLWF) spread functional and its minimisation, the finite-difference geometry
that makes it computable on a discrete k-mesh, disentanglement of composite bands, and Wannier
interpolation. The conventions here match the code and the reference Wannier90 v3.1.0; the
implementation notes in `docs/reference-notes/` cite the reference source line-by-line.

References cited throughout:

- **[MV97]** N. Marzari and D. Vanderbilt, Phys. Rev. B **56**, 12847 (1997).
- **[SMV01]** I. Souza, N. Marzari, and D. Vanderbilt, Phys. Rev. B **65**, 035109 (2001).
- **[RMP12]** N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, and D. Vanderbilt,
  Rev. Mod. Phys. **84**, 1419 (2012).

---

## 1. Wannier functions and the gauge freedom

For an isolated group of `J` bands, the Wannier functions (WFs) in the home cell are

```
|w_{nR}⟩ = (V/(2π)³) ∫_BZ dk  e^{-i k·R} Σ_m U^{(k)}_{mn} |ψ_{mk}⟩ ,
```

where `|ψ_{mk}⟩` are the Bloch states and `U^{(k)}` is a **k-dependent unitary** mixing the `J`
bands at each k. The physical content (charge density, total energy) is invariant under this
gauge choice, but the *shape and localisation* of the WFs are not. MLWFs are the choice of
`U^{(k)}` that makes the WFs maximally localised [MV97, RMP12].

## 2. The spread functional Ω = Ω_I + Ω_OD + Ω_D

Localisation is quantified by the **quadratic spread** [MV97, Eq. (11)]

```
Ω = Σ_n [ ⟨r²⟩_n − |⟨r⟩_n|² ] ,      ⟨O⟩_n ≡ ⟨w_{n0}| O |w_{n0}⟩ .
```

Marzari and Vanderbilt split `Ω` into a **gauge-invariant** part and a **gauge-dependent**
remainder [MV97, Eqs. (34)–(36)]:

```
Ω = Ω_I + Ω̃ ,        Ω̃ = Ω_OD + Ω_D .
```

- **Ω_I** (invariant): independent of `U^{(k)}`. It is a lower bound on the total spread and is
  the quantity minimised by disentanglement (§5). Physically it measures how much the WF
  subspace fails to be smooth/self-overlapping across neighbouring k-points.
- **Ω_OD** (off-diagonal): the sum over off-diagonal `|M|²`; driven to a minimum by gauge
  rotation.
- **Ω_D** (diagonal): the mismatch between the WF centre implied by the diagonal overlap phases
  and the centre itself. It vanishes exactly for a suitably symmetric problem (e.g. diamond).

Only `Ω̃ = Ω_OD + Ω_D` is minimised in the MLWF gauge step; `Ω_I` is fixed once the subspace is
fixed.

## 3. Finite-difference b-vectors and the B1 condition

On a discrete Monkhorst–Pack mesh the position expectation values are not computed from `∇_k`
directly; instead MV97 (Appendix B) use a **finite-difference** approximation built from
neighbour shells. For each k-point one selects a set of vectors `b` connecting it to nearby
k-points, with weights `w_b`, chosen so that

```
Σ_b w_b  b_α b_β = δ_{αβ}        (the "B1 condition", MV97 Eq. (B1)).
```

This is a completeness relation: it guarantees the finite-difference gradient reproduces the
continuum limit to first order. In this package the shells are found by grouping neighbour
distances, and the weights `w_b` are solved from the B1 linear system (in Å⁻² so that spreads
come out in Å²). The reference reports satisfaction of B1 as
*"Completeness relation is fully satisfied [Eq. (B1), PRB 56, 12847 (1997)]."*

The only overlap the algorithm needs is

```
M^{(k,b)}_{mn} = ⟨u_{mk} | u_{n,k+b}⟩ ,
```

read verbatim from the `.mmn` file. Everything below is expressed in terms of `M` and `w_b, b`.

### Centres and spread in finite differences

With the principal branch of the complex logarithm (`Im ln z = atan2(Im z, Re z) ∈ (−π, π]`),
the WF centre is [MV97 Eq. (31), RMP12 Eq. (28)]

```
r̄_n = −(1/N_k) Σ_{k,b} w_b  b  Im ln M^{(k,b)}_{nn} ,
```

and the mean-square extent is

```
⟨r²⟩_n = (1/N_k) Σ_{k,b} w_b [ 1 − |M^{(k,b)}_{nn}|² + (Im ln M^{(k,b)}_{nn})² ] .
```

The spread of WF `n` is `⟨r²⟩_n − |r̄_n|²`. The three components are

```
Ω_I  = (1/N_k) Σ_{k,b} w_b ( J − Σ_{mn} |M^{(k,b)}_{mn}|² )        (invariant)
Ω_OD = (1/N_k) Σ_{k,b} w_b Σ_{m≠n} |M^{(k,b)}_{mn}|²
Ω_D  = (1/N_k) Σ_{k,b} w_b ( Im ln M^{(k,b)}_{nn} + b·r̄_n )²
```

with `J = num_wann`. The branch cut is the single most delicate numerical choice: a wrong
branch shifts `Ω_D` while leaving `Ω_I` and `Ω_Total` looking plausible.

## 4. Gauge gradient and the unitary update

Write the gauge update at each k as `U^{(k)} → U^{(k)} e^{ΔW^{(k)}}` with `ΔW` **anti-Hermitian**
(so the exponential is unitary). The gradient of `Ω̃` with respect to the anti-Hermitian
generator is [MV97 Eqs. (52)–(57)]

```
G^{(k)} = (4/N_k) Σ_b w_b ( A[R^{(k,b)}] − S[T^{(k,b)}] ) ,
```

with the auxiliary matrices

```
R_{mn}  = M_{mn} · conj(M_nn)
R̃_{mn} = M_{mn} / M_nn
T_{mn}  = R̃_{mn} · q_n ,     q_n = Im ln M_nn + b·r̄_n ,
A[X] = (X − X†)/2 ,          S[X] = (X + X†)/(2i) .
```

`A[·]` is the anti-Hermitian part and `S[·]` the Hermitian-part-over-`i`; both keep `G`
anti-Hermitian. The minimisation is a **Fletcher–Reeves conjugate gradient**: the search
direction is `d = G + γ d_prev` with `γ = ‖G‖²/‖G_prev‖²` (reset to steepest descent every few
steps or if `γ` grows too large), and the step length is chosen by a **parabolic line search**
using the spread at the current point and at a trial step.

The update itself is exact unitarity, not a re-orthonormalisation: with generator `ΔW`,

```
e^{ΔW} = V diag(e^{−i λ}) V† ,   where  i ΔW = V diag(λ) V†  (λ real),
```

so `U ← U e^{ΔW}` stays unitary to machine precision. The overlaps are rotated consistently,
`M^{(k,b)} ← e^{ΔW_k}† M^{(k,b)} e^{ΔW_{k+b}}`, using the rotation at **both** k and its
neighbour `k+b`.

By default the reference disables the convergence-window check and simply runs the requested
`num_iter` iterations; this package follows the same behaviour.

## 5. Disentanglement (Souza–Marzari–Vanderbilt)

When the `J = num_wann` bands of interest are **entangled** with other bands (metals, or a
conduction manifold crossing the target group), there is no isolated set of `J` bands to
Wannierise. SMV01 solves this by first choosing, at every k, an optimal `J`-dimensional
subspace out of the larger set of `num_bands` states inside an **outer energy window**,
optionally with an inner **frozen window** whose states are kept exactly.

The subspace is chosen to **minimise Ω_I** — equivalently, to maximise the smoothness
(mutual overlap) of the subspace across neighbouring k-points [SMV01 Eqs. (11)–(20)]. Because
`Ω = Ω_I + Ω̃` and `Ω_I` is exactly the gauge-invariant piece of §2, this cleanly separates the
problem: **(1)** pick the subspace minimising `Ω_I` (disentanglement, an iterative
self-consistent eigenproblem at each k), then **(2)** run the ordinary MLWF minimisation of
`Ω̃` (§4) inside that fixed subspace. The output is a rectangular `num_bands × num_wann`
embedding `U_opt^{(k)}` plus a square initial `U^{(k)}` for step (2).

This is implemented and validated: `run_wannier` auto-selects it when `num_bands > num_wann`. On
silicon (12 → 8, frozen window) the Ω_I convergence trace matches the reference iteration by
iteration and the final spread to ~1e-8; copper (a metal, 12 → 7) matches to ~1e-7. The isolated
path (`num_bands == num_wann`) skips this step entirely.

## 6. Wannier interpolation H(R) → H(k)

Once the gauge is fixed, the Wannier-gauge Hamiltonian at each mesh point is

```
H^{(W)}(k) = U^{(k)†} diag(ε_{1k}, …, ε_{Jk}) U^{(k)} .
```

Fourier-transforming to the real-space (tight-binding) representation on the set of
Wigner–Seitz lattice vectors `{R}`:

```
H(R) = (1/N_k) Σ_k e^{−i 2π k·R} H^{(W)}(k) .
```

Because the WFs are localised, `H(R)` decays rapidly with `|R|`, so this handful of matrices is
a compact tight-binding model. To interpolate onto any k′ (e.g. a band path) one inverts the
transform:

```
H(k′) = Σ_R (1/N_deg(R)) e^{+i 2π k′·R} H(R) ,
```

and diagonalising `H(k′)` gives the interpolated band energies (ascending). Two conventions are
load-bearing and are reproduced exactly:

- **Fourier sign asymmetry.** The forward transform (k → R) carries `e^{−i2πk·R}`; the inverse
  (R → k′) carries `e^{+i2πk′·R}`.
- **Wigner–Seitz degeneracies.** `{R}` is the set of lattice vectors lying in (or on the
  boundary of) the Wigner–Seitz cell of the Born–von-Kármán supercell defined by `mp_grid`.
  Boundary vectors are shared among `N_deg(R)` equivalent images; each `H(R)` is stored
  *undivided* and the weight `1/N_deg(R)` is applied at interpolation. The set satisfies the sum
  rule `Σ_R 1/N_deg(R) = ∏ mp_grid`, which the code checks.

The reference's default `use_ws_distance` minimal-image refinement (a per-matrix-element shift
that further improves smoothness) is not yet applied here; it is on the roadmap.

---

For exact array conventions, index orders, branch-cut handling, and file:line citations to the
reference, see `docs/reference-notes/localization.md`, `interpolation.md`, and
`disentanglement.md`.
