# Validation

Every feature in this package is anchored to an oracle. This page explains the methodology
and records the headline results; per-algorithm details (conventions, exact anchors, file
formats) live in the repository's `docs/reference-notes/`.

## Methodology

**Only gauge-invariant quantities are compared.** Band energies, spread components
(Ω, Ω_I, Ω_OD, Ω_D), Wannier centres, and physical response tensors are well-defined;
individual U(k) or H(R) matrix elements are gauge-dependent and are never used as validation
targets (except where the reference-exact `:w90` optimiser makes the whole gauge trajectory
reproducible, e.g. for byte-identical file comparisons).

Three classes of oracle, in decreasing order of strength:

1. **Shipped reference benchmarks** — the golden outputs of the Wannier90 test suite. Where a
   benchmark pins a *trajectory* (e.g. exactly 10 disentanglement iterations), the trajectory
   is matched iteration by iteration, not just the fixed point.
2. **Self-generated oracles** — for features without shipped tests (converged spreads,
   ballistic transport, `-pp` for new modes), a locally built `wannier90.x`/`postw90.x` is run
   on controlled inputs and its outputs compared byte-for-byte or at file precision.
3. **Cross-package validation** — where no Fortran reference exists at all (circular injection
   current), the same tight-binding model is evaluated by an independent code (WannierBerri)
   and the results compared, with the agreement level and its limiting factor stated.

The suite has **454 tests across 35 test sets**, passing identically at 1 and 8 threads.
Tests degrade gracefully: cases needing the reference clone or optional tools skip rather
than fail.

## Core pipeline

Spread components reproduce the reference to print precision (≤ 5·10⁻¹⁰), far inside the
test-suite tolerances:

| Case | System | Bands → WF | Disentangle | Final Ω (this / reference) |
|------|--------|-----------|-------------|----------------------------|
| example01 | GaAs | 4 → 4 | no | 4.466880976 / 4.466880976 |
| example05 | diamond | 4 → 4 | no | 2.320904915 / 2.320904915 |
| example03 | silicon | 12 → 8 | frozen window | 14.499574503 / 14.499574503 |
| example04 | copper (metal) | 12 → 7 | frozen window | 4.028040058 / 4.028040058 |

For silicon the Ω_I convergence trace matches the reference **iteration by iteration**.
Against a locally built `wannier90.x`: interpolated bands to ~10⁻⁵ eV (file precision),
`_hr.dat` byte-identical, `.nnkp` byte-identical, `.chk` interchange working in both
directions (its `restart=plot` reproduces its own bands from our checkpoint).

## postw90 physics

Each module is validated against the corresponding test-suite oracle (representative
anchors; all in the suite):

| Module | Anchor |
|--------|--------|
| AHC (adaptive + Fermi scan) | Fe: (0.0334, 0.0572, 1222.1510) S/cm digit-exact; 33/33 scan values ≤ 5·10⁻⁷ |
| Kubo σ(ω) + JDOS | Fe: `kubo_S/A`/`jdos` files byte-identical (E16.8) |
| Orbital magnetisation | Fe: 0.0431 μ_B/cell = benchmark; `transl_inv_full` variant 0.0415 |
| Spin Hall (Qiao/Ryoo/tetrahedron) | Pt: 201 Fermi levels ≤ 4·10⁻⁸ rel; Ryoo and tetrahedron oracles exact |
| Shift current | GaAs: five oracles ≤ 2·10⁻⁷ abs (harness tolerance 10⁻⁶) |
| k·p | GaAs: order-0 byte-identical |
| Gyrotropic | Te: five of six oracle files byte-identical |
| DOS / pDOS / spin DOS | copper 5·10⁻⁸; example04 pdos oracle |
| BoltzWann | silicon: σ/S/κ/TDF vs oracle (see caveat in [compatibility](wannier90-compat.md)) |
| kpath / kslice | Fe/Pt: `path.kpt`, `bands.dat` and task files **byte-identical** |
| geninterp | silicon: file rows identical (Fortran G-format emulation) |
| spin moment | Fe: 3.090787 μ_B digit-exact |

The `bin/postw90.jl` driver as a whole was validated by running complete `.win`-driven jobs
and byte-comparing every output file against local `postw90.x` runs across eleven task
families.

## Advanced wannierisation

| Feature | Oracle and result |
|---------|-------------------|
| SLWF+C | Ω_C = 1.634087566 vs oracle 1.634087565; reference gradient vanishes at their optimum |
| Symmetry-adapted WFs | GaAs Ω = 10.136492662 exact; gauge symmetry-adapted to 7·10⁻¹² |
| SAWF + disentanglement | H3S: 10-iteration Ω_I trajectory exact to every printed digit; converged Ω 6.301957278 vs self-generated 6.301957261 |
| PDWF | graphene Ω = 15.803350 (3·10⁻⁷) |
| dis_spheres | LaVO₃ Ω = 7.508128 (2·10⁻¹¹) |
| Γ-only (real-orthogonal) | all four gamma oracles to every printed digit; gauges exactly real |
| Higher-order finite differences | KNbO₃ spread to 9 digits; Fe morb & Pt SHC higher-order oracles exact; `-pp` byte-identical |
| Stengel–Spaldin | state-function parity to 9–10 digits (see [compatibility](wannier90-compat.md) for why stopping points are not comparable) |
| Guiding centres / precond / select_projections | silicon/GaAs oracles |

## Beyond-parity features

| Feature | Validation |
|---------|-----------|
| Ballistic transport (`tran_bulk`) | self-generated `wannier90.x` oracle on Te: T(E) to 9·10⁻⁷ over the full scan, H00/H01 at file precision |
| Irreducible-BZ integration | wedge AHC = full-BZ to 8·10⁻⁶; Fe morb to 2·10⁻⁶ with 78/512 k-points; DOS exact; k-reduction eigenvalue multiset exact |
| Injection current | WannierBerri cross-validation on a shared TB model, ~10⁻⁴ agreement (degeneracy-regularisation limited) |
| TB-model input | H(R)/r(R) file round-trip exact; Fe AHC reproduced from `_tb.dat` alone |
| DFTK bridge | live silicon SCF: interpolated bands reproduce SCF eigenvalues to 10⁻¹² eV; SCDM and explicit projections reach the identical minimum |

## Where the details live

- `docs/reference-notes/parity-audit-2026.md` — the triaged survey of the reference feature
  surface and its status.
- `docs/reference-notes/remaining-gaps-2026-07.md` — conventions and anchors for the final
  gap-closing slate.
- One note per algorithm alongside them (kmesh, disentanglement, Berry modules, transport,
  SHC variants, …) at implementation depth: exact Fortran formats, sign conventions, and the
  quirks discovered on the way.
