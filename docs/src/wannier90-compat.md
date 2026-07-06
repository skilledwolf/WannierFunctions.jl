# Wannier90 compatibility

This page is the contract: which binaries and keywords are covered, which conventions are
followed, and every known behavioural difference from the reference implementation.

## Binaries

| Reference | Here | Notes |
|-----------|------|-------|
| `wannier90.x -pp` | `bin/wannier90.jl -pp` | `.nnkp` byte-identical (incl. spinor projections and `higher_order_n` blocks) |
| `wannier90.x` | `bin/wannier90.jl` | `.wout`, `.chk`, `_hr.dat`, `_tb.dat`, band files, cube/xsf/`_r.dat`/`.bxsf`/xyz writers, transport |
| `postw90.x` | `bin/postw90.jl` | all modules; reference-named `.dat` outputs in the reference formats (validated byte-identical per module) |
| `w90chk2chk.x` | `bin/w90chk2chk.jl` | `-export`/`-import`; data-exact (the header date is restamped) |
| — (library mode / C API) | the Julia package itself | a `wannier90` v4-API → Julia-API mapping is the intended equivalent |

## Keyword policy

The `.win` parser is **strict**: unknown keywords are an error with a did-you-mean
suggestion, checked against a catalogue generated from the reference parser's own source.
Keywords fall into three classes:

- **Supported** — consumed with reference semantics (the large majority: the wannierise,
  disentanglement, plotting, transport and all postw90 module keywords).
- **Ignored** — cosmetic/reporting controls (`num_print_cycles`, `iprint`, `timing_level`, …):
  accepted silently.
- **Recognised but unsupported** — valid wannier90 keywords whose feature is not implemented
  here; these warn once and are ignored, so reference input decks still run.

Notable keyword semantics matched to the reference (discovered from its source, not its
documentation):

- `berry_task`/`gyrotropic_task` match by **substring** (`eval_shc` enables shc; glued
  `-C-dos-D0` works); task lists split on `+`.
- The global interpolation mesh keyword is `kmesh`; module-specific `*_kmesh` override it.
- `kubo_eigval_max` defaults to `dis_froz_max + 2/3` when a frozen window exists, else
  `max(ε) + 2/3`.
- BoltzWann's TDF energy window is the **disentanglement window** ± 0.2 eV; its DOS range
  defaults to the eigenvalue range ± 0.6667 eV.
- `dis_proj_min/dis_proj_max` (PDWF): `proj_min` is the *pool* threshold (discard below),
  `proj_max` the *freeze* threshold.
- Projection blocks accept both `;` and `,` as orbital separators (`O:s,p` ≡ `O:s;p`).

## Conventions

- **Bohr radius**: CODATA2006 (0.52917720859 Å) — the reference default. With it, validated
  spreads match reference prints to ≤ 5·10⁻¹⁰.
- `use_ws_distance` defaults to **true** (as in postw90) — an ~11% effect on transport
  coefficients in the validation set; disable explicitly to compare with older runs.
- `transl_inv_full` phases and corrections act on the *expanded* minimal-image R-set (the
  reference's `operator_wigner_setup`), matching its numerics exactly.
- Γ-only mode runs the real-orthogonal reference algorithm (`wann_main_gamma`) and returns
  exactly real gauges; the model-level API can optionally minimise over the full unitary
  group instead (marginally lower Ω, complex gauge).

## Known differences and caveats

Everything on this list is deliberate and documented; nothing else is known to differ.

- **`tran_lcr` (lead–conductor–lead transport) is not implemented** — the one reference
  feature that is out of scope (no shipped oracle; dedicated NEGF codes serve this use case).
  `transport_mode = bulk` is fully supported and validated.
- **Upstream bug (reported herein): `tran_read_ht`** — the reference v4-dev binary segfaults
  on transport-from-file inputs; our `read_ht`/`transport_bulk` path works and was validated
  against the reference's full pipeline instead.
- **Upstream bug: Ryoo SHC with `shc_gamma = 2, 3`** — the reference accumulates the spin
  σ-rotation across polarisations in `get_SAA_R`/`get_SBB_R` (get_oper.F90), mixing spin
  components. We compute clean per-component values; results agree for `shc_gamma = 1`.
- **Stengel–Spaldin functional**: the SS objective surface is a near-flat valley, so the
  *stopping point* depends on the convergence criterion even between two `wannier90.x` runs
  (13.845 vs 13.312 Å² on the reference test with different `conv_tol`). We validate state
  functions (objective and gradient at shared gauges, to 9–10 digits) rather than the
  stopping point, and run the reference optimiser so trajectories track.
- **Injection current** has no Fortran reference; it is cross-validated against WannierBerri
  on a shared tight-binding model. Agreement is at the ~10⁻⁴ level, limited by the two codes'
  degenerate-state regularisation conventions.
- **BoltzWann TDF** shows ≤ 10⁻⁴-relative deviations at isolated band-crossing energies
  (degenerate-band velocity conventions); all other BoltzWann outputs are identical to the
  last printed digit.
- **kdotp orders 1–2 are gauge-dependent** (the reference's own test checks only order 0,
  which matches byte-identically); magnitudes agree to ~10⁻⁶ relative.
- `.wout`/stdout text: the CLI writes a faithful `.wout` for the wannierise stage;
  `postw90.jl` prints a compact log rather than replicating the `.wpout` prose (its `.dat`
  files are the validated deliverables).

## Feature-by-feature validation anchors

The complete list of what was validated against which oracle, with numbers, lives in
[Validation](validation.md); implementation-depth notes (conventions, oracle anchors, exact
Fortran formats) are in the repository's `docs/reference-notes/` directory — one note per
algorithm, including the parity audit (`parity-audit-2026.md`) and the gap-closing slate
(`remaining-gaps-2026-07.md`).
