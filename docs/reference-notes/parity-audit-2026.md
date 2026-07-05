# Parity audit (July 2026): WannierFunctions.jl vs the field

> **Status update (2026-07-04, after slate 2 + 3).** Everything in §1 (source-verified
> wannier90/postw90 gaps) marked *do-now* or *do-later* with a shipped oracle test is now
> **implemented and oracle-validated**, plus several *do-later* items: shift current, kdotp,
> cube/`_r.dat`/`.bxsf`/`write_hr_diag`/`write_xyz`, PDWF, dis_spheres, guiding centres +
> select_projections, preconditioned CG, TB-model input (+ 3-D tabulation/FermiSurfer), and
> **tetrahedron-method SHC**. The remaining wannier90 gaps are SLWF+C, symmetry-adapted WFs
> (`sitesym`/`.dmn`), the Stengel–Spaldin functional, and ballistic transport (skip). The
> strategic ecosystem items (DFTK bridge, irreducible-BZ symmetrisation, injection current)
> are unstarted by design. The triage below is the original survey, kept for reference.

Triaged gap list against three sources, surveyed 2026-07-04:

1. **wannier90 v4-dev** (reference Fortran at `/Users/wolft/Dev/wannier90_greenfield/reference/wannier90/`,
   source-verified with file:line pointers and test-suite oracle names).
2. **WannierBerri v1.7.0** (Feb 2026; wannier-berri.org, docs.wannier-berri.org calculator list).
3. **Wannier.jl v0.3.6** (Nov 2025; github.com/qiaojunfeng/Wannier.jl + WannierIO.jl).

**Excluded — already in progress, not gaps:** spin module (moment, spin-decomposed DOS),
kpath module (bands/curv/morb/shc colouring), kslice morb+shc tasks, Ryoo SHC (.sHu/.sIu),
gyrotropic (D0/Dw/C/K/NOA/dos), projected DOS, transl_inv.

**Already covered (verified against source/tests):** wannier90.x core (MV localization
:rcg/:w90, SMV disentanglement, Gamma-only, spinors, -pp/.nnkp, .chk both ways, SCDM,
UNK/xsf, hr/tb/band, CLI); postw90 AHC+adaptive+fermiscan, Berry curvature, Kubo
sigma/JDOS, morb, geninterp, DOS (adaptive), BoltzWann, SHC (Qiao), kslice fermi
lines+curvature.

Positioning honesty: we already beat all three packages on oracle-validated file-precision
parity with the Fortran reference — none of them reproduce wannier90 outputs byte-for-byte.
WannierBerri beats us on BZ-integration performance (irreducible-k + symmetrization +
recursive refinement) and nonlinear-response breadth. Wannier.jl beats us on
gauge-research features (parallel transport, manifold splitting) and Julia-ecosystem
integration (DFTK). The slate below is chosen to (a) close the honest gaps in the
"full wannier90/postw90" claim cheaply, (b) add the one high-value physics feature we
lack (shift current), and (c) start on the two strategic differentiators.

---

## 1. Source-verified remaining wannier90/postw90 features

### 1.1 postw90: shift current — `berry_task = sc`  ★ largest physics gap

- `src/postw90/berry.F90:285` (`eval_sc`), `berry_get_sc_klist`; Ibañez-Azpiroz, Tsirkin,
  Souza PRB 97, 245143 (2018). Needs the generalized derivative of the position operator
  (uses `get_BB_R`/`get_CC_R`-adjacent machinery plus `sc_phase_conv` phase choices).
- Keywords (defaults from `src/postw90/postw90_types.F90:170-173`):
  `sc_phase_conv = 1`, `sc_eta = 0.04` (eV), `sc_w_thr = 5.0` (in units of `sc_eta`),
  `sc_use_eta_corr = .true.` (read at `postw90_readwrite.F90:962`). Frequency grid shares
  `kubo_freq_*`/`kubo_nfreq`; requires a single Fermi energy (`berry.F90:340-344`).
- Output: 18 files `seedname-sc_abc.dat`, `a = achar(119+i)` ∈ {x,y,z},
  `(b,c)` over the 6 symmetric pairs `alpha_S/beta_S` (`berry.F90:1656-1658`); per line
  `write(file_unit,'(2E18.8E3)') ω, fac*sc_list(i,jk,ifreq)` with
  `fac = eV_seconds*π*e³/(4·ħ²·V_cell)` (`berry.F90:1644`).
- Oracles: `testpostw90_gaas_sc_xyz`, `_sc_xyz_ws`, `_sc_xyz_scphase2`,
  `_sc_xyz_scphase2_ws`, `_sc_eta_corr` — five tests covering both phase conventions,
  use_ws_distance, and the eta correction. Best oracle coverage of any missing feature.

### 1.2 postw90: k·p expansion coefficients — `berry_task = kdotp`

- `berry.F90:287,502-510,583` (`berry_get_kdotp`): quasi-degenerate (Löwdin) perturbation
  theory around `kdotp_kpoint` for `kdotp_bands`. Keywords: `kdotp_kpoint` (default 0 0 0),
  `kdotp_num_bands`/`kdotp_bands` (`postw90_readwrite.F90:970-1016`).
- Outputs `seedname-kdotp_0.dat`, `-kdotp_1.dat`, `-kdotp_2.dat`, format `(2E18.8E3)`
  (`berry.F90:1743-1755`) — order-0/1/2 coefficient matrices.
- Oracle: `testpostw90_gaas_kdotp`.

### 1.3 postw90: tetrahedron method (SHC only) — `tetrahedron.F90`

- Ghim & Park PRB 106, 075126 (2022), with Kawamura PRB 89, 094515 correction matrix
  (`tetrahedron_P_matrix_init`, `tetrahedron.F90:40-76`). Wired only into SHC:
  `berry_get_shc_tetrahedron` (`berry.F90:1076`); anything else errors.
- Keywords: `tetrahedron_method`, `tetrahedron_higher_correction` (must be `.true.` —
  plain tetrahedron errors out, `berry.F90:488`), `tetrahedron_cutoff`,
  `tetrahedron_avoid_degeneracy` (`postw90_types.F90:184-188`).
- Oracles: `testpostw90_pt_tetra_shcfermi`, `testpostw90_pt_tetra_shcfreq`.

### 1.4 wannier90.x: wannierise options

- **Guiding centres**: `guiding_centres` (default `.false.`), `num_guide_cycles`
  (default 1), `num_no_guide_iter` (default 0) — `wannier90_types.F90:178-180`,
  `wannier90_readwrite.F90:696-718`. Projection sites become phase guides for the
  branch-cut choice in the spread. Oracle: `testw90_guidingcentre_selectproj`.
- **Selective localization + constrained centres (SLWF+C)**: `slwf_num` (default
  `num_wann`; `< num_wann` switches on `selective_loc`), `slwf_constrain` (default
  `.false.`), `slwf_lambda` (default 1.0), `slwf_centres` block —
  `wannier90_readwrite.F90:745-784`, `wannier90_types.F90:161-175`. Wang et al.
  PRB 90, 165125 (2014). Oracle: `testw90_example21_As_sp`.
- **Stengel–Spaldin functional**: `use_ss_functional` (v4-dev only,
  `wannier90_readwrite.F90:700`, hooks in `wannierise.F90:334,449`).
  Oracle: `testw90_knbo3_higher_stengel_spaldin`.
- **Preconditioned CG**: `precond` keyword; oracles `testw90_precond_1/2`. We accept the
  keyword (`src/known_keywords.jl:168`) but do not implement it.
- **dis_spheres**: spherical disentanglement regions (`dis_spheres_num`,
  `dis_spheres_first_wann`). Oracle: `testw90_lavo3_dissphere`. Keyword-recognized only.
- **Projectability-disentangled WFs (PDWF)**: `dis_proj_min`/`dis_proj_max` (v4-dev;
  Qiao, Pizzi, Marzari npj Comput Mater 9, 208 (2023)) — window selection from .amn
  projectabilities instead of energies. Oracle: `testw90_graphene_pdwf`.
  Keyword-recognized only.

### 1.5 wannier90.x: output files not yet written

All recognized in `src/known_keywords.jl` but unimplemented in our `src/`:

- **Cube format** (`wannier_plot_format = cube`): `plot.F90:2075-2181`
  (`internal_cube_format`), filename `format(a,'_',i5.5,'.cube')`. Oracle:
  `testw90_cube_format`. Explicitly erroring in our `src/plot.jl:185`.
- **`seedname_r.dat`** (`write_rmn`): position matrix elements ⟨0n|r|Rm⟩ from the .mmn
  overlaps, `plot_write_rmn` (`plot.F90:301-313`). Oracle: `testw90_rmn`.
- **`write_hr_diag`**: on-site ⟨0n|H|0n⟩ dump (`plot.F90:250`,
  `wannier90_readwrite.F90:973`).
- **`write_xyz` + `translate_home_cell` + `translation_centre_frac`**: WF-centre .xyz
  with home-cell translation (`plot.F90:277`, `wannier90_readwrite.F90:1372-1381`).
- **Fermi surface plot** (`fermi_surface_plot` → `.bxsf` for XCrySDen):
  `plot_fermi_surface` (`plot.F90:1360-1544`), keywords `fermi_surface_num_points`,
  `fermi_surface_plot_format`.

### 1.6 wannier90.x: bigger subsystems

- **Symmetry-adapted WFs (sitesym)**: `sitesym.F90` (all of it: symmetrize U, gradient,
  Z-matrix, dis_extract with symmetry constraints), reads `seedname.dmn` from
  pw2wannier90. Sakuma PRB 87, 235109 (2013). Oracle: `testw90_disentanglement_sawfs`.
  Touches both disentangle and wannierise inner loops — the most invasive w90 gap.
- **Ballistic transport**: `transport.F90` (~2300 lines): `tran_bulk` (Landauer
  quantum conductance of a periodic chain), `tran_lcr` (lead-conductor-lead with 2c2
  geometry auto-sorting via integral signatures, `tran_find_integral_signatures`),
  surface Green functions (`tran_transfer`/`tran_green`), `*_htX.dat` file I/O.
  No oracle tests in the shipped test-suite.
- **Library mode**: v4 `library_interface.F90`/`c_interface.F90`/`wannier90.h` — the
  new setter-based C/Fortran API (oracles `checkpoint0_write`/`checkpoint1_read`).
- **Utilities**: `w90chk2chk.F90` (formatted↔unformatted .chk conversion for
  cross-platform transport), `w90spn2spn.F90` (same for .spn).

---

## 2. WannierBerri v1.7.0 — capabilities beyond our scope

From the docs calculator list (docs.wannier-berri.org) and README:

- **Methodology (the real differentiators):** symmetry reduction to irreducible-BZ
  k-points **plus exact tensor symmetrization of results**; recursive adaptive grid
  refinement for *all* quantities (we have w90-parity adaptive refinement for AHC only);
  Fermi-level scan at no extra cost; minimal-distance replica (= our use_ws_distance);
  FFT-based interpolation; object-oriented calculator framework; ray-based parallelism.
- **Static calculators we lack:** Ohmic conductivity (FermiSea/FermiSurf — ≈ our
  BoltzWann τ-const special case), classical Hall (σ:S/m/T), Berry-curvature dipole /
  NLAHC (≈ postw90 gyrotropic D0/Dw — **covered by our in-progress gyrotropic module**),
  nonlinear Drude, Zeeman corrections to AHC (orb/spin), eMChA, GME orb/spin
  (≈ gyrotropic K — in progress), cumDOS.
- **Dynamic calculators we lack:** **shift current** (= §1.1), **injection current**
  (circular photogalvanic; not in postw90 at all), dynamic SHC (= our Kubo SHC freq scan
  — have). SHG is *not* in the v1.7 mainline calculator list (shift + injection only).
- **Tabulation:** k-resolved Energy/Velocity/InvMass/BerryCurvature(+der,der2)/
  OrbitalMoment/Spin/SpinBerry on 3D grids, FermiSurfer-compatible output. We cover the
  2D slice case (kslice) and arbitrary-point case (geninterp) but not 3D-grid tabulation
  with a viewer-ready format.
- **Inputs:** wannier90 chk/mmn *and* directly from `_tb.dat`-style TB models, PythTB,
  TBmodels, FPLO, **k·p models**, **phonons via phonopy**. We currently require
  chk(+mmn/eig); we *write* hr/tb but do not *read* them as an interpolation source.
- **Recent:** symmetry-adapted WFs + projection search from symmetry indicators
  (overlaps §1.6 sitesym).

## 3. Wannier.jl v0.3.6 — capabilities beyond our scope

- **Parallel-transport gauge** (Gontier–Levitt–Siraj-Dine): Wannierization of isolated
  manifolds with no initial projections.
- **Automated initial projections** and **valence/conduction manifold splitting**
  (split a converged valence+conduction gauge into two separately-Wannierized groups).
- Constrained WF centres (≈ SLWF+C, §1.4).
- Real-space WF evaluation with xsf **and cube** output (cube = §1.5 gap).
- **DFTK.jl integration** (`DFTKWannierExt` in DFTK): in-memory amn/mmn/eig from a Julia
  DFT code, no file round-trip. Research-grade but strategically important — it makes an
  all-Julia DFT→Wannier pipeline possible today, and we are the more complete backend.

---

## 4. Triage table

Value = benefit to a 2026 user choosing our package. Effort = in our codebase
(S ≲ 1 day, M ≈ 2–5 days, L ≳ 1–2 weeks), given existing infrastructure.

| Feature | Who has it | Value | Effort | Verdict |
|---|---|---|---|---|
| Shift current (`berry_task=sc`) | postw90, WannierBerri | high | M | **do-now** — top physics gap, 5 oracle tests, reuses kubo freq machinery |
| Cube-format WF output | w90, Wannier.jl | high | S | **do-now** — most-requested viz format; grid data already computed for xsf; oracle exists |
| `seedname_r.dat` (`write_rmn`) | w90 | med | S | **do-now** — cheap, oracle `testw90_rmn`, feeds external optics/TB tools |
| Guiding centres | w90 | med-high | S | **do-now** — robustness of wannierisation for poor projections; oracle exists |
| `.bxsf` Fermi-surface plot | w90 | med | S | **do-now** — trivial on top of our interpolation; XCrySDen/FermiSurfer users expect it |
| `write_hr_diag` / `write_xyz` / `translate_home_cell` | w90 | low | S | **do-now** — bundle in one "output parity" sweep; closes the full-core claim |
| PDWF (`dis_proj_min/max`) | w90 v4-dev, (Wannier.jl lineage) | med | S | **do-now** — becoming the ecosystem-recommended default; tiny change to window selection; oracle exists |
| TB-model input (read `_hr.dat`/`_tb.dat` as interpolation source) | WannierBerri (System_tb, PythTB, TBmodels) | high | S | **do-now** — we already have the parsers' inverse; unlocks the whole postw90 stack for model Hamiltonians |
| DFTK.jl bridge (in-memory amn/mmn/eig, package extension) | Wannier.jl | high (strategic) | M | **do-now** — the Julia-native differentiator; we become the complete backend for an all-Julia pipeline |
| Injection current | WannierBerri | med-high | M | **do-later** — natural follow-up sharing shift-current machinery; no Fortran oracle (validate vs WannierBerri) |
| Irreducible-BZ sampling + tensor symmetrization of results | WannierBerri | high | L | **do-later** — flagship perf feature (10–50× on symmetric crystals); needs spglib (Spglib.jl) + R-matrix symmetrization; design carefully, gauge-invariant validation only |
| Adaptive recursive refinement for all quantities | WannierBerri | med | M | **do-later** — generalize our existing AHC refinement kernel; pairs with the symmetry work |
| Tetrahedron method for SHC | postw90 | med | M | **do-later** — 2 oracles; niche until users hit smearing-convergence pain; self-contained module |
| SLWF+C (selective localization, constrained centres) | w90, (Wannier.jl constrained centres) | med | M | **do-later** — real user base (embedding, defect WFs); localized change to spread/gradient; oracle exists |
| `berry_task=kdotp` | postw90 | low-med | S | **do-later** — cheap and oracle-tested but niche; nice pairing with a future k·p model reader |
| 3D tabulation grids + FermiSurfer `.frmsf` output | WannierBerri | med | S | **do-later** — small extension of geninterp/kslice machinery |
| `dis_spheres` | w90 | low-med | S | **do-later** — oracle exists; implement when a user needs k-localized disentanglement |
| Preconditioned CG (`precond`) | w90 | med | M | **do-later** — matters for large systems; our :rcg converges well on test set; oracles exist |
| Parallel-transport gauge | Wannier.jl | med | M | **do-later** — projection-free isolated-manifold wannierisation; research value, small user base |
| Automated projections + valence/conduction splitting | Wannier.jl | med | M | **do-later** — PDWF (do-now) covers most of the practical demand |
| Symmetry-adapted WFs (sitesym, `.dmn`) | w90, WannierBerri (new) | med-high | L | **do-later** — invasive (constraints inside both optimization loops); do after the WannierBerri-style symmetrization, which shares group-theory code |
| Formatted chk/spn conversion (w90chk2chk/w90spn2spn) | w90 | low | S | **do-later** — fold into io.jl when cross-platform fixtures are needed |
| Ballistic transport (`tran_bulk`/`tran_lcr`) | w90 | low | L | **skip** — community moved to kwant/TranSIESTA/NEGF codes; no oracle tests shipped; 2c2 geometry is high-maintenance |
| Library mode (v4 C API) | w90 | low | — | **skip** — a Julia package *is* the library; document a w90-v4-API → Julia-API mapping instead |
| Stengel–Spaldin functional (`use_ss_functional`) | w90 v4-dev | low | M | **skip** — undocumented dev feature; revisit if it lands in a w90 release |
| Fermi-sea vs fermi-surface formula variants | WannierBerri | low-med | M | **skip** — cross-validation nicety; our formulas are already oracle-anchored to the Fortran |
| External/internal (Wannier vs beyond) term split | WannierBerri | low | S | **skip** — diagnostic output, not physics users ask for |
| Nonlinear Drude, classical Hall, Zeeman/eMChA corrections | WannierBerri | low-med | M each | **skip** — niche; revisit on user demand after injection current |
| Second-harmonic generation | nobody mainline (not in WB v1.7 calculators) | med | L | **skip** — no reference implementation to validate against; reconsider if WB mainlines it |
| Phonon support (phonopy interface) | WannierBerri | med | L | **skip for now** — different physics domain; a Julia route would go via Phonopy.jl and deserves its own design |
| k·p model input | WannierBerri | low | S | **skip** — pairs with kdotp task if ever needed |
| Real-space WF operator evaluation | Wannier.jl | low | S | **skip** — UNK/xsf covers the use case |

---

## 5. Recommended next slate

Ordered; items 1–3 finish the current milestone honestly, 4–6 open the next one.

1. **wannier90.x output-parity sweep (all S, all oracle-tested):** cube format,
   `_r.dat`, `.bxsf`, `write_hr_diag`, `write_xyz`/`translate_home_cell`, guiding
   centres, PDWF, `dis_spheres`. Roughly a week total; after it, "full wannier90.x"
   is true without asterisks, and every item lands with a Fortran oracle test.
2. **Shift current (`berry_task=sc`).** The single highest-value physics addition:
   photogalvanics is an active field, five oracle tests exist covering both phase
   conventions + ws_distance + eta correction, and it completes the postw90 optics
   stack alongside our existing Kubo module. Extract a berry-ahc.md-style spec note
   first (generalized derivative + `sc_phase_conv` traps).
3. **TB-model input** (read `_hr.dat`/`_tb.dat` into the interpolation object). One
   small constructor that turns every postw90-side feature into a model-Hamiltonian
   tool and matches the entry path WannierBerri/PythTB users expect.
4. **DFTK.jl bridge** (package extension producing amn/mmn/eig in memory). Strategic:
   the all-Julia DFT→MLWF→transport pipeline exists in no other ecosystem with our
   level of validation.
5. **Symmetry: irreducible-BZ sampling + result symmetrization** (L). The flagship
   performance/quality feature that keeps WannierBerri ahead today; shares group-theory
   infrastructure with a later sitesym/SAWF implementation. Validate gauge-invariantly
   (bands/Ω/centres — never U(k) or H(R) elements).
6. **Injection current + tetrahedron SHC + SLWF+C** as the following physics batch,
   in whatever order user demand suggests.

Explicitly deprioritized despite being "features on a list": ballistic transport,
library mode, SHG, phonons — see table for reasons.
