# Closing the remaining gaps (2026-07-05)

The slate that closed every item still open in `parity-audit-2026.md` after the ecosystem
slate. Conventions and oracle anchors for each feature.

## 1. Symmetry-adapted disentanglement (`dis_extract_symmetry`)

- Port of `sitesym_symmetrize_zmatrix` + `sitesym_dis_extract_symmetry` (sitesym.F90:563) into
  the Z-matrix iteration ([disentangle.jl](../../src/disentangle.jl), [sitesym.jl](../../src/sitesym.jl)).
- Constrained Ω_I step at each irreducible representative: `λ = U†ZU`, `ΔU = ZU − U·λ`,
  band-by-band 2×2 generalized eigenproblem in span{u_i, Δu_i} (larger eigenvalue), then
  `symmetrize_ukirr` re-projection each sweep (≤ 50). Convergence `Σ|ΔU| < 1e-10`.
- Z symmetrised as star-sum (each distinct member once) + little-group average, ÷|G_k|.
- Bookkeeping: representative k carries star weight `nsym/|G_k|`; `Ω_I(i−1)` uses `−tr Re λ`.
  Frozen windows error (parity with the reference). Post-loop: NO rotation of U_opt to the
  subspace-H eigenbasis (only eigval_opt updated); square gauge built at representatives and
  star-propagated after `replace_d_matrix_band` (d_band := d_wann).
- **Fixed along the way**: `symmetrize_gradient!` mode 2 double-counted the identity
  (acc seeded with G *and* the loop included isym = 1). This let the H3S localisation drift
  off the symmetric manifold (8.8e-3); now 7.7e-12 over 500 iterations.
- Oracle: `testw90_disentanglement_sawfs` (H3S, 20→12, mix 0.2, 10 iterations). The full
  10-iteration Ω_I trajectory matches the benchmark to every printed digit
  (3.61187069/3.46062472 → 3.44624098/3.40892357); final Ω_I = 3.408923571 exact. Converged
  localisation Ω = 6.301957278 vs a self-generated converged wannier90.x run (num_iter=5000)
  at 6.301957261.
- `symmetrize_eps` keyword now parsed (Sitesym.eps field; H3S uses 1e-8, default 1e-3).

## 2. Higher-order finite differences (`higher_order_n`, Lihm)

- b-vector list carries every first-order b and its multiples 2b…Nb with 1D central-difference
  factors `w_mb = w_b (1/m²) Π_{j≠m} j²/(j²−m²)` (N=2: 4/3, −1/12); B1 preserved via
  `Σ_m m² fact_m = 1`. kmesh.F90's `kmesh_shell_reconstruct`.
- Read path ([kmesh.jl](../../src/kmesh.jl)): ordering-INDEPENDENT detection — each b is assigned the integer
  multiple of the shortest parallel list member; B1 least-squares runs on the multiple-1 set
  only (the min-norm solution over all shells is WRONG for parallel shells: it splits 1:4).
  Needed because shipped .mmn orderings vary (knbo3's is (sign, order, shell), not blocked).
- `-pp` path ([nnkp.jl](../../src/nnkp.jl)): first-order search as before, multiples appended in the canonical
  [block1, 2·block1, …] order, neighbours located by exact folded match (kmesh.F90:491).
  nnkpts block byte-identical to `wannier90.x -pp` on knbo3.
- postw90 side inherits automatically through `build_bvectors` (BerryModel et al.).
- `.uHu` reader now auto-detects Fortran-unformatted files (record markers contain NULs).
- Projections spec: comma is accepted as an orbital separator (`O:s,p` ≡ `O:s;p`).
- Oracles: `testw90_knbo3_higher` (Ω_I/D/OD/total match to 9 digits),
  `testpostw90_fe_morb_transl_inv_higher` (M_z −0.0617 exact),
  `testpostw90_pt_shc_higher` (fermi scan: −3.9190876 exact, 1633.0502 / 435.58483 at 1e-7).

## 3. Stengel–Spaldin functional (`use_ss_functional`)

- Single-point objective on the k-averaged diagonal overlap in uniform b-order:
  `Ω_SS = Σ_n Σ_b w_b (1 − |M̄_nn(b)|²)`, `M̄ = (1/N_k) Σ_k M_nn(k, b_ord(k))`; the reference's
  om_i + om_od + om_d(SS k-variance) split collapses exactly to this. 4-term gradient over the
  ±b pairs (`nnord`/`nnrev` maps). [ss.jl](../../src/ss.jl); runs on the :w90 optimiser; final report is the
  standard MV decomposition at the SS gauge.
- **The SS surface is a near-flat valley**: the benchmark's "converged" 13.845371018 is where
  the default Δ-criterion stops — wannier90.x itself reaches 13.312 with conv_tol = 1e-16, and
  our runs reach 13.794. Validation is therefore state-function parity, not stopping-point
  parity: objective at the shared initial gauge = w90's iter-0 spread to 10 digits
  (14.7994988992); objective at w90's converged gauge = its om_tot to 9 digits; our gradient
  norm there 3.5e-6 (their stationary point is ours); FD slope check 2e-6.
- Guiding centres are ignored in SS mode (csheet cancels in both objective and gradient).

## 4. Γ-only real-orthogonal parity (`wann_main_gamma`)

- Jacobi-sweep joint diagonalisation on the √w_b-weighted real/imag parts of the half b-set
  ([gamma.jl](../../src/gamma.jl)): per-pair closed-form angle (a11/a12/a22 stencil), M and U updated by real
  2×2 rotations; spread via `wann_omega_gamma` (atan2 centres; om_d ≡ 0 for the orthorhombic
  3-b case by the B1 identity).
- Γ-disentanglement: the complex Z-iteration runs unchanged; the optimal subspace of real Γ
  data is conjugation-closed, so the embedding is rotated to a real orthonormal basis
  (`realify_subspace`, SVD of [Re U | Im U]; guarded to 1e-8) before a real Löwdin handoff.
- Guiding centres: only affect the centre log-branch bookkeeping. Ignored for 3 half-b's
  (reference auto-disables); for larger half-sets the guided-spread reporting
  (`compute_spread(...; guides)`) unwraps the branches — essential for chains, where the
  principal branch puts Ω_D at +388 (na_chain) although the gauge itself is right.
- All four gauges come out exactly real (`maxIm(U) = 0`).
- Oracles (Ω_I and Ω to all printed digits): benzene_gamma_val 10.455472666/12.958338012,
  valcond (90→18 with frozen window!) 28.959463251/31.468492322, hexcell
  9.743320252/12.091929665(bench …660), na_chain (30→10) 36.511339464/37.505387845.

## 5. Ballistic transport (`transport_mode = bulk`)

- [transport.jl](../../src/transport.jl): López-Sancho decimation (`tran_transfer`, η = 5e-4i, ≤50 iters, 1e-7),
  surface/bulk Green functions (`tran_green`), Fisher–Lee `T(E) = Tr[Γ_L G Γ_R G†]` + GF DOS
  (`transport_bulk`), principal-layer assembly from H(R) (`transport_from_tb` =
  tran_reduce_hr + tran_cut_hr_one_dim + tran_get_ht, with home-cell-translated WF centres
  for the distance cutoff), `_htB.dat`/`_qc.dat`/`_dos.dat` read/write, `.win`-driven
  `run_transport` wired into the CLI.
- Defaults match the reference: hr_cutoff 0, dist_cutoff 1000 (clamped to L/2), mode three_dim.
- Oracle: self-generated `wannier90.x` run on tellurium (helical chains ∥ c, 12→9,
  transport along z; committed under test/data/te_tran). H00/H01 match the reference
  `_htB.dat` at file precision (3e-6), T(E) to 9e-7 over 121 energies, DOS to 3e-4 at
  GF-singular peaks. NOTE: upstream v4-dev `tran_read_ht` input path segfaults (no shipped
  test exercises it); reported oracle is the full pipeline. `tran_lcr` (2c2 auto-sort)
  remains out of scope.

## 6. Symmetrised-wedge morb + DOS

- `orbital_magnetisation_sym` (pseudovector rule, like AHC) and `density_of_states_sym`
  (scalar star weights) in [symmetry.jl](../../src/symmetry.jl). Fe (magnetic subgroup of order 8 filtered from
  Oₕ): morb matches full-BZ to 1.6e-6 with 78/512 k-points; DOS matches exactly.

## 7. DFTK live end-to-end

- [ext/WannierFunctionsDFTKExt.jl](../../ext/WannierFunctionsDFTKExt.jl) rewritten against the real DFTK API:
  `wannier_model(scfres, projections; num_wann)` uses OUR kmesh search for the b-list and
  DFTK's `overlap_Mmn_k_kpb` / `compute_amn_kpoint` for the matrix elements (requires
  `symmetries = false`).
- Live silicon LDA (Ecut 14, 4×4×4): Ω = 6.4566 Å², WFs on the four bond centres,
  interpolated bands reproduce the SCF eigenvalues on the mesh to 2.5e-12 eV.
  Example: [examples/06_dftk_end_to_end.jl](../../examples/06_dftk_end_to_end.jl); conditional testset test/dftk_e2e.jl.

## Injection-current anchors recalibrated

The WannierBerri cross-check (regenerated on the current reference data) agrees at the
~1e-4 level, limited by the two codes' degenerate-state regularisation on GaAs's exact band
degeneracies — tolerances set to rtol 3e-4 accordingly (largest component agrees at 5.8e-5).

## 8. postw90.jl drop-in binary (added after the slate)

- [bin/postw90.jl](../../bin/postw90.jl) → `postw90_main(seedname)` in [src/postw90.jl](../../src/postw90.jl): keyword → kwargs
  mapping for every module (berry ahc/morb/kubo/sc/shc/kdotp, gyrotropic, dos, kpath, kslice,
  geninterp, boltzwann, spin_moment), reference-named output files.
- New writers: `-kubo_S_ab/-kubo_A_ab/-jdos` (3E16.8/2E16.8), `-ahc/-morb-fermiscan`
  (4(F12.6,1x)), the BoltzWann set `_tdf/_elcond/_sigmas/_seebeck/_kappa/_boltzdos`
  (G18.10 rows, exact headers; Seebeck is the full 3×3 in ROW-major order), geninterp rows
  in `(I10,4G18.10)`.
- `fortran_g(x, w, d)`: Gw.d emulation — F(w−4).(d−n) + 4 blanks in the F-range, Ew.d
  otherwise (needed for byte parity of geninterp/BoltzWann).
- Keyword conventions discovered: the reference matches `berry_task`/`gyrotropic_task` by
  SUBSTRING (`eval_shc` → shc; glued `-C-dos-D0` works); task lists split on `+`; the global
  interpolation mesh keyword is `kmesh` (module-specific `*_kmesh` override); the BoltzWann
  TDF energy window is the DISENTANGLEMENT window ±0.2 eV; `kubo_eigval_max` defaults to
  dis_froz_max + 2/3 (else max eig + 2/3); boltz_dos energy range defaults to eig bounds
  ±0.6667.
- Validation: byte-identical vs local postw90.x runs on Fe kubo (4 files), Si geninterp,
  Si BoltzWann (elcond/sigmas/seebeck/kappa/boltzdos identical; tdf ≤1e-4-relative at 4
  band-crossing energies — degenerate-velocity convention), Cu dos, Fe kpath bands+curv,
  Fe kslice (3 files), Pt shc-fermiscan, Te gyrotropic (6 files), GaAs kdotp, Fe
  ahc-fermiscan. Oracle pairs committed under test/data/postw90.
- [bin/w90chk2chk.jl](../../bin/w90chk2chk.jl): `-export`/`-import` checkpoint conversion (data-exact; header
  date restamped).
