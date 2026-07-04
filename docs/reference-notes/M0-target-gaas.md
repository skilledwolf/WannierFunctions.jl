# M0/M1 validation oracle — GaAs example01 (isolated 4 valence bands)

Reference: `reference/wannier90/test-suite/tests/testw90_example01/`
Inputs shipped: `gaas.win`, `gaas.amn`, `gaas.mmn` (no `.eig` — pure localization, no disentanglement).
Reference output: `benchmark.out.default.inp=gaas.win`.

## Case parameters
- `num_wann = 4`, `num_bands = 4` (isolated valence bands; num_bands defaults to num_wann).
- `mp_grid = 2 2 2` → 8 k-points. `num_iter = 20`. `search_shells = 12`. `use_ws_distance = .false.`
- Unit cell in **bohr**; projections `As:sp3`.

## kmesh reference (validate first, in isolation)
- Nearest-neighbour shell 1: distance 0.957961 Å⁻¹, multiplicity 8.
- b-vectors chosen: **shell 1 only**, 8 nearest neighbours. Completeness (B1) fully satisfied.
- Weight per b-vector: **w_b = 0.408635 Ų** (all 8 equal).
- Example b_k (Å⁻¹): (-0.553079, 0.553079, -0.553079), (0.553079, 0.553079, 0.553079), …

## M0 — initial state (from reference `.amn`, before any iteration)
Each of 4 WFs: spread **1.11720303 Ų**. Centres (Å):
- WF1 (-0.866632, 1.973462, 1.973462)
- WF2 (-0.866632, 0.866632, 0.866632)
- WF3 (-1.973462, 1.973462, 0.866632)
- WF4 (-1.973462, 0.866632, 1.973462)
- Sum of centres+spreads: (-5.680188, 5.680188, 5.680188), total **4.46881212 Ų**
- Iteration 0 decomposition: **Ω_D = 0.0083198, Ω_OD = 0.5036294, Ω_TOT = 4.4688121**
  (⇒ Ω_I = 4.4688121 − 0.0083198 − 0.5036294 = **3.9568629**, invariant.)

## M1 — final state (after 20 MV iterations)
- Ω_I = **3.956862958** (unchanged — invariant, good check)
- Ω_D = **0.008030049**
- Ω_OD = **0.501987969**
- Ω_Total = **4.466880976**
- Each WF spread ≈ 1.11672024; centres barely move (near-optimal projection start).

## File formats confirmed from these files
- `.amn`: line1 comment; line2 `num_bands num_kpts num_wann`; then num_bands·num_wann·num_kpts lines
  `m n k re im` with **m (band) fastest, then n (wann), then k**. A[m,n,k].
- `.mmn`: line1 comment; line2 `num_bands num_kpts nntot`; then nntot·num_kpts blocks, each:
  one line `k  kb  g1 g2 g3` (kb = neighbour k index; g = reciprocal G shift), then num_bands²
  lines `re im` for M[m,n] with **m fastest**. M[m,n,b,k].
- `.eig`: lines `band_index kpt_index energy_eV` (needed for disentanglement/interpolation, not M0/M1).
