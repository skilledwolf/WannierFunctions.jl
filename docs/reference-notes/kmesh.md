# kmesh / b-vector finite-difference machinery — reference notes

Source: `reference/wannier90/src/kmesh.F90` (Wannier90 v3.1.0-greenfield tree).
Supporting: `src/utility.F90`, `src/constants.F90`, `src/types.F90`, `src/readwrite.F90`.

Goal of this module: find the neighbour-shell b-vectors `{b}` connecting each k-point to its
periodic-image neighbours, and their weights `w_b`, satisfying the completeness (B1) relation of
Marzari & Vanderbilt, PRB 56, 12847 (1997), Appendix B:

    Sum_b  w_b  b_alpha  b_beta  =  delta_{alpha beta}          (Eq. B1)

All formulas below are quoted with `file:line`. Everything here is for the **default path**
(`higher_order_n = 1`, `higher_order_nearest_shells = .false.`, `kmesh_shell_from_file = .false.`,
`skip_B1_tests = .false.`, `gamma_only = .false.`). Higher-order and gamma_only variations are noted
where they change numbers.

---

## 0. CRITICAL conventions / gotchas (read first)

These are the things most likely to cause a *silent* numerical mismatch in a Julia port.

1. **Reciprocal lattice carries a factor of 2π and is stored ROW-WISE.**
   `utility_recip_lattice_base` (utility.F90:346-351): `recip_lat = twopi * adjoint(real_lat)/volume`,
   `twopi = 2*pi` (constants.F90:56). So b-vectors are in units where BZ vectors include 2π
   (Cartesian, length units are inverse-length: Ang^-1 or Bohr^-1).
   Storage is **row-major lattice vectors**: `matmul(lmn(:,loop), recip_lattice)` computes
   `Sum_i lmn(i) * recip_lattice(i, :)`, i.e. reciprocal lattice vector `i` is
   `recip_lattice(i, 1:3)` (a *row*). Confirmed by `utility_frac_to_cart` (utility.F90:459-461):
   `cart(i) = Sum_j real_lat(j,i)*frac(j)` — real lattice vector `j` is the *row* `real_lat(j, :)`.
   In a Julia column-major port you must transpose relative to the naive reading, or the whole
   geometry is wrong.

2. **THREE DISTINCT tolerances — do NOT merge them** (the single "kmesh_tol for the SVD" framing is
   wrong):
   - `kmesh_input%tol` (user keyword `kmesh_tol`, **default 1.0e-6**, types.F90:151, readwrite.F90:936):
     used for (a) shell-distance / shell-membership comparisons (kmesh.F90:183-188, 435, 1301) and
     (b) the **B1 acceptance test** (kmesh.F90:1839,1841 in `kmesh_check_condition`; 1987,1990 in
     `kmesh_shell_fixed`; 644,649 in the higher-order per-kpoint block).
   - **SVD small-singular-value reject** — value **differs by code path**:
     `eps5 = 1e-5` in `kmesh_shell_automatic` (kmesh.F90:1526);
     `eps7 = 1e-7` in `kmesh_shell_fixed` (kmesh.F90:1956) and `kmesh_shell_from_file` (2212).
   - `eps6 = 1e-6`: parallel-shell (cosine≈1) rejection in the automatic search (kmesh.F90:1422).
   - `eps8 = 1e-8`: vector equality in `utility_compar` (utility.F90:407) and degeneracy tie in
     `internal_maxloc` (kmesh.F90:2269).
   (constants: eps5..eps8 at constants.F90:70-76.)

3. **B1 is NOT re-checked per-kpoint in the default path.** The per-kpoint B1 block
   (kmesh.F90:611-657) is guarded by `(.not. skip_B1_tests) .and. higher_order_nearest_shells`.
   With defaults (`higher_order_nearest_shells=.false.`) the **entire block is skipped**. B1 is
   enforced *only during weight determination for k-point 1* inside `kmesh_shell_automatic` /
   `kmesh_shell_fixed`. The per-kpoint loop that *does* run (kmesh.F90:586-605) only checks that
   neighbour b-vector *lengths* match between k-point 1 and k-point nkp (non-symmetric neighbour
   guard), NOT B1.

4. **Reproducible ordering.** The final order of b-vectors (and hence `.nnkp` line order and the
   `wb`/`bk` index conventions everything downstream relies on) is:
   `supercell_sort` cell order  ×  shell order  ×  (inner `nkp2`-then-cell) discovery order.
   `internal_maxloc` (kmesh.F90:2243-2277) is a *deterministic* maxloc that returns the **lowest
   index on ties** (line 2275) precisely so the cell ordering is reproducible. A Julia port that
   sorts cells differently will still satisfy B1 but produce a *different* b-vector ordering →
   mismatched `.nnkp`, `M_mn(k,b)` ordering, etc.

5. **gamma_only** keeps half the b-vectors and **doubles** their weights (`wb(na)=2*wb_local(nn)`,
   kmesh.F90:848).

6. **1-based indexing** throughout. Arrays: `bk(3, nntot, num_kpts)`, `wb(nntot)`,
   `bka(3, nnh)`, `nnlist(num_kpts, nntot)`, `nncell(3, num_kpts, nntot)`,
   `neigh(num_kpts, nnh)`. Note the **transposed** k/nn index order in `nnlist`/`nncell`/`neigh`
   (k-point is the first index) vs `bk` (k-point is the *last* index).

7. **Length units.** All internal geometry is Cartesian in whatever the lattice was given
   (`.win` `bohr`/`ang`); printed values are divided by `print_output%lenconfac` and weights
   multiplied by `lenconfac**2` for output only (kmesh.F90:715). `bohr = 0.5291772108` Å
   (constants.F90:219). The stored `wb`, `bk`, `bka` are in the internal units, NOT the printed ones.

---

## 1. Entry point and driver: `kmesh_get` (kmesh.F90:74-955)

Inputs: `kpt_latt(3,num_kpts)` (fractional k-points), `real_lattice(3,3)`, `num_kpts`, flags.
Outputs into `kmesh_info`: `nntot, nnh, nnlist, neigh, nncell, wb, wbtot, bk, bka` (+ `nnord/nninv/nnrev`).

Steps:

1. `utility_recip_lattice(real_lattice, recip_lattice, volume)` (kmesh.F90:151) — 2π, row-wise (§0.1).
   Also `utility_inverse_mat(recip_lattice, inv_lattice)` (152) — only used by higher-order path.
2. `kmesh_supercell_sort` (158) — build & distance-sort the search cells `lmn` (see §2).
3. Convert every k-point to Cartesian: `utility_frac_to_cart(kpt_latt(:,nkp), kpt_cart(:,nkp), recip_lattice)`
   (165-167).
4. **Nearest-neighbour shell distances** (172-199): loop over `search_shells`, for each find the next
   larger distance `dnn1` (from k-point 1 to all `kpt_cart(:,nkp2)+G_lmn`), counting multiplicity.
   `dnn(nlist)` = shell distance, `multi(nlist)` = shell multiplicity, `ndnntot` = number of shells
   actually found. `eta = 99999999.0` sentinel (108).
   Distance test (183-188):
   ```
   if (dist > tol .and. dist > dnn0 + tol) then
     if (dist < dnn1 - tol) then dnn1 = dist; counter = 0
     if (dist > dnn1 - tol .and. dist < dnn1 + tol) counter = counter + 1
   ```
5. **Shell selection & weights** (256-273): dispatch to exactly one of
   - `kmesh_shell_from_file` if `kmesh_shell_from_file` (§6),
   - `kmesh_shell_automatic` if `num_shells == 0` (default) (§4),
   - `kmesh_shell_fixed` if `num_shells > 0` (user `shell_list`) (§5).
   Each returns `bweight(1:num_shells)` (per-shell weights) and sets `num_shells` + `shell_list`.
6. `nntot = Sum over selected shells of multi(shell_list(loop_s))` (300-303).
   For higher_order default (`higher_order_n=1`) this multiply is a no-op (365).
7. Assemble per-b weights into `wb_local(nnx)` (399-405): each of the `multi` b-vectors in shell
   `loop_s` gets weight `bweight(loop_s)`.
8. **Build neighbour lists** `nnlist`, `nncell`, and local b-vectors `bk_local` for *every* k-point
   (423-450) — see §3.
9. Symmetry length check per k-point (586-605) (NOT B1 — see §0.3).
10. Per-kpoint B1 block (611-657) — **skipped by default** (§0.3).
11. `wbtot = Sum_nnx wb_local(nnx)` over shell-1 neighbours (668-676).
12. Build `bka` (half-set of b-directions, inversion removed) (678-701), `neigh` index array
    (735-751), copy locals into `kmesh_info%wb`, `kmesh_info%bk` (753-763).
13. gamma_only halving (767-901); permutation arrays `nnord/nninv/nnrev` for non-gamma (905-925).

---

## 2. Supercell search order: `kmesh_supercell_sort` (kmesh.F90:1199-1253)

- Module parameter `nsupcell = 5` (kmesh.F90:68). Builds all cells `lmn(1:3, .)` with
  `l,m,n ∈ [-5,5]` → `(2*5+1)^3 = 1331` cells (1229-1239). Origin `(0,0,0)` placed first.
- `dist(counter) = |matmul(lmn, recip_lattice)|` (Cartesian distance of each cell image from origin).
- Sort **descending into the tail, so final order is ascending distance** using the reproducible
  `internal_maxloc` (1241-1246): repeatedly pull the current max to the back, marking used entries
  `-1`. Result: `lmn` reordered nearest-to-origin first.
- **GOTCHA:** the sort routine uses the module constant `nsupcell = 5`, but the *neighbour search
  loops* in `kmesh_get`/`kmesh_get_bvectors` iterate `(2*search_supcell_size+1)**3` with the
  **input** `search_supcell_size` (default 5, types.F90:150). They coincide at defaults but are
  different variables — a Julia port should size the `lmn` array from `nsupcell`=5 and index the
  loops from `search_supcell_size`.

`internal_maxloc` (2243-2277): `maxloc` then, among values within `eps8` of the max, return the
**minimum index** — deterministic tie-break.

---

## 3. b-vector collection

### `kmesh_get_bvectors` (kmesh.F90:1256-1319)
For a given shell distance `shell_dist` and origin k-point `kpt`, returns the `multi` b-vectors:
```
vkpp = matmul(lmn(:,loop), recip_lattice) + kpt_cart(:,nkp2)          (1296-1298)
dist = |kpt_cart(:,kpt) - vkpp|                                        (1299-1300)
if (dist >= shell_dist - tol .and. dist <= shell_dist + tol):
    num_bvec += 1
    bvector(:, num_bvec) = vkpp(:) - kpt_cart(:, kpt)                  (1303)
```
So **`b = (k_neighbour_image) - k_origin`** (points from origin k to neighbour). Sign matters for
downstream `M_mn(k, b)`. Early-exit once `num_bvec == multi` (1306).

### Main neighbour build in `kmesh_get` (kmesh.F90:423-450)
For every k-point `nkp`, over selected shells `ndnn = shell_list(ndnnx)`, over sorted cells, over
`nkp2`:
```
vkpp = matmul(lmn(:,loop), recip_lattice) + kpt_cart(:,nkp2)          (430-432)
dist = |kpt_cart(:,nkp) - vkpp|                                       (433-434)
if (dist in [dnn(ndnn)-tol, dnn(ndnn)+tol]):
    nnx += 1
    nnlist(nkp, nnx) = nkp2                                           (438)
    nncell(1:3, nkp, nnx) = (l, m, n)                                 (439-441)
    bk_local(:, nnx, nkp) = vkpp - kpt_cart(:, nkp)                   (442)
exit shell inner loops when nnshell(nkp,ndnn) == multi(ndnn)          (445)
```
`nnlist(nkp,nn)` = index (1..num_kpts) of the home-BZ k-point that is the periodic image of `k+b`.
`nncell(:,nkp,nn) = (l,m,n)` = the reciprocal-lattice G that maps that home k-point to the actual
`k+b`. These two are exactly what is written to `.nnkp` (§7).

---

## 4. Automatic shell selection + weights: `kmesh_shell_automatic` (kmesh.F90:1322-1676)

Strategy (docstring 1326-1331): take next shell; reject if a b-vector is parallel to an existing
shell (|cos|≈1 within `eps6`, 1413-1434); add shell; solve least-squares for weights; test B1; if
not satisfied add another shell; repeat up to `search_shells`.

Default-path structures (`higher_order_n_local = 1`, so `max_shells_aux = 6`):
- **`amat`** is `6 × num_shells`, filled by `kmesh_get_amat` (1506) with rows ordered
  `(xx, yy, zz, xy, yz, zx)`. For loop_order=1, `num_of_eqs = (1+1)(1+2) = 6`; the six
  `(num_x,num_y,num_z)` monomial exponents summed over the shell's b-vectors:
  `amat(row, s) = Sum_b  b_x^num_x * b_y^num_y * b_z^num_z`  (kmesh.F90:1791-1794).
  This reproduces the explicit fixed form (see §5) for order 1.
- **Target** `target(1)=target(3)=target(6)=1`, rest 0 (1377) → RHS `(1,1,1,0,0,0)` in the
  `(xx,yy,zz,xy,yz,zx)` basis, i.e. `Sum_b w_s (b⊗b) = I`.
- **SVD**: LAPACK `dgesvd('A','A', 6, num_shells, amat, 6, singv, umat, 6, vmat, num_shells, work, 60, info)`
  (1511-1513). `dgesvd` returns **Vᵀ** in `vmat`. Small-singular-value reject: `any(|singv| < eps5)`
  (1526): if `num_shells==1` → fatal; else drop this shell (`num_shells -= 1`) and try next (1536-1537).
- **Weights (pseudoinverse)** `w = V Σ⁻¹ Uᵀ target` (1541-1552), unpacked as:
  ```
  smat(s,s) = 1/singv(s)
  tmp0 = transpose(umat);  tmp1 = tmp0·target;  tmp2 = smat·tmp1;  tmp3 = transpose(vmat)·tmp2
  bweight(1:num_shells) = tmp3
  ```
  `transpose(vmat)` = V (since `vmat` holds Vᵀ). So `bweight(loop_s)` is the weight for **all**
  b-vectors of shell `loop_s`.
- **B1 check**: `kmesh_check_condition` (1563-1566 → 1801-1845):
  `delta = Sum_s Sum_b bweight(s) b_x^nx b_y^ny b_z^nz`; for order-1 diagonal terms
  (`loop_i ∈ {1,3,6}`, i.e. xx,yy,zz) require `|delta - 1| ≤ tol`; else `|delta| ≤ tol`
  (1838-1842). If unsatisfied, add another shell; if none of `search_shells` works → fatal.

---

## 5. Fixed (user shell_list) weights: `kmesh_shell_fixed` (kmesh.F90:1848-2005)

Used when `num_shells > 0` (user gave `shell_list`). No higher-order. `amat` is
`max_shells(=6) × num_shells`, filled **explicitly** (1930-1939):
```
amat(1,s)+= b_x*b_x   amat(2,s)+= b_y*b_y   amat(3,s)+= b_z*b_z
amat(4,s)+= b_x*b_y   amat(5,s)+= b_y*b_z   amat(6,s)+= b_z*b_x
```
i.e. rows `(xx, yy, zz, xy, yz, zx)`. `target = (1,1,1,0,0,0)` (1889).
SVD `dgesvd('A','A', max_shells=6, num_shells, ...)` (1942-1943); small-sv reject `< eps7` (1956);
weights `bweight = matmul(transpose(vmat), matmul(smat, matmul(transpose(umat), target)))` (1966).
B1 check (1976-1994): loops `loop_i=1..3`, `loop_j=loop_i..3`,
`delta = Sum_s Sum_b bweight(s) b_{loop_i} b_{loop_j}`; diagonal `|delta-1|>tol` or off-diagonal
`|delta|>tol` → `b1sat=.false.` → fatal (1996-1999). Skipped if `skip_B1_tests`.

`max_shells = 6`, `num_nnmax = 12` (types.F90:128-129).

---

## 6. From-file shells: `kmesh_shell_from_file` (kmesh.F90:2008-2240)

Activated by `kmesh_shell_from_file=T`; reads `<seedname>.kshell`. Each non-comment
(`!`/`#`) line lists the integer indices (into the full ordered b-vector list) forming one shell;
number of lines = `num_shells`, tokens per line = that shell's `multi` (2107-2134).
b-vectors gathered from the full ordered list (built shell-by-shell via `kmesh_get_bvectors`,
2071-2080) into `bvec_inp(3, maxval(multi), num_shells)` (2143-2147).
Same `amat` `(xx,yy,zz,xy,yz,zx)` fill (2186-2195), `target=(1,1,1,0,0,0)` (2052), SVD with
`eps7` reject (2212), same pseudoinverse weight solve (2222). **B1 is NOT tested inside this routine**
(comment 2230); the per-kpoint B1 block back in `kmesh_get` (611+) is where it would be checked —
but note that block requires `higher_order_nearest_shells` (§0.3), so with a from-file + first-order
setup B1 is effectively not re-verified beyond the SVD residual. Neighbour lists for the from-file
path are (re)built in `kmesh_get` at 527-559 by matching `kpt_cart(:,nkp)+bvec_inp` to images.

---

## 7. Assembly of `bka`, `neigh`, and the half-set

- `nnh = nntot/2` (678).
- **`bka(3, nnh)`** (681-701): scan the `nntot` neighbours of k-point 1; add each b-vector unless its
  **inverse** is already in the set. `utility_compar(bka(:,nap), bk_local(:,nn,1), ifpos, ifneg)`
  → `ifneg==1` means `bk ≈ -bka` within `eps8` (utility.F90:404-409). So `bka` is one representative
  per ± direction. Must end with exactly `nnh` entries or fatal (698-701).
- **`neigh(nkp, na)`** (735-751): for each k-point and each `bka` direction `na`, find the neighbour
  index `nn` with `bk_local(:,nn,nkp) ≈ +bka(:,na)` (`ifpos==1`). Zero-init then fill; fatal if any
  stays 0.
- Copy locals to public arrays: `wb(loop) = wb_local(loop)` (755-757); `bk(:,loop,nkp) =
  bk_local(:,loop,nkp)` (759-763). `wb` is **per-neighbour** (length `nntot`), not per-shell.

`utility_compar` (utility.F90:388-410): `rrp=|a-b|^2`, `rrm=|a+b|^2`; `ifpos=1` iff `rrp<eps8`,
`ifneg=1` iff `rrm<eps8`. Tolerance is on the **squared** distance (`1e-8`), effectively `1e-4`
in linear distance.

### gamma_only (767-901)
Requires `num_kpts==1`. Rebuilds keeping only one of each ± pair (dedup via `utility_compar` ifneg,
834-861), sets `wb(na) = 2 * wb_local(nn)` (848), `nntot → nntot/2`, and asserts each kept `bk`
equals the corresponding `bka` (855-859).

### permutations (905-925, `kmesh_bvectors_perm` 2280-2329)
Non-gamma only. Builds `nnord(nn,ik)` (perm mapping k-point-ik's b ordering to the k=1 reference),
`nninv` (inverse), `nnrev` (index of the reversed −b). Tolerance `tol = 1e-7` on each component
(2299, 2308). Used by postw90 / spread functional; not needed to reproduce standard `M_mn`.

---

## 8. Higher-order finite differences (non-default, summary)

`higher_order_n = N > 1` (default 1). `max_shells_h = N(4N^2+15N+17)/6`, `num_nnmax_h = 2*max_shells_h`
(readwrite.F90:919-922). Two sub-modes:
- `higher_order_nearest_shells=.false.` (default when N>1): find first-order shells + weights, then
  `kmesh_shell_reconstruct` (1678-1742) adds the `2b, 3b, ..., Nb` shells with
  `dnn(order shell) = dnn(shell)*order` (1721) and rescaled weights
  `w_{order} = w_1 * prod_{j≠order} (j^2)/(j^2 - order^2) / order^2` (1726-1735). `nntot` multiplied
  by N (365). Extra b-vectors `bk_local(:,nnx2) = nn*bk_local(:,nnx)` (471), neighbour images found
  via `floor(kpt_latt+bk_latt+1e-6)` (491-506), warns "experimental".
- `higher_order_nearest_shells=.true.`: fully general multiset of monomial equations
  `num_of_eqs = (1+n)(1+2n)`; the per-kpoint B1 block (611-657) runs and enforces the generalized
  completeness for orders `1..N`.
For the standard reimplementation, `higher_order_n=1` and this whole section is inert.

---

## 9. `kmesh_write` — the `.nnkp` file (kmesh.F90:958-1128)

Writes `<seedname>.nnkp` (formatted). Structure:
- header `# File written on <date> at <time>`, `calc_only_A : <L>`.
- `begin/end real_lattice` — 3 rows `real_lattice(i,1:3)`, format `3(f12.7)` (1026-1030).
- `begin/end recip_lattice` — from `utility_recip_lattice_base` (2π), 3 rows `f12.7` (1033-1038).
- `begin/end kpoints` — `num_kpts`, then each `kpt_latt(1:3,nkp)` as `3f14.8` (fractional) (1041-1046).
- `begin/end projections` or `spinor_projections` / optional `auto_projections`.
- **`begin nnkpts`** (1100-1108): first line `nntot` (`i4`), then `num_kpts*nntot` lines, each:
  ```
  write(...,'(2i8,3x,3i4)') nkp, nnlist(nkp,nn), nncell(1:3,nkp,nn)
  ```
  i.e. `k-index  neighbour-k-index  G1 G2 G3`. This is exactly the data needed for `M_mn(k,b)`
  (M&V Eq. 25). Note: weights `wb` and b-vectors `bk` are **NOT** written to `.nnkp`; they are
  recomputed by whoever needs them (or read from `.wout`).
- `begin/end exclude_bands`.

---

## 10. What is reported in `.wout` (for isolation testing)

From the `iprint>0` write blocks in `kmesh_get`:
- "Distance to Nearest-Neighbour Shells": table of `ndnn`, `dnn(ndnn)/lenconfac`, `multi(ndnn)`
  (kmesh.F90:201-217). Units Ang^-1 or Bohr^-1.
- (`iprint>=4`) full list of b-vectors and lengths (219-254).
- "The following shells are used: ..." — the selected `shell_list` (277-297).
- "Shell # Nearest-Neighbours" table: `ndnn`, `nnshell(1,ndnn)` (414-584).
- "Completeness relation is fully satisfied [Eq. (B1) ...]" (660).
- **"b_k Vectors and Weights"** table (703-717): for `i=1..nntot`,
  `bk_local(1:3,i,1)/lenconfac` and `wb_local(i)*lenconfac**2`. Weight units Ang^2 or Bohr^2.
- **"b_k Directions"** table (718-730): `bka(1:3,i)/lenconfac` for `i=1..nnh`.
- gamma_only: reduced b_k/weights table (868-887).

**Validation recipe (kmesh in isolation):** given `real_lattice` and `kpt_latt`, reproduce
(1) shell distances `dnn` and multiplicities `multi`; (2) selected `shell_list`; (3) per-shell
`bweight` (V Σ⁻¹ Uᵀ (1,1,1,0,0,0) on the `(xx,yy,zz,xy,yz,zx)` amat); (4) the full ordered
`bk(:,:,1)` and `wb`; (5) verify `Sum_b wb(b) b⊗b = I` to `kmesh_tol=1e-6`; (6) the `.nnkp`
`nnkpts` block (`nnlist`, `nncell`). Match ordering exactly by reproducing supercell_sort +
`internal_maxloc` low-index tie-break, or the b-vector indices will permute.

---

## 11. Default parameter reference

| symbol | value | source |
|---|---|---|
| `kmesh_input%tol` (`kmesh_tol`) | 1.0e-6 | types.F90:151 |
| `search_shells` | 36 | types.F90:149 (error strings saying "default=30/36" are stale) |
| `search_supcell_size` | 5 | types.F90:150 |
| `nsupcell` (module const, sort) | 5 | kmesh.F90:68 |
| `max_shells` | 6 | types.F90:128 |
| `num_nnmax` | 12 | types.F90:129 |
| `higher_order_n` | 1 | types.F90:141 |
| `eps5 / eps6 / eps7 / eps8` | 1e-5 / 1e-6 / 1e-7 / 1e-8 | constants.F90:70-76 |
| `bohr` | 0.5291772108 Å | constants.F90:219 |
| `twopi` | 2π | constants.F90:56 |

SVD singular-value floor: **eps5 (automatic)**, **eps7 (fixed / from_file)** — not unified.
