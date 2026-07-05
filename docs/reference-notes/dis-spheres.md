# Disentanglement in spheres in k-space (`dis_spheres`)

Implementation-grade notes extracted from the Wannier90 v3.x reference source.
File:line citations are to `reference/wannier90/src/` unless otherwise noted.
Companion to `disentanglement.md` (the main SMV subspace selection) — this note
covers **only** the `dis_spheres` modification, which decides, per k-point,
*whether* the SMV disentanglement runs at all.

The feature (author-tagged `GS-...` in the source) lets you restrict the
optimally-connected-subspace optimization to k-points lying inside one or more
spheres in reciprocal space. **Everywhere else the subspace is a trivial fixed
selection of absolute band indices — no window scan, no optimization.** Read the
semantics carefully: spheres mark *where* disentanglement happens; the default
outside is the cheap fixed path, not the expensive one.

--------------------------------------------------------------------------------
## 0. Data structure and keywords

`type dis_spheres_type` (`wannier90_types.F90:153–159`):

```
integer :: first_wann = 1                 ! default 1
integer :: num        = 0                 ! default 0 (feature OFF)
real(dp), allocatable :: spheres(:, :)    ! shape (4, num): rows = (kx,ky,kz,radius)
```

Note the **column-major storage `spheres(4, num)`**: index 1..3 is the fractional
k-centre, index 4 is the radius. `spheres(1:3, i)` = centre of sphere `i`,
`spheres(4, i)` = its radius.

Keywords (parsed in `w90_wannier90_readwrite_read_disentangle`,
`wannier90_readwrite.F90:791–880`):

| Keyword | Type | Default | Meaning | Validation |
|---|---|---|---|---|
| `dis_spheres_num` | int | `0` | number of spheres; `0` disables the feature | `< 0` → fatal "cannot be negative" (`:855–857`) |
| `dis_spheres_first_wann` | int | `1` | absolute band index of the first band kept at out-of-sphere k (see §2) | `< 1` → fatal (`:844–846`); `> num_bands - num_wann + 1` → fatal (`:848–850`) |
| `dis_spheres` | block | — | `dis_spheres_num` rows, each `kx ky kz radius` | required if `num > 0` (`:868–870`); each radius must be `>= 1.0e-15` (`:872–876`) |

`first_wann` and `num` are read **unconditionally** (always present with their
type defaults). The `dis_spheres` block is read only when `num > 0`:
`allocate(spheres(4, num))` then `get_keyword_block(..., rows=num, columns=4, ...)`
(`:859–866`).

**`first_wann` upper bound** `num_bands - num_wann + 1` guarantees that the fixed
window `[first_wann, first_wann + num_wann - 1]` (§2) stays within
`[1, num_bands]`.

--------------------------------------------------------------------------------
## 1. Block format, parsing, and radius units

The block reader is `w90_readwrite_get_keyword_block`
(`readwrite.F90:3057–3231`), called with `rows = dis_spheres_num`, `columns = 4`.
Each data line is read `read(dummy,*) (r_value(i,counter), i=1,4)` — free-format,
one row per line:

```
Begin dis_spheres
  kx  ky  kz  radius
  ...
End dis_spheres
```

- `kx ky kz` — the sphere **centre, in fractional (reduced) reciprocal
  coordinates** (same convention as the `kpoints` block; i.e. crystal
  coordinates of `recip_lattice`).
- `radius` — see units below.

### Radius units: Å⁻¹ **in the 2π convention** (this is the trap)

The comparison in `dis_windows` builds a **Cartesian** k-displacement
`dk = df · recip_lattice` and tests `dot(dk,dk) < radius²` (§2). Because
`recip_lattice` is built as

```
recip_lat = twopi * inv3(real_lat) / volume        ! utility.F90:346–351
```

with `real_lat` in **Å**, the rows of `recip_lattice` are the primitive
reciprocal vectors **including the factor 2π**, in Å⁻¹. Therefore:

> **The radius is in Å⁻¹ and must be compared against distances measured with the
> 2π-including reciprocal metric.** A radius of `0.2` means `|Δk|_cart < 0.2`
> where `b_i = 2π/a_i` for an orthorhombic cell. Do **not** interpret it as
> "0.2 in fractional units" nor as "0.2 Å⁻¹ without 2π".

### Optional units line and the `bohr` landmine

`get_keyword_block` supports an optional units line (`ang` / `bohr`) as the first
line inside the block if the block has `rows+1` lines (`readwrite.F90:3188–3201`):

- `ang` (or absent) → no conversion, `lconvert = .false.`.
- `bohr` → `lconvert = .true.` and, at `:3219–3223`, **the entire block is scaled
  `r_value = r_value * bohr`** where `bohr = bohr_angstrom_internal`
  (`0.52917720859` for the default CODATA2006 build, `constants.F90:182`).

**Trap:** the `bohr` scaling multiplies **all four columns**, not just the radius —
it would also scale the fractional k-centre `kx,ky,kz`, which is almost certainly
wrong for a user. In practice the reference tests never put a units line inside
`dis_spheres` (units default to Å⁻¹), so this path is untested. A Julia port
should treat a `bohr` units line here as unsupported / suspect rather than
faithfully replicate the whole-block scaling.

--------------------------------------------------------------------------------
## 2. How spheres restrict disentanglement (`dis_windows`)

All of the sphere logic lives inside the per-k loop of `dis_windows`
(`disentangle.F90:886–1136`), specifically the block at `:1013–1035`. It runs
**after** the normal outer-window scan has set `ndimwin(nkp)` / `nfirstwin(nkp)`
and **before** the frozen (inner-window) scan.

### 2a. Normal outer-window result (unchanged, computed for every k)

For each k (`:995–1011`), scanning the ascending eigenvalues `eigval_opt(:,nkp)`:

```
imin = first band index i with win_min <= E_i <= win_max
imax = last  band index i with           E_i <= win_max
ndimwin(nkp)   = imax - imin + 1
nfirstwin(nkp) = imin
```

A fatal error is raised first if the outer window is empty at this k
(`:982–993`): `E_1 > win_max` or `E_{num_bands} < win_min`.

### 2b. Sphere membership test (`:1016–1034`)

```fortran
if (dis_spheres%num > 0) then
  dis_ok = .false.
  do i = 1, dis_spheres%num
    dk = kpt_latt(:, nkp) - dis_spheres%spheres(1:3, i)   ! fractional difference
    dk = matmul(anint(dk) - dk, recip_lattice(:, :))      ! -> Cartesian Å^-1
    if (abs(dot_product(dk, dk)) < dis_spheres%spheres(4, i)**2) then
      dis_ok = .true.
      exit
    end if
  end do
  if (.not. dis_ok) then                    ! k in NO sphere
    dis_manifold%ndimwin(nkp)   = num_wann
    dis_manifold%nfirstwin(nkp) = dis_spheres%first_wann
  end if
end if
```

Precise semantics (get every one of these right in a port):

1. **Fractional difference then minimum-image fold.**
   `df = kpt_latt(:,nkp) - centre_i` (both fractional).
   `df_folded = anint(df) - df` folds each component into `[-0.5, 0.5]`, i.e. it
   selects the **nearest periodic image** of the centre. (`anint` rounds to
   nearest integer; the overall sign is irrelevant because the vector is squared
   next.) This is essential so a sphere near a BZ boundary wraps correctly.

2. **To Cartesian via `matmul(row_vector, matrix)`.**
   `dk(j) = Σ_{i=1..3} df_folded(i) * recip_lattice(i, j)`.
   The **rows** of `recip_lattice` are the reciprocal vectors
   (`recip_lattice(1,:) = b_1`, etc.). In column-major Julia this is a transpose
   trap: `dk = df_folded' * recip_lattice`, equivalently
   `dk = recip_lattice' * df_folded` — pin the index order to `(i,j)` above.
   `recip_lattice` includes 2π (§1), so `dk` is in Å⁻¹ (2π convention).

3. **Strict inside test.** inside sphere `i` iff
   `dot(dk,dk) < spheres(4,i)^2` — **strict `<`**. A k exactly on the sphere
   surface is **excluded**; the centre itself (`dk = 0`) is included.
   `abs(dot_product(dk,dk))` is `dot(dk,dk)` (already non-negative); the `abs`
   is cosmetic.

4. **Union / early exit.** `dis_ok` becomes true if k is inside **any** sphere
   (`exit` on first hit). Spheres compose as a **union**.

### 2c. Effect of being OUT of every sphere (`:1030–1033`)

If `dis_ok == .false.`, the window result from §2a is **overwritten**:

```
ndimwin(nkp)   = num_wann
nfirstwin(nkp) = dis_spheres%first_wann
```

Consequences, made explicit:

- The subspace at this k is the **fixed set of absolute band indices**
  `[first_wann, first_wann + num_wann - 1]` — `num_wann` states, taken verbatim.
  It bypasses the energy window entirely.
- `dis_win_min` / `dis_win_max` have **no lasting effect** at an out-of-sphere k
  except for the "outer window empty" fatal check at `:982–993` — the `imin/imax`
  computed in §2a are discarded here.
- Since `ndimwin == num_wann`, the eigenvalue slim-down (`:1127–1130`) copies
  bands `first_wann … first_wann+num_wann-1` down into slots `1 … num_wann`, and
  `lwindow(j,nkp) = .true.` is set for exactly those absolute `j`
  (`disentangle.F90:203–206`, used later by writers/`postw90`).

### 2d. Effect of being INSIDE a sphere

`dis_ok == .true.` → nothing is overwritten; the k-point disentangles **as
usual** with the full outer-window `ndimwin(nkp)` from §2a. This is the only place
the SMV optimization has any rotational freedom.

### 2e. Guard shared with the normal path (`:1037–1042`)

After the sphere block, the usual check fires:
`if (ndimwin(nkp) < num_wann)` → fatal "Energy window contains fewer states than
number of target WFs". For an out-of-sphere k, `ndimwin` is exactly `num_wann`,
so it always passes.

--------------------------------------------------------------------------------
## 3. Composition with `dis_win_min/max`, frozen window, and the Z-matrix loop

### 3a. Order of operations in `dis_windows`

1. Outer-window scan → `imin/imax/ndimwin/nfirstwin` (§2a).
2. **Sphere override** (§2b–2c): out-of-sphere k gets
   `ndimwin=num_wann`, `nfirstwin=first_wann`.
3. `ndimwin < num_wann` guard (§2e).
4. **Frozen (inner-window) scan** (`:1048–1107`): loops over the **original**
   `imin..imax`, *not* the sphere-overridden range.

### 3b. Frozen + spheres is an untested / inconsistent path — flag, don't spec

The inner-window loop at `:1054` iterates `do i = imin, imax` using the **stale**
`imin/imax` computed in §2a — these are *not* updated when the sphere override
replaced `ndimwin/nfirstwin`. `ndimfroz(nkp)`, `indxfroz`, and `lfrozen` are then
built relative to that stale window, while `ndimwin/nfirstwin` describe the
`first_wann` fixed window. For an out-of-sphere k with a non-empty frozen window,
these two descriptions disagree.

The reference benchmark (`testw90_lavo3_dissphere`) has **no frozen states**
(`dis_win_min/max` only, no `dis_froz_*`), so this combination is exercised by no
test. **A Julia port should treat "`dis_spheres` together with a frozen inner
window" as unsupported / undefined** rather than replicate the stale-index
behaviour, unless a reference test is added.

### 3c. How out-of-sphere k enters the Z-matrix iteration

In `dis_extract` (`disentangle.F90:2716–2775`), the per-k update with
`ndimfroz(nkp) = 0` and `ndimwin(nkp) = num_wann`:

- `ndiff = ndimwin - ndimfroz = num_wann`, so ZHPEVX (`:2738`) diagonalizes the
  `num_wann × num_wann` Z-matrix and the update loop `:2759` takes
  `j = ndimwin-num_wann+1 … ndimwin-ndimfroz`, i.e. **all** `num_wann`
  eigenvectors. The retained subspace therefore spans the entire `ndimwin`-dim
  space: the projector is unchanged — **zero rotational freedom**.
- There is **no special short-circuit**; the fixed subspace falls out naturally
  because a full-rank eigenbasis reconstructs the same subspace every iteration.
- Its contribution to `womegai1` (`num_wann*wbtot - Σ w_j`) is therefore
  **constant across iterations**, and only in-sphere k contribute anything that
  can change. This is why, in the benchmark, `Omega_I` is bit-stable from
  iteration 1 (`Delta ~ 1e-16`).

### 3d. `Omega_I` normalization unchanged

`Omega_I = womegai / num_kpts` (`disentangle.F90:2910`) still divides by the
**total** `num_kpts`, including the frozen/fixed out-of-sphere k. Sphere
restriction changes the per-k subspaces, not the normalization.

--------------------------------------------------------------------------------
## 4. Output formatting (`.wout`)

Two blocks are emitted (both only when `dis_spheres%num > 0`):

**(a) DISENTANGLE header** (`wannier90_readwrite.F90:2103–2111`):

```
|  Number of spheres in k-space              :                 1             |
|   center n.   1 :     0.500   0.500   0.500,    radius   =   0.200         |
|  Index of first Wannier band               :                 1             |
```

Formats: count line `'(1x,a46,10x,I8,13x,a1)'`; each centre line
`'(1x,a13,I4,a2,2x,3F8.3,a15,F8.3,9x,a1)'` (centre `1:3` as `3F8.3`, radius as
`F8.3`); first-wann line `'(1x,a46,10x,I8,13x,a1)'`.

**(b) Energy Windows box** (`disentangle.F90:960–978`) is printed unconditionally
and is *not* sphere-aware — it just echoes `win_min/win_max` and the frozen
window. There is no per-k "which k are inside spheres" table in the `.wout`.

--------------------------------------------------------------------------------
## 5. Reference test: `testw90_lavo3_dissphere`

Path: `test-suite/tests/testw90_lavo3_dissphere/`. Driver
`WANNIER90_WOUT_OK`, `extract_fn = parse_wout`, benchmark
`benchmark.out.default.inp=LaVO3.win`. Jobconfig: `tests/jobconfig:230–233`.

### Input (`LaVO3.win`)

```
num_wann = 3
num_bands = 6
exclude_bands : 1-20,27-40
dis_win_min = 15.0
dis_win_max = 18.5
dis_spheres_first_wann = 1
dis_spheres_num = 1
Begin dis_spheres
  0.5  0.5  0.5  0.2
End dis_spheres
mp_grid = 6 6 6          ! 216 k-points
```

Orthorhombic cell `a=b=3.7201157722 Å`, `c=4.0549261917 Å` →
`b_1=b_2=2π/3.72 ≈ 1.689 Å⁻¹`, `b_3=2π/4.05 ≈ 1.549 Å⁻¹`. No frozen states.

### Worked geometry (single sphere, single k inside)

Sphere centre `(0.5,0.5,0.5)` fractional, radius `0.2 Å⁻¹`. On the 6×6×6 grid the
only grid point at the centre is `(0.5,0.5,0.5)` itself → `dk = 0`, inside. Its
nearest neighbours differ by `1/6` in one fractional component:
`|Δk|_cart = b_i/6 ≈ 1.689/6 ≈ 0.281 Å⁻¹` (or `1.549/6 ≈ 0.258` along z), both
`> 0.2`. So **exactly one of the 216 k-points is inside the sphere**; the other
215 get the fixed `ndimwin = num_wann = 3`, `nfirstwin = first_wann = 1`
(absolute bands 1–3 of the 6 in-window bands). This is why the disentanglement
converges instantly (`Delta ~ 1e-16` from iteration 1; benchmark lines 246–254).

### Benchmark anchor quantities (gauge-invariant — use these to validate)

Final state (`.wout`, benchmark lines 425–434):

| Quantity | Value | Abs tol | Rel tol | Label in `userconfig` |
|---|---|---|---|---|
| `Omega_I` | `7.457463597` Å² | `1.0e-6` | `1.0e-6` | `omegaI` |
| `Omega_D` | `0.000000000` Å² | `1.0e-6` | `5.0e-6` | `omegaD` |
| `Omega_OD` | `0.050664432` Å² | `1.0e-6` | `1.0e-6` | `omegaOD` |
| `Omega_Total` | `7.508128029` Å² | `1.0e-6` | `1.0e-6` | `omegaTotal` |

`Final Omega_I` printed by `dis_extract` = `7.45746360` Å² (benchmark line 257;
same value, fewer digits).

Final WF centres and spreads (all three WFs sit on the V site
`(1.860058, 1.860058, 2.027463) Å` by symmetry; benchmark lines 425–428):

| WF | centre (Å) | spread (Å²) |
|---|---|---|
| 1 | (1.860058, 1.860058, 2.027463) | 2.57756832 |
| 2 | (1.860058, 1.860058, 2.027463) | 2.57756832 |
| 3 | (1.860058, 1.860058, 2.027463) | 2.35299138 |

Sum of centres/spreads: `(5.580174, 5.580174, 6.082389)`, `7.50812803`.

Tolerances (`test-suite/tests/userconfig`, `[WANNIER90_WOUT_OK]`, lines 9–16):
`final_centres_x/y/z` abs=rel=`1.0e-5`; `final_spreads` abs=rel=`3.0e-6`;
`omegaI/omegaOD/omegaTotal` abs=rel=`1.0e-6`; `omegaD` abs=`1.0e-6` rel=`5.0e-6`.

--------------------------------------------------------------------------------
## 6. Traps summary (for the Julia port)

1. **Radius units = Å⁻¹ in the 2π convention.** `recip_lattice = 2π·inv(A)/V`
   (Å⁻¹, with 2π). Compare `dot(dk,dk) < radius²` with `dk` in that metric. Not
   fractional, not 2π-free.
2. **`matmul(anint(df)-df, recip_lattice)` — minimum-image fold + `row·matrix`.**
   `dk(j) = Σ_i (anint(df_i)-df_i)·recip_lattice(i,j)`; rows are the b-vectors.
   Transpose trap in column-major.
3. **Strict `<`.** On-surface k is excluded; centre included. Spheres compose as a
   **union** (first-hit `exit`).
4. **Out-of-sphere ⇒ fixed absolute bands `[first_wann, first_wann+num_wann-1]`,**
   bypassing the energy window. This is the *default* everywhere outside spheres,
   not an exceptional case. `dis_win_min/max` only survive as the empty-window
   fatal check at out-of-sphere k.
5. **`bohr` units line scales all 4 columns** (incl. the fractional centre) —
   untested; treat as unsupported.
6. **`dis_spheres` + frozen window uses stale `imin/imax`** for the inner-window
   scan — untested/inconsistent; treat "spheres + inner window" as unsupported.
7. **Defaults:** `dis_spheres_num = 0` (feature off), `dis_spheres_first_wann = 1`.
   `first_wann` must satisfy `1 <= first_wann <= num_bands - num_wann + 1`.
8. **Block storage is `spheres(4, num)`** (column-major): `spheres(1:3,i)` centre,
   `spheres(4,i)` radius. Each `.win` row is `kx ky kz radius`.
9. **`Omega_I` still divided by full `num_kpts`.** Fixed out-of-sphere k
   contribute a constant, iteration-invariant amount; only in-sphere k carry
   rotational freedom → instant convergence when few k are inside.
