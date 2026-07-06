# How-to guides

Task-oriented recipes. Each snippet is self-contained given a wannierised `model`/`res`
(Workflow A) or a checkpoint on disk; the `.win`-keyword route via `bin/postw90.jl` is noted
where it exists. The docstrings (see [API](api.md)) carry the full keyword lists.

## Run the drop-in binaries

```bash
julia --project=. bin/wannier90.jl -pp seedname     # .nnkp generation (byte-identical)
julia --project=. bin/wannier90.jl seedname         # wannierise; writes .wout/.chk/_hr/_tb/bands
julia --project=. bin/postw90.jl seedname           # all postw90 modules per the .win keywords
julia --project=. bin/w90chk2chk.jl -export seedname   # .chk → .chk.fmt   (-import: reverse)
```

`postw90.jl` dispatches on the same keywords as `postw90.x`: `berry`/`berry_task`
(ahc, morb, kubo, sc, shc, kdotp), `gyrotropic`, `dos`, `kpath`, `kslice`, `geninterp`,
`boltzwann`, `spin_moment`, plus the module-specific controls (`berry_kmesh`, `fermi_energy*`,
`kubo_*`, `shc_*`, `dos_*`, `boltz_*`, …) — and writes the reference-named `.dat` files in the
reference formats. See [compatibility](wannier90-compat.md) for the exact keyword semantics.

## Disentanglement recipes

```julia
# energy windows (classic)
res = run_wannier(model; win_min=-5.0, win_max=17.0, froz_max=6.4)

# PDWF: freeze by projectability, not energy (robust against intruding states)
dis = disentangle(model; win_min=εF-20, win_max=εF+10,
                  froz_proj=true, proj_min=0.02, proj_max=0.95)

# dis_spheres: disentangle only inside k-spheres (k-localised entanglement)
dis = disentangle(model; win_max=20.0, spheres=[(SVector(0.5,0.5,0.5), 0.4)],
                  sphere_first_wann=1)

# SCDM for the entangled case: energy-smeared automatic projections
A = scdm_projections(model; dir="UNK_dir", mode=:erfc, mu=εF, sigma=2.0)   # from UNK files
# (or in-memory: wannier_model(scfres; num_wann, num_bands, scdm_mode=:erfc, scdm_mu, scdm_sigma))
```

## Localisation variants

```julia
res = wannierise(model; algorithm=:rcg)                     # native CG (default)
res = wannierise(model; algorithm=:w90, num_iter=500)       # reference-exact trajectory
res = wannierise(model; algorithm=:gamma)                   # Γ-only real-orthogonal sweeps
# via .win keywords, all through run_wannier(seedname):
#   guiding_centres, precond, slwf_num/slwf_constrain/slwf_lambda/slwf_centres (SLWF+C),
#   site_symmetry + .dmn (symmetry-adapted WFs), use_ss_functional (Stengel–Spaldin),
#   gamma_only, higher_order_n (higher-order finite differences)
```

Symmetry-adapted WFs read `seedname.dmn` (pw2wannier90's `write_dmn`) and work for both the
isolated case and combined with disentanglement (the constrained Ω_I optimiser); the gauge
then satisfies `U(Rk) = d(R)·U(k)·D(R)†` to the `symmetrize_eps` tolerance.

## Berry-phase physics (postw90 parity)

All modules consume a `BerryModel` (from `seedname.chk(.fmt)` + `.eig` + `.mmn`, or a TB file,
or in-memory operators):

```julia
bm = BerryModel("seedname")                                  # honours use_ws_distance

ahc = anomalous_hall(bm; fermi_energy=12.6279, kmesh=(25,25,25))          # S/cm
scan = ahc_fermiscan(bm; fermi_energies=11.0:0.1:13.0, kmesh=(25,25,25),
                     adpt_kmesh=5, adpt_thresh=100.0)                     # adaptive refinement
kubo = optical_conductivity(bm; fermi_energy, freqs=0:0.01:7, kmesh=(25,25,25))
σS, σA = kubo_S(kubo, 1, 1, 10), kubo_A(kubo, 1, 2, 10)                  # tensor components

mm = MorbModel("seedname"; transl_inv_full=true)             # needs .uHu
M = orbital_magnetisation(mm; fermi_energy, kmesh=(25,25,25))            # μ_B/cell

sm = ShcModel("seedname")                                    # needs .spn; ShcRyooModel: .sHu/.sIu
shc = shc_fermiscan(sm; fermi_energies, kmesh=(15,15,15))                # (ħ/e)·S/cm
shc = shc_tetra(sm; kmesh=(15,15,15), fermi_energies)                    # tetrahedron method

sc = shift_current("seedname"; fermi_energy, freqs, kmesh=(25,25,25))    # σ_abc(0;ω,−ω)
η  = injection_current(bm; freqs, fermi_energy, kmesh=(12,12,12), smr_width=0.1)
kp = kdotp(bm; kpoint=[0,0.5,0], bands=[4,5])                            # Löwdin orders 0–2
gy = gyrotropic(mm; tasks=(:D0,:Dw,:C,:K,:NOA,:dos), fermi_energies, freqs,
                kmesh=(25,25,25), smr_width=0.1)
```

DOS, transport coefficients and structure plots:

```julia
es, dos = density_of_states(bm; energies=8:0.01:14, kmesh=(25,25,25))    # adaptive smearing
es, dos, up, dn = density_of_states(bm; energies, spin=SpinModel("seedname"))  # spin-decomposed
r = boltzwann(bm; kmesh=(30,30,30), relax_time=10.0, mus=[5.0], temps=[300.0],
              win=(minimum(eig), maximum(eig)))                          # σ, S, κ, TDF

kp = kpath(bm; segments=kpath_segments(win), tasks=(:bands,:curv), fermi_energy)
write_kpath("seedname", kp)                                              # reference file set
ks = kslice(bm; b1=[...], b2=[...], mesh=(200,200), tasks=(:bands,:curv), fermi_energy)
write_kslice("seedname", ks.coords, ks.bands; curv=ks.curv)
tabulate_3d(...) ; write_frmsf(...)                                      # FermiSurfer export
```

## Ballistic transport (Landauer)

For a system periodic along one axis (`transport = true`, `transport_mode = bulk` in the
`.win`, or directly):

```julia
model, win, res = run_wannier("seedname")
energies, qc, dos = run_transport(model, win, res)     # writes _qc.dat/_dos.dat (+_htB.dat)
# or from explicit principal-layer blocks:
energies, qc, dos = transport_bulk(H00, H01; win_min=-3.0, win_max=3.0, energy_step=0.01)
```

`tran_lcr` (lead–conductor–lead) is not implemented — see
[compatibility](wannier90-compat.md).

## Symmetrised (irreducible-BZ) integration

With a `.sym` file (pw2wannier90) or a point group filtered from Oₕ:

```julia
sym = read_sym("seedname.sym")                        # Cartesian rotations
reps, wts, _ = irreducible_kmesh((25,25,25), sym; kaction=:cart, lattice=bm.lattice)
ahc, n = anomalous_hall_sym(bm, sym; fermi_energy, kmesh=(25,25,25))     # pseudovector rule
M, n   = orbital_magnetisation_sym(mm, sym; fermi_energy, kmesh=(25,25,25))
es, d, n = density_of_states_sym(bm, sym; energies, kmesh=(25,25,25))
```

Typical speed-up: the full mesh collapses to its irreducible wedge (e.g. 78/512 points for
bcc Fe's magnetic group) with results equal to the full-BZ sum.

## Tight-binding model input (no chk/mmn/eig)

Interpolate the entire physics stack from a `_hr.dat` or `_tb.dat` alone — the WannierBerri
`System_tb` workflow:

```julia
bm = tb_model("seedname_tb.dat")                      # H(R) + r(R) → full BerryModel
anomalous_hall(bm; fermi_energy, kmesh=(25,25,25))
# hr-only files work for everything that doesn't need the position operator
```

`write_tb`/`write_hr` produce the same files from a finished run, so models round-trip.

## In-memory models (DFTK or any Julia DFT)

```julia
wmodel = wannier_model(scfres; num_wann=4)                        # SCDM, projection-free
wmodel = wannier_model(scfres, projections; num_wann, num_bands)  # explicit trial orbitals
# or fully manual, from your own arrays:
wmodel = wannier_model(; unit_cell, kpoints, mp_grid, num_wann, M, A, kpb, gpb, eig)
```

Everything downstream (disentanglement, localisation, `hamiltonian_operator`, the Berry
stack) treats these models identically to file-read ones.

## Checkpoint interchange

```julia
chk = read_chk("seedname.chk")            # or read_chk_fmt; write_chk / write_chk_fmt
```

Checkpoints interoperate with `wannier90.x`/`postw90.x` in both directions — their
`restart=plot` reproduces its own bands from our `.chk`, and our `BerryModel` consumes theirs.
