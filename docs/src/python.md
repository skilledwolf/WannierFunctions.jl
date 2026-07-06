# Using from Python

Python users have two good routes into the package. Neither requires writing any Julia
beyond a line or two of setup.

## Route 1: the drop-in binaries (no interop at all)

If you drive Wannier90 from Python today — `subprocess`, AiiDA, ASE, custom pipeline
scripts — the drop-in CLIs slot in unchanged, because they are ordinary executables reading
and writing the standard files:

```python
import subprocess
subprocess.run(["wannier90.jl", "-pp", "silicon"], check=True)
# ... pw2wannier90 ...
subprocess.run(["wannier90.jl", "silicon"], check=True)
subprocess.run(["postw90.jl", "silicon"], check=True)
```

(Install the launchers once with `using WannierFunctions; install_cli()` — see
[Getting started](getting-started.md).) Everything downstream keeps parsing the same
`.wout`/`.dat` files it parsed before.

## Route 2: juliacall — the library API with numpy in and out

For results as data rather than files, use [juliacall](https://juliapy.github.io/PythonCall.jl/stable/juliacall/)
(`pip install juliacall`). It installs Julia automatically if none is found, and Julia
arrays come back as zero-copy numpy views.

One-time setup (installs the package into juliacall's Julia environment):

```python
from juliacall import Main as jl
jl.seval('import Pkg; Pkg.add(url="https://github.com/skilledwolf/WannierFunctions.jl")')
```

Then the whole library API is available; this is the [Getting started](getting-started.md)
workflow, verbatim, from Python:

```python
import numpy as np
from juliacall import Main as jl
jl.seval("using WannierFunctions")

model = jl.read_model("diamond")                 # .win/.amn/.mmn/.eig
res   = jl.run_wannier(model)                    # keywords work: win_max=17.0, froz_max=6.4, ...

res.spread.Ω                                     # 2.320904915 (attribute access, unicode and all)
centres = np.asarray(res.spread.centres)        # (3, num_wann) numpy array, Cartesian Å

H = jl.hamiltonian_operator(model, res)          # interpolable H(R)
E = np.asarray(jl.bands(H, [[0.0, 0, 0], [0.5, 0, 0.5]]))   # plain lists of k-points are fine
```

`np.asarray` on a returned Julia array is a **view**, not a copy — slice it, plot it, or
`copy()` it if you need it to outlive the Julia object.

### Notes and gotchas

- **First-call latency.** Julia compiles on first use: expect a few seconds on the first
  `run_wannier`/`bands` call of a session, native speed afterwards. Long-running processes
  (Jupyter, workflow daemons) amortise this to nothing.
- **Keywords.** Julia keyword arguments are Python keyword arguments:
  `jl.run_wannier(model, win_max=17.0, froz_max=6.4)`. Unicode field names work as
  attributes (`res.spread.Ω`); if your editor makes that awkward, `getattr` or the ASCII
  aggregates (`res.spread.centres`, `.spreads`) cover most needs.
- **k-point lists.** Functions taking a list of k-points (`bands`, interpolators) accept a
  plain Python list of 3-lists as shown. For a numpy `(n, 3)` array, pass
  `[list(k) for k in kpts]` (row per k-point).
- **Matrices are column-major.** Julia arrays map to Fortran-ordered numpy arrays;
  `centres[:, n]` is the centre of WF `n` either way, but be mindful when reshaping.
- **Threading.** Set the environment variable `PYTHON_JULIACALL_THREADS=auto` *before*
  importing juliacall to enable the package's threaded k/R loops.
- **Which Julia runs?** juliacall manages its own Julia and environment by default; set
  `PYTHON_JULIACALL_EXE` to reuse a system Julia install.

The snippets above are tested against `examples/data/diamond` in the repository (Ω =
2.320904915 Å², reference-exact). Anything in the [How-to guides](howto.md) — Berry-physics
modules, transport, SCDM, TB-model input — works the same way through `jl.<function>`.

## What about a native pip package?

A thin `wannierfunctions` PyPI wrapper (juliacall under the hood, pythonic signatures) is a
possible future step once the package is registered in the Julia General registry. If that
would matter for your workflow, please open an issue.
