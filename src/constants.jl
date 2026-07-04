# Physical constants and tolerances, matched to the reference Wannier90 (src/constants.F90).

"""
Bohr radius in Ångström. The reference Wannier90 defaults to **CODATA2006** (compile flags select
2010/2018/2022); we match the default so bohr-specified cells reproduce reference output exactly.
"""
const BOHR = 0.52917720859

"2π."
const TWOPI = 2 * π

# Default numerical tolerances (mirror Wannier90 defaults).
const KMESH_TOL_DEFAULT = 1.0e-6      # shell-distance / B1 degeneracy tolerance
const CONV_TOL_DEFAULT  = 1.0e-10     # spread-convergence tolerance (Ų)
