# Physical constants and tolerances, matched to the reference Wannier90 (src/constants.F90).

"Bohr radius in Ångström (CODATA, matches Wannier90 `bohr_angstrom_internal`)."
const BOHR = 0.529177210903

"2π."
const TWOPI = 2 * π

# Default numerical tolerances (mirror Wannier90 defaults).
const KMESH_TOL_DEFAULT = 1.0e-6      # shell-distance / B1 degeneracy tolerance
const CONV_TOL_DEFAULT  = 1.0e-10     # spread-convergence tolerance (Ų)
