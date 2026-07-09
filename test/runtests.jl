# Test driver. The physics/behavioral suites live in thematic files (split from the original
# monolithic runtests.jl); common.jl provides the reference-tree locations, shared models,
# and helpers every suite uses.
include("common.jl")

include("qa.jl")                    # Aqua static package-quality checks (+ JET when present)
include("wannierise_core.jl")       # units, reference validation, disentanglement, optimizer
include("operators_io.jl")          # TBOperator, checkpoints, Γ-only, spinors, plotting, SCDM
include("postw90_berry.jl")         # Berry curvature/AHC, geninterp, Kubo, morb, DOS, SHC
include("postw90_responses.jl")     # spin, kpath, gyrotropic, shift current, tetrahedron
include("symmetry_slwf.jl")         # BZ reduction, SLWF+C, symmetry-adapted WFs
include("driver_transport.jl")      # postw90.jl driver, ballistic transport, Γ-only parity
include("extras_features.jl")       # higher-order FD, SS functional, injection, DFTK, CLI
