# Data provenance and licensing

The input files in this directory (`gaas.*`, `diamond.*`) originate from the
[Wannier90](https://github.com/wannier-developers/wannier90) test suite
(`test-suite/tests/testw90_example01` and `testw90_example05`) and are
redistributed here unmodified (diamond.win: band-path/plot keywords appended)
for self-contained examples and tests.

Wannier90 is distributed under the **GNU General Public License v2**; these
data files remain under that license. They are *inputs consumed by* — not part
of — the MIT-licensed WannierFunctions.jl source code, and can be deleted
without affecting the package (examples and validation tests then fetch or
skip, respectively).
