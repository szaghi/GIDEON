# Sources

| Source | What GIDEON takes from it |
|---|---|
| ADAM, `adam_fdv_operators_library` and its FNL backend (ff879596; Zaghi, Rossi) | the operator set, the pair form, the scaling by the step, the 3-D signatures; verbatim copies in `notes/adam-reference/` are the reference of the equivalence test, never the source of a row |
| B. Fornberg, *Generation of finite difference formulas on arbitrarily spaced grids*, Math. Comp. 51 (1988) 699–706 | the idea of deriving every finite-difference weight on any node set; GIDEON solves the moment system exactly instead of running the recursion |
| P. McCorquodale, P. Colella, *A high-order finite-volume method for conservation laws on locally refined grids*, Comm. App. Math. Comp. Sci. 6 (2011) 1–25 | the deconvolution between cell averages and point values, `<q> = q + h²/24 q'' + ...` (as ADAM cites it; not re-read for GIDEON) |
| SETH (`~/fortran/SETH`) | the conventions: compile-time rows, explicit-shape windows, SETH's face, the test and verification method, the mutation harness |

Every coefficient is derived by `scripts/gen_gideon.py`; the papers name the constructions, they are not transcribed.
