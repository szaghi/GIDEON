# Sources

| Source | What GIDEON takes from it |
|---|---|
| ADAM, `adam_fdv_operators_library` and its FNL backend (ff879596; Zaghi, Rossi) | the operator set, the pair form, the scaling by the step, the 3-D signatures; verbatim copies in `notes/adam-reference/` are the reference of the equivalence test, never the source of a row |
| B. Fornberg, *Generation of finite difference formulas on arbitrarily spaced grids*, Math. Comp. 51 (1988) 699–706 | the idea of deriving every finite-difference weight on any node set; GIDEON solves the moment system exactly instead of running the recursion |
| P. McCorquodale, P. Colella, *A high-order finite-volume method for conservation laws on locally refined grids*, Comm. App. Math. Comp. Sci. 6 (2011) 1–25 | the deconvolution between cell averages and point values, `<q> = q + h²/24 q'' + ...` (as ADAM cites it; not re-read for GIDEON) |
| SETH (`~/fortran/SETH`) | the conventions: compile-time rows, explicit-shape windows, SETH's face, the test and verification method, the mutation harness |
| S. K. Lele, *Compact finite difference schemes with spectral-like resolution*, J. Comput. Phys. 103 (1992) 16–42 | the compact (Padé) operators: centred node derivatives, staggered face values and derivatives. **Not read for GIDEON**: the rows are derived by meaning, and agree with the values the generator recalls from the paper (tier A, `LELE_RECALLED`); the paper is to be obtained for the cross-check and for its boundary schemes |
| D. V. Gaitonde, M. R. Visbal, *Padé-type higher-order boundary filters for the Navier–Stokes equations*, AIAA J. 38 (2000) 2103–2112 | the compact filters' form, the parameter $\alpha_f$. **Not read**: derived by meaning, agree with the recalled coefficients (tier A, `GV_RECALLED`) |
| S. Pirozzoli, *Conservative hybrid compact-WENO schemes for shock-turbulence interaction*, J. Comput. Phys. 178 (2002) 81–117 | the compact upwind face value (`gideon_fv_face_pbl5` is his linear scheme as recalled, tier A, `PIROZZOLI_RECALLED`) and the coupling: a flagged face's row replaced by the nonlinear value (`_masked`). **Not read**: to be obtained |
| A. Jameson, W. Schmidt, E. Turkel, *Numerical solution of the Euler equations by finite volume methods using Runge–Kutta time-stepping schemes*, AIAA Paper 81-1259 (1981) | the blend of second and fourth differences that the dissipation recipe writes with `gideon_fv_face_diss_c2` and `_c4` (the coefficients $k_2 = 1/2$, $k_4 = 1/32$ as commonly quoted; not read) |

Every coefficient is derived by `scripts/gen_gideon.py`; the papers name the constructions, they are not transcribed.
The papers of roadmap #2 marked *not read* were not available to the night that wrote P2: every value attributed to
them above is a recollection, checked against the derived rows and labelled as such, never a source.
