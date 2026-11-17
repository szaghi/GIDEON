# Audit of ADAM's FDV operators: what GIDEON must fix, keep and drop

Date: 2026-11-17. Source: ADAM `develop` at `ff879596`, working tree; read-only. Verbatim copies in
`notes/adam-reference/` (`adam_fdv_operators_library.F90`, `adam_fnl_fdv_operators_library.F90`,
`fdv_coefficients.py`, `fdv_coefficients.md`): the reference of the equivalence test, never compiled into GIDEON.

## 1. What ADAM has

| Part | Content | Production users in ADAM |
|---|---|---|
| FD centred | `FD0_CC` (point-value face interpolation), `FD1..FD6_CC` (node derivatives), orders 2–10 | PRISM curl; FLUME-MHD first derivatives |
| FV centred | `FV1_CC` (face value from averages), `FVR1..5_CC` (face derivatives from averages), flux-difference ladder | almost none |
| FV upwind (right, left) | `FV1_UR/UL`, `FVUR/FVUL1..4` ladders | none (tests only); **kept by owner decision for future use** |
| Deconvolution | `DECONV_A2P/P2A`, averages ↔ point values, order 2s+2 | none found |
| 3-D assemblies | gradient, divergence, curl, Laplacian | PRISM, FLUME-MHD |
| FNL backend | `*_dev` variants: the same numerics with a gathered-stencil signature (no non-contiguous sections inside `acc routine seq`) | PRISM FNL, FLUME FNL |

Design: run-time half-stencil `s`, tables `(S_MAX, S_MAX)` with `S_MAX = 5`, `!$acc declare copyin` device copies,
loops of run-time trip count. The pattern SETH's founding decision rejected (`SETH/notes/kernel-api-design.md` §2;
SETH's compile-time kernels measured 1.26–1.31x faster than ADAM's on the RTX 4070).

Not FDV clients although they could be: FLUME's Ducros sensor (inline 2nd-order velocity gradient,
`adam_flume_cpu_hybrid_kernels.F90`) and FLUME's viscous face fluxes (inline orders 2 and 4 with the Shu–Osher
`−Δx²G''/24` correction, ADAM `docs/applications/flume/numerics.md`, "Face derivatives").

## 2. Defects (measured with `notes/adam-reference/check_upwind_d1.F90`, gfortran 16.2, rerun 2026-11-17)

Test function: cell averages of sin 2πx, x0 = 0.3, h = 1/50 … 1/400; error against the exact cell average of the
derivative on cells j−1, j, j+1.

**A. `compute_derivative1_fv_lupwind` is wrong for s ≥ 2.** The `FV1_UL` rows are stored in ascending cell order
(`FV1_UL_S2 = [-1, 3]/2`, line 433) but read descending, `FV1_UL(m,s)*q(1-m)`: the row is mirrored twice. For s = 2,
with `q(-s:0)` and cell i at `q(0)`, the kernel computes `(-q_i/2 + 2 q_{i-1} - 3 q_{i-2}/2)/h`, a first-order
derivative at cell i−1, instead of a second-order one at cell i. Measured: first-order convergence for every s = 2…5,
matching none of cells j−1, j, j+1 at the design order (s = 5, h = 1/400: errors 0.38, 0.48, 0.57).

**B. `compute_derivative1_fv_rupwind`: the window contract and ADAM's own test disagree.** The kernel declares
`q(0:)` over `[0:1+s]` and its second call reads `q(1:1+s)` (line 1139); ADAM's `test_fdv_operators_trigonometric`
passes `sin0(i:i+s)`, one element short. The kernel then reads one element out of bounds and returns the derivative
of cell i+1 at the design order s (s = 5: 2.5e-5 → 9.3e-10 over three halvings against cell j+1; O(h) against cell
j). With `-fcheck=bounds`: `Index '2' of dimension 1 of array 'q' outside of expected range (1:0)` at line 1139.
Also: d1 upwind uses `q(0:)` windows, d2…d5 `q(1-s:)`: two conventions in one family.

**C. Unsupported orders return 0 silently.** `FD3_CC_S1`, `FD5/FD6_CC_S1,S2` are zero rows; `compute_derivative6_fv_rupwind`
returns 0 (`S_MAX = 5` cannot hold its 6-cell row). Pure, device-callable code with a wrong answer and no error.

**D. The tests print, they do not assert.** `test_fdv_operators_{trigonometric,step}` print L0/L2 errors and write
`.dat` files: no expected order, no exit status. A and B were never caught.

**E. Generator labels.** `fdv_coefficients.py` validates the tables against the committed tables, not against their
meaning, and its docstring says `FV1_UR` gives q at i+1/2 from cells i…i+s−1; the row [3, −1]/2 on cells (i, i+1) is
the value at i−1/2. The same face-labelling confusion as A.

**F. The upwind names flip between d1 and d2–d5.** ADAM's `compute_derivative1_fv_rupwind` reads the cells on the
right of the cell (i..i+s), its `compute_derivative{2..5}_fv_rupwind` those on the left (i−s..i); `lupwind` the
reverse. Found while mapping the windows for the equivalence test (`src/tests/gideon_test_adam.F90`); GIDEON names
the side the stencil lies on (`ur`: right side, `ul`: left side) for every operator.

**A, refined.** On a linear function the left-side d1 is right: both faces err by the same amount and the difference
cancels it. The defect shows from degree 2 (a quadratic: −4 instead of 0 at s = 2), which is why it converges at
first order instead of order s.

## 3. Rules GIDEON adopts so that A–F cannot recur

1. **One generator, checked by meaning.** Every row derived in exact rationals (Fornberg's weights on the node set for
   point values; the derivatives of the primitive function for cell averages) and checked as exact on polynomials up
   to its design degree, at its stated location; the compiled kernels re-checked by the tests (SETH's G1/G2).
2. **Mirrors are generated, never stored.** A left-upwind row is the right-upwind row mirrored, times (−1)ⁿ for the
   n-th derivative; a property test asserts left(mirrored q) = (−1)ⁿ right(q) (kills A).
3. **Explicit-shape windows with named bounds.** Every kernel takes `q(LO:HI)`, bounds named constants, the operator's
   location at offset 0 (node operators) or at SETH's face (cell 0 left of the face, 1 right); callers slice
   `q(i+LO:i+HI)`, so a wrong window is a bounds error at the call site (kills B).
4. **Unsupported = absent.** A combination without a kernel is a missing symbol, never a zero (kills C); the upwind
   ladders are generated through d6.
5. **Tests assert frozen bounds.** Design orders and tolerances written before the kernel (SETH's G3) (kills D).
6. **One signature for the 3-D operators**, the directional stencils (the FNL one): the FNL backend disappears; a
   host-only cube wrapper may remain for ADAM's CPU callers.
7. **Run-time order through a dispatch layer** (PRISM's `fdv_order`): one `select case` per call over the per-order
   kernels, uniform across cells.
8. **Equivalence with ADAM**: bitwise where ADAM is correct, documented deviations where it is fixed (A, B).

## 5. Equivalence, measured (2026-11-17, `gideon_test_adam`, gfortran 16.2 and nvfortran 26.9, debug and optimised)

| ADAM | GIDEON | Result |
|---|---|---|
| FD centred face, d1–d6; FV centred face, face d1–d5, cell d1–d6; deconvolution; FNL and CPU 3-D operators | `gideon_fd_*_c*`, `gideon_fv_*_c*`, `gideon_fv_point_c*`, `gideon_fd_average_c*`, `gideon_{fd,fv}_{gradient,divergence,curl,laplacian}_c*` | bitwise |
| `compute_reconstruction_r_fv_rupwind`, `compute_derivative1_fv_rupwind` (window `[0:1+s]`, cell at `q(1)`) | `gideon_fv_face_ur*`, `gideon_fv_d1_ur*` | bitwise |
| `compute_derivative{2..5}_fv_rupwind`, `compute_reconstruction_r_fv_rupwind_d{1..4}` | `gideon_fv_d*_ul*`, `gideon_fv_face_d*_ul*` (left side, audit F) | bitwise |
| `compute_derivative{2..5}_fv_lupwind`, `compute_reconstruction_r_fv_lupwind_d{1..4}` | `gideon_fv_d*_ur*`, `gideon_fv_face_d*_ur*` (right side, audit F) | within 64 eps: ADAM sums the farthest cell first |
| `compute_reconstruction_r_fv_lupwind`, `compute_derivative1_fv_lupwind`, s ≥ 2 | `gideon_fv_face_ul*`, `gideon_fv_d1_ul*` | ADAM wrong (A), GIDEON exact |
| zero rows (FD d3–d6, FV face d2–d5, cell d3–d6 at small s; upwind d6) | absent kernels; `gideon_fv_d6_u*` from 6 cells | ADAM returns 0 (C) |

Totals: 279 bitwise, 20 within tolerance, 8 defect cases reproduced, 43 ADAM zero rows.

## 4. Relation to SETH (owner decisions, 2026-11-17)

GIDEON and SETH stay separate libraries: no SETH → GIDEON library dependency. SETH's verification harness may depend on
GIDEON (test-only) to verify the hybrids that compose GIDEON's linear operators with SETH's nonlinear ones. Rows that
exist in both (SETH's `CEN4/6/8_C` = `FV1_CC_S2..S4`, `LIN*`, the KEEP-PE a_l = `FD1_CC`) are derived exactly in each;
GIDEON adopts SETH's face convention so the two cannot drift in meaning.
