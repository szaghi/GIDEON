# Compact-WENO rows

**Problem.** A compact scheme resolves smooth waves with fewer points than an explicit one, but its face values are
coupled along the whole line: at a discontinuity it rings, and the ringing reaches every face of the line. A hybrid
(Pirozzoli 2002) keeps the compact rows where the solution is smooth and takes a nonlinear (WENO, TENO) value at the
faces a sensor flags. For the solve, a flagged face's row becomes the identity and its right side the value given.

## The pieces

| Piece | Who | Here |
|---|---|---|
| the compact operator's right side and band | GIDEON | `gideon_fv_face_pbl5(q(-1:1), r)`, `FV_FACE_PBL5_LEFT` (the linear upwind-biased compact face value, flow to the right) |
| the solve with rows replaced | GIDEON | `gideon_tri_periodic_masked(n, a, b, c, mask, d, w, z)` |
| the sensor and the threshold | the caller | the second difference of the averages over the line's scale, > 1/4 |
| the nonlinear value at a flagged face | the caller (SETH) | here the first-order upwind value, the cell average on the left: a stand-in |

`mask(i)` true makes row $i$ the identity: the solve returns `d(i)` there exactly, and couples it to its unflagged
neighbours through their rows (decision D7: GIDEON owns the operator and the solve, the switch is the caller's).

## Code

Face values of the cell averages of a smooth bump and a square wave on 100 cells: the explicit upwind-biased `bl5`,
the compact `pbl5` on every face, and `pbl5` with the faces the sensor flags given the upwind value:

<<< @/examples/src/cookbook_compact_weno.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_compact_weno.txt-->
```

- **Compact everywhere.** The most accurate on the bump, and the worst at the square wave: its overshoot is larger
  than the explicit scheme's, since the band spreads the jump's error along the line.
- **Flagged.** Six faces (three at each jump) take the given value; the overshoot is gone, and the bump keeps the
  compact accuracy: the band's inverse decays geometrically (by about 0.55 per cell for `pbl5`), so the flagged faces'
  influence fades well before the bump. [A5](/verification/coupling) measures the order away from a step: 5.00 and
  7.00 for `pbl5` and `pbl7`.
- The sensor matters: one normalised by the local magnitude flags the bump's tails, where $u \approx 0$, and gives the
  bump first-order values; the line's scale (or a pressure, as in Jameson's switch) does not.

**Reference.** Roadmap issue #2, phase 4 (4.1, 4.2); the operators and the solvers: [Compact operators](/guide/compact).
Pirozzoli (2002), as recalled ([Sources](/guide/sources)).
