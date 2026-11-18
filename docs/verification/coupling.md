# A5 Compact-WENO coupling

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A5 | Accuracy | A | CPU (gfortran R4P, R8P, R16P; nvfortran; ifx); GPU (OpenACC) | no | `gideon_test_coupling` |

## Case and scope

A hybrid compact–WENO scheme replaces the compact rows of the faces a sensor flags by the nonlinear values there. This
case checks what GIDEON owns of it (the operator and the masked solve, roadmap #2, 4.1 and 4.2) on one line with a
step: the flagged faces keep their values, the scheme stays conservative, and away from the step the compact order is
kept.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Data | $u = \sin 2\pi x + 0.5\cos(6\pi x + 0.3) + H(x)$, $H = 1$ on $[1/2, 1)$: two jumps on faces of a periodic line of $N$ cells, exact cell averages | [`gideon_test_coupling.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_coupling.F90) |
| Operators | `gideon_fv_face_pbl5`, `gideon_fv_face_pbl7` and their bands; `gideon_tri_periodic_masked` | |
| Flags | the faces within 3 cells of a jump; their value: the cell average on the left (the stand-in for WENO) | |
| Given values | every flagged face returns its value, bitwise | |
| Conservation | the flux differences summed in real128 vanish within $64\,\varepsilon \sum\vert f\vert$ | |
| Order | $N = 64 \dots 131072$; the error at the faces at least 0.2 from both jumps; a plateau of three resolved sizes, both observed orders $\ge$ design − 0.25; measured in R16P (required) | |
| Device | 1024 coupled lines of random data and masks, one per thread, against the host | |

## Program

```bash
./exe-r16/gideon_test_coupling
```

## Results

<!--@include: ./_tables/coupling-builds.md-->

## Findings

- Away from the step the coupled line keeps the compact order exactly: 5.00 and 7.00 in quadruple precision. The
  band's inverse decays geometrically from a flagged face (by about 0.55 per cell downstream and 0.18 upstream for
  `pbl5`), so the step's influence at a fixed distance falls faster than any power of $h$.
- The flagged faces return their values bitwise, in the interior of the line and on its periodic border.
- Conservation is a property of the flux form, whatever the face values: the check guards that every face has one
  value, the one both its cells use.
- Two floors of this test were corrected after its first runs, each stated in its header: the resolution floor
  omitted the cell averages' own round-off ($2\max\vert P\vert/h$, which grows like $N$), and the conservation sum was
  compared with $64\varepsilon$ in the working kind although a sum of $N$ terms rounds to $N\varepsilon$.
