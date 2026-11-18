# C8 The line solvers

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C8 | Contracts | C | CPU (gfortran R4P, R8P, R16P; nvfortran; ifx); GPU (OpenACC) | no | `gideon_test_solvers` |

## Case and scope

Every compact operator ends in a solve along a line. This case checks the twenty solvers of `gideon_line_solvers`
(tridiagonal and pentadiagonal, bounded and periodic, constant band, variable coefficients, masked rows, factored) against a
dense solve that shares no code with them, on bands harder than any GIDEON operator's, and on the device against the
host (roadmap #2, gates L3 and L5).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Bands | random strictly diagonally dominant: off-diagonals in $[-1, 1]$, the diagonal of either sign, its dominance margin in $[0.25, 1.25]$; lines of 3 (5) to 64 unknowns; 10 systems per length | [`gideon_test_solvers.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_solvers.F90) |
| Reference | Gaussian elimination with partial pivoting in real128 (real64 where the compiler has no real128) | |
| Forward error | $\max\vert x - x_{ref}\vert \le 64\,\varepsilon\,\mathrm{cond}\,\max\vert x_{ref}\vert$, cond from Varah's bound $\Vert A\Vert / \min_i(\vert b_i\vert - \sum_j \vert a_{ij}\vert)$ | |
| Residual | $\max\vert Ax - d\vert \le 64\,\varepsilon\,(\Vert A\Vert \max\vert x\vert + \max\vert d\vert)$, in real128 | |
| Never read | on a bounded line, the coefficients outside it are NaN in the variable variant | |
| Variants | `_const` = `_var` on the same band, bitwise; `_masked` (30 % of the rows the identity) = `_var` on that band, bitwise, and each masked row returns exactly the value given | |
| Device | 4 × 4096 lines in one `!$acc parallel loop`, one line per thread, against the host within the forward bound | |
| Factored | on every constant band, one `_factor` and three right sides solved with it by `_solve`, each within the forward and residual bounds; on the device, one factor copied in once and 4 × 4096 lines solved with it, against the host | |

## Program

```bash
./exe/gideon_test_solvers
```

## Results

<!--@include: ./_tables/solvers-builds.md-->

On the GPU (the OpenACC build):

<!--@include: ./_tables/solvers-device.md-->

## Findings

- The solvers meet the forward and residual bounds on every band and length, in every kind; the three variants agree
  bitwise because they are one template with the coefficient written three ways.
- The periodic pentadiagonal border was first solved by Cramer's rule: a border row replaced by the identity came back
  rounded off its given value (`r1*s22/det` is not `r1`). Elimination with the first border row as pivot returns it
  exactly; the Schur complement of a dominant band is dominant, so no pivoting is needed.
- The factored solves (`_factor` once, `_solve` per line) meet the same bounds; one factor serves every line, on the
  device as on the host. Their first pentadiagonal factor read the previous rows of `fac` back inside its row loop, and
  nvfortran 26.9 at `-O2` vectorised that loop as if its iterations were independent (rows 2–4 unreduced; right at
  `-O0`, with gfortran and with ifx): this test caught it on the nvfortran builds. The previous rows are now carried in
  scalars.
- Every GIDEON band is diagonally dominant or symmetric positive definite (`gideon_fd_d1_p10`'s is the latter only):
  elimination without pivoting is stable for both. The test's bands are dominant; the positive definite one is
  exercised by [A4](./lines).
- On the device the lines are bitwise equal to the host's: the solvers' arithmetic has nothing for the compiler to
  contract differently.
