# Golden values

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| R1 | Regression | C | CPU (gfortran R4P, R8P, R16P; nvfortran); the CI matrix | no | `gideon_test_golden` |

## Case and scope

The verification cases assert exactness, orders and identities. A change of a few ulp, or an error that grows inside
a bound, passes them. The golden test pins the values themselves: every kernel of the registry on fixed windows,
compared with a file written once by a reference build.

## Setup

| Item | Value |
|---|---|
| Kernels | every kernel of the registry (`apply`) |
| Windows | 8 windows of 33 values, made by an integer linear congruential generator: multiples of $2^{-10}$ in $[-32, 32)$, exact in every kind and with every compiler (no libm, no `random_number`); $\Delta s = 0.75$, a dyadic step |
| Comparison | $\vert r - g\vert \le {\rm TOL}\cdot\sum_k\vert c_k\vert\max\vert q\vert/\Delta s^n$, TOL per kind (below): other compilers and flags contract multiply-adds differently from the writer |
| Reported | the values bitwise equal to the file, the worst difference relative to the magnitude |
| Files | `src/tests/golden/`, chosen by the build's kind: `gideon_golden_v1.txt` (R8P), `gideon_golden_r4_v1.txt` (R4P), `gideon_golden_r16_v1.txt` (R16P); each starts with the compiler and options that wrote it |

<!--@include: ./_tables/golden-count.md-->

<!--@include: ./_tables/golden-files.md-->

## A new version

A golden file is never edited. A deliberate change of results (a new kernel, a changed row) gets a new file, written
by the test itself, and the test's default path moves to it:

```bash
./exe/gideon_test_golden write src/tests/golden/gideon_golden_v2.txt   # refuses to overwrite an existing file
```

The writer opens the file with `status='new'`: an existing file stops it with an error, so v1 cannot be replaced by
accident. The new file is reviewed like code: for the kernels already in v1, its lines must equal v1's character for
character when written by the same build (then v1 can be retired), and the new lines are the values of the new
kernels. Then `GOLDEN` in `gideon_test_golden.F90` points to it, for each kind.

## Program

```bash
./exe/gideon_test_golden         # from the repository root (the path is relative)
```

## Results

<!--@include: ./_tables/golden-builds.md-->

## Findings

- The gfortran builds that wrote the files reproduce them bitwise, in all three kinds.
- The nvfortran OpenACC build (its host code, `-Kieee`) differs from the gfortran writer in the last bits of a share
  of the values, within the tolerance: the reason the comparison is a tolerance and the bitwise count is reported.
  The CI matrix runs the test with gfortran, nvfortran and ifx on every push ([Continuous integration](./ci)).
- The data are exact in binary on purpose, so that the test measures the kernels and not libm: SETH's first golden
  file, built with libm, differed between gfortran `-O0` and `-O3` before any kernel ran.
