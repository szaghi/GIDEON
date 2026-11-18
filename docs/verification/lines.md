# A4 Compact operators

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A4 | Accuracy | A | Python (100 digits); CPU (gfortran R4P, R8P, R16P; nvfortran; ifx); GPU (OpenACC) | no | `scripts/gen_gideon.py`, `gideon_test_lines` |

## Case and scope

The 29 compact operators: their rows by their meaning, their orders on periodic lines, their symbols, their right
sides and solves on the device (roadmap #2, gates L1, L2, L4, L5). A compact operator is a right-side kernel and a
band; the test drives both through a generated registry whose band constants are the library's own.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| L1 rows | every row exact on $x^j$, $j \le K$ = order + n − 1: $\vert\sum_o \alpha_o T(x^j)(o) - r\vert \le 64\,\varepsilon(\sum\vert w\vert\max\vert q\vert + \sum\vert\alpha\vert\max\vert T\vert)$; on $x^{K+1}$ the generator's exact residual (filters at $\alpha_f$ = 0, 0.3, 0.45) | [`gideon_test_lines.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_lines.F90) |
| L2 orders | $g = \sin(x + 0.4) + 0.5\cos(2x - 0.2) + 0.3\sin(3x + 1.1)$ on 8 … 65536 points of $[0, 2\pi)$, solved; a plateau of three resolved sizes, both observed orders $\ge$ design − 0.25; every operator measured in R16P (required) | |
| Round-off floor | $\varepsilon \sum\vert w\vert \max\vert q\vert / h^n \cdot F$, $F = (1 + \sum\vert\alpha\vert/s)/s$, $s$ the band's smallest symbol: the right side's round-off and the solve's backward error, through the inverse band | |
| L4 symbols | on 64 points, the solve of each mode = the kernel's symbol over the band's; centred: real for an even $n$, imaginary for an odd one; biased: the flux difference dissipative on its side; filters: a real transfer function in $[0, 1]$, the odd-even mode removed | |
| Mirrors | `pbl`'s right side on mirrored data = `pbr`'s, bitwise | |
| L5 device | 16384 right sides and 512 lines of 128 points per operator on the device against the host | |
| Registry | [`gideon_test_lines_registry.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_lines_registry.F90), generated: the bands are `gideon_tables`' constants | |

## Program

```bash
./exe/gideon_test_lines
./exe-r16/gideon_test_lines   # every order measured
```

## Results

The operators, from the generator (the bands and the residuals are exact rationals):

<!--@include: ./_tables/lines-operators.md-->

<!--@include: ./_tables/lines-generator.md-->

The compiled operators:

<!--@include: ./_tables/lines-builds.md-->

On the GPU:

<!--@include: ./_tables/lines-device.md-->

## Findings

- Every compact operator reaches its design order on a periodic line in quadruple precision; double precision
  resolves 10 of the 29 (the floor $F$ is up to 200 for a filter at $\alpha_f = 0.45$), single precision none.
- Each centred operator is non-dissipative to round-off, each biased one dissipative on its upwind side at every mode,
  each filter's transfer function lies in $[0, 1]$ and kills the odd-even mode.
- The first version of this test had the floor $1/s$: it omitted the solve's backward error, and the symbol check
  failed by up to 1.3 times on the fourth-order operators before any change of the operators. The derivation was
  corrected (a missing term, not a looser constant), stated in the program's header.
- The first version of the registry carried the bands as literals: a wrong library band passed. The registry now reads
  the library's constants, and a planted wrong band fails L1 and L2.
- On the device the right sides are bitwise equal to the host's except for the filters, whose $(R + \alpha A)\,q$ the
  device contracts to multiply-adds: within the bound, $2 \cdot 10^{-17}$ of its scale.
