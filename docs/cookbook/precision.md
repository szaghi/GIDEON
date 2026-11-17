# Single and quadruple precision

**Problem.** Build GIDEON in single precision (R4P) for a solver that runs in it, or in quadruple precision (R16P) to
see what double precision's round-off hides; know what each build verifies and which compiler flags the bitwise
checks need.

## One kind per build

The kernels compute in PENF's `R_P`, fixed when the library is compiled ([Precision](/guide/precision)):

| Kind | Macros | Library modes | Test modes | Build directories |
|---|---|---|---|---|
| R8P (default) | none | `static-{gnu,nvf,intel}`, `shared-*` | `tests-{gnu,nvf,intel}`, `-debug`, `tests-gnu-coverage`, `tests-nvf-oac` | `lib/<compiler>/{static,shared}/`, `exe/` |
| R4P | `-D_R_P_IS_R4P` | `static-{gnu,nvf,intel}-r4`, `shared-*-r4`: `libgideon_r4.a`, `.so` | `tests-gnu-r4`, `tests-gnu-debug-r4`, `tests-nvf-r4`, `tests-nvf-oac-r4`, `tests-intel-r4` | `lib/<compiler>/{static,shared}-r4/`, `exe-r4/` |
| R16P | `-DPENF_R16P -D_R_P_IS_R16P` | none: a verification kind | `tests-gnu-r16` (gfortran, host only) | `exe-r16/` |

```bash
fobis build --mode static-gnu-r4                                         # lib/gnu/static-r4/libgideon_r4.a
fobis build --mode tests-gnu-r4  && bash scripts/run_tests.sh --exe-dir exe-r4
fobis build --mode tests-gnu-r16 && bash scripts/run_tests.sh --exe-dir exe-r16
```

The macros must reach PENF and GIDEON alike (the modes compile PENF from its sources with the same flags);
`gideon_kind_parameters` stops a build where they disagree, with a division by zero in a constant expression. There
are no generic interfaces: the R4P and R8P libraries export the same names, so a program links one of them. An
OpenACC build in R16P stops at compile time (`#error`): R16P is not a device kind.

The program below prints the kind it was compiled against, checks `GIDEON_NAN`, and follows the error of
`gideon_fd_d1_c8` down to the round-off floor:

<<< @/examples/src/cookbook_precision.F90{fortran}

With the R8P library (`lib/gnu/static`):

```text
<!--@include: ../examples/output/cookbook_precision.txt-->
```

The same source, compiled with `-I lib/gnu/static-r4/mod` and linked with `lib/gnu/static-r4/libgideon_r4.a`:

```text
R_P: 32 bits, 6 digits, epsilon  1.19E-07
GIDEON_NAN is a NaN: T
gideon_fd_d1_c8 on sin(x) at x = 0.5: error, and its round-off scale eps sum|c| max|q| / h
  h = 2^-1    5.13E-06  4.95E-07
  h = 2^-3    1.19E-07  1.67E-06
  h = 2^-5    4.17E-07  4.65E-06
  h = 2^-7    5.96E-07  1.61E-05
  h = 2^-9    4.17E-06  6.18E-05
  h = 2^-11   4.15E-05  2.45E-04
```

In R8P the eighth-order error falls by $2^{16}$ per two halvings until it meets round-off at $h = 2^{-5}$; in R4P it
meets it at the second step and then grows as $h$ shrinks, inside its round-off scale. This is why the order of a kernel cannot be observed in every
kind (next section).

## What each kind verifies

| Test | R8P | R4P | R16P |
|---|---|---|---|
| exactness on monomials and the leading error (`gideon_test_tables`), 64 ε | yes | yes | yes |
| orders observed on a smooth function (`gideon_test_order`) | 96 of 243 kernels | 11 of 243 | **all 243** |
| mirrors and symmetries, dispatch, 3-D operators, host = device | yes | yes | yes (host) |
| golden values (`src/tests/golden/`) | `gideon_golden_v1.txt`, tolerance $10^{-13}$ | `gideon_golden_r4_v1.txt`, $10^{-5}$ | `gideon_golden_r16_v1.txt`, $10^{-30}$ |
| equivalence with ADAM (`gideon_test_adam`) | yes | no: ADAM's code is R8P | no |

The order test counts a kernel as measured only on a plateau of three resolved steps above ten round-off floors; a
kernel without one is "unresolved in this precision", not failed. The widest upwind stencils leave their
pre-asymptotic range only where R8P's round-off already dominates, whatever the function, so R16P is where every
compiled kernel's order is observed; the generator measures the rows themselves in 100-digit decimal arithmetic
([The generator](/guide/generator)). CI runs the R16P build in the coverage job (`fobis rule --ex
makecoverage-analysis`).

## `GIDEON_NAN`

A quiet NaN of the build's kind, defined bit by bit (`transfer` of an integer constant) in `gideon_kind_parameters`:
`0x7FC00000` in R4P, `0x7FF8000000000000` in R8P, `0x7FFF8000…0` in R16P (the low 64 bits first, little endian). It is
the result of a run-time dispatcher called with an order that has no kernel: `error stop` is not device code, and a
silent 0, ADAM's answer (audit item C), is a wrong number. Test it with `ieee_is_nan`
([A run-time order, fast](/cookbook/run-time-order)).

## The flags of the bitwise checks

The tests compare values bitwise (dispatcher against kernel, host against device, GIDEON against ADAM, golden values)
and rely on the operations being executed as written. Two compilers reorder them by default:

- **ifx: `-fp-model=precise`** on every Intel mode, library and tests (`$FP_INT` in the `fobos`). ifx's default model
  (`fast=1`) reassociates inlined arithmetic: the pair form's exact symmetries break, and 7 curl values differ from
  ADAM's by 1 ulp.
- **nvfortran: `-Kieee`** on the test modes (`$IEEE_NVF`). The release profile's `-fast` vectorises for the CPU it
  runs on; on another CPU (the CI runner) it reorders the longer sums and the bitwise equality with ADAM breaks at
  $s = 4, 5$. The nvfortran libraries and benchmarks keep `-fast` alone: a library built that way computes the same rows
  to round-off, not to the bit.

**Reference.** The kind module: [API](/guide/api#kinds); the tests: [Orders of accuracy](/verification/orders),
[Golden values](/regression/golden).
