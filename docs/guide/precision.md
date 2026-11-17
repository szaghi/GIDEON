# Precision

GIDEON computes in one real kind per build, PENF's `R_P`: R8P (double precision) by default, R4P (single precision)
with the `*-r4` modes, R16P (quadruple precision) in one test mode, for verification only. This page says what the
kind changes, what is verified in each, and where each stops; the build commands and the test flags are the recipe
[Single and quadruple precision](/cookbook/precision).

## One kind, chosen at build time

`gideon_kind_parameters` takes `R_P` from PENF, which PENF sets from a macro:

| Kind | Macros (GIDEON and PENF alike) | Modes |
|---|---|---|
| R8P | none (or `-D_R_P_IS_R8P`) | every mode without a suffix |
| R4P | `-D_R_P_IS_R4P` | `*-r4`: the libraries `libgideon_r4.a`, `libgideon_r4.so`, the tests in `exe-r4/` |
| R16P | `-DPENF_R16P -D_R_P_IS_R16P` | `tests-gnu-r16`, the tests in `exe-r16/` |

- **One kind per build.** No generic interfaces: every kernel is compiled once, for `R_P`. The R4P and R8P libraries
  export the same module and procedure names, so a program links one.
- **A compile-time check.** `KIND_CHECK = 1/merge(1, 0, R_P == R4P)` (and its R8P and R16P twins) divides by zero in a
  constant expression when the macro reached GIDEON but not PENF, or the reverse: the build stops.
- **R16P is host only.** An OpenACC compilation of `gideon_kind_parameters` in R16P stops on `#error`; the R16P
  modes are built with gfortran (`PENF_R16P` enables PENF's quadruple kind).
- **The rows are the same rationals.** Every coefficient is a literal `p._R_P/q._R_P` folded by the compiler in the
  build's kind. In R8P every numerator and denominator is exact, so each coefficient is the correctly rounded rational.
  In R4P, 15 of the 466 literals have a numerator or denominator above $2^{24}$, rounded before the division; 2
  coefficients come out one unit in the last place from the correctly rounded value. The exactness test's bound,
  64 ε of the scale, covers them.

## `GIDEON_NAN`

The one kind-dependent literal of the library: a quiet NaN of `R_P`, defined bit by bit with `transfer` (no
arithmetic, so no floating-point exception at compile or run time):

| Kind | Bits |
|---|---|
| R4P | `0x7FC00000` |
| R8P | `0x7FF8000000000000` |
| R16P | `0x7FFF8000…0` (two 64-bit words, the low one first: little endian) |

The run-time order dispatchers return it for an order without a kernel. `error stop` is not device code, and a
silent number, ADAM's 0 for an unsupported order (audit item C), is a wrong result that looks right. Test it with
`ieee_is_nan`.

## What each kind verifies

The whole suite runs in each kind, except the equivalence with ADAM, whose code is R8P only. The checks are
bitwise (the mirrors, the symmetries, the dispatchers against their kernels, the host against the device), in units
of the kind's ε (64 ε of the scale for exactness and the leading error), or per kind (the golden files: $10^{-13}$,
$10^{-5}$, $10^{-30}$ relative to the magnitude), so the same test asserts the same property in each.

The orders are the exception: an order is observed only where truncation dominates round-off over three successive
halvings. In R8P that holds for 119 of the 288 kernels, in R4P for 17; the widest one-sided stencils (7 to 9 cells)
leave their pre-asymptotic range only where R8P's round-off already dominates, whatever the function, since
truncation over round-off depends on $kh$ alone. R16P resolves all 288: that is where the order of every compiled
kernel is observed ([Orders of accuracy](/verification/orders)), and CI runs it. The generator measures the rows
themselves at 100 digits, before any kind is involved.

## Where R4P stops

- **Round-off floor.** A kernel's error cannot fall below about $\varepsilon \sum_k |c_k| \max|q| / h^n$, with
  $\varepsilon = 1.19 \cdot 10^{-7}$. The eighth-order first derivative of $\sin x$ meets it at $h = 2^{-3}$ in R4P,
  at $h = 2^{-5}$ in R8P ([the recipe's program](/cookbook/precision)); a derivative of order $n$ meets it earlier as
  $n$ grows.
- **Range.** The step enters as $h^{-n}$: the tests stop refining where $h^{-(n+1)}$ approaches R4P's range
  ($2^{127}$), and a sixth derivative at a step of $10^{-7}$ already loses the result to the floor above.
- **Speed.** Not measured for GIDEON: the benchmark ([Performance](/guide/performance)) runs in R8P only.
