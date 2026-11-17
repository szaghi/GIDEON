# C1 Rows by their meaning

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C1 | Contracts | C | Python 3 (exact rationals, 100-digit decimals) | yes | `scripts/gen_gideon.py --check` |

## Case and scope

Every coefficient of GIDEON is derived by the generator from what the row means, a datum (a point value or a cell
average) and a target (a derivative at a point, or its average over a cell), never copied from a table. This case is
the generator's own checks, run before any Fortran is written: a row that fails one of them is not emitted, and the
generator exits 1. What it cannot check, the compiled kernels, is the subject of C2 to C6.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Derivation | the moment system $\sum_o w_o L_o(x^j) = T(x^j)$, $j < N$, solved by Gauss–Jordan in `fractions.Fraction` | [`derive`](https://github.com/szaghi/GIDEON/blob/main/scripts/gen_gideon.py) |
| Exact degree, order | the largest $K$ with the row exact on $x^0 \dots x^K$ (checked up to $N + 4$); order $K + 1 - n$ | `exact_degree`, `composite_order` |
| Pair form | a centred row is symmetric or antisymmetric, exactly, or the generator stops | `check_pairs` |
| Mirrors | every left-side row solved independently on its own window, and compared with the right-side row mirrored, times $(-1)^n$; the cell derivatives' orders recomputed on the mirrored operator | `build` |
| ADAM | ADAM's committed tables (`notes/adam-reference/`) parsed into exact fractions and compared with the derived rows, never adopted | `parse_adam_tables`, `cross_check_adam` |
| Orders, free of round-off | $\vert k - T\vert$ on $g(x) = \sin(1.3x + 0.4) + 0.5\cos(0.7x - 0.2)$ at $x_0 = 0.37$, $h = 2^{-11}, 2^{-12}$, every weight exact, in 100-digit decimal arithmetic; bound: observed $\ge$ design $- 0.25$ | `check_orders` |
| Read-back | the emitted Fortran rows parsed back and compared with the derived ones | `write_sources`, `parse_emitted_tables` |

The figure shows four rows as the generator derives them: the weights by offset (exact), the window, the target.

![The windows of the upwind FV kernels with the weights of their rows](/guide/windows-upwind.png)

*The right-side rows and their mirrors: face values (a face at 0) and cell-average first derivatives (cell 0). The
left-side row is the right one read backwards, times $(-1)^n$; GIDEON stores only the right one.*

## Program

```bash
python3 scripts/gen_gideon.py --check                    # derive and verify, write nothing
python3 scripts/gen_gideon.py --check --json rows.json   # the data of these pages (scripts/run_report.sh)
```

## Results

<!--@include: ./_tables/rows-checks.md-->

By family, with the smallest margin of the 100-digit order against the design order:

<!--@include: ./_tables/rows-families.md-->

## Findings

- No mismatch with ADAM's tables: ADAM's coefficients are right where it has them. ADAM's defects (audit items A, B,
  F) are in how its kernels read the rows and their windows, not in the rows ([C5](./adam)).
- The upwind rows sit closest to the bound: at $h = 2^{-12}$ their 100-digit orders are a few thousandths below the
  design order, the pre-asymptotic term of a one-sided stencil. The compiled kernels' orders, on coarser plateaus
  where quadruple precision resolves them, are in [A1](./orders).
- The FV centred cell-average derivatives have the rows of the FD node derivatives (below).

<!--@include: ./_tables/rows-fd-fv.md-->

The reason: for a Fourier mode the cell average multiplies the data and the target by the same factor,
$\sin(\xi/2)/(\xi/2)$, so the moment systems of the two have the same solution. The kernels compute it differently,
by a flux difference of face rows (`fv`) and by pairs (`fd`), so their values agree to round-off, not bitwise in general
([A3](./spectral) draws one curve for both).
