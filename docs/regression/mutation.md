# Mutation testing

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| R2 | Regression | C | CPU (gfortran) | yes | `scripts/mutation_test.py` |

## Case and scope

A test suite is only as good as the defects it catches. The mutation run plants one defect at a time in a scratch
copy of the tree, rebuilds the tests, runs the whole suite (`scripts/run_tests.sh`) and records which test programs
fail. A mutant no test catches is a defect the suite cannot see. The mutants edit the generated Fortran, not the
generator: what is measured is the tests, and the tests never run the generator.

## Setup

Each mutant replaces the first occurrence of a text after an anchor in code (never in a comment: a mutated comment
is an equivalent mutant). They cover one coefficient of each row family (one of them the smallest coefficient of
the widest FD first derivative), the pair signs and offsets of the centred kernels, the division by the step, the
mirror of a left-side kernel and its sign, the order of a flux difference, a dispatcher's case and its NaN, the NaN
constant itself, and the wiring of the 3-D operators; for the FD face kernels of phase 2, a coefficient of each
new row type, the pair sign and step of a face derivative, the mirror read of a one-sided and of a biased kernel, and a
biased dispatcher's case; and for roadmap #2, a face-difference row, a dissipation part's sign and constant, an FV
biased row and its mirror, a filter row and pair sign, a dissipation dispatcher's case, the line solvers' pivot,
periodic column, border, factor reuse and masks, and a compact operator's band, right side, pair sign, step,
$\alpha_f$ row and term and biased mirror.

| Item | Choice | Code |
|---|---|---|
| Copy | the tracked and untracked files of the working tree and the fetched PENF, into `--workdir` (default a new temporary directory); the working tree is never modified | [`mutation_test.py`](https://github.com/szaghi/GIDEON/blob/main/scripts/mutation_test.py) |
| Build | `fobis build --mode tests-gnu --force-compile` per mutant: FoBiS does not rebuild the dependents of a changed cpp `#include`, so an incremental build could keep a stale mutant | `failing_tests` |
| Baseline | the unmutated copy must be green first | |
| Killed | at least one test program fails (`FAIL` or `XPASS` in the report of `run_tests.sh`) | |

```bash
python3 scripts/mutation_test.py --mode tests-gnu [--workdir DIR] [--only NAME ...]
python3 scripts/mutation_test.py > notes/mutation-run-<date>[-<label>].txt   # recorded: the tables pick it up
```

The tables below are those of the newest recorded run that lists exactly the current mutants
(`docs/scripts/make_tables.py`): a run of an older mutant list documents an older suite.

## Results

The recorded run of the current suite and its 53 mutants, `notes/mutation-run-2026-11-18-p2.txt`. The earlier
runs stay recorded: the first release's 22 mutants (`notes/mutation-run-2026-11-17.txt`), the same against the suite
with the spectral test (`notes/mutation-run-2026-11-17-spectral.txt`) and phase 2's 30
(`notes/mutation-run-2026-11-17-phase2.txt`), every mutant killed in each.

<!--@include: ./_tables/mutation-total.md-->

<!--@include: ./_tables/mutation-mutants.md-->

![Which test kills which mutant](/regression/mutation.png)

*Every mutant (rows) and the test programs that fail on it (filled), from the recorded run.*

<!--@include: ./_tables/mutation-tests.md-->

## Findings

- Every mutant is killed, most by four or more programs: the suite is redundant by design, each contract checked
  from more than one side.
- Only `gideon_test_dispatch` kills `dispatch-no-nan` and `dispatch-wrong-kernel-bl`: an order without a kernel, and
  the biased dispatchers (which have no ADAM counterpart to compare with), are exercised nowhere else among the 1-D
  operators. `gideon_test_golden` and `gideon_test_properties` are the only ones that kill `scale-fv-d6-ur2`: the
  sixth-derivative upwind kernels have no ADAM counterpart (audit C), the exactness test runs at $\Delta s = 1$ where a
  wrong power of the step is invisible, and in R8P the kernel has no plateau ([A1](/verification/orders)): a kernel
  that does not converge is reported unresolved, not failed (the mutant's error stays near its target while the
  floor rises past it: one resolved step, no plateau). The mutation run uses the R8P build; the R16P order test,
  rebuilt with this mutant by hand, fails on it (observed order 0, design 2).
- The mutant `dispatch-wrong-kernel` (a d2 kernel in a d1 dispatcher's case) is caught by the dispatch test and by
  the comparison with ADAM, whose run-time-order operators the test calls through the dispatchers.
- One mutant of roadmap #2's first draft was equivalent and was replaced: dropping the guard before a solver's read of
  a coefficient outside the line (`if (i < n-1) u2(i) = f(i)/bp`) changes only factor entries that back substitution
  never uses, so no result can change and no test of results can kill it. The guard keeps a corner of the band from
  being read at all; its replacement, `penta-resolve-pivot`, mutates the factor's reuse for the border columns.
- The P2 mutants are killed by the new programs. `gideon_test_solvers` alone kills the periodic-column and the two
  masked mutants: the compact operators' tests use the constant-band solvers, and the coupling test's masked band is
  tridiagonal with a unit diagonal, on which `merge(1, b, mask)` equals `b` and the pentadiagonal mutant is absent.
  Every compact mutant is killed by `gideon_test_lines` and its golden test, the biased ones also by the coupling test.
- `gideon_vv_advection` kills no mutant on its own: it is a demonstration with frozen physical bounds, and every
  defect it would see is already caught by the contracts.
