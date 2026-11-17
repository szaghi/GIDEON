# Quick Start

```fortran
use gideon
real(R_P) :: u(-2:102), du
integer   :: i
call gideon_fd_d1_c6(u(i-3:i+3), 0.01_R_P, du)        ! d1 at node i, centred, order 6
call gideon_fd_d1_c(6, u(i-3:i+3), 0.01_R_P, du)      ! the same, the order chosen at run time
```

Names, windows and conventions: [Conventions](./conventions); every kernel: [Operators](./operators).

## Tests

```bash
fobis fetch
fobis build --mode tests-gnu && bash scripts/run_tests.sh
fobis build --mode tests-gnu-r4 && bash scripts/run_tests.sh --exe-dir exe-r4   # single precision
fobis rule --ex makecoverage-analysis   # what CI runs: coverage into docs/public/coverage.json
```
