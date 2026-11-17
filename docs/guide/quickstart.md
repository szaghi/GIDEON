# Quick Start

```fortran
use gideon, only : I4P, R_P   ! the working kinds, re-exported from PENF
```

The operators are being rebuilt from ADAM's; this page follows them.

## Tests

```bash
fobis fetch
fobis build --mode tests-gnu && bash scripts/run_tests.sh
fobis build --mode tests-gnu-r4 && bash scripts/run_tests.sh --exe-dir exe-r4   # single precision
fobis rule --ex makecoverage-analysis   # what CI runs: coverage into docs/public/coverage.json
```
