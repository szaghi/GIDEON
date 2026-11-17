# Regression

The verification cases ([Verification](/verification/)) assert what must hold: exactness, orders, identities,
equivalence with ADAM. They leave room: a value may move by a few ulp, or a wrong but still exact row of a kernel
outside a test's reach may pass them. Regression closes that room, and checks the checks:

- **[Golden values](./golden)** pin the values of every kernel on fixed data, one file per kind, written once and
  never overwritten.
- **[Mutation testing](./mutation)** plants one defect at a time in a copy of the tree and records which tests catch
  it: a defect no test catches is a hole in the suite.
- **[Continuous integration](./ci)** runs the suite on every push: gfortran with coverage and the quadruple-precision
  orders, and the matrix of gfortran, nvfortran and ifx, each with the floating-point flags that keep the bitwise
  contracts.

```bash
fobis build --mode tests-gnu && bash scripts/run_tests.sh   # golden values: gideon_test_golden, in the suite
python3 scripts/mutation_test.py --mode tests-gnu           # every mutant: rebuild, run, record
fobis rule --ex makecoverage-analysis                       # what CI runs: coverage, then the R16P orders
```
