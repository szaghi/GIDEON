# C6 Host and device

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C6 | Contracts | C | GPU (nvfortran OpenACC, `tests-nvf-oac`); CPU builds check the device path on the host | no | `gideon_test_device` |

## Case and scope

Every kernel is `!$acc routine seq` and `!$omp declare target`, with compile-time coefficients: nothing is set up and
nothing is copied to the device. This case runs every kernel inside an `!$acc parallel loop` on the GPU and compares
each value with the host's. The device and the host may contract multiply-adds differently (FMA), so the comparison
is a tolerance, and the number of bitwise-equal values is reported.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Data | 65536 windows of random numbers in $[-1, 1)$, a fixed seed, $\Delta s = 0.37$ | [`gideon_test_device.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_device.F90) |
| Device call | `apply_dev`, the registry's device-callable dispatcher (one `select case` per call, as a caller with a run-time choice would do), in one `!$acc parallel loop` per kernel | [`gideon_test_registry.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_registry.F90) |
| Host call | `apply`, the same kernels on the host | |
| Asserted | $\vert r_{dev} - r_{host}\vert \le 64\,\varepsilon \sum_k\vert c_k\vert \max\vert q\vert/\Delta s^n$ | |
| Without OpenACC | the loop runs on the host: `apply_dev` against `apply`, which must agree bitwise | |

## Program

```bash
fobis build --mode tests-nvf-oac --build-dir exe/oac/     # GPU of compute capability $GPU_CC (fobos varsets)
CUDA_VISIBLE_DEVICES=1 ./exe/oac/gideon_test_device
scripts/run_report.sh --gpu                               # the same, into the report
```

## Results

<!--@include: ./_tables/device-builds.md-->

The OpenACC run:

<!--@include: ./_tables/device-hardware.md-->

## Findings

- Every kernel computes on the device within the contraction tolerance; no failure.
- On the GPU a share of the values (in the table) differs from the host's in the last bits, by a few ε of the scale
  at most: the device code and the host code contract multiply-adds differently. On the host, the three kinds of the
  gfortran build agree bitwise with the host path of the same kernels.
- R16P is host only: an OpenACC build of `gideon_kind_parameters` in quadruple precision stops on `#error`
  ([Precision](/guide/precision)).
