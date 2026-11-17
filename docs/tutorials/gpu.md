# T6 On the GPU

**Goal.** The order-6 Laplacian and gradient of T5 at every node, inside one `!$acc parallel loop`, built with
nvfortran `-acc` for an NVIDIA GPU; the same program built without OpenACC on the host, printing the same lines.

**Concepts.** Every GIDEON operator is `pure`, `!$acc routine seq` and `!$omp declare target`, with its coefficients
compile-time constants ([Conventions](/guide/conventions)): a kernel is called from inside a compute region like any
device routine, and nothing of GIDEON's is copied to the device, only your data. T5 for the gather;
[C6 Host and device](/verification/device) for the verification of every kernel on the device.

## The program

<<< @/examples/src/gpu.F90{fortran}

- **The data region.** `copyin(phi, dxyz)` and `copyout(lap, grad)`: the field and the steps go in, the results come
  out. The rows are parameters folded into the device code.
- **The loop.** `!$acc parallel loop collapse(3) private(qx, qy, qz)`: each iteration gathers its three stencils into
  its own private arrays and calls two kernels, `routine seq` calls inside the parallel loop. The private arrays are
  the gather of T5 moved onto the device: no strided section of `phi` reaches a kernel.
- **The results.** `lap(i,j,k)` and the contiguous column `grad(:,i,j,k)` are passed to the kernels' scalar `lap`
  and explicit-shape `g(3)`.
- **The host reference.** The same loop without directives fills `lap_h` and `grad_h`. Device and host code are
  compiled apart, and a compiler may contract their multiply-adds differently, so the program compares them within
  64 ulps of the operator's scale, not bitwise. The output does not say where the loop ran: built either way, the
  program prints the same lines.

## Build and run

On the host (gfortran ignores the `!$acc` directives without `-fopenacc`):

```bash
fobis build --mode examples-gnu && ./exe-examples/gpu
```

On the GPU the library must be compiled with `-acc` too: `!$acc routine seq` gives each kernel a device version only
in an OpenACC compile. Against `lib/nvf/static` (built without `-acc`) nvfortran stops with
`NVFORTRAN-S-1061-Procedures called in a compute region must have acc routine information - gideon_fd_laplacian_c6`.
FoBiS builds the library and the program in one go; `-f /dev/null` sets the project's `fobos` aside, whose modes fix
their compiler and flags:

```bash
fobis build -f /dev/null --compiler nvfortran \
            --cflags "-c -acc -gpu=cc89 -Minfo=accel" --lflags "-acc -gpu=cc89" \
            --src src/lib/ --src src/third_party/PENF/src/lib/ --src docs/examples/src/ \
            --target docs/examples/src/gpu.F90 --build-dir exe-examples/acc/
CUDA_VISIBLE_DEVICES=1 NVCOMPILER_ACC_NOTIFY=1 ./exe-examples/acc/gpu
```

`cc89` is the compute capability of an RTX 40 series (Ada); `cc80` for an A100, `cc90` for an H100.
`CUDA_VISIBLE_DEVICES` picks the GPU, `NVCOMPILER_ACC_NOTIFY=1` reports each kernel launch on the standard error.

## Output

<<< @/examples/output/gpu.txt{text}

This is the host build's output, the one `scripts/run_examples.sh` checks. Built as above with nvfortran 26.9 and run
on an RTX 4070, the program printed the same lines, compared with `diff`, and the launch report on the standard error
was (paths shortened):

```text
launch CUDA kernel  file=.../docs/examples/src/gpu.F90 function=gpu line=30 device=0 threadid=1 num_gangs=256 num_workers=1 vector_length=128 grid=256 block=128
launch CUDA kernel  file=.../docs/examples/src/gpu.F90 function=gpu line=30 device=0 threadid=1 num_gangs=2048 num_workers=1 vector_length=128 grid=2048 block=128
```

One launch per grid, $32^3$ and $64^3$ iterations in gangs of 128. `-Minfo=accel` reports the loop as
`!$acc loop gang, vector(128) collapse(3)` and `Local memory used for qx,qz,qy`: the private stencils.

## What to notice

- **The kernels need nothing on the device but your data**: no table to copy, no setup call. The compiler calls the
  `routine seq` versions, whose coefficients are constants folded into the device code.
- **Device and host agree.** In this build they agree bit for bit, every value of both grids (measured with
  nvfortran 26.9, default flags and `-fast`). The program's tolerance is of the kind [C6](/verification/device)
  verifies for every kernel, a few ulps of the operator's scale; the bitwise equality here is an observation of one
  compiler on one program, not a property to rely on.
- **The accuracy is T5's**: order 5.99, the same errors to the digits printed.

## Exercises

1. The order at run time on the device: `gideon_fd_laplacian_c(order, qx(-s:s), qy(-s:s), qz(-s:s), dxyz, lap(i,j,k))`
   with `order` a variable. The `select case` takes the same branch in every iteration: no divergence, and 1.093
   times the kernel's time on the RTX 4070 ([Performance](/guide/performance#the-run-time-dispatcher-budget-missed)).
2. Time the parallel loop with `system_clock` (not `cpu_time`, which counts the host's processor time). Add `async`
   to the loop and the timing measures the launch only, unless `!$acc wait` precedes the second reading.
3. Count the values where `lap == lap_h` with `-fast` and with `-gpu=cc89,nofma -Mnofma`.
4. OpenMP offload: the kernels also carry `!$omp declare target`. Write the loop as `!$omp target teams distribute
   parallel do collapse(3) private(qx, qy, qz)` with a `map` for the arrays, and compile the library and the program
   with `-mp=gpu`. GIDEON's CI does not build this path.
