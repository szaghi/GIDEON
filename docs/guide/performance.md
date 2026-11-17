# Performance

Program `src/bench/gideon_bench_adam.F90` (modes `bench-nvf-oac`, `bench-nvf`, `bench-gnu`), raw outputs in
`notes/bench-data/`. Method (SETH's): N = 2^24 points of a 1-D array, one kernel at every point inside one loop
(`!$acc parallel loop` on the GPU, a plain loop on one CPU core pinned with `taskset -c 2`), 2 warm-up passes, the
median of 20 timed passes (`system_clock` around the loop, `!$acc wait`), four runs of the program for the spread.
ADAM's operators are compiled verbatim (ff879596) in the same build.

Hardware: NVIDIA GeForce RTX 4070 (GPU1, driver 610.88), nvfortran 26.9 `-acc -gpu=cc89`; Intel Core i9-13900KF, one
core, gfortran 16.2 `-O3` and nvfortran 26.9 (FoBiS release profiles). The RTX 4070 runs FP64 at 1/64 of its FP32
rate: a ratio measured here is not a claim about a data-centre GPU (PERF5), where these kernels are bound by memory,
not by FP64 throughput.

## Ratio GIDEON / ADAM (median of the four runs' medians; range in brackets)

| Case | RTX 4070 | one core, gfortran | one core, nvfortran |
|---|---|---|---|
| fd d1 c2 | 0.380 (0.375–0.384) | 0.423 (0.423–0.430) | 0.363 (0.355–0.373) |
| fd d1 c4 | 0.465 (0.452–0.482) | 0.428 (0.409–0.442) | 0.361 (0.359–0.363) |
| fd d1 c6 | 0.532 (0.522–0.538) | 0.425 (0.381–0.431) | 0.397 (0.385–0.415) |
| fd d1 c10 | 0.625 (0.617–0.633) | 0.433 (0.405–0.458) | 0.379 (0.374–0.382) |
| fd d2 c6 | 0.605 (0.595–0.615) | 0.475 (0.452–0.479) | 0.401 (0.399–0.408) |
| fv face c6 | 0.344 (0.337–0.357) | 0.436 (0.381–0.438) | 0.417 (0.416–0.421) |
| fv d1 c6 | 0.792 (0.785–0.797) | 0.331 (0.304–0.333) | 0.173 (0.171–0.174) |
| fv d2 c6 | 0.819 (0.807–0.825) | 0.342 (0.314–0.344) | 0.179 (0.173–0.182) |
| fv face ur3 | 0.345 (0.343–0.356) | 0.530 (0.509–0.558) | 0.553 (0.541–0.560) |

Absolute (run 1): on the GPU GIDEON 0.034–0.093 ns per point against ADAM's 0.099–0.128; on a gfortran core 0.94–1.53
against 1.82–4.87.

**PERF1 (GIDEON <= 1.0 x ADAM): met, every case, every device.** The gain is the compile-time row: ADAM's kernels read
a `(S_MAX, S_MAX)` table (a device copy on the GPU) through a run-time `s`, with assumed-shape windows; GIDEON's rows
are folded constants and the loops unrolled. The FV flux differences gain most on the CPU (0.17–0.34): ADAM calls
the face reconstruction twice through the assumed-shape interface.

## The run-time dispatcher: budget missed

`gideon_fd_d1_c(6, ...)` against `gideon_fd_d1_c6`: 1.093 (1.093–1.104) on the GPU, 1.183 (1.134–1.220) with gfortran,
1.115 (1.102–1.136) with nvfortran, against a budget of 1.05 set before coding. The `select case` is evaluated per
call and is not hoisted out of the loop even with a loop-invariant order; on the CPU it also blocks the vectorisation
of the loop. Recorded, not edited. What it means for a caller: an order known at compile time should call the kernel;
a run-time order (PRISM's `fdv_order`) pays 9–22 % over the kernel, and still about 0.45–0.65 of ADAM's run-time-`s`
path (1.093 x 0.532 on the GPU, 1.183 x 0.425 with gfortran). A faster dispatch would hoist the choice out of the
loop: one loop per order in the caller, or a procedure pointer chosen once (not device code). Left to the owner.
