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

![Time of GIDEON's kernel over ADAM's for nine cases on the RTX 4070, one core with gfortran and one core with nvfortran, all below 1; right: the run-time dispatcher fd_d1_c(6) over the kernel fd_d1_c6, 1.09 to 1.18](/guide/performance.png)

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
The first, as code: [A run-time order, fast](/cookbook/run-time-order).

## Compact operators: the band factored once, budget met

Program `src/bench/gideon_bench_compact.F90` (modes `bench-compact-nvf-oac`, `bench-compact-nvf`, `bench-compact-gnu`),
raw outputs `notes/bench-data/compact-*`. 2^24 points as 65536 periodic lines of 256, a line contiguous (`q(i,l)`, the
layout the line solvers take). The explicit operator and a compact right side: one kernel call per point; the solve:
one line per call, one line per thread on the GPU. Same method and hardware as above. Budget, set before the program
was written (roadmap #2, L6): `gideon_fd_d1_p6` (right side + `gideon_tri_periodic_const`) <= 4 x `gideon_fd_d1_c6` on
one core, <= 8 x on the GPU.

| Ratio to the explicit operator of equal order (median of four; range) | one core, gfortran | one core, nvfortran | RTX 4070 |
|---|---|---|---|
| `fd_d1_p6`, right side only | 0.93 (0.90–0.94) | 0.95 (0.87–0.96) | 0.88 (0.86–0.90) |
| `fd_d1_p6`, solve only | 5.73 (5.55–5.77) | 6.25 (4.70–6.48) | 11.59 (10.76–11.92) |
| **`fd_d1_p6`, both** (budget 4 / 8) | **6.71** (6.64–6.76) | **7.24** (5.64–7.42) | **12.54** (11.45–12.86) |
| `fd_d1_p8`, both (pentadiagonal; against `c8`) | 10.89 (10.85–10.93) | 9.76 (6.32–11.94) | 31.81 (31.69–32.21) |
| `fd_filter_p8`, both (against `fd_filter_c8`) | 6.72 (6.70–6.73) | 7.56 (7.20–7.59) | 18.17 (17.82–18.80) |

Absolute (gfortran, run 1): `c6` 1.70 ns per point, `p6` 11.5 (its right side 1.5, its solve 9.4); on the GPU `c6` 0.070,
`p6` 0.80.

The right side of a compact operator costs what an explicit kernel costs; the excess is the solve. On a core, each row
of the periodic solve carries dependent divisions (two in the factor, one more in the resolve of the border column)
over three sweeps of the line; on the GPU, one thread per line reads addresses 256 values apart, so the loads of a warp
are not coalesced. Per unit of accuracy the comparison is different (`p6` is about 15 times more accurate than `c6` on
the same grid: [T7](/tutorials/compact-derivative)), but the budget was set per point, and with `_const` it is missed.

For a constant band the factor depends only on the band and the line length, so it is computed once
(`gideon_{tri,penta}_periodic_factor`) and every line solved with it by multiplications alone
(`gideon_{tri,penta}_periodic_solve`; roadmap #2, decision 4a). Measured in the same program, the same runs, cases
9–12 (raw outputs `notes/bench-data/compact2-*`, where cases 1–8 were measured again):

| Ratio to the explicit operator of equal order (median of four; range) | one core, gfortran | one core, nvfortran | RTX 4070 |
|---|---|---|---|
| `fd_d1_p6`, factored solve only | 1.96 (1.93–1.99) | 1.70 (1.47–1.86) | 3.95 (3.75–4.07) |
| **`fd_d1_p6`, both, factored** (budget 4 / 8) | **2.96** (2.91–3.05) | **2.71** (2.38–2.92) | **4.81** (4.54–4.95) |
| `fd_d1_p8`, both, factored (against `c8`) | 2.81 (2.77–2.83) | 2.87 (2.10–3.17) | 3.36 (3.32–3.40) |
| `fd_filter_p8`, both, factored (against `fd_filter_c8`) | 2.96 (2.75–3.13) | 2.75 (2.39–3.07) | 7.39 (7.32–7.42) |

**Budget met, on both**, with the factored solve. On the GPU the factor also removes the per-line scratch arrays:
one factor, read by every thread, stays in the cache, and the pentadiagonal `p8` drops from 32 to 3.4 times its
explicit twin. The loads of the line itself are still uncoalesced (one line per thread, a line contiguous); a layout
with the line index fastest is the next measurement, with no budget pending. Use `_const` for a single line or a band
that changes from line to line, `_factor` + `_solve` whenever lines share the band and the length.
