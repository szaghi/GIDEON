# Replace ADAM's FDV calls

**Problem.** Move a caller of ADAM's `adam_fdv_operators_library` (and of its FNL backend,
`adam_fnl_fdv_operators_library`) to GIDEON: the call that replaces each routine, the window to pass, and what
changes in the results.

The reference is ADAM `develop` at `ff879596`, copied verbatim in `notes/adam-reference/` and compiled into the test
`gideon_test_adam`, which calls every ADAM routine at every half-stencil $s = 1 \dots 5$ next to its GIDEON
replacement. Its result (gfortran 16.2 and nvfortran 26.9, debug and optimised): **279 values bitwise equal, 20 within
64 ε, 8 ADAM defects reproduced, 43 ADAM zero rows** (`notes/adam-fdv-audit.md`, §5;
[Equivalence with ADAM](/verification/adam)).

## The index map

ADAM's kernels declare an assumed-shape window with the node or cell at index 1, `q(1-s:)`; GIDEON's an explicit-shape
window with it at offset 0, `q(-S:S)`. The **caller's section does not change** for the centred operators: both are
called with `u(i-s:i+s)` for node or cell `i`, and the face kernels of both with `u(i+1-s:i+s)` for the face
`i+1/2`. What changes:

1. **The half-stencil becomes an order.** ADAM's `s` is a run-time argument; GIDEON's order is in the kernel's name,
   or the first argument of the dispatcher. The order of a window of half-width $s$ depends on the operator:

   | ADAM routine | GIDEON order |
   |---|---|
   | `compute_derivativeN_fd_centered`, `compute_derivativeN_fv_centered` | $2s$ for N = 1, 2; $2s-2$ for 3, 4; $2s-4$ for 5, 6 |
   | `compute_reconstruction_r_fd_centered`, `compute_reconstruction_r_fv_centered` | $2s$ |
   | `compute_reconstruction_r_fv_centered_dN` | $2s$ for N = 1; $2s-2$ for 2, 3; $2s-4$ for 4, 5 |
   | `compute_avg_to_point`, `compute_point_to_avg` | $2s+2$ |
   | `compute_reconstruction_r_fv_{r,l}upwind`, `compute_derivative1_fv_{r,l}upwind` | $s$ |
   | `compute_derivativeN_fv_{r,l}upwind`, N = 2…5 | $s+1-N$ |
   | `compute_reconstruction_r_fv_{r,l}upwind_dN` | $s-N$ |
   | the 3-D operators | $2s$ |

2. **The window is exact.** An explicit-shape dummy reads exactly `HI − LO + 1` elements: pass the window of the
   kernel, no more, no less. A section of the wrong length is not detected at run time, not even by `-fcheck=all`:
   the kernel reads the elements in order and returns the result of a shifted point. The window of each kernel is the
   caller's responsibility ([Conventions](/guide/conventions#windows)); a checked build catches only a section that
   leaves the caller's array. ADAM's audit item B was caught the other way: `-fcheck=bounds` stopped inside the kernel,
   because an assumed-shape dummy carries the section's length, which an explicit-shape one does not.
3. **Fixed-size buffers.** An ADAM caller that passes a whole buffer declared with `S_MAX` bounds (as
   `compute_reconstruction_r_fd_centered_dev` expects, `q(1-FDV_S_MAX:FDV_S_MAX)`) passes the section of the order
   instead: `buf(1-s:s)`.
4. **The upwind sides are named by the stencil.** ADAM's `rupwind` d1 and face value read the cells on the right; its
   `rupwind` d2–d5 and face derivatives those on the left; `lupwind` the reverse (audit item F). GIDEON's `_ur` is
   always the right side, `_ul` the left.

## Before and after

ADAM's calls (not compiled here; their equivalence is `gideon_test_adam`'s) and their GIDEON replacements, for node or
cell `i`, face `i+1/2`, half-stencil `s`:

| ADAM | GIDEON |
|---|---|
| `compute_derivativeN_fd_centered(s, ds, u(i-s:i+s), d)` | `gideon_fd_dN_c(p, u(i-s:i+s), ds, d)`, p from the table |
| `compute_reconstruction_r_fd_centered(s, u(i+1-s:i+s), f)` and its `_dev` | `gideon_fd_face_c(2*s, u(i+1-s:i+s), f)` |
| `compute_derivativeN_fv_centered(s, ds, a(i-s:i+s), d)` | `gideon_fv_dN_c(p, a(i-s:i+s), ds, d)` |
| `compute_reconstruction_r_fv_centered(s, a(i+1-s:i+s), f)` | `gideon_fv_face_c(2*s, a(i+1-s:i+s), f)` |
| `compute_reconstruction_r_fv_centered_dN(s, a(i+1-s:i+s), f)` | `gideon_fv_face_dN_c(p, a(i+1-s:i+s), 1._R_P, f)`: ADAM's is unscaled ($h^N$ times the derivative); `ds = 1` keeps it, the step gives the derivative |
| `compute_avg_to_point(s, a(i-s:i+s), v)` | `gideon_fv_point_c(2*s+2, a(i-s:i+s), v)` |
| `compute_point_to_avg(s, u(i-s:i+s), v)` | `gideon_fd_average_c(2*s+2, u(i-s:i+s), v)` |
| `compute_reconstruction_r_fv_rupwind(s, a(i:i+s), f)` | `gideon_fv_face_ur(s, a(i+1:i+s), f)`: ADAM's first element is unread |
| `compute_reconstruction_r_fv_lupwind(s, a(i+1-s:i), f)` | `gideon_fv_face_ul(s, a(i+1-s:i), f)` (audit A) |
| `compute_derivative1_fv_rupwind(s, ds, a(i-1:i+s), d)` | `gideon_fv_d1_ur(s, a(i:i+s), ds, d)`: ADAM's first element is unread |
| `compute_derivative1_fv_lupwind(s, ds, a(i-s:i), d)` | `gideon_fv_d1_ul(s, a(i-s:i), ds, d)` (audit A) |
| `compute_derivativeN_fv_rupwind(s, ds, a(i-s:i), d)`, N = 2…5 | `gideon_fv_dN_ul(s+1-N, a(i-s:i), ds, d)`: the left side (F) |
| `compute_derivativeN_fv_lupwind(s, ds, a(i:i+s), d)`, N = 2…5 | `gideon_fv_dN_ur(s+1-N, a(i:i+s), ds, d)`: the right side (F) |
| `compute_reconstruction_r_fv_rupwind_dN(s, a(i+1-s:i), f)` | `gideon_fv_face_dN_ul(s-N, a(i+1-s:i), 1._R_P, f)`: left side, unscaled |
| `compute_reconstruction_r_fv_lupwind_dN(s, a(i+1:i+s), f)` | `gideon_fv_face_dN_ur(s-N, a(i+1:i+s), 1._R_P, f)`: right side, unscaled |
| `compute_derivative6_fv_{r,l}upwind` | `gideon_fv_d6_ur`, `gideon_fv_d6_ul` (ADAM returned 0: audit C) |
| `compute_gradient_{fd,fv}_centered_dev(s, dxyz, qsx, qsy, qsz, g)` | `gideon_{fd,fv}_gradient_c(2*s, qx, qy, qz, dxyz, g)` |
| `compute_divergence_{fd,fv}_centered_dev(s, dxyz, qsx, qsy, qsz, div)` | `gideon_{fd,fv}_divergence_c(2*s, vx_x, vy_y, vz_z, dxyz, div)` |
| `compute_curl_{fd,fv}_centered_dev(s, dxyz, qsx_y, qsx_z, qsy_x, qsy_z, qsz_x, qsz_y, c)` | `gideon_{fd,fv}_curl_c(2*s, vy_x, vz_x, vx_y, vz_y, vx_z, vy_z, dxyz, c)`: the same stencils in the same order |
| `compute_laplacian_fv_centered_dev(s, dxyz, qsx, qsy, qsz, lap)` | `gideon_fv_laplacian_c(2*s, qx, qy, qz, dxyz, lap)` |
| `compute_laplacian_fd_centered_dev(s, dxyz, q, lap)`, a cube | `gideon_fd_laplacian_c(2*s, qx, qy, qz, dxyz, lap)`: gather the three axes |
| `compute_{gradient,divergence,curl,laplacian}_{fd,fv}_centered(s, dxyz, q, r)`, a cube | the same GIDEON operator on the stencils through the point, e.g. `f(i-s:i+s,j,k)`, `f(i,j-s:j+s,k)`, `f(i,j,k-s:k+s)` |
| the tables `FD0_CC(m,s)`, `FD1_CC(m,s)`, `FV1_CC(m,s)` | the rows `FD_FACE_S<s>(m)`, `FD_D1_S<s>(m)`, `FV_FACE_S<s>(m)` of `gideon_tables`: compile-time constants, no `!$acc declare copyin` |

With `s` a compile-time constant, call the kernel of that order instead of the dispatcher
(`gideon_fd_d1_c6` for `s = 3`): the dispatcher costs 9 to 18 % more ([A run-time order, fast](/cookbook/run-time-order)).
On the host, a non-contiguous section (`f(i,j-s:j+s,k)`) passed to an explicit-shape dummy is copied by the compiler;
in device code, gather the stencils into local arrays first, as ADAM's FNL callers already do.

The program below makes each replacement on the data of a cubic, where every operator in it is exact, and prints the
value next to the exact one ($s = 3$, a run-time integer, so the dispatchers are called):

<<< @/examples/src/cookbook_from_adam.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_from_adam.txt-->
```

## What changes numerically

| ADAM | Change after the migration |
|---|---|
| FD and FV centred (face, face derivatives, cell derivatives), deconvolution, the 3-D operators (FNL and cube signatures) | none: **bitwise** (the same rows, the pair form, the step divided inside in ADAM's way) |
| `compute_reconstruction_r_fv_rupwind`, `compute_derivative1_fv_rupwind` | none: bitwise |
| `compute_derivative{2..5}_fv_rupwind`, `compute_reconstruction_r_fv_rupwind_d{1..4}` (the left side) | none: bitwise |
| `compute_derivative{2..5}_fv_lupwind`, `compute_reconstruction_r_fv_lupwind_d{1..4}` (the right side) | round-off, within 64 ε of the scale: ADAM sums the farthest cell first |
| `compute_reconstruction_r_fv_lupwind`, `compute_derivative1_fv_lupwind`, $s \ge 2$ | **fixed** (audit A): ADAM's results were wrong, GIDEON's are exact to their order |
| a zero row: centred d3, d4 (FD and FV cell) and face d2, d3 at $s = 1$; centred d5, d6 and face d4, d5 at $s \le 2$; upwind d2–d5 at $s < N$; upwind d6 at every $s$ (43 calls of the test) | **no kernel**: the call has no GIDEON equivalent, the dispatcher returns `GIDEON_NAN` (audit C) |

The audit items behind the changes (`notes/adam-fdv-audit.md`, measured with `notes/adam-reference/check_upwind_d1.F90`):

- **A.** `FV1_UL` is stored in ascending cell order and read descending, `FV1_UL(m,s)*q(1-m)`: the row is mirrored
  twice. For $s \ge 2$ the left-side face value misses a linear function and the left-side d1 a quadratic one, and
  `compute_derivative1_fv_lupwind` converges at first order for every $s = 2 \dots 5$. GIDEON generates every
  left-side row as the mirror of the right one, and a property test asserts it bitwise.
- **B.** `compute_derivative1_fv_rupwind` declares `q(0:)` over `[0:1+s]` with the cell at `q(1)`; ADAM's own test
  passed `u(i:i+s)`, one element short, so the kernel read one element past the end and returned the derivative of
  cell $i+1$. A caller that passed `u(i:i+s)` meaning cell $i$ gets cell $i$ from `gideon_fv_d1_ur(s, u(i:i+s), ...)`:
  its results change, and are now those it meant.
- **C.** Unsupported orders returned 0 from pure, device-callable code. GIDEON has no kernel for them (a missing
  symbol at compile time); its dispatchers return `GIDEON_NAN`, a quiet NaN, for an order without a kernel. Check the
  order once, before the loop ([A run-time order, fast](/cookbook/run-time-order)).
- **F.** ADAM's `rupwind` and `lupwind` swap sides between d1 and d2–d5. A migrated call names the side its stencil
  lies on: check each `rupwind`/`lupwind` call against the table above, not by its name.
