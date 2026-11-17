# Conventions

Every GIDEON operator is one `pure` subroutine, `!$acc routine seq` and `!$omp declare target`, that reads a window
of data by **offset** and returns one number (or a vector, for the 3-D operators). The coefficients are compile-time
constants (`gideon_tables`), so nothing is set up and nothing is copied to a device.

## Names

```
gideon_<data>_<op>_<bias><order>
```

| Part | Values | Meaning |
|---|---|---|
| `data` | `fd`, `fv` | point values (finite differences) or cell averages (finite volumes) |
| `op` | `face`, `dN`, `face_dN`, `point`, `average`, `gradient`, `divergence`, `curl`, `laplacian` | the value at a face; the N-th derivative at the node (`fd`) or averaged over the cell (`fv`); the N-th derivative at a face; the point value at a cell centre from averages; the cell average from point values; the 3-D operators |
| `bias` | `c`, `ur`, `ul` | centred; upwind with the stencil on the **right** side; upwind with the stencil on the **left** side |
| `order` | integer | the order of accuracy, measured by the generator (not a stencil half-width) |

`gideon_fd_d1_c6` is the sixth-order centred first derivative at a node from point values; `gideon_fv_face_ul3` the
third-order value at a face from the three cells on its left (the upwind value when the flow goes to the right).

The upwind bias is named by **where the stencil lies**, never by the face it reconstructs. ADAM named its upwind
operators the other way, inconsistently between the first derivative and the others (audit item F in
`notes/adam-fdv-audit.md`).

## Windows

All in units of the grid step; the window is an explicit-shape dummy `q(LO:HI)`, and the caller passes the section
`q(i+LO:i+HI)`:

| Geometry | Data at offset `o` | The result at |
|---|---|---|
| `fd`, node operators | the point `x = o` | node 0 |
| `fd`, face operators | the point `x = o − 1/2` | the face `x = 0`, between nodes 0 and 1 |
| `fv`, face operators | the cell `[o − 1, o]` | the face `x = 0`, between cells 0 and 1 (SETH's face) |
| `fv`, cell operators | the cell `[o − 1/2, o + 1/2]` | cell 0 (the average over it) |

| Operator | Window |
|---|---|
| centred node or cell, half-stencil S | `q(-S:S)` |
| centred face | `q(1-S:S)` |
| upwind face, right side (S cells) | `q(1:S)` |
| upwind face, left side | `q(1-S:0)` |
| upwind cell, right side | `q(0:S)` |
| upwind cell, left side | `q(-S:0)` |

The window bounds are the caller's responsibility. The explicit-shape dummy names them in the interface, but a
section of the wrong length is not an error when its bounds are run-time expressions: the dummy is associated with
the section's elements in order, so a section two elements short (`u(i-3:i+3)` for `gideon_fd_d1_c8`, whose window is
`q(-4:4)`) is read from `u(i-3)` as offset −4 and returns the derivative at node `i+1`, silently, even with gfortran
16.2's `-fcheck=all`. Only a section of constant wrong length fails, at compile time. What a checked build does catch
is a section that leaves the caller's array (`-fcheck=bounds`, `-C`, `-check all`): a halo too thin at the last point
of a loop. (ADAM's audit item B, a read past the end of an assumed-shape window, was caught inside the kernel: an
assumed-shape dummy carries the section's length, an explicit-shape one does not.) Write the section as
`u(i+LO:i+HI)` with the kernel's own `LO` and `HI`, as listed in [Operators](/guide/operators#every-kernel); for a
dispatcher, with the bounds of its window as expressions of the order ([API](/guide/api#run-time-order-dispatchers)).
The ghost cells a loop needs follow from the window: [Ghost cells for a window](/cookbook/ghosts).

The windows drawn, offsets on the axis, the data each kernel reads and the point of its result:

![The windows of gideon_fd_d1_c4, q(-2:2) with node 0 at offset 0, and of gideon_fd_face_c4, q(-1:2) with the points at x = o - 1/2 around the face x = 0, each weight above its offset](/guide/windows-fd.png)

![The windows of gideon_fv_face_c4, q(-1:2) with cell o on [o - 1, o] and the face between cells 0 and 1, and of gideon_fv_d1_c4, q(-2:2) with cell o on [o - 1/2, o + 1/2] and the average over cell 0](/guide/windows-fv.png)

![The windows of the upwind FV operators of three cells: gideon_fv_face_ur3, q(1:3), and its mirror gideon_fv_face_ul3, q(-2:0), at the face x = 0; gideon_fv_d1_ur3, q(0:3), and gideon_fv_d1_ul3, q(-3:0), over cell 0](/guide/windows-upwind.png)

## Arithmetic

- **Pair form** (decision D3): a centred row multiplies pairs, `c_m (q(m) ± q(−m))` for a node or cell and
  `c_m (q(m) ± q(1−m))` for a face, as ADAM does. Fewer flops, the (anti)symmetry exact, and ADAM's values bitwise.
- **The step inside** (decision D4): the kernels divide by the step as ADAM does, by repeated divisions for `fd`
  (`d/ds/ds`) and by one product for `fv` (`d/(ds*ds)`).
- **Flux difference**: an `fv` cell derivative is `(g(+1/2) − g(−1/2))/ds^n`, `g` the face row of the derivative
  below; the face value is the flux a conservative scheme would difference.
- **Mirrors**: a left-side kernel reads the right-side row mirrored, times (−1)ⁿ; no left row is stored, and the
  identity is bitwise (audit item A).

## Run-time order

`gideon_<data>_<op>_<bias>(order, q, [ds,] r)` dispatches to the kernel of that order: one `select case` per call. The
window bounds are expressions of `order` (`q(-(order/2):order/2)` for `gideon_fd_d1_c`), so the whole window goes to
the kernel. An order without a kernel returns `GIDEON_NAN`, a quiet NaN: `error stop` is not device code, and ADAM's
silent 0 is what audit item C fixes.

## Kinds

The kernels compute in PENF's `R_P`, chosen when GIDEON (and PENF, compiled by the same mode) is built:
`-D_R_P_IS_R4P` or nothing (R8P). One precision per build. R16P (`-DPENF_R16P -D_R_P_IS_R16P`, mode `tests-gnu-r16`)
is a host-only verification kind: there the orders of every kernel are observed, which double precision's round-off
hides for the widest stencils. An OpenACC build refuses it. The builds and what each verifies:
[Precision](/guide/precision).
