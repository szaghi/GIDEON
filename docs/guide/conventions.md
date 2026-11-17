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

A wrong window is a bounds error at the call site in a checked build (`-fcheck=bounds`), never a silent shift (audit
item B).

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
hides for the widest stencils. An OpenACC build refuses it.
