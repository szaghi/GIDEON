# C3 Mirrors and symmetries

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C3 | Contracts | C | CPU (gfortran R4P, R8P, R16P; nvfortran); the CI matrix | no | `gideon_test_properties` |

## Case and scope

Exactness (C2) checks each kernel on polynomials. This case checks identities between kernels and invariances of
each kernel on random data, bitwise: they hold in floating point because of how the kernels are written, and a
change of the arithmetic (an operand order, a sign, a step) breaks them before it changes an order. The mirror
identity is the rule that keeps ADAM's audit item A from recurring: a left-side kernel reads the right-side row
backwards, and no left row is stored.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Data | 20 trials of uniform random numbers in $[-1/2, 1/2)$ on the widest window, a fixed seed | [`gideon_test_properties.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_properties.F90) |
| Mirror | $k_{ul}(\tilde q) = (-1)^n\,k_{ur}(q)$, $\tilde q_o = q_{1-o}$ (face windows) or $q_{-o}$ (cell windows) | |
| Symmetry | a centred kernel on data of the parity that makes its result odd returns exactly 0: the pair form $c_m(q_m \pm q_{-m})$ (node, cell) or $c_m(q_m \pm q_{1-m})$ (face) cancels term by term | |
| Scaling | $k(8q) = 8\,k(q)$ and $k(q, 1.5) = k(q, 0.75)/2^n$: powers of two are exact in binary | |
| Constant data | a derivative in pair or flux-difference form returns exactly 0; a value returns the constant within 64 ε Σ\|c\| | |
| Registry | which kernel mirrors which (`kmirror`), which return 0 on a constant (`kzero_const`), all generated | [`gideon_test_registry.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_registry.F90) |

## Program

```bash
./exe/gideon_test_properties
```

## Results

<!--@include: ./_tables/properties-counts.md-->

<!--@include: ./_tables/properties-builds.md-->

The dispatch columns are the run-time dispatchers of [C4](./dispatch), run on the same builds.

## Findings

- Every identity holds bitwise in every kind and with both compilers, on the GPU build's host side too: the mirror is
  in the arithmetic, not only in the rows.
- The bitwise mirror needs the left-side kernel to add the same products in the same order as its right twin. The
  generator emits both from one template, so a hand edit of one side breaks C3 at once (the mutants
  `mirror-fv-face-ul3` and `mirror-sign-fv-face-d1-ul2` of [Mutation testing](/regression/mutation)).
- The properties do not test accuracy: a kernel with a wrong but symmetric row passes C3 and fails C2 and A1.
