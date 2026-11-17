---
layout: home

hero:
  name: GIDEON
  text: Generated Interpolations and Derivatives, Exact to Order N
  tagline: Finite-difference and finite-volume stencil operators for HPC, in modern Fortran, on every CPU and GPU
  actions:
    - theme: brand
      text: Guide
      link: /guide/
    - theme: alt
      text: Tutorials
      link: /tutorials/
    - theme: alt
      text: Verification
      link: /verification/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/GIDEON

features:
  - title: Exact coefficients, any order
    details: Every row derived in exact rational arithmetic from what it must compute, cross-checked against ADAM's 93 tables, its order measured at 100 digits. Nothing is typed by hand.
    link: /guide/generator
  - title: Every operator with its mirror
    details: 243 kernels — derivatives d1–d6, face values, face derivatives and cell averages, centred or upwind on stencils of 1 to 9 cells, on point values and on cell averages; left-side rows generated as mirrors of the right ones.
    link: /guide/operators
  - title: Called from inside your loops
    details: One pure, device-callable kernel per operator and order, its row a compile-time constant; 0.17–0.82 of ADAM's cost on a GPU and a CPU core. Run-time order dispatchers and the 3-D gradient, divergence, curl and Laplacian on top.
    link: /guide/performance
  - title: Verified, then frozen
    details: Exactness to the leading error constant in every kind, orders observed in quadruple precision, bitwise equivalence with ADAM, golden values, the GPU against the host, 22 of 22 mutants killed.
    link: /verification/
  - title: Learn by compiling
    details: Tutorials whose code is compiled and run by the test suite, the output on the page the program's own.
    link: /tutorials/
  - title: Recipes
    details: Choosing an operator, sizing the ghost cells, replacing ADAM's FDV calls, converting averages and point values, a fast run-time order, other precisions, a new operator.
    link: /cookbook/
---
