---
layout: home

hero:
  name: GIDEON
  text: Generated Interpolations and Derivatives, Exact to Order N
  tagline: Finite-difference and finite-volume stencil operators for HPC, in modern Fortran
  actions:
    - theme: brand
      text: Guide
      link: /guide/
    - theme: alt
      text: API Reference
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/GIDEON

features:
  - title: Exact coefficients, any order
    details: Every stencil row generated in exact rational arithmetic and checked by its meaning, exactness on polynomials at the stated point, not against another table.
  - title: Every operator with its mirror
    details: Derivatives, interpolations and extrapolations on point values and cell averages, centred and upwind; the left-upwind rows generated as mirrors of the right ones, never stored twice.
  - title: CPU and GPU
    details: One pure, device-callable kernel per operator and order, coefficients as compile-time constants; tested with gfortran, nvfortran and ifx, OpenACC with nvfortran.
---
