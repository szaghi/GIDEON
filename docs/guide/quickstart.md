# Quick Start

## Build and test

```bash
git clone https://github.com/szaghi/GIDEON && cd GIDEON
fobis fetch                                            # PENF
fobis build --mode tests-gnu && bash scripts/run_tests.sh
fobis build                                            # the library: lib/gnu/static/libgideon.a
```

Every mode, fpm and the installer: [Installation](/guide/installation).

## A first derivative

`use gideon` gives every operator. Each is a `pure` subroutine that takes the window of data around one point as an
explicit-shape section, by offset, and returns one number: here the first derivative at node `i` from the seven
nodes `i-3 … i+3`, with the order fixed in the name (`gideon_fd_d1_c6`) or chosen at run time (`gideon_fd_d1_c`):

<<< @/examples/src/guide_quickstart.F90{fortran}

```bash
gfortran -I lib/gnu/static/mod guide_quickstart.F90 lib/gnu/static/libgideon.a
```

```text
<!--@include: ../examples/output/guide_quickstart.txt-->
```

## Where next

- the names, the windows and the face convention: [Conventions](/guide/conventions);
- every kernel, its window and its order: [Operators](/guide/operators);
- which operator for which data, and which order: [Choose an operator](/cookbook/choose);
- the halo a loop needs: [Ghost cells for a window](/cookbook/ghosts);
- step-by-step programs (FD, FV, upwind, 3-D, GPU): [Tutorials](/tutorials/);
- coming from ADAM: [Replace ADAM's FDV calls](/cookbook/from-adam).
