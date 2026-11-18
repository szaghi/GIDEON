# Cookbook

Recipes: a problem, the code that solves it, what it prints, and where the reference is. Every recipe's code is a
compiled program in `docs/examples/src/` and its output, below the code, is that program's output, checked
(`fobis build --mode examples-gnu`, `scripts/run_examples.sh`).

| Recipe | The question it answers |
|---|---|
| [Choose an operator and its order](/cookbook/choose) | FD or FV data, centred or upwind, the derivative or its face value, which order for which window |
| [Ghost cells for a window](/cookbook/ghosts) | how many ghost cells each operator needs on each side, by order, derived from its window |
| [Replace ADAM's FDV calls](/cookbook/from-adam) | the GIDEON call for each ADAM routine, the index map, what changes numerically |
| [Averages and point values](/cookbook/deconvolution) | when a finite-volume solver must convert averages to point values and back, and at which order |
| [Face values from point values](/cookbook/face-from-points) | a viscous flux at the faces, an upwind-biased face value, the value at a boundary face, from point values |
| [A run-time order, fast](/cookbook/run-time-order) | the dispatcher's cost, and how to pay it once instead of at every point |
| [Single and quadruple precision](/cookbook/precision) | the R4P and R16P builds, what each verifies, `GIDEON_NAN` |
| [Dissipation dialled by a sensor](/cookbook/dissipation) | the face differences, the upwind-central split, a sensor between the central and the biased scheme, Jameson's blend |
| [Filtering a central scheme](/cookbook/filters) | the odd-even mode removed by an explicit or a compact filter, the smooth waves kept to its order |
| [Compact-WENO rows](/cookbook/compact-weno) | compact face values where smooth, given values where a sensor flags: one masked solve |
| [Add an operator](/cookbook/new-operator) | a new family through the generator: derive, cross-check, emit, test, the golden-file rule |

The conventions every recipe relies on (names, windows by offset, the face between offsets 0 and 1) are in
[Conventions](/guide/conventions); every kernel, with its window, is in [Operators](/guide/operators).
