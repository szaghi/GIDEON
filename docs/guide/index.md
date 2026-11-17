# About

GIDEON (**G**enerated **I**nterpolations and **D**erivatives, **E**xact to **O**rder **N**) is the library of the
linear stencil operators of finite-difference and finite-volume schemes: derivatives of any order, interpolations and
extrapolations to faces, on point values and on cell averages, centred and upwind, and the conversions between
averages and point values.

It is the single source of truth for the FDV operators of ADAM (`adam_fdv_operators_library` and its FNL backend),
rebuilt with the conventions of SETH: one `pure`, device-callable kernel per operator and order, coefficients as
compile-time constants, generated in exact rational arithmetic.

The name comes from Gideon's fleece (Judges 6:36–40): a test, then the same test with the condition reversed, a wet
fleece on dry ground and a dry fleece on wet ground. A test and its mirrored control: verification first, and every
upwind operator generated as the mirror of its twin.

## Status

The repository is being seeded: the build, the tests and the CI are in place, and the operators are being rebuilt
from ADAM's (`notes/adam-fdv-audit.md`).
