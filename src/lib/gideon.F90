!< GIDEON, Generated Interpolations and Derivatives, Exact to Order N: the library's entry point.
module gideon
!< GIDEON, Generated Interpolations and Derivatives, Exact to Order N: the library's entry point.
!<
!< Seed of the repository: it re-exports the working kinds, so a client needs `use gideon` only. The finite-difference
!< and finite-volume stencil operators (the single source of truth for ADAM's FDV operators) are added behind this
!< module, each a `pure`, device-callable kernel whose coefficients are compile-time constants generated in exact
!< rational arithmetic.
use penf, only : I4P, I8P, R4P, R8P
use gideon_kind_parameters
implicit none
private
public :: I4P, I8P, R4P, R8P, R_P
endmodule gideon
