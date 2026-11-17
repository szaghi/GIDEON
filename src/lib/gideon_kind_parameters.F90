!< GIDEON, the working real kind of the operators.
module gideon_kind_parameters
!< GIDEON, the working real kind of the operators.
!<
!< Every operator computes in PENF's R_P, chosen when GIDEON (and PENF, compiled from its sources by the same build
!< mode) is built: -D_R_P_IS_R4P, -D_R_P_IS_R8P, or nothing (R8P, PENF's default). One precision per build, no
!< generic interfaces: the operators are compiled once, for that kind. A literal whose meaning depends on the kind
!< lives here, one set per kind behind the same macros, never in an operator (SETH's convention).
!<
!< A compile-time check ties the macro to PENF's R_P: a build that defines it for GIDEON but not for PENF (or the
!< reverse) stops with a division by zero in a constant expression.
use penf, only : R_P, R4P, R8P
implicit none
private
public :: R_P

#if defined _R_P_IS_R4P
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R4P) !< Stops the build if PENF's R_P is not R4P.
#elif defined _R_P_IS_R16P
#error "GIDEON: R16P is not supported (no device code, no caller)"
#else
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R8P) !< Stops the build if PENF's R_P is not R8P.
#endif
endmodule gideon_kind_parameters
