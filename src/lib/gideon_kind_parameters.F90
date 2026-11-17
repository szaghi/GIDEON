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
use penf, only : I4P, I8P, R_P, R4P, R8P, R16P
implicit none
private
public :: R_P
public :: GIDEON_NAN

#if defined _R_P_IS_R4P
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R4P) !< Stops the build if PENF's R_P is not R4P.
real(R_P), parameter :: GIDEON_NAN = transfer(2143289344_I4P, 1._R_P) !< Quiet NaN (0x7FC00000): the result of a
                                                                      !< run-time dispatch to an order that has no
                                                                      !< kernel (error stop is not device code).
#elif defined _R_P_IS_R16P
! R16P: a host-only verification kind (gfortran, with PENF_R16P): the orders of the widest stencils are measured where
! double precision's round-off hides them (gideon_test_order). Never a device build.
#if defined _OPENACC
#error "GIDEON: R16P is a host-only verification kind, not a device one"
#endif
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R16P) !< Stops the build if PENF's R_P is not R16P.
real(R_P), parameter :: GIDEON_NAN = transfer([0_I8P, 9223231299366420480_I8P], 1._R_P) !< Quiet NaN (0x7FFF8000...0,
                                                                                         !< little endian: the low
                                                                                         !< 64 bits first).
#else
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R8P) !< Stops the build if PENF's R_P is not R8P.
real(R_P), parameter :: GIDEON_NAN = transfer(9221120237041090560_I8P, 1._R_P) !< Quiet NaN (0x7FF8000000000000): the
                                                                               !< result of a run-time dispatch to an
                                                                               !< order that has no kernel (error stop
                                                                               !< is not device code).
#endif
endmodule gideon_kind_parameters
