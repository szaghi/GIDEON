!< GIDEON test: the working kinds reach a client through `use gideon` (the PENF dependency is wired).
program gideon_test_kinds
!< GIDEON test: the working kinds reach a client through `use gideon` (the PENF dependency is wired).
use gideon, only : I4P, I8P, R4P, R8P, R_P
implicit none

if (storage_size(1_I4P) /= 32)  error stop 'I4P is not a 32-bit integer'
if (storage_size(1_I8P) /= 64)  error stop 'I8P is not a 64-bit integer'
if (storage_size(1._R4P) /= 32) error stop 'R4P is not a 32-bit real'
if (storage_size(1._R8P) /= 64) error stop 'R8P is not a 64-bit real'
#if defined _R_P_IS_R4P
if (R_P /= R4P) error stop 'R_P is not R4P in a single-precision build'
#else
if (R_P /= R8P) error stop 'R_P is not R8P in a double-precision build'
#endif
print '(A)', 'gideon_test_kinds: ok'
endprogram gideon_test_kinds
