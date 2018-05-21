module precision
!-------------------------------------------------------------------------------
!
! Purpose:
!	Define the precision to use for floating point and integer operations
!	throughout the model.
!
! Author: CCM Core group
!
!-------------------------------------------------------------------------------
  integer, parameter :: r4 = selected_real_kind(5)
! integer, parameter :: r8 = selected_real_kind(12)
  integer, parameter :: r8 = selected_real_kind(5)
  integer, parameter :: i8 = selected_int_kind(5)
!#if (defined DOUBLE_PRECISION)
!    integer, parameter :: r8 = selected_real_kind(12)
!#else
!    integer, parameter :: r8 = selected_real_kind(5)
!#endif
end module precision
