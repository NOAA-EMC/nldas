!===============================================================================
! CVS: $Id: shr_vmath_mod.F90,v 1.1.1.1 2003/02/06 16:10:56 jgottsch Exp $
! CVS: $Source: /YUKON/LDASCVS/LDAS/SRC/CLM2SRC/csm_share/shr_vmath_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
! PURPOSE: 
!   provides a uniform, platform-independent API for vector math functions
!===============================================================================

module shr_vmath_mod

   !----------------------------------------------------------------------------
   ! routines that evaluate various math functions for vector arguments
   ! intended to provide platform independent access to vendor optimized code
   !----------------------------------------------------------------------------

   use shr_kind_mod

   implicit none

   private
   public :: shr_vmath_sqrt, &
      shr_vmath_exp, shr_vmath_log, &
      shr_vmath_sin, shr_vmath_cos

   contains

!===============================================================================

subroutine shr_vmath_sqrt(X, Y, n)

   !----- arguments ---
   integer(shr_kind_in),intent(in)  ::   n  ! vector length
   real   (shr_kind_r8),intent(in)  :: X(n) ! input vector argument
   real   (shr_kind_r8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: sqrt for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = sqrt(X)
#else

#if (defined AIX)
   call vsqrt(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vsqrt(X, Y, n)
#endif

#if (defined OSF1)
   call vsqrt(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = sqrt(X)
#endif
#endif

end subroutine shr_vmath_sqrt

!===============================================================================

subroutine shr_vmath_exp(X, Y, n)

   !----- arguments ---
   integer(shr_kind_in),intent(in)  ::   n  ! vector length
   real   (shr_kind_r8),intent(in)  :: X(n) ! input vector argument
   real   (shr_kind_r8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: exp for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = exp(X)
#else

#if (defined AIX)
   call vexp(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vexp(X, Y, n)
#endif

#if (defined OSF1)
   call vexp(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = exp(X)
#endif
#endif

end subroutine shr_vmath_exp

!===============================================================================

subroutine shr_vmath_log(X, Y, n)

   !----- arguments ---
   integer(shr_kind_in),intent(in)  ::   n  ! vector length
   real   (shr_kind_r8),intent(in)  :: X(n) ! input vector argument
   real   (shr_kind_r8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: log for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = log(X)
#else

#if (defined AIX)
   call vlog(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vlog(X, Y, n)
#endif

#if (defined OSF1)
   call vlog(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = log(X)
#endif
#endif

end subroutine shr_vmath_log

!===============================================================================

subroutine shr_vmath_sin(X, Y, n)

   !----- arguments ---
   integer(shr_kind_in),intent(in)  ::   n  ! vector length
   real   (shr_kind_r8),intent(in)  :: X(n) ! input vector argument
   real   (shr_kind_r8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: sin for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = sin(X)
#else

#if (defined AIX)
   call vsin(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vsin(X, Y, n)
#endif

#if (defined OSF1)
   call vsin(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = sin(X)
#endif
#endif

end subroutine shr_vmath_sin

!===============================================================================

subroutine shr_vmath_cos(X, Y, n)

   !----- arguments ---
   integer(shr_kind_in),intent(in)  ::   n  ! vector length
   real   (shr_kind_r8),intent(in)  :: X(n) ! input vector argument
   real   (shr_kind_r8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: cos for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = cos(X)
#else

#if (defined AIX)
   call vcos(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vcos(X, Y, n)
#endif

#if (defined OSF1)
   call vcos(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = cos(X)
#endif
#endif

end subroutine shr_vmath_cos

!===============================================================================

end module shr_vmath_mod
