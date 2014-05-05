!-------------------------------------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
! Released October 2005
!
! See SOFTWARE DISTRIBUTION POLICY for software distribution policies
!
! The LIS source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See COPYRIGHT.TXT for copyright details.
!
!-------------------------------------------------------------------------
module infnan
!-------------------------------------------------------------------------
! Purpose:
!
!	Set parameters for the floating point flags "inf" Infinity
!	and "nan" not-a-number. As well as "bigint" the point
!	at which integers start to overflow. These values are used
!	to initialize arrays with as a way to detect if arrays
!	are being used before being set.
!
! Author: CCM Core group
!-------------------------------------------------------------------------
  use precision
#if (defined DOUBLE_PRECISION)
  real(r8), parameter :: inf = O'777600000000000000000'
  real(r8), parameter :: nan = O'777677777777777777777'
#else

#if (defined SYSLINUX)
  real(r8), parameter :: inf = Z'7F800000'
  real(r8), parameter :: nan = Z'7FC00000'
#else
!  real(r8), parameter :: inf = O'17740000000'
!  real(r8), parameter :: nan = O'17757777777'
  real(r8), parameter :: inf = 0
  real(r8), parameter :: nan = 0
#endif


#endif
!  integer, parameter  :: bigint = 100000000
  integer, parameter  :: bigint = 100
end module infnan

 
