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
#include "misc.h"
!BOP
! !ROUTINE: check_error
! 
! !DESCRIPTION:
! Error check; Program exits in case of error
! 
! !INTERFACE:
subroutine check_error(ierr,msg,iam)
!EOP
  implicit none
  
  integer :: ierr, iam
  character*40 :: msg
!BOC  
  if ( ierr /= 0 ) then
     print*,'ERR: ',msg,' Stopping.',' (',iam,')'
     call endrun
  endif
!EOC
end subroutine check_error

!BOP
! !ROUTINE: lis_check_error
! 
! !DESCRIPTION:
! Error check; Program exits in case of error
! 
! !INTERFACE:
subroutine lis_check_error(ierr,msg)
!EOP
  implicit none
  
  integer :: ierr
  character(len=*) :: msg
!BOC  
  if ( ierr /= 0 ) then
     call lis_log_msg(msg)
     call endrun
  endif
!EOC
end subroutine lis_check_error

!BOP
! !ROUTINE: check_nc
! 
! !DESCRIPTION:
! Checks status from netcdf calls.
! 
! !INTERFACE:
#if ( defined USE_NETCDF )
subroutine check_nc(status)
!EOP
   use netcdf

   implicit none
  
   integer, intent(in) :: status
!BOC  
   if ( status /= nf90_noerr ) then
      call lis_log_msg(trim(nf90_strerror(status)))
   endif
!EOC
end subroutine check_nc
#endif

