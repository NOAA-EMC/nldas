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
! !ROUTINE: clm2_dynsetup.F90
!
! !DESCRIPTION:
!  
!  Updates the time dependent CLM variables
!
! !REVISION HISTORY:
! 15 Apr 2002: Sujay Kumar   Initial Specification
! 
! !INTERFACE:
subroutine clm2_dynsetup()
! !USES:
  use lisdrv_module, only: lis,tile 
  use spmdMod, only : masterproc, npes
  use clm2pardef_module
  use clm2_laitable, only : readlaitable
!EOP
  implicit none

  integer :: ier
  integer :: t, n
!=== End Variable List ===================================================
!BOC
#if ( ! defined OPENDAP )
  if ( masterproc ) then 
#endif
     call clm2lairead
!     call readlaitable
#if ( ! defined OPENDAP )
  endif
#if (defined SPMD)
!  call MPI_BCAST(lis%p%laiflag,1,MPI_INTEGER,0, & 
!       MPI_COMM_WORLD,ier)
!  call MPI_BCAST(lis%p%saiflag,1,MPI_INTEGER,0, & 
!       MPI_COMM_WORLD,ier)
!  if(npes > 1 .and. lis%p%laiflag==1 .and. &
!       lis%p%saiflag ==1 ) then
!<notes>
! For now always do a scatter because clm2lairead temporally interpolates
! lai/sai each time-step
!</notes>
  if ( npes > 1 ) then
!     call clm2_scatter
     call clm2_scatterlai
  endif
#endif 
!EOC  
#endif
end subroutine clm2_dynsetup

