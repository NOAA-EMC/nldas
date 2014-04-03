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
! !ROUTINE: noah_dynsetup.F90
!
! !DESCRIPTION:
!  
!  Updates the time dependent NOAH variables
!
! !REVISION HISTORY:
! 15 Apr 2002: Sujay Kumar   Initial Specification
! 
! !INTERFACE:
subroutine noah_dynsetup()
! !USES:
  use lisdrv_module, only: lis,tile 
  use noah_varder
  use spmdMod, only : masterproc, npes
  use noahpardef_module
!EOP
  IMPLICIT NONE
!BOC
  integer :: t, n,ier
#if ( ! defined OPENDAP )
  if ( npes > 1 ) then
     call noah_gather
  endif
  if(masterproc) then 
#endif
     call noah_gfrac()
     call noah_alb()  
#if ( ! defined OPENDAP )
  endif
#if (defined SPMD)
  call MPI_BCAST(noahdrv%noah_gflag,1,MPI_INTEGER,0, &
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(noahdrv%noah_aflag,1,MPI_INTEGER,0, &
       MPI_COMM_WORLD,ier)
  if( npes > 1 .and. ( noahdrv%noah_gflag==1 .or. &
       noahdrv%noah_aflag ==1 ) ) then 
     call noah_scatter
  endif
#endif 
#endif
!EOC  
end subroutine noah_dynsetup

