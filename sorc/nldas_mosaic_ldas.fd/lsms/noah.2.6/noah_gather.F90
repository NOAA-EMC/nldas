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
! !ROUTINE: noah_gather.F90
!
! !DESCRIPTION:
!  Gathers noah tiles
!
! !REVISION HISTORY:
! 
!  Apr 2003 ; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine noah_gather()
! !USES:
  use tile_spmdMod
  use noah_varder
  use noahpardef_module
!EOP
  IMPLICIT NONE

  integer ierr
!BOC
#if  (defined SPMD)
  call MPI_GATHERV(noah(1:di_array(iam)),di_array(iam), &
       MPI_NOAH_STRUCT,noah,di_array,displs,MPI_NOAH_STRUCT, &
       0,MPI_COMM_WORLD, ierr)
#endif
!EOC  
end subroutine noah_gather


