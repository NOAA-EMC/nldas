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
!
! !ROUTINE: clm2_gather
!
! !DESCRIPTION:
!  Gathers CLM tiles
!
! !REVISION HISTORY:
! 30 Jan 2003: Sujay Kumar Initial Specification
!
! !INTERFACE:
subroutine clm2_gather()
! !USES:
  use clm_varder
  use tile_spmdMod
  use clm2pardef_module
!EOP
  implicit none
  integer ierr
!BOC
#if  (defined SPMD)
  call MPI_GATHERV(clm(1:di_array(iam)),di_array(iam), & 
       MPI_CLM_STRUCT,clm,di_array,displs,MPI_CLM_STRUCT, & 
       0,MPI_COMM_WORLD, ierr)
#endif
!EOC  
end subroutine clm2_gather



