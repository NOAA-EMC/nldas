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
! !ROUTINE : hyssib$_-$scatter.F90
!
! !DESCRIPTION:
!  Distributes HY-SSiB tiles on to compute nodes
!
! !REVISION HISTORY:
!    Apr 2003: Sujay Kumar, Initial Code
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      subroutine hyssib_scatter()
! !USES:
      use tile_spmdMod
      use hyssib_varder
      use hyssibpardef_module
!EOP
      implicit none

      integer :: t
!=== End Variable List ===================================================
      integer ierr
#if (defined SPMD)
      call MPI_SCATTERV(hyssib,di_array,displs, & 
          MPI_HYSSIB_STRUCT,hyssib,di_array(iam),MPI_HYSSIB_STRUCT, & 
          0,MPI_COMM_WORLD,ierr)
#endif
      end subroutine hyssib_scatter

