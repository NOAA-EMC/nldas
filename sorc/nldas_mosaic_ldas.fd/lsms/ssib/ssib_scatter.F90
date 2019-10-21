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
! !ROUTINE: ssib_scatter.F90
!
! !DESCRIPTION:
!  Distributes SSiB tiles on to compute nodes
!
! !REVISION HISTORY:
!    Apr 2003: Sujay Kumar, Initial Code
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
!  6 May 2004: David Mocko, made compatible with SiB-lings
!
! !INTERFACE:
      subroutine ssib_scatter
! !USES:
      use tile_spmdMod
      use ssib_varder
      use ssibpardef_module
!EOP
      implicit none

      integer :: t
      integer ierr
!BOC
!=== End Variable List ===================================================
#if (defined SPMD)
      call MPI_SCATTERV(ssib,di_array,displs, & 
          MPI_SSIB_STRUCT,ssib,di_array(iam),MPI_SSIB_STRUCT, & 
          0,MPI_COMM_WORLD,ierr)
#endif
      return
!EOC
      end subroutine ssib_scatter

