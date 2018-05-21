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
! !MODULE: hyssibpardef_module.F90
!
! !DESCRIPTION:
! This module contains routines that defines MPI derived data types
! for HY-SSiB LSM
!
! !REVISION HISTORY:
! 06 Oct 2003: Sujay Kumar, Initial Specification
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      module hyssibpardef_module
! !USES:
      use hyssib_module
      use hyssibdrv_module
      use spmdMod
!EOP
      implicit none

! !ARGUMENTS:
#if (defined SPMD)
      integer:: MPI_HYSSIB_STRUCT  !MPI derived type for hyssib$_-$module
      integer :: MPI_HYSSIBDRV_STRUCT !MPI derived type for hyssibdrv$_-$module
!EOP
      integer, parameter :: hyssib_ntypes = 2
      integer, dimension(hyssib_ntypes) :: hyssib_blkcnts =(/6,91/)
      integer, dimension(hyssib_ntypes) :: hyssib_types = & 
              (/MPI_INTEGER, MPI_REAL/)
      integer, dimension(hyssib_ntypes) :: hyssib_displs

      integer, parameter :: hyssibdrv_ntypes = 4
      integer, dimension(hyssibdrv_ntypes) :: hyssibdrv_blkcnts =(/10,340,1,3/)
      integer, dimension(hyssibdrv_ntypes) :: hyssibdrv_types = & 
              (/MPI_INTEGER, MPI_CHARACTER,MPI_REAL8,MPI_REAL/)
      integer, dimension(hyssibdrv_ntypes) :: hyssibdrv_displs
#endif
      contains
!BOPI
!
! !DESCRIPTION:
! 
! Routine that defines MPI derived data types for HY-SSiB
!
! !INTERFACE:
      subroutine def_hyssibpar_struct()
!EOPI
      integer:: t,l, ierr
      type(hyssibdec) :: hyssib
      type(hyssibdrvdec) :: hyssibdrv
#if (defined SPMD)
      call MPI_ADDRESS(hyssib%ts, hyssib_displs(1),ierr)
      call MPI_ADDRESS(hyssib%vegp(1), hyssib_displs(2),ierr)

      do l=hyssib_ntypes, 1, -1
         hyssib_displs(l) = hyssib_displs(l)-hyssib_displs(1)
      enddo
      call MPI_TYPE_STRUCT(hyssib_ntypes, hyssib_blkcnts, hyssib_displs, & 
          hyssib_types, MPI_HYSSIB_STRUCT, ierr)
      call MPI_TYPE_COMMIT(MPI_HYSSIB_STRUCT, ierr)

      call MPI_ADDRESS(hyssibdrv%hyssibopen, hyssibdrv_displs(1),ierr)
      call MPI_ADDRESS(hyssibdrv%hyssib_rfile, hyssibdrv_displs(2),ierr)
      call MPI_ADDRESS(hyssibdrv%hyssib_gfractime, hyssibdrv_displs(3),ierr)
      call MPI_ADDRESS(hyssibdrv%hyssib_ism, hyssibdrv_displs(4),ierr)

      do l=hyssibdrv_ntypes, 1, -1
         hyssibdrv_displs(l) = hyssibdrv_displs(l)-hyssibdrv_displs(1)
      enddo
      call MPI_TYPE_STRUCT(hyssibdrv_ntypes, hyssibdrv_blkcnts, hyssibdrv_displs, & 
          hyssibdrv_types, MPI_HYSSIBDRV_STRUCT, ierr)
      call MPI_TYPE_COMMIT(MPI_HYSSIBDRV_STRUCT, ierr)
#endif
      end subroutine def_hyssibpar_struct
      end module

