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
!BOP
!
! !MODULE: mospardef_module.F90
! 
! !DESCRIPTION: 
!
! This module contains routines that defines MPI derived data types
! for Mosaic LSM
!
! !REVISION HISTORY:
! 
! 15 Oct 2003; Sujay Kumar  Initial Specification 
!
#include "misc.h"
! !INTERFACE:
      module mospardef_module
! !USES:
      use mos_module
      use mosdrv_module
      use spmdMod
!EOP
      implicit none
! !ARGUMENTS:
#if (defined SPMD)
      integer:: MPI_MOS_STRUCT  !MPI derived type for noah$_-$module
      integer :: MPI_MOSDRV_STRUCT !MPI derived type for noahdrv$_-$module
!EOP
      integer, parameter :: mos_ntypes = 2
      integer, dimension(mos_ntypes) :: mos_blkcnts =(/4,93/)
      integer, dimension(mos_ntypes) :: mos_types = & 
              (/MPI_INTEGER, MPI_REAL/)
      integer, dimension(mos_ntypes) :: mos_displs

      integer, parameter :: mosdrv_ntypes = 3
      integer, dimension(mosdrv_ntypes) :: mosdrv_blkcnts =(/9,400,3/)
      integer, dimension(mosdrv_ntypes) :: mosdrv_types = & 
              (/MPI_INTEGER, MPI_CHARACTER,MPI_REAL/)
      integer, dimension(mosdrv_ntypes) :: mosdrv_displs
#endif
      contains
!BOPI
!
! !DESCRIPTION:
! 
! Routine that defines MPI derived data types for Mosaic
!
! !INTERFACE:
      subroutine def_mospar_struct()
!EOPI
      integer:: t,l, ierr
      type(mosdec)::mos
      type(mosdrvdec) :: mosdrv
#if (defined SPMD)
      call MPI_ADDRESS(mos%ts, mos_displs(1),ierr)
      call MPI_ADDRESS(mos%vegp(1), mos_displs(2),ierr)

      do l=mos_ntypes, 1, -1
         mos_displs(l) = mos_displs(l)-mos_displs(1)
      enddo
      call MPI_TYPE_STRUCT(mos_ntypes, mos_blkcnts, mos_displs, & 
          mos_types, MPI_MOS_STRUCT, ierr)
      call MPI_TYPE_COMMIT(MPI_MOS_STRUCT, ierr)

      call MPI_ADDRESS(mosdrv%mosopen, mosdrv_displs(1),ierr)
      call MPI_ADDRESS(mosdrv%mos_rfile, mosdrv_displs(2),ierr)
      call MPI_ADDRESS(mosdrv%mos_ism, mosdrv_displs(3),ierr)

      do l=mosdrv_ntypes, 1, -1
         mosdrv_displs(l) = mosdrv_displs(l)-mosdrv_displs(1)
      enddo
      call MPI_TYPE_STRUCT(mosdrv_ntypes, mosdrv_blkcnts, mosdrv_displs, & 
          mosdrv_types, MPI_MOSDRV_STRUCT, ierr)
      call MPI_TYPE_COMMIT(MPI_MOSDRV_STRUCT, ierr)
#endif
      end subroutine def_mospar_struct
      end module 
