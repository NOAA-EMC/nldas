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
! !MODULE: ssibpardef_module.F90
! 
! !DESCRIPTION: 
!
! This module contains routines that defines MPI derived data types
! for SSiB LSM
!
! !REVISION HISTORY:
! 06 Oct 2003: Sujay Kumar, Initial Specification 
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
! 22 May 2004: David Mocko, made compatible with SiB-lings
!
#include "misc.h"
! !INTERFACE:
      module ssibpardef_module
! !USES:
      use ssib_module
      use ssibdrv_module
      use spmdMod
!EOP
      implicit none
! !ARGUMENTS:
#if (defined SPMD)
      integer:: MPI_SSIB_STRUCT  !MPI derived type for ssib$_-$module
      integer :: MPI_SSIBDRV_STRUCT !MPI derived type for ssibdrv$_-$module
!EOP
      integer, parameter :: ssib_ntypes = 2
      integer, dimension(ssib_ntypes) :: ssib_blkcnts =(/9,131/)
      integer, dimension(ssib_ntypes) :: ssib_types = & 
              (/MPI_INTEGER, MPI_REAL/)
      integer, dimension(ssib_ntypes) :: ssib_displs

      integer, parameter :: ssibdrv_ntypes = 4
      integer, dimension(ssibdrv_ntypes) :: ssibdrv_blkcnts =(/13,340,1,3/)
      integer, dimension(ssibdrv_ntypes) :: ssibdrv_types = & 
              (/MPI_INTEGER, MPI_CHARACTER,MPI_REAL8,MPI_REAL/)
      integer, dimension(ssibdrv_ntypes) :: ssibdrv_displs
#endif
      contains
!BOPI
!
! !DESCRIPTION:
! 
! Routine that defines MPI derived data types for SSiB
!
! !INTERFACE:
      subroutine def_ssibpar_struct
!EOPI
      integer:: t,l, ierr
      type(ssibdec)::ssib
      type(ssibdrvdec) :: ssibdrv
#if (defined SPMD)
      call MPI_ADDRESS(ssib%ts, ssib_displs(1),ierr)
      call MPI_ADDRESS(ssib%vegp(1), ssib_displs(2),ierr)

      do l=ssib_ntypes, 1, -1
         ssib_displs(l) = ssib_displs(l)-ssib_displs(1)
      enddo
      call MPI_TYPE_STRUCT(ssib_ntypes, ssib_blkcnts, ssib_displs, & 
          ssib_types, MPI_SSIB_STRUCT, ierr)
      call MPI_TYPE_COMMIT(MPI_SSIB_STRUCT, ierr)

      call MPI_ADDRESS(ssibdrv%ssibopen, ssibdrv_displs(1),ierr)
      call MPI_ADDRESS(ssibdrv%ssib_rfile, ssibdrv_displs(2),ierr)
      call MPI_ADDRESS(ssibdrv%ssib_gfractime, ssibdrv_displs(3),ierr)
      call MPI_ADDRESS(ssibdrv%ssib_ism, ssibdrv_displs(4),ierr)

      do l=ssibdrv_ntypes, 1, -1
         ssibdrv_displs(l) = ssibdrv_displs(l)-ssibdrv_displs(1)
      enddo
      call MPI_TYPE_STRUCT(ssibdrv_ntypes, ssibdrv_blkcnts, ssibdrv_displs, & 
          ssibdrv_types, MPI_SSIBDRV_STRUCT, ierr)
      call MPI_TYPE_COMMIT(MPI_SSIBDRV_STRUCT, ierr)
#endif
    end subroutine def_ssibpar_struct
  end module ssibpardef_module
      
