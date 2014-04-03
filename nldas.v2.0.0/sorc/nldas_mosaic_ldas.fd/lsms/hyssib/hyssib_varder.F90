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
! !ROUTINE : hyssib$_-$varder.F90
!
! !DESCRIPTION:
!  Module for 1-D HY-SSiB land model driver variable initialization
!
! !REVISION HISTORY:
!    Apr 2003: Sujay Kumar, Initial Code
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      module hyssib_varder
!EOP
      use hyssib_module
      use tile_spmdMod
      use hyssibpardef_module
      use hyssibdrv_module

      type(hyssibdrvdec) :: hyssibdrv
      type(hyssibdec), allocatable :: hyssib(:)
      SAVE
      contains

      subroutine hyssib_varder_ini(nch)

      integer :: nch
      if (masterproc) then
         call readhyssibcrd(hyssibdrv)
      endif
      call def_hyssibpar_struct
#if (defined SPMD)
      call MPI_BCAST(hyssibdrv, 1, MPI_HYSSIBDRV_STRUCT, 0, & 
                     MPI_COMM_WORLD, ierr)
#endif 
      if (masterproc) then
         allocate(hyssib(nch))
      else
         allocate(hyssib(di_array(iam)))
      endif
      end subroutine hyssib_varder_ini
      end module hyssib_varder

