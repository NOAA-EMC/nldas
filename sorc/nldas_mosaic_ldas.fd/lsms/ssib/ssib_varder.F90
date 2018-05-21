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
! !ROUTINE: ssib_varder.F90
!
! !DESCRIPTION:
!  Module for 1-D SSiB land model driver variable initialization
!
! !REVISION HISTORY:
!    Apr 2003: Sujay Kumar, Initial Code
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
!  6 May 2004: David Mocko, made compatible with SiB-lings
!
! !INTERFACE:
      module ssib_varder
! !USES:
      use ssib_module
      use tile_spmdMod
      use ssibpardef_module
      use ssibdrv_module
!EOP
      type(ssibdrvdec) :: ssibdrv
      type(ssibdec), allocatable :: ssib(:)
      SAVE
      contains
!BOP
!
! !ROUTINE: ssib_varder_ini
!
! !DESCRIPTION:
!  Reads in runtime ssib parameters, allocates memory for variables
!
! !INTERFACE:
      subroutine ssib_varder_ini(nch)
! !USES:
#if ( defined OPENDAP )
      use opendap_module
#endif
!EOP
      integer :: nch
!BOC
      if (masterproc) then
         call readssibcrd(ssibdrv)
      endif
      call def_ssibpar_struct
#if (defined SPMD)
      call MPI_BCAST(ssibdrv, 1, MPI_SSIBDRV_STRUCT, 0, & 
                     MPI_COMM_WORLD, ierr)
#endif
      if (masterproc) then
         allocate(ssib(nch))
      else
         allocate(ssib(di_array(iam)))
      endif
      return
      end subroutine ssib_varder_ini
!EOC
      end module ssib_varder

