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
! !MODULE: clm2pardef_module.F90
! 
! !DESCRIPTION: 
!
! This module contains routines that defines MPI derived data types
! for CLM LSM
!
! !REVISION HISTORY:
! 
! 06 Oct 2003; Sujay Kumar  Initial Specification 
!
! !INTERFACE:
module clm2pardef_module
! !USES:
  use clmtype
  use clm2drv_module
  use spmdMod
!EOP
  implicit none
! !ARGUMENTS:
#if (defined SPMD)
  integer:: MPI_CLM_STRUCT   !MPI derived type for clm$_-$module
  integer:: MPI_CLMDRV_STRUCT   !MPI derived type for clm2drv$_-$module
!EOP

  integer, parameter :: clm_ntypes = 3
  integer, dimension(clm_ntypes) :: clm_blkcnts =(/25,3,363/)
  integer, dimension(clm_ntypes) :: clm_types = & 
       (/MPI_INTEGER, MPI_LOGICAL, MPI_REAL/)
  integer, dimension(clm_ntypes) :: clm_displs
  
  integer, parameter :: clmdrv_ntypes = 3
  integer, dimension(clmdrv_ntypes) :: clmdrv_blkcnts =(/31,120,4/)
  integer, dimension(clmdrv_ntypes) :: clmdrv_types = & 
       (/MPI_INTEGER, MPI_CHARACTER, MPI_REAL/)
  integer, dimension(clmdrv_ntypes) :: clmdrv_displs
#endif
  
  
contains
!BOP
! !ROUTINE: def_clmpar_struct
!
! !DESCRIPTION:
! 
! Routine that defines MPI derived data types for CLM
!
! !INTERFACE:
  subroutine def_clmpar_struct()
!EOP
#if (defined SPMD)
    integer:: t,l, ierr
    type(clm1d):: clm
    type(clmdrvdec) :: clmdrv
    
    call MPI_ADDRESS(clm%nstep, clm_displs(1),ierr)
    call MPI_ADDRESS(clm%lakpoi, clm_displs(2),ierr)
    call MPI_ADDRESS(clm%lat, clm_displs(3),ierr)
    do l=clm_ntypes, 1, -1
       clm_displs(l) = clm_displs(l)-clm_displs(1)
    enddo
    call MPI_TYPE_STRUCT(clm_ntypes, clm_blkcnts, clm_displs, & 
         clm_types, MPI_CLM_STRUCT, ierr)
    call MPI_TYPE_COMMIT(MPI_CLM_STRUCT, ierr)
    
    call MPI_ADDRESS(clmdrv%numout, clmdrv_displs(1),ierr)
    call MPI_ADDRESS(clmdrv%clm2_rfile, clmdrv_displs(2),ierr)
    call MPI_ADDRESS(clmdrv%clm2_ism, clmdrv_displs(3),ierr)
    do l=clmdrv_ntypes, 1, -1
       clmdrv_displs(l) = clmdrv_displs(l)-clmdrv_displs(1)
    enddo
    call MPI_TYPE_STRUCT(clmdrv_ntypes, clmdrv_blkcnts, clmdrv_displs, & 
         clmdrv_types, MPI_CLMDRV_STRUCT, ierr)
    call MPI_TYPE_COMMIT(MPI_CLMDRV_STRUCT, ierr)
#endif 
  end subroutine def_clmpar_struct
end module clm2pardef_module
