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
! !MODULE: vicpardef_module.F90
! 
! !DESCRIPTION: 
!
! This module contains routines that defines MPI derived data types
! for VIC LSM
!
! !REVISION HISTORY:
! 
! 14 Oct 2003; Sujay Kumar  Initial Specification 
!
! !INTERFACE:
module vicpardef_module
! !USES:
  use vicdrv_module
  use spmdMod
!EOP
  implicit none
! !ARGUMENTS:
#if (defined SPMD)  
  integer :: MPI_VICDRV_STRUCT !MPI derived type for vicdrv$_-$module
!EOP
  
  integer, parameter :: vicdrv_ntypes = 3
  integer, dimension(vicdrv_ntypes) :: vicdrv_blkcnts =(/9,400,2/)
  integer, dimension(vicdrv_ntypes) :: vicdrv_types = & 
       (/MPI_INTEGER, MPI_CHARACTER, MPI_REAL/)
  integer, dimension(vicdrv_ntypes) :: vicdrv_displs
#endif  
contains
!BOP
! !ROUTINE: def_vicpar_struct
!
! !DESCRIPTION:
! 
! Routine that defines MPI derived data types for VIC
!
! !INTERFACE:
  subroutine def_vicpar_struct()
!EOP
    integer:: t,l, ierr
    type(vicdrvdec) :: vicdrv
#if (defined SPMD)    
    call MPI_ADDRESS(vicdrv%vic_nlayer, vicdrv_displs(1),ierr)
    call MPI_ADDRESS(vicdrv%vic_sfile, vicdrv_displs(2),ierr)
    call MPI_ADDRESS(vicdrv%writeintvic, vicdrv_displs(3),ierr)
    
    do l=vicdrv_ntypes, 1, -1
       vicdrv_displs(l) = vicdrv_displs(l)-vicdrv_displs(1)
    enddo
    call MPI_TYPE_STRUCT(vicdrv_ntypes, vicdrv_blkcnts, vicdrv_displs, & 
         vicdrv_types, MPI_VICDRV_STRUCT, ierr)
    call MPI_TYPE_COMMIT(MPI_VICDRV_STRUCT, ierr)
#endif
  end subroutine def_vicpar_struct
end module vicpardef_module
