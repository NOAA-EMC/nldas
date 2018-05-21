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
! !ROUTINE: clm2_scatter
!
! !DESCRIPTION:
!  Distributes clm2 tiles to compute nodes
!
! !REVISION HISTORY:
!
! Apr 2003 ; Sujay Kumar, Initial Code
!
! !INTERFACE:
subroutine clm2_scatter()
! !USES:
  use clm_varder
  use tile_spmdMod
  use clm2pardef_module
!EOP
  implicit none
  integer ierr
#if (defined SPMD)
  call MPI_SCATTERV(clm,di_array,displs, & 
       MPI_CLM_STRUCT,clm,di_array(iam),MPI_CLM_STRUCT, & 
       0,MPI_COMM_WORLD,ierr)
#endif  
end subroutine clm2_scatter

!BOP
!
! !ROUTINE: clm2_scatterlai
!
! !DESCRIPTION:
!  Distributes clm2 lai tiles to compute nodes
!
! !REVISION HISTORY:
!
! 17 Sep 2004 ; James Geiger, Initial Code
!
! !INTERFACE:
subroutine clm2_scatterlai()
! !USES:
   use clm_varder
   use tile_spmdMod
   use clm2pardef_module
   use lisdrv_module, only : lis
!EOP
   implicit none
 
   integer :: ierr
   real, allocatable, dimension(:) :: temp_array

#if ( defined SPMD )

   if ( masterproc ) then
      allocate(temp_array(lis%d%glbnch))
   else
      allocate(temp_array(di_array(iam)))
   endif

   ! scatter tlai
   if ( masterproc ) then
      temp_array = clm(:)%tlai
   endif

   call MPI_SCATTERV(temp_array,di_array,displs,    & 
        MPI_REAL,temp_array,di_array(iam),MPI_REAL, & 
        0,MPI_COMM_WORLD,ierr)

   if ( .not. masterproc ) then
      clm(:)%tlai = temp_array
   endif

   ! scatter tsai
   if ( masterproc ) then
      temp_array = clm(:)%tsai
   endif

   call MPI_SCATTERV(temp_array,di_array,displs,    & 
        MPI_REAL,temp_array,di_array(iam),MPI_REAL, & 
        0,MPI_COMM_WORLD,ierr)

   if ( .not. masterproc ) then
      clm(:)%tsai = temp_array
   endif

   ! scatter htop
   if ( masterproc ) then
      temp_array = clm(:)%htop
   endif

   call MPI_SCATTERV(temp_array,di_array,displs,    & 
        MPI_REAL,temp_array,di_array(iam),MPI_REAL, & 
        0,MPI_COMM_WORLD,ierr)

   if ( .not. masterproc ) then
      clm(:)%htop = temp_array
   endif

   ! scatter hbot
   if ( masterproc ) then
      temp_array = clm(:)%hbot
   endif

   call MPI_SCATTERV(temp_array,di_array,displs,    & 
        MPI_REAL,temp_array,di_array(iam),MPI_REAL, & 
        0,MPI_COMM_WORLD,ierr)

   if ( .not. masterproc ) then
      clm(:)%hbot = temp_array
   endif

   deallocate(temp_array)

#endif  

end subroutine clm2_scatterlai
