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
! !ROUTINE: mos_varder.F90
!
! !DESCRIPTION:
!  Module for 1-D MOSAIC land model driver variable initialization
!
! !REVISION HISTORY:
! Jun 2003; Jon Gottschalck, Initial Code
!
! !INTERFACE:
    module mos_varder
! !USES:        
    use mos_module
    use mosdrv_module
    use mospardef_module
    use tile_spmdMod
!EOP       
    type(mosdec), allocatable :: mos(:)
    type(mosdrvdec) :: mosdrv
    SAVE
contains
!BOP  
!
! !ROUTINE: mos_varder_ini
!
! !DESCRIPTION:
! Reads in runtime mos parameters, allocates memory for variables
!
! !INTERFACE:
       subroutine mos_varder_ini(nch)
! !USES:
#if ( defined OPENDAP )
  use opendap_module
#endif
!EOP  
  integer :: nch
!BOC
  if(masterproc) then
     call readmoscrd(mosdrv)
  endif
  call def_mospar_struct
#if (defined SPMD)
  call MPI_BCAST(mosdrv,1,MPI_MOSDRV_STRUCT, 0, &
       MPI_COMM_WORLD, ierr)
#endif 
  if(masterproc) then 
    allocate(mos(nch))
  else
    allocate(mos(di_array(iam)))
  endif
  end subroutine mos_varder_ini
!EOC
end module mos_varder



