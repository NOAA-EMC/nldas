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
! !MODULE: vic_varder.F90
!
! !DESCRIPTION:
! Module for 1-D VIC land model driver initialization. 
!
! !REVISION HISTORY:
! 
! 21 Nov 2003; Sujay Kumar : Initial Specification
! 
! !INTERFACE:
module vic_varder
! !USES:  
  use tile_spmdMod
  use vicdrv_module
  use vicpardef_module
!EOP  
  type(vicdrvdec) :: vicdrv
  SAVE
contains
!BOP
! 
! !ROUTINE: vic_varder_ini
! 
! !DESCRIPTION:   
! Reads in runtime VIC parameters, allocates memory for variables
!
! !INTERFACE:
  subroutine vic_varder_ini(nch)
!EOP
    integer :: nch
    integer :: ier
!BOC
    if(masterproc) then
       call readviccrd(vicdrv)
    endif
    call def_vicpar_struct
#if (defined SPMD)
    call MPI_BCAST(vicdrv, 1, MPI_VICDRV_STRUCT, 0, & 
         MPI_COMM_WORLD, ier)
#endif
#if ( defined OPENDAP )
!    call opendap_init_c_struct
#endif
    call vic_allocate(nch, di_array,displs,vicdrv%vic_snowband)
!EOC
  end subroutine vic_varder_ini
end module vic_varder
