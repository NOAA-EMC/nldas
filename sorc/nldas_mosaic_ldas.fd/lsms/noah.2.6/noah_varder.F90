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
! !ROUTINE: noah_varder.F90
!
! !DESCRIPTION:
!  Module for 1-D NOAH land model driver variable initialization
!
! !REVISION HISTORY:
! Apr 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:
module noah_varder
! !USES:        
  use noah_module
  use tile_spmdMod
  use noahpardef_module
  use noahdrv_module
!EOP  
  type(noahdrvdec) :: noahdrv
  type(noahdec), allocatable :: noah(:)
  SAVE
contains
!BOP
! 
! !ROUTINE: noah_varder_ini
! 
! !DESCRIPTION:        
! Reads in runtime noah parameters, allocates memory for variables
! 
! !INTERFACE:
  subroutine noah_varder_ini(nch)
! !USES:
#if ( defined OPENDAP )
    use opendap_module
#endif
!EOP
    integer :: nch
!BOC
    if(masterproc) then
       call readnoahcrd(noahdrv)
    endif
    call def_noahpar_struct
#if (defined SPMD)
    call MPI_BCAST(noahdrv, 1, MPI_NOAHDRV_STRUCT, 0, & 
         MPI_COMM_WORLD, ierr)
#endif

#if ( defined OPENDAP )
    noahdrv%noah_albfile = trim(opendap_data_prefix)//'/'// &
                           trim(adjustl(ciam))//'/'//noahdrv%noah_albfile
    noahdrv%noah_mgfile  = trim(opendap_data_prefix)//'/'// &
                           trim(adjustl(ciam))//'/'//noahdrv%noah_mgfile
    noahdrv%noah_mxsnal  = trim(opendap_data_prefix)//'/'// &
                           trim(adjustl(ciam))//'/'//noahdrv%noah_mxsnal
    noahdrv%noah_tbot    = trim(opendap_data_prefix)//'/'// &
                           trim(adjustl(ciam))//'/'//noahdrv%noah_tbot
#endif
    if(masterproc) then 
       allocate(noah(nch))
    else
       allocate(noah(di_array(iam)))
    endif
  end subroutine noah_varder_ini
!EOC
end module noah_varder



