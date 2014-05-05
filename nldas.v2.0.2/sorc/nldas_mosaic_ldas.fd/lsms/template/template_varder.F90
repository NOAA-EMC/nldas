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
! !ROUTINE: template_varder.F90
!
! !DESCRIPTION:
!  Module for 1-D land model driver variable initialization
!
! !REVISION HISTORY:
! 21 July 2004; Sujay Kumar, Initial Code
!
! !INTERFACE:
module template_varder
! !USES:        
  use tile_spmdMod
  use template_module
  use templatedrv_module
!EOP  
  type(templatedec), allocatable :: template(:)
  type(templatedrvdec) :: templatedrv
  SAVE
contains
!BOP
! 
! !ROUTINE: template_varder_ini
! 
! !DESCRIPTION:        
! Reads in runtime template parameters, allocates memory for variables
! 
! !INTERFACE:
  subroutine template_varder_ini(nch)
! !USES:
    integer :: nch
!BOC
    if(masterproc) then 
       call readtemplatecrd(templatedrv)
    endif

    if(masterproc) then 
       allocate(template(nch))
    else
       allocate(template(di_array(iam)))
    endif
  end subroutine template_varder_ini
!EOC
end module template_varder



