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
! !ROUTINE: clm2_restart.F90
!
! !DESCRIPTION:
!  This program reads the restart files for CLM
!
! !REVISION HISTORY:
! 20 Jan 2003; Sujay Kumar Initial Specification
! 
! !INTERFACE:
subroutine clm2_restart()
! !USES:
  use spmdMod, only : masterproc,npes
  use lisdrv_module, only : lis
  use restFileMod   , only : restrd
!EOP
!  integer :: rw
!BOC  
  if(masterproc) then 
     if(lis%o%startcode.eq.1) then 
        print*, 'Reading restart files..'
        call restrd()
     endif
  endif
#if ( defined SPMD )
  if ( ( lis%o%startcode == 1 ) .and. ( npes >  1 ) ) then
     call clm2_scatter()
  endif
#endif
!EOC  
end subroutine clm2_restart
