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
!BOP
!
! !ROUTINE: vic_binout.F90
!
! !DESCRIPTION:  
!  LIS VIC data writer: Writes vic output in binary format
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine vic_binout(ftn)
! !USES:
  use lisdrv_module, only : lis
  use drv_output_mod, only : drv_writevar_bin
  use vic_varder
 
  implicit none
  
  integer :: ftn
  real    :: tmp(lis%d%glbnch)
!EOP
  integer :: t, c, r, nvars
!BOC
  nvars = 30
  if ( lis%o%wfor == 1 ) then
     nvars = 38
  endif

  do t = 1, nvars 
     call get_vicvar(t, lis%d%glbnch, tmp)
     call drv_writevar_bin(ftn,tmp) !Net shortwave radiation (surface) (W/m2)
  enddo
  
!EOC
 end subroutine vic_binout
 
