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
! !ROUTINE: template_output.F90
!
! !DESCRIPTION:
!  
!  Template for calling the output routines. 
!
! !REVISION HISTORY:
! 21 Jul 2004: Sujay Kumar   Initial Specification
! 
! !INTERFACE:
subroutine template_output()
! !USES:
  use lisdrv_module, only : lis, tile, glbgindex
  use spmdMod, only : masterproc,npes
  use template_varder, only : templatedrv
!EOP
  implicit none
!BOC
  if(mod(lis%t%gmt, templatedrv%writeint).eq.0)then
     if(masterproc) then 
        call template_out()
     endif
  endif
!EOC  
end subroutine template_output

