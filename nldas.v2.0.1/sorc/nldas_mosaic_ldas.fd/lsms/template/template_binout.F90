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
! !ROUTINE: template_binout.F90
!
! !DESCRIPTION:  
!  LIS data writer: Writes template output in binary format
!
! !REVISION HISTORY:
! 21 Jul 2004; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine template_binout(ftn)
! !USES:
  use lisdrv_module, only : lis
  use drv_output_mod, only : drv_writevar_bin
  use template_varder
 
  implicit none
  real :: rainf(lis%d%glbnch)
  real :: snowf(lis%d%glbnch)
  integer :: ftn
!EOP
  integer :: t,c,r
!BOC
  do t=1,lis%d%glbnch
     if(template(t)%forcing(1) < 273.15) then
        rainf(t) = 0.0
        snowf(t) = template(t)%forcing(8)
     else
        rainf(t) = template(t)%forcing(8)
        snowf(t) = 0.0
     endif
  enddo
   if(lis%o%wfor.eq.1) then
      call drv_writevar_bin(ftn, sqrt(template%forcing(5)*template%forcing(5)+ & 
           template%forcing(6)*template%forcing(6)))
      call drv_writevar_bin(ftn,rainf)
      call drv_writevar_bin(ftn,snowf)
      call drv_writevar_bin(ftn,template%forcing(1))
      call drv_writevar_bin(ftn,template%forcing(2))
      call drv_writevar_bin(ftn,template%forcing(7))
      call drv_writevar_bin(ftn,template%forcing(3))
      call drv_writevar_bin(ftn,template%forcing(4))
   endif

!EOC
 end subroutine template_binout
 
