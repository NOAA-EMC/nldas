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
! !ROUTINE: template_f2t.F90
!
! !DESCRIPTION:
!  
!  template for calling the forcing transfer routines.
!
! !REVISION HISTORY:
! 21 Jul 2004: Sujay Kumar   Initial Specification
! 
! !INTERFACE:
subroutine template_f2t(t, forcing)
! !USES:
  use template_varder
  use lisdrv_module, only : lis
!EOP
  implicit none
  real :: forcing(10)
  integer :: t, f    ! Loop counters
!BOC
  do f=1,lis%f%nforce
     template(t)%forcing(f)=forcing(f)
  enddo

!EOC  
end subroutine template_f2t

