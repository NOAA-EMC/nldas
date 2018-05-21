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
! !ROUTINE: noah_totinit.F90
!
! !DESCRIPTION:
!  Initialize NOAH output arrays
!
! !REVISION HISTORY:
! 
!  14 Jun 2002  Sujay Kumar  Initial Specification
!
! !INTERFACE:
subroutine noah_totinit()
! !USES:
  use noah_varder      ! NOAH LSM module  
  use tile_spmdMod
  use lisdrv_module, only : lis
!EOP
  IMPLICIT NONE

!=== End Variable List ===================================================
  integer t, i
!BOC
  do t = 1, di_array(iam)
     if(mod(lis%t%gmt,noahdrv%writeintn).eq.0)then
        noah(t)%soilm_prev=noah(t)%smc(1)*1000.0*0.1+ &
             noah(t)%smc(2)*1000.0*0.3 + & 
             noah(t)%smc(3)*1000.0*0.6 + & 
             noah(t)%smc(4)*1000.0
        noah(t)%swe_prev =  noah(t)%sneqv*1000.0
     endif
  enddo
  do t = 1, di_array(iam)
     noah(t)%swnet = 0
     noah(t)%lwnet = 0
     noah(t)%qle = 0
     noah(t)%qh = 0
     noah(t)%qg = 0
     noah(t)%snowf = 0
     noah(t)%rainf = 0
     noah(t)%evap = 0
     noah(t)%qs = 0
     noah(t)%qsb = 0
     noah(t)%qsm = 0
     noah(t)%ecanop = 0
     noah(t)%tveg = 0
     noah(t)%esoil = 0
     noah(t)%count = 0
  enddo
!EOC  
end subroutine noah_totinit

