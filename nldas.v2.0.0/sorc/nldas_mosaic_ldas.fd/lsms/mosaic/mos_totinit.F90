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
! !ROUTINE: mos_totinit.F90
!
! !DESCRIPTION:
!  Initialize Mosaic output arrays
!
! !REVISION HISTORY:
!
!  1 Aug 2003  Sujay Kumar  Initial Specification
!
! !INTERFACE:
      subroutine mos_totinit()
! !USES:
      use mos_varder      ! Mosaic LSM module
      use tile_spmdMod
      use lisdrv_module, only : lis
!EOP      
  IMPLICIT NONE

!=== End Variable List ===================================================
  integer t
!BOC 
  do t = 1, di_array(iam)
     if(mod(lis%t%gmt,mosdrv%writeintm).eq.0)then
        mos(t)%soilm_prev=mos(t)%water1 + & 
             mos(t)%water2 + & 
             mos(t)%water3
        mos(t)%swe_prev =  mos(t)%snow
     endif
  enddo
  do t = 1, di_array(iam)
     
    mos(t)%swnet = 0
    mos(t)%lwnet = 0
    mos(t)%qle = 0
    mos(t)%qh = 0
    mos(t)%qg = 0
    mos(t)%rainf = 0
    mos(t)%snowf = 0
    mos(t)%snohf = 0
    mos(t)%evap = 0
    mos(t)%qs = 0
    mos(t)%qsm = 0
    mos(t)%qsb = 0
    mos(t)%sbsno = 0
    mos(t)%swe = 0
    mos(t)%soilmoist1 = 0
    mos(t)%soilmoist2 = 0
    mos(t)%soilmoist3 = 0
    mos(t)%soilwet = 0
    mos(t)%ecanop = 0
    mos(t)%tveg = 0
    mos(t)%esoil = 0
    mos(t)%rootmoist = 0
    mos(t)%canopint = 0
    mos(t)%count = 0
    mos(t)%lwrad = 0
    mos(t)%swrad = 0
  enddo
!EOC
end subroutine mos_totinit
