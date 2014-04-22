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
! !ROUTINE: mos_setup.F90
!
! !DESCRIPTION:
!
! Complete the setup routines for mosaic
!
! !REVISION HISTORY:
! 4 Nov. 1999: Paul Houser; Initial Code
! 3 Jun  2003: Jon Gottschalck; Modified code
!
! !INTERFACE:

 subroutine mos_setup()
! !USES:
      use lisdrv_module, only : lis, tile 
      use mos_varder
      use sibalb_module
      use spmdMod, only : masterproc, npes
!EOP
      IMPLICIT NONE
      integer :: t, n
!=== End Variable List ===================================================
!BOC
#if ( ! defined OPENDAP )
  if ( masterproc ) then
#endif  
      call setmosp()
      call mapsib2umd()
      call mos_coldstart()
      do t=1,lis%d%nch
         mos(t)%swnet = 0
         mos(t)%lwnet = 0
         mos(t)%qle = 0
         mos(t)%qh = 0
         mos(t)%qg = 0
         mos(t)%swrad = 0
         mos(t)%lwrad = 0
         mos(t)%rainf = 0
         mos(t)%snowf = 0
         mos(t)%evap = 0
         mos(t)%qs = 0
         mos(t)%qsm = 0
         mos(t)%qsb = 0
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
         mos(t)%soilm_prev = 0
         mos(t)%swe_prev = 0
         mos(t)%count = 0
         mos(t)%dtcanal = 0
         mos(t)%soilmr=0
	mos(t)%soilm1=0
	mos(t)%soilmtot=0
	mos(t)%mstavr=0
	mos(t)%soilv1=0
	mos(t)%soilv2=0
        mos(t)%soilv3=0
        mos(t)%soilv4=0
        mos(t)%soilv5=0
        mos(t)%soilv6=0
        mos(t)%soilv7=0
        mos(t)%soilv8=0
	mos(t)%green=0
	mos(t)%laiout=0
	mos(t)%snwfrcout=0
	mos(t)%acond=0
	mos(t)%ccond=0
	mos(t)%albedo=0
	mos(t)%avgsurft=0
	mos(t)%snohf=0
	mos(t)%sbsno=0
	mos(t)%snod=0
	mos(t)%snoc=0

      enddo

#if ( ! defined OPENDAP )
   endif
   if ( npes > 1 ) then
      call mos_scatter
   endif
#endif
!EOC
 end subroutine mos_setup

