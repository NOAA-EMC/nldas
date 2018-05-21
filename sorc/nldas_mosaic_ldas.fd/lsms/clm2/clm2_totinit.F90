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
! !ROUTINE: clm2_totinit.F90
!
! !DESCRIPTION:
!  Initialize CLM output arrays
!
! !REVISION HISTORY:
! 
!  14 Jun 2002  Sujay Kumar  Initial Specification
!
! !INTERFACE:
subroutine clm2_totinit()
! !USES:
  use clm_varder
  use tile_spmdMod
!EOP
  implicit none
  integer t, m
  real :: soilm(di_array(iam),1:nlevsoi)
  real :: soilmtc(di_array(iam))
       
!=== End Variable List ===================================================
!BOC
  do t = 1, di_array(iam)
     clm(t)%totfsa=0.              ! solar absorbed solar radiation [W/m2]
     clm(t)%toteflx_lwrad_net=0.   ! net longwave radiation [W/m2]
     clm(t)%toteflx_lh_tot=0.      ! total latent heat flux [W/m2]
     clm(t)%toteflx_sh_tot=0.      ! total sensible heat flux [W/m2]      
     clm(t)%toteflx_soil_grnd=0.   ! ground heat flux [W/m2]
     clm(t)%totqflx_snomelt=0.     ! snowmelt heat flux [W/m2]
     clm(t)%totrain=0.             ! accumulation of rain [mm]
     clm(t)%totsnow=0.             ! accumulation of snow [mm]
     clm(t)%totqflx_evap=0.        ! total evaporation [mm]
     clm(t)%totqflx_surf=0.        ! surface runoff [mm]
     clm(t)%totqflx_drain=0.       ! subsurface runoff [mm]
     clm(t)%totqflx_ecanop=0.      ! interception evaporation [W/m2]
     clm(t)%totqflx_tran_veg=0.    
     clm(t)%totqflx_evap_grnd=0.
     clm(t)%totqflx_sub_snow=0.
     clm(t)%count=0
   enddo
   soilmtc = 0.0
   do m=1,nlevsoi 
      do t=1,di_array(iam)
         soilm(t,m)=clm(t)%h2osoi_liq(m)+clm(t)%h2osoi_ice(m)
      enddo
   enddo
   do m=1,nlevsoi
      do t=1,di_array(iam)
         soilmtc(t)=soilmtc(t)+soilm(t,m)
      enddo
   enddo
   do t=1,di_array(iam)
      clm(t)%soilmtc_prev = soilmtc(t)
      clm(t)%h2osno_prev = clm(t)%h2osno
   enddo
!EOC   
 end subroutine clm2_totinit

