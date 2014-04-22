!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! clm2_initout.f: 
!
! DESCRIPTION:
!  Initialize CLM2 output arrays for use in LDAS
!
! REVISION HISTORY:
! 30 Jun 2002: Jon Gottschalck; Initial Code
!=========================================================================

       SUBROUTINE clm2_initout ()

! Declare modules and data structures

       use clm_varder
       implicit none

!=== End Variable List ===================================================

      clm%totfsa=0.              ! solar absorbed solar radiation [W/m2]
      clm%toteflx_lwrad_net=0.   ! net longwave radiation [W/m2]
      clm%toteflx_lh_tot=0.      ! total latent heat flux [W/m2]
      clm%toteflx_sh_tot=0.      ! total sensible heat flux [W/m2]      
      clm%toteflx_soil_grnd=0.   ! ground heat flux [W/m2]
      clm%totqflx_snomelt=0.     ! snowmelt heat flux [W/m2]
      clm%totsolisbd=0.          ! total downward surface shortwave radiation [W/m2]
      clm%totforc_lwrad=0.       ! atmospheric infrared (longwave) radiation [W/m2]
      clm%totrain=0.             ! accumulation of rain [mm]
      clm%totsnow=0.             ! accumulation of snow [mm]
      clm%totqflx_evap=0.        ! total evaporation [mm]
      clm%totqflx_surf=0.        ! surface runoff [mm]
      clm%totqflx_drain=0.       ! subsurface runoff [mm]
      clm%totqflx_ecanop=0.      ! interception evaporation [W/m2]
      clm%totqflx_tran_veg=0.    
      clm%totqflx_evap_grnd=0.
      clm%totqflx_sub_snow=0.
 
      clm%count=0
      clm%acond=0.

      END SUBROUTINE clm2_initout




