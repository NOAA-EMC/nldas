!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://www.clm.org?
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! clm1_initout.f: 
!
! DESCRIPTION:
!  Initialize CLM output arrays for use in LDAS
!
! REVISION HISTORY:
! 29 Oct. 1999: Jon Radakovich; Initial Code
! 27 Sep. 2000: Brian Cosgrove; Revisions to enable CLM to 
!               output ALMA/LDAS variables
!=========================================================================

       SUBROUTINE clm1_initout (drv,clm1)

! Declare modules and data structures

       use drv_module      ! clm1 driver variables 
       use clm1type         ! 1-D ClM variables
       implicit none
       type (drvdec) ,intent(inout)     :: drv
       type (clm11d)  ,intent(inout)     :: clm1 (drv%nch)

!=== End Variable List ===================================================

      clm1%totfsa=0.                   
      clm1%toteflx_lwrad_net=0.
      clm1%toteflx_lh_tot=0.
      clm1%toteflx_sh_tot=0.
      clm1%toteflx_soil_grnd=0.
      clm1%totqflx_snomelt=0.
      clm1%totsolisbd=0.
      clm1%totforc_lwrad=0.
      clm1%totsnow=0.
      clm1%totrain=0.
      clm1%totqflx_evap=0.
      clm1%totqflx_surf=0.
      clm1%totqflx_drain=0.
      clm1%totqflx_ecanop=0.
      clm1%totqflx_tran_veg=0.
      clm1%totqflx_evap_grnd=0.
      clm1%totqflx_sub_snow=0.
 
      clm1%count=0 

      END SUBROUTINE clm1_initout




