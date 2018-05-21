!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale 
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project 
!  S                                      A  This is the GSFC-LDAS code. 
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov 
!                                                                  
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers         
!
!  NOAH LAND-SURFACE MODEL, UNCOUPLED 1-D COLUMN: VERSION 2.5 OCT 2001
!
!=========================================================================
! noah_module.f: 
!
! DESCRIPTION:
!  Module for 1-D NOAH land model driver variable specification.
!
! REVISION HISTORY:
! 28 Apr, 2002: K. Arsenault added NOAH LSM 2.5 code to LDAS. 
! 04 Nov, 2002: K. Arsenault; Included new bottom temperature fields
! 24 Jun, 2002: K. Arsenault; Added VEGMP1,VEGMP2,ALBSF1,ALBSF2  
!=========================================================================     

      MODULE noah_module

      implicit none
      public noahdec

      type noahdec

!=== NOAH Parameters ======================================================

      INTEGER :: ts                 !Timestep (seconds)
      INTEGER :: maxt               !Maximum tiles per grid
      INTEGER :: SIBVEG             !UMD to SIB Vegetation Class Index value
      INTEGER :: NSLAY              !Number of NOAH soil layers (4)
      INTEGER, pointer :: ZOBSOIL(:) !Zobler Soil Classes (LDAS%NCH)

      REAL, pointer :: VEGP(:)      !Static vegetation parameter values, dim(NOAH_NVEGP)
      REAL, pointer :: VEGIP(:)     !Interpolated Green Fraction from monthly parameters 
      REAL, pointer :: VEGMP1(:)    !Month 1 Green Fraction Values
      REAL, pointer :: VEGMP2(:)    !Month 2 Green Fraction Values
      REAL, pointer :: SOILP(:)     !Static soil parameter values, dim(NOAH_NSOILP)
      REAL, pointer :: ALBSF(:)     !Quarterly Snow-Free Albedo dataset
      REAL, pointer :: ALBSF1(:)    !Date 1 Snow-Free albedo values
      REAL, pointer :: ALBSF2(:)    !Date 2 Snow-Free albedo values
      REAL, pointer :: MXSNALB(:)   !Maximum snow albedo dataset
      REAL, pointer :: TEMPBOT(:)   !Bottom boundary temperature

!=== NOAH-State Variables ================================================

      REAL :: T1                 !NOAH Skin Temperature (K)
      REAL :: CMC                !NOAH Canopy Water Content 
      REAL :: SNOWH              !NOAH Actual Snow depth (m) 
      REAL :: SNEQV              !NOAH Water Equivalent Snow Depth (m)
      REAL, pointer :: STC(:)    !NOAH Soil Temperaure (4 layers)
      REAL, pointer :: SMC(:)    !NOAH Soil (4 layers)
      REAL, pointer :: SH2O(:)   !NOAH Liquid-only soil moisture (4 layers)
      REAL :: CH                 !NOAH Heat/moisture exchange coef.
      REAL :: CM                 !NOAH Momentum exchange coef.

!=== NOAH Output Variables ===============================================

      INTEGER :: COUNT             !NOAH Output Counter
      REAL, pointer :: RETURN(:)   !NOAH Output Array
      REAL, pointer :: TOTRET(:)   !NOAH Output Averaging Array

!=== End Variable List ===================================================

      end type

      end module noah_module
