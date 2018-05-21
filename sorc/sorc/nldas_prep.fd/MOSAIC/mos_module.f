!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! mos_module.f: 
!
! DESCRIPTION:
!  Module for MOASIC tile space variable specification.
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial code
!  11 Feb 2002: Jon Gottschalck; Added AVHRR derived variables
!=========================================================================

      MODULE mos_module 

      IMPLICIT NONE
      public mosdec

      type mosdec

!=== LDAS-MOSAIC Parameters ===============================================
      REAL, pointer :: VEGP(:)       !Static vegetation parameter values, dim(MOS_NVEGP)
      REAL, pointer :: VEGMP(:,:)    !Monthly vegetation parameter values dim(MOS_NMVEGP,12)
      REAL, pointer :: VEGIP(:)      !Interpolated from monthly parameters (MOS_NMVEGP)
      REAL, pointer :: SOILP(:)      !Static soil parameter values, dim(NOS_NSOILP)
      INTEGER :: NSLAY               !Number of MOSAIC soil layers

      REAL, pointer :: LAI_T1_F      !AVHRR Derived LAI Time 1
      REAL, pointer :: LAI_T2_F      !AVHRR Derived LAI Time 2
      REAL, pointer :: DSAI_T1       !AVHRR Derived DSAI Time 1
      REAL, pointer :: DSAI_T2       !AVHRR Derived DSAI Time 2
      REAL, pointer :: GREEN1        !AVHRR Derived GREENNESS Time 1
      REAL, pointer :: GREEN2        !AVHRR Derived GREENNESS Time 2
      REAL, pointer :: LAI           !AVHRR Derived LAI
      REAL, pointer :: GREEN         !AVHRR Derived GREENNESS
      REAL, pointer :: DSAI          !AVHRR Derived DSAI

!=== LDAS-MOSAIC States ===============================================
      REAL :: CT                     !MOSAIC Canopy/Soil Temperature 
      REAL :: QA                     !MOSAIC Canopy Humidity
      REAL :: ICS                    !MOSAIC Interception Canopy Storage
      REAL :: SNOW                   !MOSAIC Snow Depth
      REAL :: SoT                    !MOSAIC Deep Soil Temperaure
      REAL, pointer :: SoWet(:)      !MOSAIC Soil Wetness (3 layers)
!=== Analysis and Bias Correction Variables ===========================
      REAL :: DTCANAL                !MOSAIC Change in Temperature based on Analysis

!=== LDAS-MOSAIC OUTPUT States ============================================
      INTEGER :: COUNT                !MOSAIC Output Counter
      REAL, pointer :: RETURN(:)      !MOSAIC Output Array
      REAL, pointer :: TOTRET(:)      !MOSAIC Output Averaging Array

      end type
      end module mos_module





















