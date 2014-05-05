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
! cat_module.f: 
!
! DESCRIPTION:
!  Module for CATCHMENT space variable specification.
!
! REVISION HISTORY:
!  4 April 2000: Jeffrey Walker; Initial code
!=========================================================================

      MODULE cat_module 

      IMPLICIT NONE
      public catdec

      type catdec

!=== LDAS-CATCHMENT Variables ============================================
      INTEGER :: NFF     !Pointer array to go from North America to LDAS

!=== LDAS-CATCHMENT Input ================================================
      REAL :: LATT       !Lattitude
      REAL :: LONN       !Longitude

!=== LDAS-CATCHMENT Parameters ===========================================
      INTEGER :: VEGCLS  !Vegetation class
      REAL :: TSA(2)     !Catchment topography 
      REAL :: TSB(2)     !Catchment topography 
      REAL :: ARW(4)     !Catchment topography 
      REAL :: ARA(4)     !Catchment topography 
      REAL :: ARS(3)     !Catchment topography 
      REAL :: BF(3)      !Catchment topography 
      REAL :: ZDPTH(3)   !Soil depth
      REAL :: CDCR(2)    !Moisture storage in catchment
      REAL :: POROS      !Soil porosity
      REAL :: WPWET      !Wilting point wetness
      REAL :: BEE        !Clapp and Hornberger b
      REAL :: PSIS       !Saturated soil matric suction
      REAL :: COND       !Saturated soil hydraulic conductivity
      REAL :: GNU        !Vertical conductivity decay factor
      REAL :: VGWMAX     !Moisture storage in root zone
      REAL :: GREEN(12)  !Greenness fraction
      REAL :: LAI(12)    !Leaf area index
      REAL :: ZOL(12)    ! ?
      REAL :: ALBEDO(12) !Snowfree albedo ?
      REAL :: VGD(12)    !Height of canopy leaves ?

!=== LDAS-CATCHMENT Forcing ==============================================
      REAL :: TMP2M      !2m air temperature
      REAL :: DEW2M      !2m dew point temperature
      REAL :: SFCPRS     !Surface pressure
      REAL :: CPCP       !Convective precipitation
      REAL :: TPTP       !Total precipitation
      REAL :: LWDN       !Longwave downward radiation
      REAL :: SWDN       !Shortwave downward radiation
      REAL :: WND        !Windspeed

!=== LDAS-CATCHMENT States ===============================================
      REAL :: TC(3)      !Canopy temperature
      REAL :: TSURF      !Land surface temperature
      REAL :: TSNOW      !Snow surface temperature
      REAL :: QA(3)      !Specific humidity in the canopy
      REAL :: INT        !Inteception storage
      REAL :: GHT(6)     !Ground heat content
      REAL :: WESN(3)    !Snow water equivalent
      REAL :: HTSN(3)    !Snow heat content
      REAL :: SNDZ(3)    !Snow depth
      REAL :: CATDEF     !Catchment defecit
      REAL :: RZEXC      !Root zone excess
      REAL :: SRFEXC     !Surface excess

!=== LDAS-CATCHMENT Output ===============================================
      REAL :: EVAP       !Evaporation
      REAL :: EINT       !Sensible heat flux
      REAL :: EVEG       !Transpiration
      REAL :: ESOI       !Bare soil evaporation
      REAL :: ESNO       !Snowpack evaporation
      REAL :: SHFLX      !Sensible heat flux
      REAL :: LHFLX      !Latent heat flux
      REAL :: GHFLX      !Ground heat flux
      REAL :: LWUP       !Longwave upward radiation
      REAL :: INFIL      !Infiltration
      REAL :: AR(3)      !Catchment moisture regime area fractions
      REAL :: RZEQ       !Root zone equilibrium moisture storage
      REAL :: SRFMC      !Surface moisture content
      REAL :: RZMC       !Root zone moisture content
      REAL :: COLMC      !Column moisture content
      REAL :: RUNOFF     !Total runoff
      REAL :: RUNSRF     !Surface runoff
      REAL :: BFLOW      !Base flow
      REAL :: SMELT      !Snow melt
      REAL :: SNOWTERM   ! ?
      REAL :: ASNOW      !Fractional coverage of snow
      REAL :: TS(6)      !Soil temperaure
      REAL :: MODALB     !Albedo accoounting for snow
      REAL :: ASSV(11)   !Assimilation variables for covariance forecasting

!=== LDAS-CATCHMENT Output Averages ======================================
      INTEGER :: COUNT   !Number of Timesteps
      REAL :: SUMTP      !Total Precipitation
      REAL :: SUMSNOW    !Total Snowfall
      REAL :: SUMQA      !Canopy Humidity
      REAL :: SUMTC      !Canopy Temperature
      REAL :: SUMTS6     !Deep Soil Temperature
      REAL :: SUMSWUP    !Outgoing Shortwave Radiation
      REAL :: SUMLWUP    !Outgoing Longwave Radiation
      REAL :: SUMSHFLX   !Sensible Heat Flux
      REAL :: SUMGHFLX   !Ground Heat Flux
      REAL :: SUMLHFLX   !Latent Heat Flux
      REAL :: SUMCATDEF  !Catchment Deficit
      REAL :: SUMRZEXC   !Root Zone Excess
      REAL :: SUMSRFEXC  !Surface Excess
      REAL :: SUMSRFMC   !Surface Soil Moisture Content
      REAL :: SUMRZMC    !Root Zone Soil Moisture Content
      REAL :: SUMCOLMC   !Column Average Soil Moisture Content
      REAL :: SUMAR1     !Saturated Fraction
      REAL :: SUMAR3     !Wilting Fraction
      REAL :: SUMINT     !Interception Depth
      REAL :: SUMEINT    !Interception Loss
      REAL :: SUMEVAP    !Evaporation
      REAL :: SUMEVEG    !Transpiration
      REAL :: SUMESOI    !Evaporation From Bare Soil
      REAL :: SUMSNDZ    !Snow Depth
      REAL :: SUMSMELT   !Snowmelt
      REAL :: SUMINFIL   !Infiltration
      REAL :: SUMBFLOW   !Baseflow
      REAL :: SUMRUNSRF  !Overland Flow
      REAL :: SUMRUNOFF  !Total Runoff

      end type
      end module cat_module





















