
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C NOAH LAND-SURFACE MODEL, UNCOUPLED 1-D COLUMN: VERSION 2.5 OCT 2001
C
C THIS MAIN PROGRAM AND ITS FAMILY OF SUBROUTINES COMPRISE VERSION 2.5
C OF THE PUBLIC RELEASE OF THE UNCOUPLED 1-D COLUMN VERSION OF THE
C "NOAH" LAND-SURFACE MODEL (LSM). THE NOAH LSM IS A DESCENDANT OF AN
C EARLIER GENERATION OF THE OREGON STATE UNIVERSITY (OSU) LSM, BUT IT
C INCLUDES SUBSTANTIAL PHYSICS EXTENSIONS AND RECODING ACCOMPLISHED
C ALONG THE WAY BY NCEP, HL (NWS), AFGWC, AND AFGL/AFPL/AFRL.  HENCE
C THE ACRONYM "NOAH" DENOTES N-NCEP, O-OSU, A-AIR FORCE, H-HYDRO LAB.
C ----------------------------------------------------------------------
C FOR DOCUMENTATION OF THIS CODE AND INSTRUCTIONS ON ITS EXECUTION AND
C INPUT/OUTPUT FILES, SEE "NOAH LSM USER'S GUIDE" IN FILE README_2.5
C IN THE SAME PUBLIC SERVER DIRECTORY AS THIS SOURCE CODE.
C ----------------------------------------------------------------------
C PROGRAM HISTORY LOG
C VERSION 1.0  --  01 MAR 1999
C VERSION 1.1  --  08 MAR 1999
C VERSION 2.0  --  27 JUL 1999
C VERSION 2.1  --  23 OCT 2000
C VERSION 2.2  --  07 FEB 2001
C VERSION 2.3  --  07 MAY 2001 = operational Eta implementation
C VERSION 2.4  --  27 JUN 2001 = ops Eta with NO physics changes
C VERSION 2.5  --  18 OCT 2001
C LDAS VERSION --  28 APR 2002 = NOAH Main added to LDAS Driver
C                                (NASA GSFC)
C VERSION 2.5.1--  28 MAY 2002 = Updated changes in NOAH LSM along
C                                with correction to SOILRZ and SOIL1M.  
C VERSION 2.5.2--  12 JUN 2002 = Updated to latest NOAL LSM Version 
C v2.5.2 Update--  03 SEP 2002 = Assigned RC and CCOND to be 0.0 when
C                                vegetative greenness fraction is 0.0 
C v2.5.2 Update--  04 NOV 2002 = Added TBOT as 2-D field (all res.)
C v2.5.2 Update--  07 JAN 2003 = Corrected T1 variable  
C VERSION 2.6  --  24 JUN 2003 = Updated Noah LSM v2.6, and corrected
C                                output variable, NOAH%RETURN(8)  
C
C   Physics changes:
C     in SUBROUTINE SFLX change CSOIL from 1.26E+6 to 2.00E+6
C     in SUBROUTINE SFLX change ZBOT from -3.0 to -8.0
C     Replaced de-bugged SUBROUTINE TDFCND
C     Removed SUBROUTINE REDPRM and moved the parameters to other
C      locations throughout noah_main and noah_physics subroutines. 
C VERSION 2.5.2  --  31 MAY 2002
C     fix in FUNCTION DEVAP related to FX calculation
C VERSION 2.6    --  Includes changes to certain parameters and 
C                     snow-soil physics
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NOAH_MAIN(NTL,LDAS,TILE,NOAH)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      USE noah_module      ! NOAH tile variables

      IMPLICIT NONE

      type (ldasdec) LDAS
      type (tiledec) TILE
      type (noahdec) NOAH

C=== Local Variables =====================================================

      INTEGER, PARAMETER:: NSOLD=20   ! Maximum number of soil layers
      INTEGER K          ! Loop integer

      INTEGER ICE         ! flag to invoke sea-ice physics (0=land)
      INTEGER NROOT       ! Number of root layers, dep. on veg. type 
      INTEGER NSOIL       ! Number of soil layers (used in noah_main.f)
      INTEGER NSLAY       ! Number of soil layers (from LDAS driver)     
      INTEGER SOILTYP     ! Soil type (integer index)
      INTEGER VEGTYP      ! Vegetation type (integer index)
      INTEGER SIBVEG      ! SIB Vegetation class type (read into noah_main.f)
      INTEGER IJ          ! Elements of array by number of soil layers
      INTEGER NTL         ! Current tile NOAH model is called on, can be removed
      INTEGER SLOPETYP    ! Class of surface slope (integer index)
      REAL SLOPE          ! Slope Estimate of linear reservoir coefficient  

      REAL BETA           ! Ratio of actual/potential EVAP (dimensionless) 
      REAL DRIP           ! Excess canopy moisture (m)
      REAL EC             ! Canopy evaporation (m s-1)
      REAL EDIR           ! Direct soil evaporation (m s-1)
      REAL ET(NSOLD)      ! Plant transp. from each root/soil layer(m s-1)
      REAL ETT            ! Accum plant transpiration (m s-1)
      REAL ESNOW          ! Sublimation from snowpack (kg m-2 s-1) (or dep.)
      REAL F              ! Net flux (TOT downward radiation)
      REAL FXEXP          ! Bare soil evaporation exponent used in DEVAP
      REAL FLX1           ! Precip-snow Surface  (W M-2)
      REAL FLX2           ! Freezing rain latent heat flux (W M-2)
      REAL FLX3           ! Phase-change heat flux from snowmelt (W M-2)
      REAL DEW            ! Dewfall amount  (m s-1)
      REAL RIB            ! Bulk Richardson Number
      REAL RUNOFF1        ! Ground surface runoff (m/s)
      REAL RUNOFF2        ! Underground runoff (m/s)
      REAL RUNOFF3        ! Runoff within soil layers (m/s)
      REAL SIGMA          ! Stefan-Boltzmann constant (5.672E-8 W m-2 T-4)
      REAL Q1             ! Effective mixing ratio at grnd sfc (kg kg-1)
      REAL AET            ! Actual evapotranspirative energy (J m-2 s-1)
      REAL ALB            ! Quarterly Snow-free albedo (fraction)
      REAL ALBEDO         ! Total Sfc. albedo fraction (incl. snow-max alb)
      REAL CH             ! Sfc. exchange coef. for heat and moisture
      REAL CM             ! Sfc. exhcange coef. for momentum
      REAL CMC            ! initial canopy water content (m)
      REAL CCOND          ! Canopy Conductance (m s-1) 
      REAL DQSDT          ! Slope of sat. specific hum. curve for Penman
      REAL CZIL           ! (Zilintinkevich coef)
      REAL REFKDT         ! Scalar surface runoff parameter 
      REAL REFDK          ! A reference value for KDT
      REAL KDT            ! Based on REFKDT, DKSAT, and REFDK 
      REAL FRZK           ! Ice content threshold in soil 
      REAL FRZX           ! Adjust FRZK Parameter to actual soil type 
      REAL FRZFACT        ! Used in ice content in soil 
      REAL DQSDT2         ! Slope of sat. specific hum. curve for Penman
      REAL DT             ! time step length (sec) of physical integration
      REAL EMISS          ! Longwave "emission" = epsilon
      REAL ESAT           ! Saturation vapor pressure for water (Pa)
      REAL QSAT           ! Saturation Specific humidity (Kg/Kg)
      REAL E              ! Saturation vapor pressure

      REAL ETA            ! Actual latent heat flux (W m-2) (NEG: up from sfc)
      REAL EVP            ! Total Evaporation (kg/m2) (=ETA/2.501E6)
      REAL ETP            ! Final Potential Evapotransp. (W m-2)
      REAL FUP            ! Upward ground LW radiation (W m-2)
      REAL SHTFLX         ! Sensible heat flux (W m-2)
      REAL SOLDN          ! Solar downward radiation (W m-2; positive,
                          !   not net shortwave)
      REAL LWDN           ! Downward Longwave radiation (W/m2)
      REAL PRCP           ! Precipitation rate conversion (kg m-2 s-1)
      REAL CPCP           ! Convective Precipitation (kg m-2 s-1)
      REAL ASNOW          ! The snowfall rate (kg m-2 s-1)
      REAL PTU            ! Phota Thermal Unit (dep. on air temp and rad)

      REAL Q2             ! Mixing ratio at 1st middle level above skin
      REAL Q2SAT          ! Sat. Mix. ratio at 1st middle level above skin
      REAL RES            ! Energy Balance equation residual (W m-2)
      REAL RNET           ! net radiation (W/m2)
      REAL SFCSPD         ! Wind speed, sqrt(u*u+v*v)
      REAL UWIND          ! U-Wind component (m s-1)
      REAL VWIND          ! V-Wind component (m s-1)
      REAL SFCPRS         ! Surface pressure (Pascals)
      REAL SFCTMP         ! Surface temperature (K)
      REAL SHDFAC         ! 12-month green vegetation fraction
      REAL SHDMIN         ! Fixed minimum green veg fraction

C      REAL MLAI(13)       ! Monthly leaf area index values

      REAL SMC(NSOLD)     ! Initial vol. soil moisture of each layer
      REAL SNOMLT         ! Snow melt (m) (water equivalent)
      REAL SNOALB         ! Maximum albedo expected over deep snow
      REAL MSTAVRZ        ! Avail root zone soil moisture (unitless fraction)
      REAL MSTAVTOT       ! Total Column Soil Moisture Availability (%)
      REAL SOILM          ! Total soil column water content (m->kg/m2)
      REAL SOILRZ         ! Root zone soil column water content (kg/m2)
      REAL SOIL1M         ! Top 1-m soil column water content (kg/m2)
      REAL STC(NSOLD)     ! Initial soil temperature (K) of each soil layer
      REAL GFLX           ! Soil heat flux (W m-2)
      REAL T1             ! Initial skin temperature (K) 
      REAL T14            ! Ground sfc. temp. to the 4th power (K^+4)
      REAL T1V            ! Virtual temperature at ground (sub 1)
      REAL T2V            ! Virtual temp. at 1st mid. lev. above grnd (2)
C LW    REAL TAK            ! Air temperature in K
      REAL TBOT           ! Annually-fixed, soil-temp condition at ZBOT
      REAL TH2            ! Potential temperature (K) at hgt z above grnd
      REAL TH2V           ! Virtual potential temperature (K) at z
      REAL Z              ! Height (meters) above ground in atmos. forcing
      REAL Z0             ! Roughness length parameter

      REAL PCPDAY         ! Daily accumulation of precipitation (mm)
      REAL PCPSUM         ! Daily accumulation of precipitation (mm)
      REAL ETADAY         ! Daily accumulation of actual evapotranp. (mm)
      REAL ETSUM          ! Daily accumulation of actual evapotranp. (mm)
      REAL RUNOFFSUM      ! Accum. daily values in mm for runoff (mm)
      REAL RUNOFFDAY      ! Accum. daily values in mm for runoff
      REAL SMCNW          ! SLDPTH(ij)*1000.*SMC(ij) + SMCNW

      REAL SH2O(NSOLD)    ! Initial vol. liq. soil moisture of each layer
      REAL SLDPTH(4)      ! Thickness values for NSOIL layers (meters)
      REAL SNOWH          ! Snow depth (meters)
      REAL SNEQV          ! Snow water-equivalent (m) above snowdepth
      REAL SNCOVR         ! FRACTIONAL SNOW COVER (UNITLESS FRACTION, 0-1)

      REAL RSMIN          ! MINIMUM CANOPY RESISTANCE (S M-1)
      REAL RGL            ! FROM SOLAR RAD TERM OF CANOPY RESISTANCE FUNCTION
      REAL HS             ! USED IN VAPOR PRESS. DEF. TERM OF CAN. RES. FUNCTION
      REAL SNUP           ! THRESHOLD SNOW DEPTH (IN WATER EQUIVALENT M) THAT
                          !  IMPLIES 100% SNOW COVER
      REAL XLAI           ! LEAF AREA INDEX (DIMENSIONLESS)
      REAL RC             ! CANOPY RESISTANCE (S M-1)
      REAL PC             ! PLANT COEFFICIENT (UNITLESS FRACTION, 0-1) WHERE
                          !   PC*ETP = ACTUAL TRANSP
      REAL RCS            ! INCOMING SOLAR RC FACTOR (DIMENSIONLESS)
      REAL RCT            ! AIR TEMPERATURE RC FACTOR (DIMENSIONLESS)
      REAL RCQ            ! ATMOS VAPOR PRESSURE DEFICIT RC FACTOR (DIMENSIONLESS)
      REAL RCSOIL         ! SOIL MOISTURE RC FACTOR (DIMENSIONLESS)

      REAL SMCMAX         ! MAX SOIL MOISTURE CONTENT (POROSITY)
      REAL SMCREF         ! REFERENCE SOIL MOISTURE (ONSET OF SOIL MOISTURE
                          !  STRESS IN TRANSPIRATION) - Volumetric
      REAL SMCREF1
      REAL SMCWLT         ! WILTING Point (VOLUMETRIC)
      REAL SMCWLT1
      REAL SMCDRY         ! DRY SOIL MOIST THRESHOLD - DIRECT EVAP 
                          !  FRM TOP LAYER ENDS (VOLUMETRIC)
      REAL PSISAT         ! SATURATED SOIL POTENTIAL
      REAL DKSAT          ! SATURATED SOIL HYDRAULIC CONDUCTIVITY
      REAL BEXP           ! THE 'B' PARAMETER
      REAL DWSAT          ! SATURATED SOIL DIFFUSIVITY
      REAL F1             ! USED TO COMPUTE SOIL DIFFUSIVITY/CONDUCTIVITY
      REAL QUARTZ         ! SOIL QUARTZ CONTENT
      REAL SMLOW          ! Soil moisture wilt and reference parameter
      REAL SMHIGH         ! Soil moisture wilt and reference parameter

      REAL, PARAMETER:: R = 287.04
      REAL, PARAMETER:: CP = 1004.5

      REAL, PARAMETER:: T0 = 273.15       ! Freezing point in Kelvin (273.15 K)
      REAL, PARAMETER:: LVH2O = 2.501000E+6 ! Latent heat for evapo for water  

      REAL, PARAMETER:: EPS = 0.622 ! Water/(dry air) molec mass ratio (epsilon)

!=== End Variable Definition ==========================================

!       print *, ' ** TILE #: ',NTL
!       print *, 'Lat/Lon of Tile: ',TILE%LAT,TILE%LON

!=== Convert LDAS Timestep varname to NOAH timestep varname (DT) (sec)

      DT = LDAS%TS

!=== Bottom Temperature Field

      TBOT = NOAH%TEMPBOT(1)

! THE FOLLOWING SECTIONS BREAKS DOWN THE THREE PARAMETER ARRAYS
! INTO THEIR INDIVIDUAL VARIABLE NAMES FOR LATER USAGE

!=== VEGETATION CLASS

       VEGTYP = TILE%VEGT 

C=== STATIC VEGETATION PARAMETERS
C ----------------------------------------------------------------------
C VEGETATION PARAMETERS ARE DEPENDENT ON VEGETATION TYPE (INDEX)
C   SHDFAC: VEGETATION GREENNESS FRACTION
C   RSMIN:  MIMIMUM STOMATAL RESISTANCE
C   RGL:    PARAMETER USED IN SOLAR RAD TERM OF
C           CANOPY RESISTANCE FUNCTION
C   HS:     PARAMETER USED IN VAPOR PRESSURE DEFICIT TERM OF
C           CANOPY RESISTANCE FUNCTION
C   SNUP:   THRESHOLD SNOW DEPTH (IN WATER EQUIVALENT M) THAT
C           IMPLIES 100% SNOW COVER
C   Z0:     Roughness Length
C   XLAI:   Leaf Area Index
C ----------------------------------------------------------------------
C SSIB VEGETATION TYPES (DORMAN AND SELLERS, 1989; JAM)
C  1:  BROADLEAF-EVERGREEN TREES  (TROPICAL FOREST)
C  2:  BROADLEAF-DECIDUOUS TREES
C  3:  BROADLEAF AND NEEDLELEAF TREES (MIXED FOREST)
C  4:  NEEDLELEAF-EVERGREEN TREES
C  5:  NEEDLELEAF-DECIDUOUS TREES (LARCH)
C  6:  BROADLEAF TREES WITH GROUNDCOVER (SAVANNA)
C  7:  GROUNDCOVER ONLY (PERENNIAL)
C  8:  BROADLEAF SHRUBS WITH PERENNIAL GROUNDCOVER
C  9:  BROADLEAF SHRUBS WITH BARE SOIL
C 10:  DWARF TREES AND SHRUBS WITH GROUNDCOVER (TUNDRA)
C 11:  BARE SOIL
C 12:  CULTIVATIONS (THE SAME PARAMETERS AS FOR TYPE 7)
C 13:  GLACIAL (THE SAME PARAMETERS AS FOR TYPE 11)
C ----------------------------------------------------------------------

       NROOT = NOAH%VEGP(1)
       RSMIN = NOAH%VEGP(2)
       RGL   = NOAH%VEGP(3)
       HS    = NOAH%VEGP(4)
       SNUP  = NOAH%VEGP(5)
       Z0    = NOAH%VEGP(6)
       XLAI  = NOAH%VEGP(7)
 
C=== MONTHLY VEGETATION PARAMETERS

       SHDFAC = NOAH%VEGIP(1)

!    Minimum greenness fraction
       SHDMIN=0.0

!    If ground surface is bare:
       IF (VEGTYP .EQ. 11) SHDFAC = 0.0

C=== STATIC SOIL PARAMETERS
C ----------------------------------------------------------------------
C  SOIL PARAMETERS ARE DEPENDENT ON SOIL TYPE (INDEX)
C    SMCMAX: MAX SOIL MOISTURE CONTENT (POROSITY)
C    SMCREF: REFERENCE SOIL MOISTURE (ONSET OF SOIL MOISTURE
C             STRESS IN TRANSPIRATION)
C    SMCWLT: WILTING PT SOIL MOISTURE CONTENT
C    SMCDRY: AIR DRY SOIL MOIST CONTENT LIMITS
C    PSISAT: SATURATED SOIL POTENTIAL
C    DKSAT:  SATURATED SOIL HYDRAULIC CONDUCTIVITY
C    BEXP:   THE 'B' PARAMETER
C    DWSAT:  SATURATED SOIL DIFFUSIVITY
C    F1:     USED TO COMPUTE SOIL DIFFUSIVITY/CONDUCTIVITY
C    QUARTZ: SOIL QUARTZ CONTENT
C ----------------------------------------------------------------------
C SOIL TYPES   ZOBLER (1986)      COSBY ET AL (1984) (quartz cont.(1))
C  1        COARSE            LOAMY SAND         (0.82)
C  2        MEDIUM            SILTY CLAY LOAM    (0.10)
C  3        FINE              LIGHT CLAY         (0.25)
C  4        COARSE-MEDIUM     SANDY LOAM         (0.60)
C  5        COARSE-FINE       SANDY CLAY         (0.52)
C  6        MEDIUM-FINE       CLAY LOAM          (0.35)
C  7        COARSE-MED-FINE   SANDY CLAY LOAM    (0.60)
C  8        ORGANIC           LOAM               (0.40)
C  9        GLACIAL LAND ICE  LOAMY SAND         (NA using 0.82)
C ----------------------------------------------------------------------
       
       NSOIL  =  4             ! 4 soil layers in NOAH 

       SLDPTH(1) = 0.1         ! Soil layer thicknesses (m)
       SLDPTH(2) = 0.3
       SLDPTH(3) = 0.6
       SLDPTH(4) = 1.0

       SOILTYP = NOAH%ZOBSOIL(1)    ! Zobler Soil Class Value
!       print * ,' Zobler Class for tile: ',NOAH%ZOBSOIL(1)

       SMCMAX =  NOAH%SOILP(1)
       PSISAT =  NOAH%SOILP(2)
       DKSAT  =  NOAH%SOILP(3)
       BEXP   =  NOAH%SOILP(4)
       QUARTZ =  NOAH%SOILP(5)

!   The following 5 parameters are just given here for reference
!    and to force static storage allocation.

!       SMCREF =  NOAH%SOILP(6)
!       SMCWLT =  NOAH%SOILP(7)
!       SMCDRY =  NOAH%SOILP(8)
!       DWSAT  =  NOAH%SOILP(9)
!       F1     =  NOAH%SOILP(10)

!-- Here is where the above five parameters are actually derived.
!    SET TWO SOIL MOISTURE WILT, SOIL MOISTURE REFERENCE PARAMETERS
       SMLOW = 0.5
C     changed in 2.6 from 3 to 6 on June 2nd 2003
!       SMHIGH = 3.0
       SMHIGH = 6.0

       DWSAT  = BEXP * DKSAT * (PSISAT/SMCMAX)
       F1     = ALOG10(PSISAT) + BEXP*ALOG10(SMCMAX) + 2.0
       SMCREF1 = SMCMAX*(5.79E-9/DKSAT)**(1.0/(2.0*BEXP+3.0))
       SMCREF = SMCREF1 + (SMCMAX-SMCREF1) / SMHIGH
       SMCWLT1 = SMCMAX * (200.0/PSISAT)**(-1.0/BEXP)
       SMCWLT = SMCWLT1 - SMLOW * SMCWLT1
!    Current version SMCDRY values equate to SMCWLT
       SMCDRY = SMCWLT

C ----------------------------------------------------------------------
C KDT IS DEFINED BY REFERENCE REFKDT AND DKSAT; REFDK=2.E-6 IS THE SAT.
C DK. VALUE FOR THE SOIL TYPE 2
C ----------------------------------------------------------------------
       REFDK=2.0E-6
       REFKDT=3.0
       KDT = REFKDT * DKSAT/REFDK

C ----------------------------------------------------------------------
C FROZEN GROUND PARAMETER, FRZK, DEFINITION: ICE CONTENT THRESHOLD ABOVE
C WHICH FROZEN SOIL IS IMPERMEABLE REFERENCE VALUE OF THIS PARAMETER FOR
C THE LIGHT CLAY SOIL (TYPE=3) FRZK = 0.15 M.
C ----------------------------------------------------------------------
      FRZK=0.15
C ----------------------------------------------------------------------
C TO ADJUST FRZK PARAMETER TO ACTUAL SOIL TYPE: FRZK * FRZFACT
C ----------------------------------------------------------------------
      FRZFACT = (SMCMAX / SMCREF) * (0.412 / 0.468)
      FRZX = FRZK * FRZFACT

C=== SLOPE TYPE
C ----------------------------------------------------------------------
C CLASS PARAMETER 'SLOPETYP' WAS INCLUDED TO ESTIMATE LINEAR RESERVOIR
C COEFFICIENT 'SLOPE' TO THE BASEFLOW RUNOFF OUT OF THE BOTTOM LAYER.
C LOWEST CLASS (SLOPETYP=0) MEANS HIGHEST SLOPE PARAMETER = 1.
C DEFINITION OF SLOPETYP FROM 'ZOBLER' SLOPE TYPE:
C SLOPE CLASS  PERCENT SLOPE
C 1            0-8
C 2            8-30
C 3            > 30
C 4            0-30
C 5            0-8 & > 30
C 6            8-30 & > 30
C 7            0-8, 8-30, > 30
C 8            GLACIAL ICE
C 9            OCEAN/SEA
C ----------------------------------------------------------------------

!--  SLOPETYP = 3
      SLOPE = 1.0

C=== MONTHLY (QUARTERLY, for now) ALBEDO (SNOW-FREE)

      ALB =  NOAH%ALBSF(1)

!    Maximum Albedo over very Deep Snow

      SNOALB = NOAH%MXSNALB(1)

C=== THE FOLLOWING BREAKS DOWN THE FORCING VARIABLES

      SFCTMP = TILE%FORCING(1)
      Q2     = TILE%FORCING(2)
      SOLDN  = TILE%FORCING(3)
      LWDN   = TILE%FORCING(4)
      UWIND  = (TILE%FORCING(5))*(TILE%FORCING(5))
      VWIND  = (TILE%FORCING(6))*(TILE%FORCING(6))
       SFCSPD = SQRT( UWIND + VWIND )
      SFCPRS = TILE%FORCING(7)
      PRCP   = TILE%FORCING(8)
      CPCP   = TILE%FORCING(9)   !Convective Precipitation (kg/m2sec)

C-- Height of observations (this needs to be modified)
      Z = 6.0       ! Height of observations (m)

C-- Prevent Numerical Instability for Wind Speed
      if(SFCSPD.le.0.01) SFCSPD=0.01

C-- Prevent Numerical Instability with HUMIDITY

      IF (Q2 .LT. 0.1E-5)  Q2 = 0.1E-5

CCC Calculate Saturation Specific Humidity (Kg/Kg) and
C    Saturation vapor pressure for water (Pa) based on Specific
C    Humidity(Kg/Kg), Temperature(K), and Pressure (Pa)
C
C      FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989: 'A
C      SHORT COURSE IN CLOUD PHYSICS', PERGAMON PRESS, 3rd ED.
C        Pablo J. Grunmann, 3/6/98.
C
C      QSAT  = Saturation Specific humidity (Kg/Kg)
C      ESAT  = Saturation vapor pressure for water (Pa)
C      EPS   = Water/(dry air) molecular mass ratio (epsilon)
C      E     = Saturation vapor pressure
C
C-- Function E(SFCTMP) = Sat. vapor pressure (in Pascal) at
C                   temperature T (uses Clausius-Clapeyron).

       ESAT = E(SFCTMP)

C-- CALCULATE SATURATION MIXING RATIO (PABLO GRUNMANN, 05/28/98)

       Q2SAT = 0.622 * ESAT /(SFCPRS - (1.-0.622)*ESAT)
       IF (Q2 .GE.  Q2SAT)  Q2 = Q2SAT*0.99

C-- CALCULATE SLOPE OF SAT. SPECIFIC HUMIDITY CURVE FOR PENMAN: DQSDT2

       DQSDT2 = DQSDT (SFCTMP, SFCPRS)

C-- CALC VIRTUAL TEMPS AND POTENTIAL TEMPS AT GRND (SUB 1) AND AT
C    THE 1ST MDL LVL ABV THE GRND (SUB 2). EXPON IS CP DIVD BY R.

       TH2 = SFCTMP + ( 0.0098 * Z )
       T2V = SFCTMP * (1.0 + 0.61 * Q2 )

       T1V  =  NOAH%T1 * (1.0 + 0.61 * Q2 )
       TH2V = TH2 * (1.0 + 0.61 * Q2 )

C==  OPTIONAL SUBROUTINE:  Calculate LW Radiation (Down)  =================
!      CALL OBTLWDN(SFCTMP,LWDN)

C--  Phota Thermal Unit (PTU)

!       PTU    = 0.0    ! Initial Value
       PTU    = 0.10

C--  Initialize SOILM for 1st timestep water balance
!       SOILM = 0.0
C--  Initialize ROOT ZONE COLUMN SOIL MOISTURE IN METERS (SOILRZ)
        SOILRZ = 0.0
C--  Initialize TOP 1-METER COLUMN SOIL MOISTURE IN METERS (SOIL1M)
        SOIL1M = 0.0

C==  CALCULATE CH (EXCHANGE COEFFICIENT)  =================================
C
C   CH IS THE SFC EXCHANGE COEFFICIENT FOR HEAT/MOISTURE
C   CM IS THE SFC EXCHANGE COEFFICIENT FOR MOMENTUM
C
C IMPORTANT NOTE: TO CALCULATE THE SFC EXCHANGE COEF (CH) FOR HEAT AND
C                 MOISTURE, SUBROUTINE SFCDIF BELOW CAN:
C
C    A)  BE CALLED HERE FROM THE DRIVER, THUS CH IS INPUT TO SFLX
C           (AS IS TYPICAL IN A COUPLED ATMOSPHERE/LAND MODEL), OR
C  * B)  BE CALLED INTERNALLY IN ROUTINE SFLX (THUS CH IS OUTPUT FROM SFLX),
C          BEFORE THE CALL TO ROUTINE "PENMAN"
C
C   OPTION B IS THE DEFAULT HERE.  THAT IS, IN THE UNCOUPLED, OFF-LINE LSM
C   REPRESENTED HEREIN BY THIS DRIVER, WE CALL SFCDIF LATER IN ROUTINE SFLX.
C
C   THE ROUTINE SFCDIF REPRESENTS THE SO-CALLED "SURFACE LAYER" OR THE
C   "CONSTANT FLUX LAYER" (THE LOWEST 20-100 M OF THE ATMOSPHERE).
C   HENCE ROUTINE SFCDIF EMBODIES THE "ATMOSPHERIC AERODYNAMIC RESISTANCE".
C
C   TO ENABLE THE FLEXIBILITY OF EITHER OPTION A OR B, WE PASS
C   THE ARGUMENTS "CH", "CM", AND "SFCSPD" (WIND SPEED:JUST CALCULATED ABOVE)
C   TO ROUTINE SFLX TO SUPPORT OPTION B -- THAT IS, FOR INPUT TO THE CALL TO
C   ROUTINE SFCDIF THEREIN.  IN OPTION A, THE ARGUMENTS "SFCSPD" AND "CM"
C   ARE NEITHER NEEDED IN ROUTINE SFLX, NOR ALTERED BY ROUTINE SFLX.
C
C     IF ONE CHOOSES OPTION A, THEN ONE MUST
C      1 - ACTIVATE (UNCOMMENT) THE CALL TO SFCDIF BELOW,
C      2 - ACTIVATE (UNCOMMENT) THE ASSIGNMENT OF "Z0" AND "CZIL" NEXT BELOW
C      3 - DE-ACTIVATE (COMMENT OUT) THE CALL TO SFCDIF IN ROUTINE SFLX.
C
C   Z0 and CZIL:
C
C   THE ROUGHNESS LENGTH PARAMETERS "Z0" AND "CZIL" MUST BE SET HERE IN THE
C   DRIVER TO SUPPORT THE "OPTION-A", I.E. THE CALL TO SFCDIF BELOW.  IN SO
C   DOING, THE "Z0" AND "CZIL" ASSIGNED HERE MUST CORRESPOND TO THEIR
C   VALUES, CALLED BY SFLX JUST BEFORE CALL SFCDIF.
C   THUS THE VALUE OF "Z0" ASSIGNED HERE MUST CORRESPOND TO THAT ASSIGNED
C   FOR THE CHOSEN VEG CLASS THAT WAS ALREADY INPUT EARLIER IN THE DRIVER.
C
C   BECAUSE OF THE IMPLICIT ITERATIVE NATURE OF THE "PAULSON" SURFACE-LAYER
C   SCHEME USED IN ROUTINE SFCDIF, CH AND CM ARE CO-DEPENDENT.  SIMILARLY,
C   THE IMPLICIT NATURE OF THE SFCDIF SCHEME ALSO REQUIRES THAT FOR EITHER
C   OPTION A OR B, CH AND CM MUST BE INITIALIZED EARLIER IN THE DRIVER BEFORE
C   THE START OF THE TIME-STEP LOOP, AS WELL AS BE CARRIED FORWARD FROM
C   TIME STEP TO TIME STEP AS "STATE VARIABLES", BECAUSE THE VALUES OF
C   CH AND CM FROM A PREVIOUS TIME STEP REPRESENT THE FIRST-GUESS VALUES FOR
C   THE CALL TO SFCDIF IN THE PRESENT TIME STEP.
C
C   SOME USERS MAY CHOOSE TO EXECUTE AN ENTIRELY DIFFERENT SCHEME IN PLACE OF
C   ROUTINE SFCDIF HERE, E.G. AN EXPLICIT SCHEME SUCH AS LOUIS (1979) THAT
C   EMPLOYS NO ITERATION AND HAS NO REQUIREMENT TO CARRY CH AND CM FORWARD
C   AS STATE VARIABLES FROM TIME STEP TO TIME STEP.  IN THAT CASE, IN
C   OPTION A, THE ROUTINE SHOULD BE CALLED HERE IN THE DRIVER AFTER ALL
C   NECESSARY INPUT ARGUMENTS FOR IT ARE DEFINED AT THIS POINT, OR CALLED IN
C   ROUTINE SFLX, AT THE POINT SFCDIF IS CALLED.
C
C      CALL SFCDIF ( Z, Z0, T1V, TH2V, SFCSPD,CZIL, CM, CH )
C
C ---------------------------------------------------------------------|
C INITIALIZE CH, CM (NOTE: initial these before time loop)
c      CH=1.E-4
c      CM=1.E-4
C 1998 May 22 0030 (Julian= 142) typical values initialization
!      CH=  0.0150022404
!      CM=  0.0205970779
! **NOTE: TRYING THESE VALUES AS TEST! 
C ---------------------------------------------------------------------|


C=== MAIN CALL TO LAND-SURFACE PHYSICS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

       CALL SFLX (NTL,
     C  ICE,DT,Z,NSOIL,SLDPTH,
     F  LWDN,SOLDN,SFCPRS,PRCP,SFCTMP,Q2,SFCSPD,
     I  TH2,Q2SAT,DQSDT2,
     S  SLOPE,SHDFAC,SHDMIN,PTU,ALB,SNOALB,
     S  RSMIN,RGL,HS,SNUP,Z0,XLAI,NROOT,
     S  PSISAT,BEXP,DKSAT,SMCMAX,QUARTZ,DWSAT,
     S  SMCWLT,SMCREF,SMCDRY,F1,KDT,FRZX,FRZFACT,TBOT,
     H  NOAH%CMC,NOAH%T1,NOAH%STC,NOAH%SMC,NOAH%SH2O,
     H  NOAH%SNOWH,NOAH%SNEQV,ALBEDO,NOAH%CH,NOAH%CM,
     O  EVP,ETA,SHTFLX,
     O  EC,EDIR,ET,ETT,ESNOW,DRIP,DEW,
     O  BETA,ETP,GFLX,
     O  FLX1,FLX2,FLX3,
     O  SNOMLT,SNCOVR,
     O  RUNOFF1,RUNOFF2,RUNOFF3,
     O  RC,PC,RCS,RCT,RCQ,RCSOIL,
     D  MSTAVRZ,MSTAVTOT,SOILM,
     M  LDAS,TILE)


C******************************************************************

C   CALCULATE UPWARD LONGWAVE RAD USING UPDATED SKIN TEMPERATURE

      T14 = NOAH%T1*NOAH%T1*NOAH%T1*NOAH%T1 
      FUP = 5.67E-8 * T14

C   CALCULATE RESIDUAL OF ALL SURFACE ENERGY BALANCE EQN TERMS.

!      GFLX = -GFLX
!      F = SOLDN*(1.0-ALBEDO) + LWDN
!      RES    = F - SHTFLX - GFLX - ETA - FUP - FLX1 - FLX2 - FLX3
!         ENDIF

c      PRINT*,'  --------------------------------------'
c      PRINT*,'  State Variables '
c      PRINT*,'  --------------------------------------'
c      WRITE(*,*) NOAH%T1,' T1...Skin temperature (K)'
c      WRITE(*,*)(NOAH%STC(IJ), IJ=1,NSOIL),' STC'
c      WRITE(*,*)(NOAH%SMC(IJ), IJ=1,NSOIL),' SMC'
c      WRITE(*,*)(NOAH%SH2O(IJ), IJ=1,NSOIL),' SH2O'
c      WRITE(*,*) NOAH%CMC,' CMC...Canopy water content (m)'
c      WRITE(*,*) NOAH%SNOWH,' SNOWH...Actual snow depth (m)'
c      WRITE(*,*) NOAH%SNEQV,' SNEQV...Water equiv snow depth (m)'
c      WRITE(*,*) 'CH= ',NOAH%CH,'   CM= ',NOAH%CM
c      PRINT*,'  --------------------------------------'

C=== Collect the output variables into NOAH%RETURN

       NOAH%RETURN(1)=SFCTMP
       NOAH%RETURN(2)=Q2
       NOAH%RETURN(3)=SFCPRS
       IF (SFCTMP .LE. T0) THEN
         NOAH%RETURN(4)=PRCP*LDAS%TS        ! SNOWFALL
       ELSE
         NOAH%RETURN(5)=PRCP*LDAS%TS        ! RAINFALL
       ENDIF
       NOAH%RETURN(6)=SOLDN
       NOAH%RETURN(7)=LWDN
       NOAH%RETURN(8)=(-1.0)*SOLDN*(1.0-ALBEDO)
       NOAH%RETURN(9)=(5.67E-8)*(NOAH%T1**4.0)-LWDN
       NOAH%RETURN(10)=TILE%FORCING(5)       ! UWIND
       NOAH%RETURN(11)=TILE%FORCING(6)       ! VWIND
       NOAH%RETURN(12)=CPCP*LDAS%TS

       NOAH%RETURN(13)=NOAH%CMC*1000.0
       NOAH%RETURN(14)=NOAH%T1
       NOAH%RETURN(15)=NOAH%STC(1)
       NOAH%RETURN(16)=NOAH%STC(2)
       NOAH%RETURN(17)=NOAH%STC(3)
       NOAH%RETURN(18)=NOAH%STC(4)
       NOAH%RETURN(19)=NOAH%SMC(1)*1000.0*SLDPTH(1)
       NOAH%RETURN(20)=NOAH%SMC(2)*1000.0*SLDPTH(2)
       NOAH%RETURN(21)=NOAH%SMC(3)*1000.0*SLDPTH(3)
       NOAH%RETURN(22)=NOAH%SMC(4)*1000.0*SLDPTH(4)
       NOAH%RETURN(23)=NOAH%SH2O(1)*1000.0*SLDPTH(1)
       NOAH%RETURN(24)=NOAH%SH2O(2)*1000.0*SLDPTH(2)
       NOAH%RETURN(25)=NOAH%SH2O(3)*1000.0*SLDPTH(3)
       NOAH%RETURN(26)=NOAH%SH2O(4)*1000.0*SLDPTH(4)
       NOAH%RETURN(27)=NOAH%SNOWH
       NOAH%RETURN(28)=NOAH%SNEQV*1000.0
!     CH is considered here as aerodynamic cond. (ACOND)
       NOAH%RETURN(29)=NOAH%CH      
       NOAH%RETURN(30)=NOAH%CM
       NOAH%RETURN(31)=ALBEDO
       NOAH%RETURN(32)=SHDFAC

       NOAH%RETURN(33)=ETA
       NOAH%RETURN(34)=SHTFLX
       NOAH%RETURN(35)=EC*1000.0*2.50100E+6
       NOAH%RETURN(36)=EDIR*1000.0*2.50100E+6
       NOAH%RETURN(37)=ETT*1000.0*2.50100E+6
       NOAH%RETURN(38)=ESNOW*2.50100E+6
       NOAH%RETURN(39)=EVP*LDAS%TS 
       NOAH%RETURN(40)=DRIP
       NOAH%RETURN(41)=DEW
       NOAH%RETURN(42)=BETA
       NOAH%RETURN(43)=ETP
       NOAH%RETURN(44)=GFLX
       NOAH%RETURN(45)=FLX1+FLX2+FLX3
       NOAH%RETURN(46)=SNOMLT*1000.0
       NOAH%RETURN(47)=SNCOVR
       NOAH%RETURN(48)=RUNOFF1*1000.0*LDAS%TS
       NOAH%RETURN(49)=RUNOFF2*1000.0*LDAS%TS
       IF(SHDFAC.EQ.0.0) THEN        
        NOAH%RETURN(50)=0.0        ! Canopy Resistence (RC)
        NOAH%RETURN(51)=0.0        ! Canopy conductance (CCOND)
       ELSE
        NOAH%RETURN(50)=RC
        NOAH%RETURN(51)=1.0/RC     ! CCOND
       ENDIF
       NOAH%RETURN(52)=MSTAVRZ
       NOAH%RETURN(53)=MSTAVTOT
       NOAH%RETURN(54)=SOILM*1000.0

C ROOT ZONE COLUMN SOIL MOISTURE IN METERS (SOILRZ)
       DO K = 1,NROOT
         SOILRZ = SOILRZ+(NOAH%SMC(K)*SLDPTH(K)*1000.0)
       END DO
       NOAH%RETURN(55)=SOILRZ

C TOP 1-METER COLUMN SOIL MOISTURE IN METERS (SOIL1M)
      DO K = 1,3
        SOIL1M = SOIL1M+(NOAH%SMC(K)*SLDPTH(K)*1000.0)
      END DO
      NOAH%RETURN(56)=SOIL1M 


C>>>  END OF NOAH_MAIN <<<<<<<<
      RETURN
      END


CC*** NOAH FUNCTIONS ****************************************************

      FUNCTION DQS (T)

      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  PURPOSE:  TO CALCULATE VALUES OF VAPOR PRESSURE (E)
CC            AND P * DQS/DT (P TIMES CHG IN SAT MXG RATIO WITH RESPECT
CC            TO THE CHG IN TEMP) IN SUBSTITUTION TO THE LOOK-UP TABLES.
CC
CC            FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989.
CC                         ADDED BY PABLO J. GRUNMANN, 6/30/97.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL DESDT
      REAL DQS
CK    REAL ESD
      REAL LW
      REAL T
      REAL ES

CK      REAL, PARAMETER:: CP = 1005.
CK      REAL, PARAMETER:: CV = 718.
CK      REAL, PARAMETER:: CVV = 1410.
      REAL, PARAMETER:: CPV = 1870.
      REAL, PARAMETER:: RV = 461.5
      REAL, PARAMETER:: CW = 4187.
      REAL, PARAMETER:: EPS = 0.622
      REAL, PARAMETER:: ESO = 611.2
      REAL, PARAMETER:: TO = 273.15
      REAL, PARAMETER:: LVH2O = 2.501000E+6


C     ABOUT THE PARAMETERS:
C
C     EPS ---------- WATER - DRY AIR MOLECULAR MASS RATIO, EPSILON
C
C   VALUES FOR SPECIFIC HEAT CAPACITY AND INDIVIDUAL GAS CONSTANTS
C   IN [JOULES/(KG*KELVIN)] UNITS.
C
C     DRY AIR:
C             CP, CV
C     WATER VAPOR:
C                 CVV = 1410.
C                 CPV = 1870.
C                 RV  =  461.5
C     LIQUID WATER:
C                  CW = 4187.
C
C     ESO = ES(T=273.15 K) = SAT. VAPOR PRESSURE (IN PASCAL) AT T=TO
C      TO = 273.15
C
C     SAT. MIXING  RATIO: QS ~= EPS*ES/P
C     CLAUSIUS-CLAPEYRON: DES/DT = L*ES/(RV*T^2)
C     @QS/@T =  (EPS/P)*DES/DT

          LW = LVH2O - ( CW - CPV ) * ( T - TO )
          ES = ESO*EXP (LW*(1/TO - 1/T)/RV)
          DESDT = LW*ES/(RV*T*T)

C    FOR INSERTION IN DQSDT FUNCTION:
C    DQSDT = DQS/P , WHERE DQS = EPS*DESDT

          DQS = EPS*DESDT

          RETURN
          END

C----------------------------------------------------------------------

      FUNCTION DQSDT ( SFCTMP, SFCPRS )

      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    PURPOSE:  TO RETRIEVE THE APPROPRIATE VALUE OF DQSDT (THE CHANGE
CC    =======   OF THE SATURATION MIXING RATIO WITH RESPECT TO THE
CC              CHANGE IN TEMPERATURE) FROM:
CC
CC              FORMULAS INTRODUCED IN FUNCTION DQS
CC                                  (MODIFIED BY PABLO GRUNMANN, 7/9/97).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL SFCTMP
      REAL SFCPRS
      REAL DQS
      REAL DQSDT

      IF ((SFCTMP .GE. 173.0) .AND. (SFCTMP  .LE.  373.0)) THEN

C  IF THE INPUT SFC AIR TEMP IS BTWN 173 K AND 373 K, USE
C   FUNCTION DQS TO DETERMINE THE SLOPE OF SAT.MIX RATIO FUNCTION

        DQSDT = DQS (SFCTMP) / SFCPRS

      ELSE

C  OTHERWISE, SET DQSDT EQUAL TO ZERO

        DQSDT = 0.0

      END IF

      RETURN
      END

C---------------------------------------------------------------

      FUNCTION E(T)

      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  PURPOSE:  TO CALCULATE VALUES OF SAT. VAPOR PRESSURE (E)
CC            FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989.
CC                         ADDED BY PABLO J. GRUNMANN, 7/9/97.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL LW
      REAL T
      REAL E

CK      REAL, PARAMETER:: EPS = 0.622 
CK      REAL, PARAMETER:: CP = 1005.
CK      REAL, PARAMETER:: CV = 718.
CK      REAL, PARAMETER:: CVV = 1410.
      REAL, PARAMETER:: CPV = 1870.
      REAL, PARAMETER:: RV = 461.5
      REAL, PARAMETER:: CW = 4187.
      REAL, PARAMETER:: ESO = 611.2
      REAL, PARAMETER:: TO = 273.15
      REAL, PARAMETER:: LVH2O = 2.501000E+6

C   ABOUT THE PARAMETERS:
C
C    EPS --- WATER - DRY AIR MOLECULAR MASS RATIO, EPSILON
C
C    VALUES FOR SPECIFIC HEAT CAPACITY AND INDIVIDUAL GAS CONSTANTS
C    IN [JOULES/(KG*KELVIN)] UNITS.
C
C     DRY AIR:
C             CP, CV
C     WATER VAPOR:
C             CVV = 1410.
C             CPV = 1870.
C             RV  =  461.5
C     LIQUID WATER:
C             CW = 4187.
C
C     ESO = ES(TO) = SAT. VAPOR PRESSURE (IN PASCAL) AT T=TO
C      TO = 273.15
C
C     CLAUSIUS-CLAPEYRON: DES/DT = L*ES/(RV*T^2)

          LW = LVH2O - ( CW - CPV ) * ( T - TO )
          E = ESO*EXP (LW*(1/TO - 1/T)/RV)

          RETURN
          END


CCCC 1. DRIVER SUBROUTINE ==> SUBROUTINE OBTLWDN CCCCCCCCCCCCCCCCC

C      SUBROUTINE OBTLWDN(SFCTMP,LWDN)

C                      RADIATION
C
C The following step (OBTENTION OF LWDN) is used if
C user wants to calculate longwave downward.
C
C OBTENTION OF LWDN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C COMPUTATION OF LWDN (INCOMING LW RADIATION) FROM TAIR AND Q:
C
C...............LWDN = EMISS*SIGMA*(TAK)^4.
C
C WHERE:   TAK = AIR TEMP IN KELVIN
C        EMISS = 0.7  OR  (IDSO AND JACKSON, 1969):
C
C        EMISS = (1 - 0.261 EXP(-7.77*10^(-4)X(273-TAK)^2)
C
C      NEED STEFAN-BOLTZMANN CONSTANT, SIGMA
C         SIGMA = 5.672 * 10^-8  W M^-2 T^-4
C
C           SIGMA = 5.672E-8
C           TAK = SFCTMP
C           EMISS = 1 - 0.261*EXP((-7.77E-4)*(273-TAK)^2.)
C
C           LWDN = EMISS*SIGMA*TAK^4.
C
C        RETURN
C        END

CCCCC  END OF DRIVER SUBROUTINES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


