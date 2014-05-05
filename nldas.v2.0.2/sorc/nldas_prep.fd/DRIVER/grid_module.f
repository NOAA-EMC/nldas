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
! grid_module.f: 
!
! DESCRIPTION:
!  LDAS non-model-specific grid variables only.
!
!  LDAS%FORCING() ARRAY:
!  1. T 2m    Temperature interpolated to 2 metres [K]
!  2. q 2m    Instantaneous specific humidity interpolated to 2 metres[kg/kg]
!  3. radswg  Downward shortwave flux at the ground [W/m^2]
!  4. lwgdwn  Downward longwave radiation at the ground [W/m^2]
!  5. u 10m   Instantaneous zonal wind interpolated to 10 metres [m/s]
!  6. v 10m   Instantaneous meridional wind interpolated to 10 metres[m/s]
!  7. ps      Instantaneous Surface Pressure [Pa]
!  8. preacc  Total precipitation [mm/s]
!  9. precon  Convective precipatation [mm/s]
! 10. albedo  Surface albedo (0-1)
! 11. sfctyp  Surface types
! 12. gwet    Soil Moisture (percent of field capacity)
! 13. snow    Snow depth (mm water equivalent)
! 14. T 10m   Temperature interpolated to 10 metres [K]
! 15. q 10m   Instantaneous specific humidity interpolated to 10 metres[kg/kg]
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial code
!  11 Apr 2000: Brian Cosgrove; Added Forcing Mask variables
!  23 Feb 2001: Urszula Jambor; Added GEOS & GDAS forcing variables
!  27 Feb 2001: Brian Cosgrove; Added Catchment forcing data variables
!  23 Mar 2001: Jon Radakovich; Added variables for PSAS assimilation
!  04 Sep 2001: Brian Cosgrove; Added variabes for humidity, precip,par
!               brightness temp,precip mask, removed awips2ldas and 
!               pinker2ldas variables, GRIB interp. package used now
!  15 Oct 2001: Jesse Meng; Revised doc block with forcing array definition
!  15 Oct 2001: Jesse Meng; Added oblwdata1 and oblwdata2
!  30 Jul 2002: Jon Gottschalck; Chnaged name of global observed precip array
!  20 Nov 2002: Jon Radakovich; Added CLM specific variables for assimilation
!               and bias correction
!=========================================================================

      MODULE grid_module 

      IMPLICIT NONE
      public griddec

      type griddec

!=== LDAS Non-Model-Specific GRID Variables ==============================
      REAL    :: MASK             !Real Land/water mask (1=land,0=water) FOR MODELING
      REAL    :: FMASK            !Real Land/water mask (1=US,2=Mexico,
                                  ! 3=Canada,0=water) FOR FORCING DATA
      INTEGER :: IMASK            !Integer Land/Water Mask (1=land,0=water) FOR MODELING
      INTEGER :: FIMASK           !Integer Land/Water Mask (1=US,2=Mexico,
                                  ! 3=Canada,0=water) FOR FORCING DATA
      REAL    :: LAT              !Latitude of grid point
      REAL    :: LON              !Longitude of grid point
      REAL, pointer :: FGRD(:)    !Fraction of vegetation class in grid     
      INTEGER, pointer :: PVEG(:) !Predominance of vegetation class in grid
      INTEGER :: SOILT            !Soil Type (may need DIM in future)
      REAL    :: PRECIPWEIGHT     !RAW weight applied to precip sources
                                  ! when reading in NLDAS precip data

 
!=== Forcing variables ===================================================
      REAL, pointer :: NCEPDATA1(:)    !Past Forcing Data
      REAL, pointer :: NCEPDATA2(:)    !Future Forcing Data
      REAL, pointer :: ETADATA1(:)     !Past Forcing Data from ETA sources
      REAL, pointer :: ETADATA2(:)     !Future Forcing Data from ETA sources
      REAL, pointer :: GLBDATA1(:)     !Past Global Forcing Data from GDAS or GEOS sources
      REAL, pointer :: GLBDATA2(:)     !Future Global Forcing Data from GDAS or GEOS sources
      REAL, pointer :: CATCHDATA1(:)   !Past Forcing Data from CATCHMENT sources
      REAL, pointer :: CATCHDATA2(:)   !Future Forcing Data from CATCHMENT sources
      REAL, pointer :: PRECIPDATA1(:)  !Past Forcing Data from ETA sources
      REAL, pointer :: PRECIPDATA2(:)  !Future Forcing Data from ETA sources
      REAL :: ETASW                    !Zenith interpolated ETA SW
      REAL :: CATCHSW                  !Zenith interpolated ETA SW  
      REAL :: OBSW                     !Observed shortwave from pinker or nesdis
      REAL :: OBSBT                    !Observed brightness temperature
      REAL :: PAR                      !Observed PAR from Pinker
      REAL :: PINKDATA1                !Past Observed Radiation Data
      REAL :: PINKDATA2                !Future Observed Radiation Data
      REAL :: PARDATA1                 !Past Observed PAR Data
      REAL :: PARDATA2                 !Future Observed PAR Data
      REAL :: BRTTMPDATA1              !Past Observed BRT TEMP Data
      REAL :: BRTTMPDATA2              !Future Observed BRT TEMP Data
      REAL :: NESDATA1                 !Past Observed Radiation Data
      REAL :: NESDATA2                 !Future Observed Radiation Data
      REAL :: OBSWDATA1                !Past Observed Radiation Data (global)
      REAL :: OBSWDATA2                !Future Observed Radiation Data (global)
      REAL :: OBLWDATA1                !Past Observed Radiation Data (global)
      REAL :: OBLWDATA2                !Future Observed Radiation Data (global)
      REAL :: MERGEPRECIP              !Current MERGED StageIV/CPC Precip
      REAL :: S4PRECIP                 !Current STAGEIV Precipitation Data
      REAL :: OBSPRECIP                !Forcing Data from observed precip sources
      REAL :: TEMPDATA1                !Past 0-30mb Temp Data from ETA6 sources
      REAL :: TEMPDATA2                !Future 0-30mb Temp Data from ETA6 sources
      REAL :: PRESDATA1                !Past 0-30mb Pres Data from ETA6 sources
      REAL :: PRESDATA2                !Future 0-30mb Pres Data from ETA6 sources
      REAL :: QDATA1                   !Past 0-30mb Q Data from ETA6 sources
      REAL :: QDATA2                   !Future 0-30mb Q Data from ETA6 sources
      REAL :: RADBCDATA1               !Bias Correction Ratio Past Data
      REAL :: RADBCDATA2               !Bias Correction Ratio Future Data

      REAL, pointer :: FORCING(:)      !Interpolated LDAS Forcing Array
      REAL, pointer :: PRECIP(:)
      REAL, pointer :: ETAPRECIP(:)

!=== Output arrays for forcing variables
      REAL :: TOTALTMP,TOTALSPFH,TOTALDSWRF     !Temperature, Specific Humidity and Solar Forcing
      REAL :: TOTALDLWRF,TOTALUGRD,TOTALVGRD    !Longwave forcing, u winds and v winds
      REAL :: TOTALPRES,TOTALAPCP,TOTALACPCP    !Surface pressure, total and convective precipitation
      REAL :: TOTALOAPCP,TOTALODSWRF            !Observed Total Precipitation, Observed Solar Forcing
      REAL :: TOTALPAR,TOTALSAPCP               !Observed PAR,Stageiv precip
      REAL :: TOTALOBRTTMP                      !Observed Brightness Temp

      INTEGER :: countfor

!=== Variables for temperature assimilation and bias correction
      REAL :: MOSBETAK(5),MOSDELT,FBIAS,CLMDELT,CLMBETAK(5)
      REAL :: TOVSTS                !TOVS skin temperature data

      end type
      end module grid_module

