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
!
! noah_out.f
!
! DESCRIPTION:  
!  LDAS NOAH data writer:  HDF, Binary, GRIB 
!
! REVISION HISTORY:
!  4 Nov. 1999: Jon Radakovich; Initial Code
! 28 Apr. 2002: Kristi Arsenault; Added NOAH LSM to LDAS
!======================================================================

      SUBROUTINE noah_out (ldas, tile, grid, noah)

!=== Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use tile_module      ! LDAS non-model-specific tile variables
      use grid_module      ! LDAS non-model-specific grid variables
      use noah_module      ! NOAH-specific variables

      implicit none 
      type (ldasdec) LDAS
      type (tiledec) TILE(LDAS%NCH)
      type (griddec) GRID(LDAS%NC,LDAS%NR)
      type (noahdec) NOAH(LDAS%NCH)

!=== Local Variables ==================================================
      INTEGER :: T,C,R,M,I,N           
      CHARACTER*80 MKFYRMO,FILENMT,FILENMG,CDIR,NAMET,NAMEG,FILENGB
      CHARACTER*80 MKFYRMO2
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13),FTIMEC(4)
      CHARACTER*1  FYRMODIR(26),FSUBFT(80)
      CHARACTER*1  FSUBFG(80),FTIMEB(10),FSUBGB(9)

!=== Variables used for writing output in HDF format

      INTEGER,PARAMETER :: NVARSG=56,NVARST=56,KMG=1
      CHARACTER*80 :: TITLEG,TITLET,SOURCE,CONTACT
      CHARACTER*12 :: LEVUNITS
      CHARACTER*80 :: VNAME(NVARSG),VTITLEG(NVARSG),VUNITSG(NVARSG)
      CHARACTER*80 :: VNAMET(NVARST),VTITLET(NVARST),VUNITST(NVARST)
      INTEGER      :: KMVARG(NVARSG),KMVART(NVARST),TIMINC
      REAL         :: LAT(LDAS%NR),LON(LDAS%NC)
      REAL         :: LEVSG(KMG),LEVST(LDAS%MAXT)
      REAL         :: VALID_RANGEG(2,NVARSG),PACKING_RANGEG(2,NVARSG)
      REAL         :: VALID_RANGET(2,NVARST),PACKING_RANGET(2,NVARST)
      REAL,DIMENSION(LDAS%NC,LDAS%NR,LDAS%MAXT) :: TEMPVARTS
      INTEGER      :: PREC,KBEGT,KOUNTT

      DATA TITLEG /"NOAH Grid-space Output based on LDAS Forcing"/
      DATA TITLET /"NOAH Tile-space Output based on LDAS Forcing"/
      DATA SOURCE /"NASA GSFC/Data Assimilation Office"/
      DATA CONTACT /"houser@dao.gsfc.nasa.gov"/
      DATA LEVUNITS /"sigma_level"/

      DATA VTITLEG /"Two Meter Air Temperature (K)",
     &  "Two Meter Specific Humidity (kg/kg)",
     &  "Surface Pressure (Pa)",
     &  "Snowfall (frozen precipitation) (kg/m2)", 
     &  "Rainfall (unfrozen precipitation) (kg/m2)",
     &  "Downward Surface Shortwave Radiation (W/m2)",
     &  "Downward Surface Longwave Radiation (W/m2)",
     &  "Net shortwave radiation (surface) (W/m2)",
     &  "Net longwave radiation (surface) (W/m2)",
     &  "Ten Meter U Wind (m/s)",
     &  "Ten Meter V Wind (m/s)",
     &  "Convective Precipitation (kg/m2)",
     &  "Canopy water Content (kg/m2)",
     &  "Skin Temperature (K)",
     &  "Layer 1 Soil Temperature (K)",
     &  "Layer 2 Soil Temperature (K)",
     &  "Layer 3 Soil Temperature (K)",
     &  "Layer 4 Soil Temperature (K)",
     &  "Layer 1 Vol. Soil Moisture = liq+frzn (kg/m2)",
     &  "Layer 2 Vol. Soil Moisture = liq+frzn (kg/m2)",
     &  "Layer 3 Vol. Soil Moisture = liq+frzn (kg/m2)",
     &  "Layer 4 Vol. Soil Moisture = liq+frzn (kg/m2)",
     &  "Layer 1 Vol. Liquid Soil Moisture (kg/m2)",
     &  "Layer 2 Vol. Liquid Soil Moisture (kg/m2)",
     &  "Layer 3 Vol. Liquid Soil Moisture (kg/m2)",
     &  "Layer 4 Vol. Liquid Soil Moisture (kg/m2)",
     &  "Snow Depth (m)",
     &  "Water Equivalent Snow Depth (kg/m2)",
     &  "Surface Exchange Coef. for Heat and Moisture",
     &  "Surface Exchange Coef. for Momentum",
     &  "Surface Albedo fraction (%)",
     &  "Vegetation Greenness (%)",
     &  "Latent Heat Flux (W/m2)",
     &  "Sensible Heat Flux (W/m2)",
     &  "Canopy evaporation (W/m2)",
     &  "Direct Soil Evaporation (W/m2)",
     &  "Transpiration (W/m2)",
     &  "Sublimation from Snowpack (W/m2)",
     &  "Total Evaporation (kg/m2)",
     &  "Excess Canopy Moisture (m)",
     &  "Dewfall Amount (m s-1)",
     &  "Ratio of actual/potential EVAP (dimensionless)",
     &  "Final Potential Evapotranspiration (W/m2)",
     &  "Ground heat flux (W/m2)",
     &  "Snow phase-change heatflux (W/m2)",
     &  "Snow melt (kg/m2) (water equivalent)",
     &  "Fractional Snow Cover (Unitless fraction)",
     &  "Storm surface runoff (kg/m2)",
     &  "Baseflow-groundwater runoff (kg/m2)",
     &  "Canopy resistence (s m-1)",
     &  "Canopy conductance (m s-1)",
     &  "Avail Soil Moist in Root Zone (%)",
     &  "Avail Soil Moist in Total Col (%)",
     &  "Total Col Soil Moist Content (kg/m2)",
     &  "Root Zone Col Soil Moist Content (kg/m2)",
     &  "Top 1-Meter Col Soil Moist Content (kg/m2)" /

      DATA VTITLET /"Two Meter Air Temperature (K)",
     &  "Two Meter Specific Humidity (kg/kg)",
     &  "Surface Pressure (Pa)",
     &  "Snowfall (frozen precipitation) (kg/m2)",
     &  "Rainfall (unfrozen precipitation) (kg/m2)",
     &  "Downward Surface Shortwave Radiation (W/m2)",
     &  "Downward Surface Longwave Radiation (W/m2)",
     &  "Net shortwave radiation (surface) (W/m2)",
     &  "Net longwave radiation (surface) (W/m2)",
     &  "Ten Meter U Wind (m/s)",
     &  "Ten Meter V Wind (m/s)",
     &  "Convective Precipitation (kg/m2)",
     &  "Canopy water Content (kg/m2)",
     &  "Skin Temperature (K)",
     &  "Layer 1 Soil Temperature (K)",
     &  "Layer 2 Soil Temperature (K)",
     &  "Layer 3 Soil Temperature (K)",
     &  "Layer 4 Soil Temperature (K)",
     &  "Layer 1 Vol. Soil Moisture = liq+frzn (kg/m2)",
     &  "Layer 2 Vol. Soil Moisture = liq+frzn (kg/m2)",
     &  "Layer 3 Vol. Soil Moisture = liq+frzn (kg/m2)",
     &  "Layer 4 Vol. Soil Moisture = liq+frzn (kg/m2)",
     &  "Layer 1 Vol. Liquid Soil Moisture (kg/m2)",
     &  "Layer 2 Vol. Liquid Soil Moisture (kg/m2)",
     &  "Layer 3 Vol. Liquid Soil Moisture (kg/m2)",
     &  "Layer 4 Vol. Liquid Soil Moisture (kg/m2)",
     &  "Snow Depth (m)",
     &  "Water Equivalent Snow Depth (kg/m2)",
     &  "Surface Exchange Coef. for Heat and Moisture",
     &  "Surface Exchange Coef. for Momentum",
     &  "Surface Albedo fraction (%)",
     &  "Vegetation Greenness (%)",
     &  "Latent Heat Flux (W/m2)",
     &  "Sensible Heat Flux (W/m2)",
     &  "Canopy evaporation (W/m2)",
     &  "Direct Soil Evaporation (W/m2)",
     &  "Transpiration (W/m2)",
     &  "Sublimation from Snowpack (W/m2)",
     &  "Total Evaporation (kg/m2)",
     &  "Excess Canopy Moisture (m)",
     &  "Dewfall Amount (m s-1)",
     &  "Ratio of actual/potential EVAP (dimensionless)",
     &  "Final Potential Evapotranspiration (W/m2)",
     &  "Ground heat flux (W/m2)",
     &  "Snow phase-change heatflux (W/m2)",
     &  "Snow melt (kg/m2) (water equivalent)",
     &  "Fractional Snow Cover (Unitless fraction)",
     &  "Storm surface runoff (kg/m2)",
     &  "Baseflow-groundwater runoff (kg/m2)",
     &  "Canopy resistence (s m-1)",
     &  "Canopy conductance (m s-1)",
     &  "Avail Soil Moist in Root Zone (%)",
     &  "Avail Soil Moist in Total Col (%)",
     &  "Total Col Soil Moist Content (kg/m2)",
     &  "Root Zone Col Soil Moist Content (kg/m2)",
     &  "Top 1-Meter Col Soil Moist Content (kg/m2)" /


      DATA VUNITSG / "K",
     &  "kg/kg",
     &  "Pa",
     &  "kg/m2",
     &  "kg/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "m/s",
     &  "m/s",
     &  "kg/m2",
     &  "kg/m2",
     &  "K",
     &  "K",
     &  "K",
     &  "K",
     &  "K",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "m",
     &  "kg/m2",
     &  "unitless",
     &  "unitless",
     &  "%",
     &  "%",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "kg/m2",
     &  "m",
     &  "m/s",
     &  "unitless",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "kg/m2",
     &  "%",
     &  "kg/m2",
     &  "kg/m2",
     &  "s/m",
     &  "m/s", 
     &  "%" ,
     &  "%" ,
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2" /


      DATA VUNITST / "K",
     &  "kg/kg",
     &  "Pa",
     &  "kg/m2",
     &  "kg/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "m/s",
     &  "m/s",
     &  "kg/m2",
     &  "kg/m2",
     &  "K",
     &  "K",
     &  "K",
     &  "K",
     &  "K",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2",
     &  "m",
     &  "kg/m2",
     &  "unitless",
     &  "unitless",
     &  "%",
     &  "%",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "kg/m2",
     &  "m",
     &  "m/s",
     &  "unitless",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "kg/m2",
     &  "%",
     &  "kg/m2",
     &  "kg/m2",
     &  "s/m",
     &  "m/s",
     &  "%" ,
     &  "%" ,
     &  "kg/m2",
     &  "kg/m2",
     &  "kg/m2" /


      DATA VNAME / "TMP",
     & "SPFH","PRES","ASNOW","ARAIN",
     & "DSWRF","DLWRF","NSWRS","NLWRS",
     & "UGRD","VGRD","ACPCP",
     & "CWAT","AVSFT",
     & "TSOIL1","TSOIL2","TSOIL3","TSOIL4",
     & "SOILM1","SOILM2","SOILM3","SOILM4",
     & "LSOIL1","LSOIL2","LSOIL3","LSOIL4",
     & "SNOD","WEASD","CH","CM",
     & "ALBDO","VEG",
     & "LHTFL","SHTFL",
     & "EVCW","EVBS","TRANS",
     & "SBSNO","EVP",
     & "DRIP","DEW","BETA",
     & "PEVPR","GFLUX",
     & "SNOWHF","SNOM","SNOC",
     & "SSRUN","BGRUN",
     & "RC","CCOND",
     & "MSTAVRZ","MSTAVC","SOILMT","SOILRZ","SOILMT1M" /

      DATA VNAMET /"TMP",
     & "SPFH","PRES","ASNOW","ARAIN",
     & "DSWRF","DLWRF","NSWRS","NLWRS",
     & "UGRD","VGRD","ACPCP",
     & "CWAT","AVSFT",
     & "TSOIL1","TSOIL2","TSOIL3","TSOIL4",
     & "SOILM1","SOILM2","SOILM3","SOILM4",
     & "LSOIL1","LSOIL2","LSOIL3","LSOIL4",
     & "SNOD","WEASD","CH","CM",
     & "ALBDO","VEG",
     & "LHTFL","SHTFL",
     & "EVCW","EVBS","TRANS",
     & "SBSNO","EVP",
     & "DRIP","DEW","BETA",
     & "PEVPR","GFLUX",
     & "SNOWHF","SNOM","SNOC",
     & "SSRUN","BGRUN",
     & "RC","CCOND",
     & "MSTAVRZ","MSTAVC","SOILMT","SOILRZ","SOILMT1M" /

      CHARACTER*40 FILE
      CHARACTER*80 NAME

!=== End Variable List ===================================================

      NOAH%COUNT=NOAH%COUNT+1

!=== Total arrays hold a running total of each output variable for time
!=== averaging, between output writes
      DO N=1,LDAS%NOAH_NRET
       DO T=1,LDAS%NCH
        NOAH(T)%TOTRET(N)=NOAH(T)%TOTRET(N)+NOAH(T)%RETURN(N)
       ENDDO
      ENDDO

!=== Test to see if output writing interval has been reached
      IF(MOD(LDAS%GMT,LDAS%WRITEINTN).EQ.0)THEN
       LDAS%NUMOUTNH=LDAS%NUMOUTNH+1    !Counts number of output times
!=== Generate directory structure and file names for NOAH output ==================
 91    FORMAT(A7,I3,A1)
 92    FORMAT(80A1)
 93    FORMAT(A80)
 94    FORMAT(I4,I2,I2)
 95    FORMAT(8A1)
 96    FORMAT(A40)
 97    FORMAT(A4,I3,A4)
 98    FORMAT(A4,I3,A6,I4,A1,I4,I2,I2)
100    FORMAT(A9)
101    FORMAT(A9)
102    FORMAT(A3)
103    FORMAT(A80,A12,A26,A4,A26,A35)
104    FORMAT(A183)
105    FORMAT(A3,A80)
106    FORMAT(A83)
107    FORMAT(40A1)
108    FORMAT(A13)
109    FORMAT(I4,I2,I2,I2)
110    FORMAT(10A1)
111    FORMAT(I4)

       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
       WRITE(90,94,REC=1)LDAS%YR,LDAS%MO,LDAS%DA
       READ(90,95,REC=1)FTIME
       DO I=1,8
        IF(FTIME(I).EQ.(' '))FTIME(I)='0'
       ENDDO

       WRITE(90,111,REC=1)LDAS%YR
       READ(90,95,REC=1)FTIMEC
       DO I=1,4
        IF(FTIMEC(I).EQ.(' '))FTIMEC(I)='0'
       ENDDO

       WRITE(90,91,REC=1)'/LDAS.E',LDAS%EXPCODE,'.'
       READ(90,92,REC=1) (FNAME(I),I=1,11)
       DO I=1,11
        IF(FNAME(I).EQ.(' '))FNAME(I)='0'
       ENDDO

       WRITE(90,96,REC=1) LDAS%ODIR
       READ(90,107,REC=1) (FBASE(I),I=1,40)
       C=0
       DO I=1,40
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO

       WRITE(90,98,REC=1)'/EXP',LDAS%EXPCODE,'/NOAH/',
     &  LDAS%YR,'/',LDAS%YR,LDAS%MO,LDAS%DA
       READ(90,92,REC=1) (FYRMODIR(I),I=1,26)
       DO I=1,26
        IF(FYRMODIR(I).EQ.(' '))FYRMODIR(I)='0'
       ENDDO

       WRITE(90,100,REC=1)'mkdir -p '
       READ(90,92,REC=1)(FMKDIR(I),I=1,9)

       WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &  (FYRMODIR(I),I=1,26)
       READ(90,93,REC=1)MKFYRMO

       WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &  '/',FTIMEC,'/',FTIME
       READ(90,93,REC=1)MKFYRMO2
        CLOSE(90)

!=== Make the directories for the NOAH output data files
       CALL SYSTEM(MKFYRMO)
c       CALL SYSTEM(MKFYRMO2)

!=== Generate file name for BINARY output
       IF(LDAS%WBIN.EQ.1)THEN
       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
        WRITE(90,109,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
        READ(90,110,REC=1)FTIMEB
        DO I=1,10
         IF(FTIMEB(I).EQ.(' '))FTIMEB(I)='0'
        ENDDO
        WRITE(90,101,REC=1)'.NOAHgbin'
        READ(90,92,REC=1) (FSUBGB(I),I=1,9)

        WRITE(90,92,REC=1)(FBASE(I),I=1,C),(FYRMODIR(I),I=1,26),
     &       (FNAME(I),I=1,11),(FTIMEB(I),I=1,10),(FSUBGB(I),I=1,9)
        READ(90,93,REC=1)FILENGB
        CLOSE(90)
       ENDIF

!=== Open HDF daily output file
       IF (LDAS%NUMOUTNH.EQ.1.AND.LDAS%WHDF.EQ.1) THEN

        OPEN(90,FILE='temp',FORM='FORMATTED',
     &    ACCESS='DIRECT',RECL=80)
        WRITE(90,101,REC=1)'.NOAHgrid'
        READ(90,92,REC=1) (FSUBFG(I),I=1,9)

        WRITE(90,92,REC=1)(FBASE(I),I=1,C),(FYRMODIR(I),I=1,26),
     &       (FNAME(I),I=1,11),(FTIME(I),I=1,8),(FSUBFG(I),I=1,9)
        READ(90,93,REC=1)FILENMG

        IF(LDAS%WTIL.EQ.1)THEN
          WRITE(90,101,REC=1)'.NOAHtile'
          READ(90,92,REC=1) (FSUBFT(I),I=1,8)
          WRITE(90,92,REC=1)(FBASE(I),I=1,C),(FNAME(I),I=1,29),
     &                      (FTIME(I),I=1,8),(FSUBFT(I),I=1,8 )
          READ(90,93,REC=1)FILENMT
        ENDIF

!=== Set up HDF Input Parameters ===============================================

        PREC = 0     !Data precision (0=32 bit, 1=64 bit)
        TIMINC= 10000*INT(LDAS%WRITEINTN)

        DO C=1,2
         DO R=1,NVARSG
          VALID_RANGEG(C,R)=LDAS%UDEF      !Do not use packing algorithm
          PACKING_RANGEG(C,R)=LDAS%UDEF
         ENDDO
        ENDDO

        DO C=1,2
         DO R=1,NVARST
          VALID_RANGET(C,R)=LDAS%UDEF      !Do not use packing algorithm
          PACKING_RANGET(C,R)=LDAS%UDEF
         ENDDO
        ENDDO

        DO C=1,LDAS%NC
         DO R=1,LDAS%NR
          LAT(R)=GRID(C,R)%LAT
          LON(C)=GRID(C,R)%LON
         ENDDO
        ENDDO

        DO M=1,NVARSG
         KMVARG(M)=0         !Number of levels for each variable in grid space; 0 for 2D variables
        ENDDO

        LEVSG=1.             !Vertical level units for grid space

!=== Create HDF file for grid space output

 
!=== Create HDF file for tile space output

        IF(LDAS%WTIL.EQ.1)THEN      !Set up tile space HDF input parameters
         DO M=1,NVARST
           KMVART(M)=LDAS%MAXT      !Number of levels for each variable in tile space
         ENDDO

         DO M=1,LDAS%MAXT
           LEVST(M)=M               !Vertical level units for tile space
         ENDDO


        ENDIF
        CLOSE(90)

       ENDIF       !End HDF/Grib daily output file setup

!=== Open statistical output file
      IF(LDAS%NOAHopen.EQ.0)THEN
       FILE='NOAHstats.dat'
       CALL OPENFILE(NAME,LDAS%ODIR,LDAS%EXPCODE,FILE)
       IF(LDAS%STARTCODE.EQ.1)THEN
        OPEN(65,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN',
     1   POSITION='APPEND')
       ELSE
        OPEN(65,FILE=NAME,FORM='FORMATTED',STATUS='REPLACE')
       ENDIF
       LDAS%NOAHopen=1
      ENDIF

       WRITE(65,996)'       Statistical Summary of NOAH Output for:  ',
     & LDAS%MO,'/',LDAS%DA,'/',LDAS%YR,LDAS%HR,':',LDAS%MN,':',LDAS%SS
996    FORMAT(A47,I2,A1,I2,A1,I4,1X,I2,A1,I2,A1,I2)
       WRITE(65,*)
       WRITE(65,997)
997    FORMAT(T27,'Mean',T41,'StDev',T56,'Min',T70,'Max')

!=== Write output in HDF and binary (if WBIN=1) format
!--- and GRIB if WGRB==1

        IF(LDAS%WGRB.EQ.1)THEN
          CALL GRIBOUTNOAH (LDAS,GRID,TILE,NOAH,FBASE,FYRMODIR)
        ENDIF

       IF(LDAS%WBIN.EQ.1) OPEN(58,file=FILENGB,FORM='UNFORMATTED')

       IF(LDAS%WBIN.EQ.1.OR.LDAS%WHDF.EQ.1) THEN

       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(1),1,0,LDAS,GRID,TILE,NOAH)                        !Two Meter Air Temperature 
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(2),2,0,LDAS,GRID,TILE,NOAH)                        !Two Meter Specific Humidity
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(3),3,0,LDAS,GRID,TILE,NOAH)                        !Surface Pressure
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(4),4,2,LDAS,GRID,TILE,NOAH)                        !Frozen Precipitation
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(5),5,2,LDAS,GRID,TILE,NOAH)                        !Liquid Precipitation
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(6),6,1,LDAS,GRID,TILE,NOAH)                        !Downward Surface Shortwave Radiation
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(7),7,1,LDAS,GRID,TILE,NOAH)                        !Downward Surface Longwave Radiation
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(8),8,1,LDAS,GRID,TILE,NOAH)                        !Net shortwave radiation (surface)
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(9),9,1,LDAS,GRID,TILE,NOAH)                        !Net longwave radiation (surface)
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(10),10,0,LDAS,GRID,TILE,NOAH)                      !Ten Meter U Wind
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(11),11,0,LDAS,GRID,TILE,NOAH)                      !Ten Meter V Wind
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(12),12,2,LDAS,GRID,TILE,NOAH)                      !Convective Precipitation
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(13),13,0,LDAS,GRID,TILE,NOAH)                      !Canopy Water Content
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(14),14,0,LDAS,GRID,TILE,NOAH)                      !Skin Temperature
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(15),15,0,LDAS,GRID,TILE,NOAH)                      !Layer 1 Soil Temperature
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(16),16,0,LDAS,GRID,TILE,NOAH)                      !Layer 2 Soil Temperature
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(17),17,0,LDAS,GRID,TILE,NOAH)                      !Layer 3 Soil Temperature
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(18),18,0,LDAS,GRID,TILE,NOAH)                      !Layer 4 Soil Temperature
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(19),19,0,LDAS,GRID,TILE,NOAH)                      !Layer 1 Vol Soil Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(20),20,0,LDAS,GRID,TILE,NOAH)                      !Layer 2 Vol Soil Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(21),21,0,LDAS,GRID,TILE,NOAH)                      !Layer 3 Vol Soil Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(22),22,0,LDAS,GRID,TILE,NOAH)                      !Layer 4 Vol Soil Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(23),23,0,LDAS,GRID,TILE,NOAH)                      !Layer 1 Vol Liquid Soil Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(24),24,0,LDAS,GRID,TILE,NOAH)                      !Layer 2 Vol Liquid Soil Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(25),25,0,LDAS,GRID,TILE,NOAH)                      !Layer 3 Vol Liquid Soil Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(26),26,0,LDAS,GRID,TILE,NOAH)                      !Layer 4 Vol Liquid Soil Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(27),27,0,LDAS,GRID,TILE,NOAH)                      !Snow Depth
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(28),28,0,LDAS,GRID,TILE,NOAH)                      !Water Equivalent Snow Depth
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(29),29,0,LDAS,GRID,TILE,NOAH)                      !Sfc Exchange Coef for Heat and Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(30),30,0,LDAS,GRID,TILE,NOAH)                      !Sfc Exchange Coef for Momentum
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(31),31,0,LDAS,GRID,TILE,NOAH)                      !Surface Albedo Fraction
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(32),32,0,LDAS,GRID,TILE,NOAH)                      !Vegetation Greenness Fraction
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(33),33,1,LDAS,GRID,TILE,NOAH)                      !Latent Heat Flux
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(34),34,1,LDAS,GRID,TILE,NOAH)                      !Sensible Heat Flux
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(35),35,1,LDAS,GRID,TILE,NOAH)                      !Canopy Evaporation
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(36),36,1,LDAS,GRID,TILE,NOAH)                      !Bare Soil Evaporation
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(37),37,1,LDAS,GRID,TILE,NOAH)                      !Transpiration
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(38),38,1,LDAS,GRID,TILE,NOAH)                      !Sublimation from Snowpack
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(39),39,2,LDAS,GRID,TILE,NOAH)                      !Total Evaporation
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(40),40,0,LDAS,GRID,TILE,NOAH)                      !Excess Canopy Moisture
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(41),41,0,LDAS,GRID,TILE,NOAH)                      !Dewfall Amount
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(42),42,1,LDAS,GRID,TILE,NOAH)                      !Ratio of Actual/Potential EVAP
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(43),43,1,LDAS,GRID,TILE,NOAH)                      !Final Potential Evapotranspiration
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(44),44,1,LDAS,GRID,TILE,NOAH)                      !Ground Heat Flux
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(45),45,1,LDAS,GRID,TILE,NOAH)                      !Snow phase-change heat flux
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(46),46,2,LDAS,GRID,TILE,NOAH)                      !Snowmelt (Water Equivalent)
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(47),47,0,LDAS,GRID,TILE,NOAH)                      !Fractional Snow Cover
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(48),48,2,LDAS,GRID,TILE,NOAH)                      !Ground Surface Runoff
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(49),49,2,LDAS,GRID,TILE,NOAH)                      !Baseflow Runoff
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(50),50,0,LDAS,GRID,TILE,NOAH)                      !Canopy Resistence
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(51),51,0,LDAS,GRID,TILE,NOAH)                      !Canopy Conductance
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(52),52,0,LDAS,GRID,TILE,NOAH)                      !Avail Root Zone SMC
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(53),53,0,LDAS,GRID,TILE,NOAH)                      !Avail Tot Col SMC
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(54),54,0,LDAS,GRID,TILE,NOAH)                      !Total Column SMC
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(55),55,0,LDAS,GRID,TILE,NOAH)                      !Root Zone Col SMC
       CALL NOAHWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(56),56,0,LDAS,GRID,TILE,NOAH)                      !Top 1-m Col SMC

        ENDIF


      DO M=1,LDAS%NOAH_NRET
       DO T=1,LDAS%NCH
         NOAH(T)%TOTRET(M)=0.0     !Reset totalizing Array
       ENDDO
      ENDDO

!=== Write out vegetation type and grid fraction for tile space output
!=== They were excluded from the output totalizing array
!=== Initialize temporary tile array
       IF(LDAS%WTIL.EQ.1.AND.LDAS%WHDF.EQ.1)THEN

         DO M=1,LDAS%MAXT
          DO R=1,LDAS%NR
           DO C=1,LDAS%NC
             TEMPVARTS(C,R,M)=LDAS%UDEF
           ENDDO
          ENDDO
         ENDDO

         KBEGT=1              !First level to write in tile space
         KOUNTT=LDAS%MAXT     !Number of levels to write in tile space

        DO T=1,LDAS%NCH
         TEMPVARTS(TILE(T)%COL,TILE(T)%ROW,TILE(T)%PVEG)=
     &    FLOAT(TILE(T)%VEGT)
        ENDDO

        DO T=1,LDAS%NCH
         TEMPVARTS(TILE(T)%COL,TILE(T)%ROW,TILE(T)%PVEG)=TILE(T)%FGRD
        ENDDO

       ENDIF

       NOAH%COUNT=0  !Reset Counters

       IF(LDAS%WHDF.EQ.1)THEN
        IF((24-LDAS%GMT).LE.LDAS%WRITEINTN.OR.
     &      LDAS%ENDTIME.EQ.1)THEN
         LDAS%NUMOUTNH=0     !Reset output time counter

        ENDIF
      ENDIF

       IF(LDAS%WBIN.EQ.1) CLOSE(58)

       WRITE(65,*)
       WRITE(65,*)
      ENDIF

      END SUBROUTINE NOAH_OUT

!*** SUBROUTINE NOAHWRT ********************************************

      SUBROUTINE NOAHWRT(FIDG,RCG,FIDT,RCT,VNAME,VAR,AVGOPT,LDAS,GRID,
     &   TILE,NOAH)

! Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use tile_module      ! LDAS non-model-specific tile variables
      use grid_module      ! LDAS non-model-specific grid variables
      use noah_module
      implicit none
      type (ldasdec) ldas
      type (tiledec) tile(ldas%nch)
      type (noahdec) noah(ldas%nch)
      type (griddec) grid(ldas%nc,ldas%nr)

!=== Local variables =====================================================
      INTEGER :: T,C,R,M,I
      INTEGER :: AVGOPT    !0=Instaneuos Output, 1=Time Average Output,2=Totalize
      INTEGER :: VAR       !RETURN Variable to be written
      REAL,DIMENSION(LDAS%NC,LDAS%NR,LDAS%MAXT) :: TEMPVARTS

!=== HDF output variables
      INTEGER :: FIDG,RCG,FIDT,RCT,KBEGG,KOUNTG
      INTEGER :: KBEGT,KOUNTT
      CHARACTER*80 :: VNAME

!=== Statistcal output variables
      REAL :: VMEAN,VSTDEV,VMIN,VMAX

      REAL :: Ttmp(LDAS%NCH)
      REAL :: Gtmp(LDAS%NC,LDAS%NR)

!=== End variable Definition ===============================================

      Ttmp=0.0
      Gtmp=0.0

!=== Parameters for HDF subroutine

      KBEGG=0              !First level to write in grid space; 0 for 2D variables
      KOUNTG=1             !Number of levels to write in grid space
      KBEGT=1              !First level to write in tile space
      KOUNTT=LDAS%MAXT     !Number of levels to write in tile space
      IF(AVGOPT.EQ.0)THEN  !Write Instantenous Output
       DO T=1,LDAS%NCH
        Ttmp(T)=NOAH(T)%RETURN(VAR)  !Prepare TILE Output
       ENDDO
       DO T=1,LDAS%NCH
        Gtmp(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   Gtmp(TILE(T)%COL,TILE(T)%ROW)+NOAH(T)%RETURN(VAR)*TILE(T)%FGRD

       ENDDO
      ENDIF
      IF(AVGOPT.EQ.1)THEN  !Write Time-Averaged Output
       DO T=1,LDAS%NCH
        Ttmp(T)=NOAH(T)%TOTRET(VAR)/FLOAT(NOAH(T)%COUNT)  !Prepare TILE Output
       ENDDO
       DO T=1,LDAS%NCH
        Gtmp(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1    Gtmp(TILE(T)%COL,TILE(T)%ROW) +
     2    NOAH(T)%TOTRET(VAR)*TILE(T)%FGRD/FLOAT(NOAH(T)%COUNT)
       ENDDO
      ENDIF
      IF(AVGOPT.EQ.2)THEN  !Write Totalized Output
       DO T=1,LDAS%NCH
        Ttmp(T)=NOAH(T)%TOTRET(VAR)  !Prepare TILE Output
       ENDDO
       DO T=1,LDAS%NCH
        Gtmp(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1  Gtmp(TILE(T)%COL,TILE(T)%ROW) + NOAH(T)%TOTRET(VAR)*TILE(T)%FGRD
       ENDDO
      ENDIF

      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
c       IF( ((VAR.GE.54).AND.(VAR.LE.59)).OR.
c     &  (VAR.EQ.37).OR.(VAR.EQ.38) ) THEN
c          IF(GRID(C,R)%FIMASK.EQ.0)Gtmp(C,R)=LDAS%UDEF  !Set Ocean to undefined
c       ELSE
          IF(GRID(C,R)%IMASK.EQ.0)Gtmp(C,R)=LDAS%UDEF  !Set Water to undefined
c       ENDIF
       ENDDO
      ENDDO
      IF(LDAS%WTIL.EQ.1)THEN
!=== Prepare tile space output for the HDF subroutines
       DO M=1,LDAS%MAXT
        DO R=1,LDAS%NR
         DO C=1,LDAS%NC
          TEMPVARTS(C,R,M)=LDAS%UDEF
         ENDDO
        ENDDO
       ENDDO

       DO T=1,LDAS%NCH
        TEMPVARTS(TILE(T)%COL,TILE(T)%ROW,TILE(T)%PVEG)=Ttmp(T)
       ENDDO
      ENDIF

!=== Write binary grid space output
      IF(LDAS%WBIN.EQ.1) WRITE(58)GTMP

!=== Write HDF grid space output
      IF(LDAS%WHDF.EQ.1)THEN

!=== Write HDF tile space output
       IF(LDAS%WTIL.EQ.1)THEN
       ENDIF
      ENDIF

!=== Write statistical output
      CALL STATS(Ttmp,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,VMAX)
      IF(VAR.EQ.3.OR.VAR.EQ.23.OR.VAR.EQ.27.OR.VAR.EQ.30.OR.VAR.EQ.31
     & .OR.VAR.EQ.32.OR.VAR.EQ.33)THEN
       WRITE(65,998)VNAME,VMEAN,VSTDEV,VMIN,VMAX
998    FORMAT(1X,A18,4E14.3)
      ELSE
       WRITE(65,999)VNAME,VMEAN,VSTDEV,VMIN,VMAX
999    FORMAT(1X,A18,4F14.3)
      ENDIF
      END

