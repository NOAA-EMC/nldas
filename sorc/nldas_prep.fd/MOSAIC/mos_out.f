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
! mos_out.f: 
!
! DESCRIPTION:
!  LDAS MOSAIC data writer.
!
! REVISION HISTORY:
! 4 Nov. 1999: Jon Radakovich; Initial Code
! 6 Nov. 1999: Paul Houser; Revision for Moasic Writing
!22 Aug. 2000: Brian Cosgrove; Modified code for output of
!              standard LDAS output variables.  Replaced
!              old output variables with new calls and
!              variables
!07 Sep. 2000: Brian Cosgrove; changed code so that downward lw
!              and sw values are output as averged values and not
!              instantaneous
!16 Nov. 2000: Brian Cosgrove; changed code to allow for output of
!              GRIB files using modified subroutine from NCEP.  Old
!              references to GRIB output based on LATS4D were
!              removed
!05 Sep. 2001: Brian Cosgrove; Added Close/Open statements, altered
!              output directory structure to match LDAS standards
!20 Sep. 2001: Urszula Jambor; Altered hdf and binary output directory 
!              structure to match NLDAS standards
!==========================================================================

      SUBROUTINE MOS_OUT (LDAS,TILE,GRID,MOS)
 
! Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use tile_module      ! LDAS non-model-specific tile variables
      use grid_module      ! LDAS non-model-specific grid variables
      use mos_module
      IMPLICIT NONE
      TYPE (LDASDEC) LDAS
      TYPE (TILEDEC) TILE(LDAS%NCH)
      TYPE (MOSDEC)  MOS(LDAS%NCH)
      TYPE (GRIDDEC) GRID(LDAS%NC,LDAS%NR)

!=== Local variables =====================================================
      INTEGER :: T,C,R,M,I
      
      CHARACTER*80 MKFYRMO,FILENMT,FILENMG,CDIR,NAMET,NAMEG,FILENGB
      CHARACTER*80 MKFYRMO2
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13),FTIMEC(4)
      CHARACTER*1  FYRMODIR(25),FSUBFT(80)
      CHARACTER*1  FSUBFG(80),FTIMEB(10),FSUBGB(8)

!=== Variables used for writing output in HDF format
      INTEGER,PARAMETER :: NVARSG=43,NVARST=38,KMG=1
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
     
      DATA TITLEG /"MOSAIC Grid-space Output based on LDAS Forcing"/
      DATA TITLET /"MOSAIC Tile-space Output based on LDAS Forcing"/
      DATA SOURCE /"NASA GSFC/Data Assimilation Office"/
      DATA CONTACT /"houser@dao.gsfc.nasa.gov"/
      DATA LEVUNITS /"sigma_level"/

      DATA VTITLEG /"Net Surface Shortwave Radiation",
     &  "Net Surface Longwave Radiation",
     &  "Latent Heat Flux (W/m2)",
     &  "Sensible Heat Flux (W/m2)",
     &  "Ground Heat Flux (W/m2)",
     &  "Snow Phase Change Heat Flux (W/m2)",
     &  "Downward Surface Shortwave Radiation (W/m2)",
     &  "Downward Surface Longwave Radiation (W/m2)",
     &  "Snowfall (Kg/m2)",
     &  "Rainfall (Kg/m2)",
     &  "Total Evaporation (Kg/m2)",
     &  "Surface Runoff (Kg/m2)",
     &  "Subsurface Runoff (Kg/m2)",
     &  "Snowmelt (Kg/m2)",
     &  "Average Surface Temperature (K)",
     &  "Surface Albedo, All Wavelengths (%)",
     &  "Snowpack Water Equivalent (Kg/m2)",
     &  "Plant Canopy Surface Water Storage (Kg/m2)",
     &  "Layer 3 Soil Temperature (K)",
     &  "Total Column Soil Moisture (Kg/m2)",
     &  "Root Zone Soil Moisture (Kg/m2)",
     &  "Top 1-meter Soil Moisture (Kg/m2)",
     &  "Layer 1 Soil Moisture (Kg/m2)",
     &  "Layer 2 Soil Moisture (Kg/m2)",
     &  "Layer 3 Soil Moisture (Kg/m2)", 
     &  "Total Soil Column Wetness (%)",
     &  "Root Zone Wetness (%)",
     &  "Canopy Surface Water Evaporation (W/m2)",
     &  "Canopy Transpiration (W/m2)",
     &  "Bare Soil Evaporation (W/m2)",
     &  "Snow Evaporation (W/m2)",
c     &  "Potential Evaporation (W/m2)",
     &  "Aerodynamic Conductance (M/S)",
     &  "Canopy Conductance (M/S)",
     &  "Vegetation Greenness (%)",
     &  "Leaf Area Index",
     &  "Snow Depth (M)",
     &  "Snow Cover (%)",
     &  "Two Meter Temperature (K)",
     &  "Two Meter Humidity (Kg/Kg)",
     &  "Ten Meter U Wind (m/s)",
     &  "Ten Meter V Wind (m/s)",
     &  "Surface Pressure (mb)",
     &  "Convective Precipitation (Kg/m2)" /



      DATA VTITLET /"Canopy Temperature","Deep Soil Temperature",
     & "Canopy Humidity","Interception Depth","Snow Depth",
     & "Swet1 (Relative)","Swet2 (Relative)","Swet3 (Relative)",
     & "Swet Column Average (Relative)","Swet1 (mm)","Swet2 (mm)",
     & "Swet3 (mm)","Swet Column Average (mm)","Swet1 (m3/m3)",
     & "Swet2 (m3/m3)","Swet3 (m3/m3)","Swet Column Average (m3/m3)",
     & "Net Shortwave Radiation","Outgoing Longwave Radiation Flux",
     & "Sensible Heat Flux",
     & "Latent Heat Flux","Ground Heat Flux","Evaporation",
     & "Interception Loss","Evaporation Rate from Bare Soil",
     & "Transpiration Rate","Total Precipitation","Snowfall Rate",
     & "Rate of Snowmelt","Infiltration of Rainwater Rate",
     & "Water Diffusion across Root Zone Bottom","Runoff",
     & "Overland Flow","Leaf Area Index",
     & "Greenness Fraction of Vegetation",
     & "Relative Soil Moisture","Vegetation Type of Tile",
     & "Fraction of Grid Covered by Tile" /



      DATA VUNITSG / "watts per meter squared",
     &  "watts per meter squared",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2", 
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "K", 
     &  "%",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "K",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "Kg/m2",
     &  "%",
     &  "%",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
     &  "W/m2",
c     &  "W/m2",
     &  "M/S",
     &  "M/S",
     &  "%",
     &  "Unitless",
     &  "M",
     &  "%",
     &  "K",
     &  "Kg/Kg",
     &  "M/S",
     &  "M/S",
     &  "mb",
     &  "Kg/m2" /

      DATA VUNITST /"Kelvin","Kelvin","Kilograms per kilogram",
     & "Millimeters","Meters","Relative","Relative","Relative",
     & "Relative","Millimeters","Millimeters","Millimeters",
     & "Millimeters","Volumetric","Volumetric","Volumetric",
     & "Volumetric","Watts per meter squared",
     & "Watts per meter squared","Watts per meter squared",
     & "Watts per meter squared","Watts per meter squared",
     & "Watts per meter squared","Watts per meter squared",
     & "Watts per meter squared","Watts per meter squared",
     & "Kilograms per meter squared per second",
     & "Kilograms per meter squared per second",
     & "Kilograms per meter squared per second",
     & "Kilograms per meter squared per second",
     & "Kilograms per meter squared per second",
     & "Kilograms per meter squared per second",
     & "Kilograms per meter squared per second",
     & "Unitless", "Unitless","Unitless","Type","Fraction" /

      DATA VNAME / "NSWRS",
     &  "NLWRS",
     & "LHTFL","SHTFL","GFLUX",
     & "SNOHF",
     &  "DSWRF","DLWRF","ASNOW",
     & "ARAIN","EVP","SSRUN","BGRUN","SNOM",
     & "AVSFT", 
     & "ALBDO","WEASD","CWAT",
     &  "SOILT3",
     &  "SOILMC",
     &  "SOILMR","SOILMT1",
     &  "SOILM1",
     & "SOILM2","SOILM3","MSTAVC",
     &  "MSTAVR","EVCW", 
     &  "TRANS",
     &  "EVBS",
     &  "SBSNO",
c     &  "PEVPR",
     &  "ACOND","CCOND","VEG",
     &  "LAI","SNOD","SNOC","TMP2M","HUMID","UWIND","VWIND",
     &  "SFCPRS","CPCP" /

      DATA VNAMET /"CT","TD","CH","INT","SNO","S1r","S2r","S3r","Sr",
     & "S1d","S2d","S3d","Sd","S1v","S2v","S3v","Sv","SWnet","LWout",
     & "H","Le","G","E","Eint","Ebs","Etrans","TP","SF","SM","INF",
     & "Ddif","RO","Rsurf","LAI","GREEN","SMrel","VEGT","FGRD" /

      CHARACTER*40 FILE
      CHARACTER*80 NAME
!=== End Variable List ===================================================
      MOS%COUNT=MOS%COUNT+1

!=== Total arrays hold a running total of each output variable for time
!=== averaging, between output writes
      DO M=1,LDAS%MOS_NRET
       DO T=1,LDAS%NCH
        MOS(T)%TOTRET(M)=MOS(T)%TOTRET(M)+MOS(T)%RETURN(M)
       ENDDO
      ENDDO

!=== Test to see if output writing interval has been reached
      IF(MOD(LDAS%GMT,LDAS%WRITEINTM).EQ.0)THEN
       LDAS%NUMOUTM=LDAS%NUMOUTM+1    !Counts number of output times
!=== Generate directory structure and file names for Mosaic output ==================
 91    FORMAT(A7,I3,A1)
 92    FORMAT(80A1)
 93    FORMAT(A80)
 94    FORMAT(I4,I2,I2)
 95    FORMAT(8A1)
 96    FORMAT(A40)
 97    FORMAT(A4,I3,A4) 
 98    FORMAT(A4,I3,A5,I4,A1,I4,I2,I2)
100    FORMAT(A9)
101    FORMAT(A8)
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

       WRITE(90,98,REC=1)'/EXP',LDAS%EXPCODE,'/MOS/',
     &  LDAS%YR,'/',LDAS%YR,LDAS%MO,LDAS%DA
       READ(90,92,REC=1) (FYRMODIR(I),I=1,25)
       DO I=1,25
        IF(FYRMODIR(I).EQ.(' '))FYRMODIR(I)='0'
       ENDDO

       WRITE(90,100,REC=1)'mkdir -p '
       READ(90,92,REC=1)(FMKDIR(I),I=1,9)

       WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &  (FYRMODIR(I),I=1,25)
       READ(90,93,REC=1)MKFYRMO

       WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &  '/',FTIMEC,'/',FTIME
       READ(90,93,REC=1)MKFYRMO2
        CLOSE(90)

!=== Make the directories for the Mosaic output data files     
       CALL SYSTEM(MKFYRMO)
c       CALL SYSTEM(MKFYRMO2)

!=== Generate file name for binary output 
       IF(LDAS%WBIN.EQ.1)THEN
       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
        WRITE(90,109,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
        READ(90,110,REC=1)FTIMEB
        DO I=1,10
         IF(FTIMEB(I).EQ.(' '))FTIMEB(I)='0'
        ENDDO
        WRITE(90,101,REC=1)'.MOSgbin'
        READ(90,92,REC=1) (FSUBGB(I),I=1,8)

        WRITE(90,92,REC=1)(FBASE(I),I=1,C),(FYRMODIR(I),I=1,25),
     &       (FNAME(I),I=1,11),(FTIMEB(I),I=1,10),(FSUBGB(I),I=1,8 ) 
        READ(90,93,REC=1)FILENGB
        CLOSE(90)
       ENDIF

!=== Open HDF daily output file 
       IF (LDAS%NUMOUTM.EQ.1.AND.LDAS%WHDF.EQ.1) THEN
       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
        WRITE(90,101,REC=1)'.MOSgrid'
        READ(90,92,REC=1) (FSUBFG(I),I=1,8) 

        WRITE(90,92,REC=1)(FBASE(I),I=1,C),(FYRMODIR(I),I=1,25), 
     &       (FNAME(I),I=1,11),(FTIME(I),I=1,8),(FSUBFG(I),I=1,8 ) 
        READ(90,93,REC=1)FILENMG
        
        IF(LDAS%WTIL.EQ.1)THEN
         WRITE(90,101,REC=1)'.MOStile'
         READ(90,92,REC=1) (FSUBFT(I),I=1,8)
         WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FNAME(I),I=1,29),
     &                     (FTIME(I),I=1,8),(FSUBFT(I),I=1,8 ) 
         READ(90,93,REC=1)FILENMT
        ENDIF
         
!=== Set up HDF Input Parameters ===============================================
        PREC = 0     !Data precision (0=32 bit, 1=64 bit) 
        TIMINC= 10000*INT(LDAS%WRITEINTM)
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
         KMVARG(M)=0               !Number of levels for each variable in grid space; 0 for 2D variables
        ENDDO

        LEVSG=1.                   !Vertical level units for grid space

!=== Create HDF file for grid space output

        IF(LDAS%WTIL.EQ.1)THEN      !Set up tile space HDF input parameters    
         DO M=1,NVARST
          KMVART(M)=LDAS%MAXT       !Number of levels for each variable in tile space
         ENDDO

         DO M=1,LDAS%MAXT
          LEVST(M)=M                 !Vertical level units for tile space    
         ENDDO      
!=== Create HDF file for tile space output
        ENDIF
        CLOSE(90)
       ENDIF       !End HDF/Grib daily output file setup
!=== Open statistical output file
      IF(LDAS%MOSopen.EQ.0)THEN
       FILE='MOSstats.dat'
       CALL OPENFILE(NAME,LDAS%ODIR,LDAS%EXPCODE,FILE)
       IF(LDAS%STARTCODE.EQ.1)THEN
        OPEN(65,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN',
     1  POSITION='APPEND')
       ELSE
        OPEN(65,FILE=NAME,FORM='FORMATTED',STATUS='REPLACE')
       ENDIF
       LDAS%MOSopen=1
      ENDIF

       WRITE(65,996)'       Statistical Summary of MOS Output for:  ',
     & LDAS%MO,'/',LDAS%DA,'/',LDAS%YR,LDAS%HR,':',LDAS%MN,':',LDAS%SS
996    FORMAT(A47,I2,A1,I2,A1,I4,1X,I2,A1,I2,A1,I2)
       WRITE(65,*)
       WRITE(65,997)
997    FORMAT(T27,'Mean',T41,'StDev',T56,'Min',T70,'Max')

!=== Write output in HDF and binary (if WBIN=1) format
!---and GRIB if WGRB==1

        IF(LDAS%WGRB.EQ.1)THEN
        CALL GRIBOUT (LDAS,GRID,TILE,MOS,FBASE,FYRMODIR)
        ENDIF

       IF(LDAS%WBIN.EQ.1) OPEN(58,file=FILENGB,FORM='UNFORMATTED')

       IF (LDAS%WBIN.EQ.1.OR.LDAS%WHDF.EQ.1) THEN

       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(1),18,1,LDAS,GRID,TILE,MOS)                           !Net Surface Shortwave Radiation
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(2),19,1,LDAS,GRID,TILE,MOS)                           !Net Surface Longwave Radiation
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(3),21,1,LDAS,GRID,TILE,MOS)                           !Latent Heat Flux
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(4),20,1,LDAS,GRID,TILE,MOS)                           !Sensible Heat Flux
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM, 
     &  VNAME(5),22,1,LDAS,GRID,TILE,MOS)                           !Ground Heat Flux 
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(6),45,1,LDAS,GRID,TILE,MOS)                           !Snow Phase Change Heat Flux
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(7),37,1,LDAS,GRID,TILE,MOS)                           !Downward Surface Shortwave Radiation
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(8),38,1,LDAS,GRID,TILE,MOS)                           !Downward Surface Longwave Radiation
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(9),28,2,LDAS,GRID,TILE,MOS)                           !Snowfall
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(10),27,2,LDAS,GRID,TILE,MOS)                          !Total Rainfall
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(11),23,2,LDAS,GRID,TILE,MOS)                          !Total Evaporation
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(12),33,2,LDAS,GRID,TILE,MOS)                          !Surface Runoff
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(13),46,2,LDAS,GRID,TILE,MOS)                          !Subsurface Runoff
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(14),29,2,LDAS,GRID,TILE,MOS)                          !Snowmelt
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(15),1,0,LDAS,GRID,TILE,MOS)                           !Average Surface Temperature
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(16),39,0,LDAS,GRID,TILE,MOS)                          !Surface Albedo, All Wavelengths
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(17),5,0,LDAS,GRID,TILE,MOS)                           !Snowpack Water Equivalent
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(18),4,0,LDAS,GRID,TILE,MOS)                           !Plant Canopy Surface Water Storage
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(19),2,0,LDAS,GRID,TILE,MOS)                           !Layer 3 Soil Temperature
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(20),13,0,LDAS,GRID,TILE,MOS)                          !Total Column Soil Moisture
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(21),52,0,LDAS,GRID,TILE,MOS)                          !Root Zone Soil Moisture
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(22),51,0,LDAS,GRID,TILE,MOS)                          !Top 1 meter Soil Moisture
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(23),10,0,LDAS,GRID,TILE,MOS)                          !Layer 1 Soil Moisture
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(24),11,0,LDAS,GRID,TILE,MOS)                          !Layer 2 Soil Moisture
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(25),12,0,LDAS,GRID,TILE,MOS)                          !Layer 3 Soil Moisture
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(26),17,0,LDAS,GRID,TILE,MOS)                           !Total Soil Column Wetness
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(27),53,0,LDAS,GRID,TILE,MOS)                           !Root Zone Soil Wetness
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(28),24,1,LDAS,GRID,TILE,MOS)                          !Canopy Surface Water Evaporation
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(29),26,1,LDAS,GRID,TILE,MOS)                          !Canopy Transpiration
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(30),25,1,LDAS,GRID,TILE,MOS)                          !Bare Soil Evaporation
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(31),40,1,LDAS,GRID,TILE,MOS)                          !Snow Evaporation
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(32),43,0,LDAS,GRID,TILE,MOS)                          !Aerodynamic Conductance
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(33),44,0,LDAS,GRID,TILE,MOS)                          !Canopy Conductance
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(34),35,0,LDAS,GRID,TILE,MOS)                          !Vegetation Greenness
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(35),34,0,LDAS,GRID,TILE,MOS)                          !Leaf Area Index
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(36),49,0,LDAS,GRID,TILE,MOS)                          !Snow Depth
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(37),50,0,LDAS,GRID,TILE,MOS)                          !Snow Cover
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(38),54,0,LDAS,GRID,TILE,MOS)                          !Two Meter Temp
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(39),55,0,LDAS,GRID,TILE,MOS)                          !Two Meter Humidity
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(40),56,0,LDAS,GRID,TILE,MOS)                          !Ten Meter U Wind
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(41),57,0,LDAS,GRID,TILE,MOS)                          !Ten Meter V Wind
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(42),58,0,LDAS,GRID,TILE,MOS)                          !Surface Pressure
       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
     &  VNAME(43),59,2,LDAS,GRID,TILE,MOS)                          !Convective Precipitation

	ENDIF


      DO M=1,LDAS%MOS_NRET
      DO T=1,LDAS%NCH
       MOS(T)%TOTRET(M)=0.0   !Reset totalizing Array
      ENDDO
      ENDDO
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(48),47,0,LDAS,GRID,TILE,MOS)                          !Surface Momentum Roughness Length
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(49),47,0,LDAS,GRID,TILE,MOS)                          !Surface Heat Roughness Length
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(39),42,0,LDAS,GRID,TILE,MOS)                          !Snow Albedo
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(41),47,1,LDAS,GRID,TILE,MOS)                          !Open Water Evaporation
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(32),48,1,LDAS,GRID,TILE,MOS)                          !Potential Evaporation
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(33),47,0,LDAS,GRID,TILE,MOS)                          !Layer 1 Liquid Soil Moisture
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(34),47,0,LDAS,GRID,TILE,MOS)                          !Layer 2 Liquid Soil Moisture
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(35),47,0,LDAS,GRID,TILE,MOS)                          !Layer 3 Liquid Soil Moisture
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(15),47,0,LDAS,GRID,TILE,MOS)                          !Groundwater Recharge
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(16),47,0,LDAS,GRID,TILE,MOS)                          !Flood Plain Recharge
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(17),47,0,LDAS,GRID,TILE,MOS)                          !Snow Temperature (Depth Avg)
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(18),47,0,LDAS,GRID,TILE,MOS)                          !Vegetation Canopy Temperature
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(19),47,0,LDAS,GRID,TILE,MOS)                          !Bare Soil Surface Temperature
c       CALL MOSWRT(LDAS%FIDGM,LDAS%RCGM,LDAS%FIDTM,LDAS%RCTM,
c     &  VNAME(21),47,0,LDAS,GRID,TILE,MOS)                          !Effective Radiative Surface Temperature


	
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

       MOS%COUNT=0  !Reset Counters

       IF(LDAS%WHDF.EQ.1)THEN  
        IF((24-LDAS%GMT).LE.LDAS%WRITEINTM.OR.
     &      LDAS%ENDTIME.EQ.1)THEN
         LDAS%NUMOUTM=0     !Reset output time counter

        ENDIF
      ENDIF
 
       IF(LDAS%WBIN.EQ.1) CLOSE(58)

       WRITE(65,*)
       WRITE(65,*)
      ENDIF

      END SUBROUTINE MOS_OUT


      SUBROUTINE MOSWRT(FIDG,RCG,FIDT,RCT,VNAME,VAR,AVGOPT,LDAS,GRID,
     &   TILE,MOS)

! Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use tile_module      ! LDAS non-model-specific tile variables
      use grid_module      ! LDAS non-model-specific grid variables
      use mos_module
      implicit none
      type (ldasdec) ldas
      type (tiledec) tile(ldas%nch)
      type (mosdec)  mos(ldas%nch) 
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
        Ttmp(T)=MOS(T)%RETURN(VAR)  !Prepare TILE Output
       ENDDO
       DO T=1,LDAS%NCH
        Gtmp(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1   Gtmp(TILE(T)%COL,TILE(T)%ROW)+MOS(T)%RETURN(VAR)*TILE(T)%FGRD

       ENDDO
      ENDIF
      IF(AVGOPT.EQ.1)THEN  !Write Time-Averaged Output
       DO T=1,LDAS%NCH
        Ttmp(T)=MOS(T)%TOTRET(VAR)/FLOAT(MOS(T)%COUNT)  !Prepare TILE Output
       ENDDO
       DO T=1,LDAS%NCH
        Gtmp(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1    Gtmp(TILE(T)%COL,TILE(T)%ROW) + 
     2    MOS(T)%TOTRET(VAR)*TILE(T)%FGRD/FLOAT(MOS(T)%COUNT)
       ENDDO
      ENDIF
      IF(AVGOPT.EQ.2)THEN  !Write Totalized Output
       DO T=1,LDAS%NCH
        Ttmp(T)=MOS(T)%TOTRET(VAR)  !Prepare TILE Output
       ENDDO
       DO T=1,LDAS%NCH
        Gtmp(TILE(T)%COL,TILE(T)%ROW) =   !Transfer Tile to Grid
     1  Gtmp(TILE(T)%COL,TILE(T)%ROW) + MOS(T)%TOTRET(VAR)*TILE(T)%FGRD
       ENDDO
      ENDIF

      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
c	IF( ((VAR.GE.54).AND.(VAR.LE.59)).OR.
c     &  (VAR.EQ.37).OR.(VAR.EQ.38) ) THEN
c          IF(GRID(C,R)%FIMASK.EQ.0)Gtmp(C,R)=LDAS%UDEF  !Set Ocean to undefined
c	ELSE
          IF(GRID(C,R)%IMASK.EQ.0)Gtmp(C,R)=LDAS%UDEF  !Set Water to undefined
c	ENDIF
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
      IF(LDAS%WHDF.EQ.1)THEN 
!=== Write HDF grid space output

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

















