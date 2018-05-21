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
! cat_out.f: 
!
! DESCRIPTION:
!  LDAS Catchment data writer.
!
! REVISION HISTORY:
! 4 Apr 2000: Jeffrey Walker; Initial Code
!11 Aug 2000: Jeff/Brian; Multiplied snow depth by fractional coverage
!=========================================================================

      SUBROUTINE CAT_OUT (LDAS,GRID,CAT)
 
! Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use grid_module      ! LDAS non-model-specific grid variables
      use cat_module
      IMPLICIT NONE
      type (ldasdec) ldas
      type (griddec) grid(ldas%nc,ldas%nr)
      type (catdec)  cat(ldas%ncatm)

!=== Local variables =====================================================
      INTEGER :: T,C,R,M,I,N,j
      
      REAL, pointer :: TP(:,:)
      REAL, pointer :: SNOW(:,:)
      REAL, pointer :: QA(:,:)
      REAL, pointer :: TC(:,:)
      REAL, pointer :: TS6(:,:)
      REAL, pointer :: SWUP(:,:)
      REAL, pointer :: LWUP(:,:)
      REAL, pointer :: SHFLX(:,:)
      REAL, pointer :: GHFLX(:,:)
      REAL, pointer :: LHFLX(:,:)
      REAL, pointer :: CATDEF(:,:)
      REAL, pointer :: RZEXC(:,:)
      REAL, pointer :: SRFEXC(:,:)
      REAL, pointer :: SRFMC(:,:)
      REAL, pointer :: RZMC(:,:)
      REAL, pointer :: COLMC(:,:)
      REAL, pointer :: AR1(:,:)
      REAL, pointer :: AR3(:,:)
      REAL, pointer :: INTC(:,:)
      REAL, pointer :: EINT(:,:)
      REAL, pointer :: EVAP(:,:)
      REAL, pointer :: EVEG(:,:)
      REAL, pointer :: ESOI(:,:)
      REAL, pointer :: SNDZ(:,:)
      REAL, pointer :: SMELT(:,:)
      REAL, pointer :: INFIL(:,:)
      REAL, pointer :: BFLOW(:,:)
      REAL, pointer :: RUNSRF(:,:)
      REAL, pointer :: RUNOFF(:,:)

      REAL :: VMEAN,VSTDEV,VMIN,VMAX

      CHARACTER*80 MKFYRMO,FILENMG,CDIR,NAMEG,FILENGB
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13)
      CHARACTER*1  FYRMODIR(80),FSUBFT(80)
      CHARACTER*1  FSUBFG(80),FTIMEB(10),FSUBGB(8)

!=== Variables used for writing output in HDF format
      INTEGER,PARAMETER :: NVARSG=29,KMG=1
      CHARACTER*80 :: TITLEG,TITLET,SOURCE,CONTACT
      CHARACTER*12 :: LEVUNITS
      CHARACTER*80 :: VNAMEG(NVARSG),VTITLEG(NVARSG),VUNITSG(NVARSG)
      INTEGER      :: KMVARG(NVARSG),TIMINC
      REAL         :: LAT(LDAS%NR),LON(LDAS%NC)
      REAL         :: LEVSG(KMG)
      REAL         :: VALID_RANGEG(2,NVARSG),PACKING_RANGEG(2,NVARSG)
      INTEGER      :: PREC,KBEG,KOUNT
     
      DATA TITLEG /"Catchment Grid-space Output based on LDAS Forcing"/
      DATA SOURCE /"NASA GSFC/Data Assimilation Office"/
      DATA CONTACT /"houser@dao.gsfc.nasa.gov"/
      DATA LEVUNITS /"sigma_level"/

      DATA VTITLEG /"Canopy Temperature","Deep Soil Temperature",
     & "Canopy Humidity","Interception Depth","Snow Depth",
     & "Surface Excess","Root Zone Excess","Catchment Deficit",
     & "Vol Surface Soil Moisture","Vol Root Zone Soil Moisture",
     & "Vol Column Average Soil Moisture","Saturated Fraction",
     & "Wilting Fraction","Outgoing Shortwave Radiation",
     & "Outgoing Longwave Radiation Flux","Sensible Heat Flux",
     & "Latent Heat Flux","Ground Heat Flux","Interception Loss",
     & "Evaporation","Bare Soil Evaporation",
     & "Transpiration","Total Precipitation","Total Snowfall",
     & "Total Snowmelt","Total Infiltration",
     & "Base Flow","Total Runoff","Overland Flow"/

      DATA VUNITSG /"Kelvin","Kelvin","Kilograms per kilogram",
     & "Millimeters","Millimeters","Millimeters","Millimeters",
     & "Millimeters","Volumetric","Volumetric","Volumetric",
     & "Unitless","Unitless","Watts per meter squared",
     & "Watts per meter squared","Watts per meter squared",
     & "Watts per meter squared","Watts per meter squared",
     & "Millimeters","Millimeters","Millimeters",
     & "Millimeters","Millimeters","Millimeters",
     & "Millimeters","Millimeters","Millimeters",
     & "Millimeters","Millimeters"/

      DATA VNAMEG /"CT","TS6","QA","INTC","SNDZ","SRFEXC","RZEXC",
     & "CATDEF","SRFMC","RZMC","COLMC","AR1","AR3","SWUP","LWUP",
     & "SHFLX","LHFLX","GHFLX","EINT","EVAP","ESOI","EVEG","TP","SNOW",
     & "SMELT","INFIL","BFLOW","RUNOFF","RSURF"/

      CHARACTER*40 FILE
      CHARACTER*80 NAME

!=== End Variable List ===================================================
  
      CAT%COUNT=CAT%COUNT+1

!=== Total arrays hold a running total of each output variable for time
!=== averaging, between output writes
      CAT%SUMTP=CAT%SUMTP+CAT%TPTP
      DO N=1,LDAS%NCATM
       IF(CAT(N)%TMP2M.LT.273.16) CAT(N)%SUMSNOW=CAT(N)%SUMSNOW+
     &                            CAT(N)%TPTP
      END DO
      CAT%SUMQA=CAT%SUMQA+CAT%QA(1)*CAT%AR(1)+CAT%QA(2)*CAT%AR(2)+
     &          CAT%QA(3)*CAT%AR(3)
      CAT%SUMTC=CAT%SUMTC+CAT%TC(1)*CAT%AR(1)+CAT%TC(2)*CAT%AR(2)+
     &          CAT%TC(3)*CAT%AR(3)
      CAT%SUMTS6=CAT%SUMTS6+CAT%TS(6)
      CAT%SUMSWUP=CAT%SUMSWUP+CAT%MODALB*CAT%SWDN
      CAT%SUMLWUP=CAT%SUMLWUP+CAT%LWUP
      CAT%SUMSHFLX=CAT%SUMSHFLX+CAT%SHFLX
      CAT%SUMGHFLX=CAT%SUMGHFLX+CAT%GHFLX
      CAT%SUMLHFLX=CAT%SUMLHFLX+CAT%LHFLX
      CAT%SUMCATDEF=CAT%SUMCATDEF+CAT%CATDEF
      CAT%SUMRZEXC=CAT%SUMRZEXC+CAT%RZEXC
      CAT%SUMSRFEXC=CAT%SUMSRFEXC+CAT%SRFEXC
      CAT%SUMSRFMC=CAT%SUMSRFMC+CAT%SRFMC
      CAT%SUMRZMC=CAT%SUMRZMC+CAT%RZMC
      CAT%SUMCOLMC=CAT%SUMCOLMC+CAT%COLMC
      CAT%SUMAR1=CAT%SUMAR1+CAT%AR(1)
      CAT%SUMAR3=CAT%SUMAR3+CAT%AR(3)
      CAT%SUMINT=CAT%SUMINT+CAT%INT
      CAT%SUMEINT=CAT%SUMEINT+CAT%EINT
      CAT%SUMEVAP=CAT%SUMEVAP+CAT%EVAP
      CAT%SUMEVEG=CAT%SUMEVEG+CAT%EVEG
      CAT%SUMESOI=CAT%SUMESOI+CAT%ESOI
      CAT%SUMSNDZ=CAT%SUMSNOW+(CAT%SNDZ(1)+CAT%SNDZ(2)+
     &  CAT%SNDZ(3))*CAT%ASNOW
      CAT%SUMSMELT=CAT%SMELT+CAT%SMELT
      CAT%SUMINFIL=CAT%SUMINFIL+CAT%INFIL
      CAT%SUMBFLOW=CAT%SUMBFLOW+CAT%BFLOW
      CAT%SUMRUNSRF=CAT%SUMRUNSRF+CAT%RUNSRF
      CAT%SUMRUNOFF=CAT%SUMRUNOFF+CAT%RUNOFF

!=== Test to see if output writing interval has been reached
      IF(MOD(LDAS%GMT,LDAS%WRITEINTCT).EQ.0)THEN
       LDAS%NUMOUTCT=LDAS%NUMOUTCT+1    !Counts number of output times

!=== Perform time averaging for the total arrays
       CAT%SUMQA=CAT%SUMQA/FLOAT(CAT%COUNT)
       CAT%SUMTC=CAT%SUMTC/FLOAT(CAT%COUNT)
       CAT%SUMTS6=CAT%SUMTS6/FLOAT(CAT%COUNT)+273.16
       CAT%SUMSWUP=CAT%SUMSWUP/FLOAT(CAT%COUNT)
       CAT%SUMLWUP=CAT%SUMSWUP/FLOAT(CAT%COUNT)
       CAT%SUMSHFLX=CAT%SUMSHFLX/FLOAT(CAT%COUNT)
       CAT%SUMGHFLX=CAT%SUMGHFLX/FLOAT(CAT%COUNT) 
       CAT%SUMLHFLX=CAT%SUMLHFLX/FLOAT(CAT%COUNT) 
       CAT%SUMCATDEF=CAT%SUMCATDEF/FLOAT(CAT%COUNT)
       CAT%SUMRZEXC=CAT%SUMRZEXC/FLOAT(CAT%COUNT)
       CAT%SUMSRFEXC=CAT%SUMSRFEXC/FLOAT(CAT%COUNT)
       CAT%SUMSRFMC=CAT%SUMSRFMC/FLOAT(CAT%COUNT)
       CAT%SUMRZMC=CAT%SUMRZMC/FLOAT(CAT%COUNT)
       CAT%SUMCOLMC=CAT%SUMCOLMC/FLOAT(CAT%COUNT)
       CAT%SUMAR1=CAT%SUMAR1/FLOAT(CAT%COUNT)
       CAT%SUMAR3=CAT%SUMAR3/FLOAT(CAT%COUNT)
       CAT%SUMINT=CAT%SUMINT/FLOAT(CAT%COUNT)
       CAT%SUMSNDZ=CAT%SUMSNDZ/FLOAT(CAT%COUNT)

!=== Allocate space for transformed output
       ALLOCATE (TP(LDAS%NC,LDAS%NR))
       ALLOCATE (SNOW(LDAS%NC,LDAS%NR))
       ALLOCATE (QA(LDAS%NC,LDAS%NR))
       ALLOCATE (TC(LDAS%NC,LDAS%NR))
       ALLOCATE (TS6(LDAS%NC,LDAS%NR))
       ALLOCATE (SWUP(LDAS%NC,LDAS%NR))
       ALLOCATE (LWUP(LDAS%NC,LDAS%NR))
       ALLOCATE (SHFLX(LDAS%NC,LDAS%NR))
       ALLOCATE (GHFLX(LDAS%NC,LDAS%NR))
       ALLOCATE (LHFLX(LDAS%NC,LDAS%NR))
       ALLOCATE (CATDEF(LDAS%NC,LDAS%NR))
       ALLOCATE (RZEXC(LDAS%NC,LDAS%NR))
       ALLOCATE (SRFEXC(LDAS%NC,LDAS%NR))
       ALLOCATE (SRFMC(LDAS%NC,LDAS%NR))
       ALLOCATE (RZMC(LDAS%NC,LDAS%NR))
       ALLOCATE (COLMC(LDAS%NC,LDAS%NR))
       ALLOCATE (AR1(LDAS%NC,LDAS%NR))
       ALLOCATE (AR3(LDAS%NC,LDAS%NR))
       ALLOCATE (INTC(LDAS%NC,LDAS%NR))
       ALLOCATE (EINT(LDAS%NC,LDAS%NR))
       ALLOCATE (EVAP(LDAS%NC,LDAS%NR))
       ALLOCATE (EVEG(LDAS%NC,LDAS%NR))
       ALLOCATE (ESOI(LDAS%NC,LDAS%NR))
       ALLOCATE (SNDZ(LDAS%NC,LDAS%NR))
       ALLOCATE (SMELT(LDAS%NC,LDAS%NR))
       ALLOCATE (INFIL(LDAS%NC,LDAS%NR))
       ALLOCATE (BFLOW(LDAS%NC,LDAS%NR))
       ALLOCATE (RUNSRF(LDAS%NC,LDAS%NR))
       ALLOCATE (RUNOFF(LDAS%NC,LDAS%NR))

!=== Transfer from catchment space to LDAS grid space
       CALL TRNSFM(2,LDAS,TP,CAT%SUMTP,CAT%NFF)
       CALL TRNSFM(2,LDAS,SNOW,CAT%SUMSNOW,CAT%NFF)
       CALL TRNSFM(2,LDAS,QA,CAT%SUMQA,CAT%NFF)
       CALL TRNSFM(2,LDAS,TC,CAT%SUMTC,CAT%NFF)
       CALL TRNSFM(2,LDAS,TS6,CAT%SUMTS6,CAT%NFF)
       CALL TRNSFM(2,LDAS,SWUP,CAT%SUMSWUP,CAT%NFF)
       CALL TRNSFM(2,LDAS,LWUP,CAT%SUMLWUP,CAT%NFF)
       CALL TRNSFM(2,LDAS,SHFLX,CAT%SUMSHFLX,CAT%NFF)
       CALL TRNSFM(2,LDAS,GHFLX,CAT%SUMGHFLX,CAT%NFF)
       CALL TRNSFM(2,LDAS,LHFLX,CAT%SUMLHFLX,CAT%NFF)
       CALL TRNSFM(2,LDAS,CATDEF,CAT%SUMCATDEF,CAT%NFF)
       CALL TRNSFM(2,LDAS,SRFMC,CAT%SUMSRFMC,CAT%NFF)
       CALL TRNSFM(2,LDAS,RZEXC,CAT%SUMRZEXC,CAT%NFF)
       CALL TRNSFM(2,LDAS,SRFEXC,CAT%SUMSRFEXC,CAT%NFF)
       CALL TRNSFM(2,LDAS,RZMC,CAT%SUMRZMC,CAT%NFF)
       CALL TRNSFM(2,LDAS,COLMC,CAT%SUMCOLMC,CAT%NFF)
       CALL TRNSFM(2,LDAS,AR1,CAT%SUMAR1,CAT%NFF)
       CALL TRNSFM(2,LDAS,AR3,CAT%SUMAR3,CAT%NFF)
       CALL TRNSFM(2,LDAS,INTC,CAT%SUMINT,CAT%NFF)
       CALL TRNSFM(2,LDAS,EINT,CAT%SUMEINT,CAT%NFF)
       CALL TRNSFM(2,LDAS,EVAP,CAT%SUMEVAP,CAT%NFF)
       CALL TRNSFM(2,LDAS,EVEG,CAT%SUMEVEG,CAT%NFF)
       CALL TRNSFM(2,LDAS,ESOI,CAT%SUMESOI,CAT%NFF)
       CALL TRNSFM(2,LDAS,SNDZ,CAT%SUMSNDZ,CAT%NFF)
       CALL TRNSFM(2,LDAS,SMELT,CAT%SUMSMELT,CAT%NFF)
       CALL TRNSFM(2,LDAS,INFIL,CAT%SUMINFIL,CAT%NFF)
       CALL TRNSFM(2,LDAS,BFLOW,CAT%SUMBFLOW,CAT%NFF)
       CALL TRNSFM(2,LDAS,RUNSRF,CAT%SUMRUNSRF,CAT%NFF)
       CALL TRNSFM(2,LDAS,RUNOFF,CAT%SUMRUNOFF,CAT%NFF)

!=== Generate directory structure and file names for Mosaic output =======
 91    FORMAT(A4,I3,A5,I4,I2,A7,I3,A1)
 92    FORMAT(80A1)
 93    FORMAT(A80)
 94    FORMAT(I4,I2,I2)
 95    FORMAT(8A1)
 96    FORMAT(A40)
 97    FORMAT(A4,I3,A4) 
 98    FORMAT(A4,I3,A5,I4,I2)
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
       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
       OPEN(89,FILE='temp2',FORM='FORMATTED',ACCESS='DIRECT',RECL=200)
       WRITE(90,94,REC=1)LDAS%YR,LDAS%MO,LDAS%DA
       READ(90,95,REC=1)FTIME
       DO I=1,8
        IF(FTIME(I).EQ.(' '))FTIME(I)='0'
       ENDDO
       WRITE(90,91,REC=1)'/EXP',LDAS%EXPCODE,'/CAT/',LDAS%YR,
     &  LDAS%MO,'/LDAS.E',LDAS%EXPCODE,'.'
       READ(90,92,REC=1) (FNAME(I),I=1,29)
       DO I=1,29
        IF(FNAME(I).EQ.(' '))FNAME(I)='0'
       ENDDO

       WRITE(90,96,REC=1) LDAS%ODIR                       
       READ(90,107,REC=1) (FBASE(I),I=1,40)
       C=0
       DO I=1,40
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO

       WRITE(90,98,REC=1)'/EXP',LDAS%EXPCODE,'/CAT/',
     &  LDAS%YR,LDAS%MO
       READ(90,92,REC=1) (FYRMODIR(I),I=1,18)
       DO I=1,18
        IF(FYRMODIR(I).EQ.(' '))FYRMODIR(I)='0'
       ENDDO

       WRITE(90,100,REC=1)'mkdir -p '
       READ(90,92,REC=1)(FMKDIR(I),I=1,9)

       WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &  (FYRMODIR(I),I=1,18)
       READ(90,93,REC=1)MKFYRMO

!=== Make the directories for the Catchment output data files     
       CALL SYSTEM(MKFYRMO)

!=== Generate file name for binary output 
       IF(LDAS%WBIN.EQ.1)THEN
        WRITE(90,109,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
        READ(90,110,REC=1)FTIMEB
        DO I=1,10
         IF(FTIMEB(I).EQ.(' '))FTIMEB(I)='0'
        ENDDO

        WRITE(90,101,REC=1)'.CATgbin'
        READ(90,92,REC=1) (FSUBGB(I),I=1,8)

        WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FNAME(I),I=1,29),
     &                    (FTIMEB(I),I=1,10),(FSUBGB(I),I=1,8 ) 
        READ(90,93,REC=1)FILENGB
       ENDIF

!=== Open HDF/GRB daily output file     
       IF((LDAS%NUMOUTCT.EQ.1.AND.LDAS%WGRB.EQ.0.AND.LDAS%WHDF.EQ.1).OR.
     &  (LDAS%WGRB.EQ.1.AND.LDAS%GMT.EQ.0.))THEN 
        WRITE(90,101,REC=1)'.CATgrid'
        READ(90,92,REC=1) (FSUBFG(I),I=1,8) 

        WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FNAME(I),I=1,29),
     &                    (FTIME(I),I=1,8),(FSUBFG(I),I=1,8 ) 
        READ(90,93,REC=1)FILENMG
        
        IF(LDAS%WGRB.EQ.1)THEN
         WRITE(90,102,REC=1)'cd '
         READ(90,92,REC=1) (FCD(I),I=1,3)
  
         WRITE(90,92,REC=1)(FCD(I),I=1,3),(FBASE(I),I=1,C),
     &   (FYRMODIR(I),I=1,18)
         READ(90,93,REC=1)CDIR

         WRITE(90,108,REC=1)'/latsdiag.dat'
         READ(90,92,REC=1)(FLATS(I),I=1,13)

         WRITE(90,92,REC=1)(FBASE(I),I=1,C),
     &   (FYRMODIR(I),I=1,18),(FLATS(I),I=1,13)
         READ(90,93,REC=1)LDAS%LATSDIAGCT
                    
         WRITE(90,92,REC=1)(FNAME(I),I=20,29),(FTIME(I),I=1,8),
     &   (FSUBFG(I),I=1,8)
         READ(90,93,REC=1)NAMEG

         WRITE(89,103,REC=1)CDIR,'; lats4d -i ',NAMEG,' -o ',NAMEG,
     &   ' -format grads_grib >> latsdiag.dat'
         READ(89,104,REC=1)LDAS%RLATSGCT


        ENDIF

        IF(LDAS%WHDF.EQ.0.AND.LDAS%WGRB.EQ.1)THEN
         WRITE(89,105,REC=1)'rm ',FILENMG
         READ(89,106,REC=1)LDAS%RMHDFGCT
        ENDIF
        
        CLOSE(90)
        CLOSE(89) 
         
!=== Set up HDF Input Parameters ===============================================
        PREC = 0     !Data precision (0=32 bit, 1=64 bit) 
        TIMINC= 10000*INT(LDAS%WRITEINTCT)
        DO C=1,2
         DO R=1,NVARSG 
          VALID_RANGEG(C,R)=LDAS%UDEF      !Do not use packing algorithm
          PACKING_RANGEG(C,R)=LDAS%UDEF
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


       ENDIF       !End HDF/Grib daily output file setup

!=== Write output in HDF and binary (if WBIN=1) format
       IF(LDAS%WHDF.EQ.1.OR.LDAS%WGRB.EQ.1)THEN
        KBEG=0     ! For 2D fields
        KOUNT=1
 
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(1),KBEG,KOUNT,           ! Canopy Temperature
     &   TC,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(2),KBEG,KOUNT,           ! Deep Soil Temperature
     &   TS6,LDAS,GRID) 
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(3),KBEG,KOUNT,           ! Canopy Humidity
     &   QA,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(4),KBEG,KOUNT,           ! Interception Depth
     &   INTC,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(5),KBEG,KOUNT,           ! Snow Depth
     &   SNDZ,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(6),KBEG,KOUNT,           ! Surface Excess
     &   SRFEXC,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(7),KBEG,KOUNT,           ! Root Zone Excess
     &   RZEXC,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(8),KBEG,KOUNT,           ! Catchment Deficit
     &   CATDEF,LDAS,GRID) 
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(9),KBEG,KOUNT,           ! Vol Surface Soil Moisture
     &   SRFMC,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(10),KBEG,KOUNT,          ! Vol Root Zone Soil Moisture
     &   RZMC,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(11),KBEG,KOUNT,          ! Vol Column Average Soil Moisture
     &   COLMC,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(12),KBEG,KOUNT,          ! Saturated Fraction
     &   AR1,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(13),KBEG,KOUNT,          ! Wilting Fraction
     &   AR3,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(14),KBEG,KOUNT,          ! Outgoing Shortwave Radiation
     &   SWUP,LDAS,GRID) 
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(15),KBEG,KOUNT,          ! Outgoing Longwave Radiation Flux
     &   LWUP,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(16),KBEG,KOUNT,          ! Sensible Heat Flux
     &   SHFLX,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(17),KBEG,KOUNT,          ! Latent Heat Flux
     &   LHFLX,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(18),KBEG,KOUNT,          ! Ground Heat Flux
     &   GHFLX,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(19),KBEG,KOUNT,          ! Interception Loss
     &   EINT,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(20),KBEG,KOUNT,          ! Evaporation
     &   EVAP,LDAS,GRID) 
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(21),KBEG,KOUNT,          ! Bare Soil Evaporation
     &   ESOI,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(22),KBEG,KOUNT,          ! Transpiration
     &   EVEG,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(23),KBEG,KOUNT,          ! Total Precipitation
     &   TP,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(24),KBEG,KOUNT,          ! Total Snowfall
     &   SNOW,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(25),KBEG,KOUNT,          ! Total Snowmelt
     &   SMELT,LDAS,GRID) 
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(26),KBEG,KOUNT,          ! Total Infiltration
     &   INFIL,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(27),KBEG,KOUNT,          ! Base Flow
     &   BFLOW,LDAS,GRID)
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(28),KBEG,KOUNT,          ! Total Runoff
     &   RUNSRF,LDAS,GRID)                            
        CALL GSHDF(LDAS%FIDGCT,LDAS%RCGCT,VNAMEG(29),KBEG,KOUNT,          ! Overland Flow
     &   RUNOFF,LDAS,GRID)

        IF((24-LDAS%GMT).LE.LDAS%WRITEINTCT.OR.
     &     LDAS%ENDTIME.EQ.1)THEN
         LDAS%NUMOUTCT=0    !Reset output time counter
         IF(LDAS%WHDF.EQ.1.AND.LDAS%WGRB.EQ.1)THEN
          OPEN(77,FILE=LDAS%LATSDIAGCT,FORM='FORMATTED')
          CALL SYSTEM(LDAS%RLATSGCT)
          CLOSE(77)
         ENDIF

         IF(LDAS%WHDF.EQ.0.AND.LDAS%WGRB.EQ.1)THEN
          OPEN(77,FILE=LDAS%LATSDIAGCT,FORM='FORMATTED')
          CALL SYSTEM(LDAS%RLATSGCT)
          CALL SYSTEM(LDAS%RMHDFGCT)
          CLOSE(77)
         ENDIF
        ENDIF

       ENDIF

!=== Write grid output in binary
       IF(LDAS%WBIN.EQ.1)THEN
        OPEN(58,file=FILENGB,FORM='UNFORMATTED')
        WRITE(58)TC
        WRITE(58)TS6
        WRITE(58)QA
        WRITE(58)INTC
        WRITE(58)SNDZ
        WRITE(58)SRFEXC
        WRITE(58)RZEXC
        WRITE(58)CATDEF
        WRITE(58)SRFMC
        WRITE(58)RZMC
        WRITE(58)COLMC
        WRITE(58)AR1
        WRITE(58)AR3
        WRITE(58)SWUP
        WRITE(58)LWUP
        WRITE(58)SHFLX
        WRITE(58)LHFLX
        WRITE(58)GHFLX
        WRITE(58)EINT
        WRITE(58)EVAP
        WRITE(58)ESOI
        WRITE(58)EVEG
        WRITE(58)TP
        WRITE(58)SNOW
        WRITE(58)SMELT
        WRITE(58)INFIL
        WRITE(58)BFLOW
        WRITE(58)RUNSRF
        WRITE(58)RUNOFF
        CLOSE(58)
       ENDIF


!=== Write statistical output
       IF(LDAS%CATopen.EQ.0)THEN
        FILE='CATstats.dat'
        CALL OPENFILE(NAME,LDAS%ODIR,LDAS%EXPCODE,FILE)
        IF(LDAS%STARTCODE.EQ.1)THEN
         OPEN(60,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN',
     1    POSITION='APPEND')
        ELSE
         OPEN(60,FILE=NAME,FORM='FORMATTED',STATUS='REPLACE')
        ENDIF
        LDAS%CATopen=1
       ENDIF
       WRITE(60,996)'       Statistical Summary of CAT Output for:  ',
     & LDAS%MO,'/',LDAS%DA,'/',LDAS%YR,LDAS%HR,':',LDAS%MN,':',LDAS%SS
996    FORMAT(A47,I2,A1,I2,A1,I4,1X,I2,A1,I2,A1,I2)
       WRITE(60,*)
       WRITE(60,997)
997    FORMAT(T27,'Mean',T41,'StDev',T56,'Min',T70,'Max') 
       
       CALL STATSG(LDAS%UDEF,TC,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(60,999)'TC (K):           ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,TS6,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(60,999)'TS6 (K):          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,QA,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(60,999)'QA (kg/kg):       ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,INTC,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'INTC (mm):         ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,SNDZ,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'SNDZ (mm):        ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,SRFEXC,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'SRFEXC (mm):      ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,RZEXC,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'RZEXC (mm):       ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,CATDEF,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'CATDEF (mm):      ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,SRFMC,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'SRFMC (Vol):      ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,RZMC,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'RZMC (Vol):       ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,COLMC,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'COLMC (Vol):      ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,AR1,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(60,999)'AR1 (-):          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,AR3,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(60,999)'AR3 (-):          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,SWUP,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'SWUP (W/m2):      ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,LWUP,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'LWUP (W/m2):      ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,SHFLX,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'SHFLX (W/m2):     ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,LHFLX,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'LHFLX (W/m2):     ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,GHFLX,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'GHFLX (W/m2):     ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,EINT,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'EINT (mm):        ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,EVAP,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'EVAP (mm):        ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,ESOI,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'ESOI (mm):        ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,EVEG,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'EVEG (mm):        ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,TP,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,VMAX)
       WRITE(60,999)'TP (mm):          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,SNOW,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'SNOW (mm):        ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,SMELT,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'SMELT (mm):       ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,INFIL,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'INFIL (mm):       ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,BFLOW,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'BFLOW (mm):       ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,RUNSRF,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'RUNSRF (mm):      ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATSG(LDAS%UDEF,RUNOFF,LDAS%NC,LDAS%NR,VMEAN,VSTDEV,VMIN,
     &            VMAX)
       WRITE(60,999)'RUNOFF (mm):      ',VMEAN,VSTDEV,VMIN,VMAX

999    FORMAT (1X,A18,4F14.3)
998    FORMAT (1X,A18,4E14.3)
       WRITE(60,*)
       WRITE(60,*)
      
!=== Deallocate space for transformed output
       DEALLOCATE (TP)
       DEALLOCATE (SNOW)
       DEALLOCATE (QA)
       DEALLOCATE (TC)
       DEALLOCATE (TS6)
       DEALLOCATE (SWUP)
       DEALLOCATE (LWUP)
       DEALLOCATE (SHFLX)
       DEALLOCATE (GHFLX)
       DEALLOCATE (LHFLX)
       DEALLOCATE (CATDEF)
       DEALLOCATE (RZEXC)
       DEALLOCATE (SRFEXC)
       DEALLOCATE (SRFMC)
       DEALLOCATE (RZMC)
       DEALLOCATE (COLMC)
       DEALLOCATE (AR1)
       DEALLOCATE (AR3)
       DEALLOCATE (INTC)
       DEALLOCATE (EINT)
       DEALLOCATE (EVAP)
       DEALLOCATE (EVEG)
       DEALLOCATE (ESOI)
       DEALLOCATE (SNDZ)
       DEALLOCATE (SMELT)
       DEALLOCATE (INFIL)
       DEALLOCATE (BFLOW)
       DEALLOCATE (RUNSRF)
       DEALLOCATE (RUNOFF)

!=== Reinitialize output variables
       CALL CAT_INITOUT(LDAS,CAT)

      ENDIF
  
      END SUBROUTINE CAT_OUT

      subroutine gshdf(fid,rc,vname,kbeg,kount,var,
     &  ldas,grid)

! Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use grid_module      ! LDAS non-model-specific grid variables
      implicit none
      type (ldasdec) ldas
      type (griddec) grid(ldas%nc,ldas%nr)

!=== Local variables =====================================================
      integer :: c,r
      real,dimension(ldas%nc,ldas%nr):: var

!=== HDF output variables
      integer :: fid,rc,kbeg,kount
      character*80 :: vname

!=== End variable Definition ===============================================

!=== Prepare grid space output for the HDF subroutines
      do c=1,ldas%nc
       do r=1,ldas%nr       
        if(grid(c,r)%imask.eq.0) var(c,r)=ldas%udef
       enddo
      enddo 

!=== Write HDF tile space output                      
     
      end subroutine gshdf
