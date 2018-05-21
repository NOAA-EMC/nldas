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
! for_out.f: 
!
! DESCRIPTION:
!  LDAS forcing data writer.
!
! REVISION HISTORY:
!  4 Nov. 1999: Jon Radakovich; Initial Code
! 11 Apr. 2000: Brian Cosgrove; Changed code so uses forcing land/sea 
!               mask (with inland water filled in)
! 22 Aug. 2000: Brian Cosgrove; Changed code so that output variables
!               are instantaneous values
! 14 Dec. 2000: Brian Cosgrove; Changed code so that standardized LDAS
!               forcing variables are output.  Also added call to new
!               grib output subroutine (GRIBOUT_FORCE), used 
!               when writing GRIB output.
! 17 Sep. 2001: Brian Cosgrove; Changed output directory structure
! 15 Feb. 2002: Urszula Jambor; Added call to gribout_glbforce to 
!               provide grib output option for 9 global forcing 
!               variables.  Important: DOES NOT APPLY TO HDF OUTPUT.
!=========================================================================
      SUBROUTINE FOR_OUT (ldas,tile,grid)
 
! Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use tile_module      ! LDAS non-model-specific tile variables
      use grid_module      ! LDAS non-model-specific grid variables
      implicit none
      type (ldasdec) ldas
      type (tiledec) tile(ldas%nch)
      type (griddec) grid(ldas%nc,ldas%nr)



!=== Local variables =====================================================
      INTEGER :: T,C,R,M,I               
      REAL :: T_TMP(LDAS%NCH),T_SPFH(LDAS%NCH) !Tile space time averages
      REAL :: T_DSWRF(LDAS%NCH),T_DLWRF(LDAS%NCH)
      REAL :: T_UGRD(LDAS%NCH),T_VGRD(LDAS%NCH)
      REAL :: T_PRES(LDAS%NCH),T_APCP(LDAS%NCH)
      REAL :: T_ACPCP(LDAS%NCH),T_PAR(LDAS%NCH)
      REAL :: T_ODSWRF(LDAS%NCH),T_OAPCP(LDAS%NCH)
      REAL :: T_SAPCP(LDAS%NCH),T_BRTTMP(LDAS%NCH)

!=== Statistcal output variables
      REAL :: VMEAN,VSTDEV,VMIN,VMAX

      CHARACTER*80 FILENFT,FILENFG,MKFYRMO,CDIR,NAMET,NAMEG,FILENFGB
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FSUBFG(80),FCD(3),FRM(3),FLATS(13)
      CHARACTER*1  FTIMEC(4)
      CHARACTER*1  FYRMODIR(80),FSUBFT(80),FTIMEB(10),FSUBFGB(8)
      

!=== Variables used for writing output in HDF format
      INTEGER,PARAMETER :: NVARS=14,KMG=1
      CHARACTER*80 :: TITLEG,TITLET,SOURCE,CONTACT
      CHARACTER*8 :: LEVUNITS
      CHARACTER*80 :: VNAME(NVARS),VTITLE(NVARS),VUNITS(NVARS)
      INTEGER :: KMVARG(NVARS),KMVART(NVARS),KBEG,TIMINC
      REAL :: LAT(LDAS%NR),LON(LDAS%NC)
      REAL :: LEVSG(KMG),LEVST(LDAS%MAXT)
      REAL :: VALID_RANGE(2,NVARS),PACKING_RANGE(2,NVARS)
      INTEGER :: PREC,KOUNT
       
      DATA TITLEG /"LDAS Grid-space Eta Forcing"/
      DATA TITLET /"LDAS Tile-space Eta Forcing"/
      DATA SOURCE /"NASA GSFC/Data Assimilation Office"/
      DATA CONTACT /"houser@dao.gsfc.nasa.gov"/
      DATA LEVUNITS /"layers"/


      DATA VTITLE /"Two Meter Temperature",
     & "Two Meter Specific Humidity","Atmospheric Surface Pressure",
     & "Ten Meter Wind in u-direction",
     & "Ten Meter Wind in v-direction","Downward Shortwave Radiation",
     & "Downward Longwave Radiation",
     & "Total Precipitation",
     & "Convective Precipitation","Observed PAR",
     & "Observed Downward Shortwave Radiation",
     & "StageIV Precipitation",
     & "Observed Total Precipitation",
     & "Brightness Temperature" /
                   

      DATA VUNITS /"Kelvin","Kilgrams per kilogram",
     & "Millibars","Meters per second","Meters per second",
     &  "Watts per meter squared",
     &  "Watts per meter squared",
     & "Millimeters","Millimeters","Watts per meter squared",
     & "Watts per meter squared","Millimeters","Millimeters",
     &  "Kelvin" /

      DATA VNAME /"TMP","SPFH","PRES","UGRD","VGRD","DSWRF",
     & "DLWRF","APCP","ACPCP","PAR","ODSWRF","SAPCP","OAPCP",
     & "BRTTMP"/

      CHARACTER*40 FILE
      CHARACTER*80 NAME
   
!=== End Variable List =========================================================
        GRID%COUNTFOR=GRID%COUNTFOR+1
    
!===Transfer forcing data into temporary arrays
      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
        IF (LDAS%DOMAIN == 1) THEN
           GRID(C,R)%TOTALAPCP=GRID(C,R)%ETADATA2(4)/3.0
           GRID(C,R)%TOTALACPCP=GRID(C,R)%ETADATA2(5)/3.0
        ENDIF
        GRID(C,R)%TOTALTMP=GRID(C,R)%FORCING(14)
        GRID(C,R)%TOTALSPFH=GRID(C,R)%FORCING(15)
        GRID(C,R)%TOTALDSWRF=GRID(C,R)%ETASW
        GRID(C,R)%TOTALDLWRF=GRID(C,R)%FORCING(21)
        GRID(C,R)%TOTALUGRD=GRID(C,R)%FORCING(1)
        GRID(C,R)%TOTALVGRD=GRID(C,R)%FORCING(2)
        GRID(C,R)%TOTALPRES=GRID(C,R)%FORCING(3)/100.0
        GRID(C,R)%TOTALPAR=GRID(C,R)%PAR
        GRID(C,R)%TOTALSAPCP=GRID(C,R)%S4PRECIP
        GRID(C,R)%TOTALOAPCP=GRID(C,R)%FORCING(4)*3600.0
        GRID(C,R)%TOTALODSWRF=GRID(C,R)%OBSW
        GRID(C,R)%TOTALOBRTTMP=GRID(C,R)%OBSBT
       ENDDO
      ENDDO
!=== Test to see if output writing interval has been reached
      IF(MOD(LDAS%GMT,LDAS%WRITEINTF).EQ.0)THEN
       LDAS%NUMOUTF=LDAS%NUMOUTF+1    !Counts number of output times
!=== Initialize tile space average arrays 
       T_TMP=0.
       T_SPFH=0.
       T_DSWRF=0.
       T_DLWRF=0.
       T_UGRD=0.
       T_VGRD=0.
       T_PRES=0.
       T_APCP=0.
       T_ACPCP=0.
       T_PAR=0.
       T_ODSWRF=0.
       T_SAPCP=0.
       T_OAPCP=0.
       T_BRTTMP=0.	

!=== Find tile space averages


!=== Mask out water points
       DO C=1,LDAS%NC
        DO R=1,LDAS%NR
         IF(GRID(C,R)%FIMASK.EQ.0)THEN
          GRID(C,R)%TOTALTMP=LDAS%UDEF  
          GRID(C,R)%TOTALSPFH=LDAS%UDEF
          GRID(C,R)%TOTALDSWRF=LDAS%UDEF
          GRID(C,R)%TOTALDLWRF=LDAS%UDEF
          GRID(C,R)%TOTALUGRD=LDAS%UDEF
          GRID(C,R)%TOTALVGRD=LDAS%UDEF 
          GRID(C,R)%TOTALPRES=LDAS%UDEF
          GRID(C,R)%TOTALAPCP=LDAS%UDEF 
          GRID(C,R)%TOTALPAR=LDAS%UDEF
          GRID(C,R)%TOTALODSWRF=LDAS%UDEF
          GRID(C,R)%TOTALOAPCP=LDAS%UDEF
          GRID(C,R)%TOTALSAPCP=LDAS%UDEF
          GRID(C,R)%TOTALACPCP=LDAS%UDEF
          GRID(C,R)%TOTALOBRTTMP=LDAS%UDEF
         ENDIF
        ENDDO
       ENDDO
       DO T=1,LDAS%NCH
        T_TMP(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALTMP
        T_SPFH(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALSPFH
        T_DSWRF(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALDSWRF
        T_DLWRF(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALDLWRF
        T_UGRD(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALUGRD
        T_VGRD(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALVGRD
        T_PRES(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALPRES
        T_APCP(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALAPCP
        T_ACPCP(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALACPCP
        T_PAR(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALPAR
        T_OAPCP(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALOAPCP
        T_SAPCP(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALSAPCP
        T_ODSWRF(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALODSWRF
        T_BRTTMP(T)=GRID(TILE(T)%COL,TILE(T)%ROW)%TOTALOBRTTMP
       ENDDO

!=== Open statistical output file
      FILE='FORstats.dat'
      CALL OPENFILE(NAME,LDAS%ODIR,LDAS%EXPCODE,FILE)
      IF(LDAS%STARTCODE.EQ.1)THEN
       OPEN(66,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN',
     1 POSITION='APPEND')
!      ELSE
!       OPEN(66,FILE=NAME,FORM='FORMATTED',STATUS='REPLACE')
!      OPEN(66,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN')
      ENDIF

       WRITE(66,996)'      Statistical Summary of Forcing Data for:  ',
     & LDAS%MO,'/',LDAS%DA,'/',LDAS%YR,LDAS%HR,':',LDAS%MN,':',LDAS%SS
996    FORMAT(A48,I2,A1,I2,A1,I4,1X,I2,A1,I2,A1,I2)
       WRITE(66,*)
       WRITE(66,997)
997    FORMAT(T27,'Mean',T41,'StDev',T56,'Min',T70,'Max')
               
       CALL STATS(T_TMP,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,999)'T(K)              ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_SPFH,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,998)'Q(kg/kg)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_DSWRF,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,999)'SW(W/m2)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_DLWRF,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,999)'LW(W/m2)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_UGRD,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,999)'U(m/s)            ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_VGRD,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,999)'V(m/s)            ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_PRES,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,999)'PS(Pa)            ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_APCP,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,998)'TP(mm)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_ACPCP,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,998)'CP(mm)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_PAR,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,999)'PAR(W/m2)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_ODSWRF,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,999)'OSW(W/m2)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_SAPCP,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,998)'STP(mm)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_OAPCP,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,998)'OTP(mm)          ',VMEAN,VSTDEV,VMIN,VMAX
       CALL STATS(T_BRTTMP,LDAS%UDEF,LDAS%NCH,VMEAN,VSTDEV,VMIN,
     &  VMAX)
       WRITE(66,998)'OBT(K)          ',VMEAN,VSTDEV,VMIN,VMAX

       WRITE(66,*)
       WRITE(66,*)
998    FORMAT(1X,A18,4E14.3)
999    FORMAT(1X,A18,4F14.3)
!=== Generate directory structure and file names for forcing data ======
 91    FORMAT(A4,I3,A5,I4,I2,A7,I3,A1)
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

       WRITE(90,91,REC=1)'/EXP',LDAS%EXPCODE,'/FOR/',LDAS%YR,
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

       WRITE(90,98,REC=1)'/EXP',LDAS%EXPCODE,'/FOR/',
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
        CLOSE(90)

!=== Make the directories for the Mosaic output data files
       CALL SYSTEM(MKFYRMO)

!=== Generate file name for binary output
       IF(LDAS%WBIN.EQ.1)THEN
        WRITE(90,109,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
        READ(90,110,REC=1)FTIMEB
        DO I=1,10
         IF(FTIMEB(I).EQ.(' '))FTIMEB(I)='0'
        ENDDO

        WRITE(90,101,REC=1)'.FORgbin'
        READ(90,92,REC=1) (FSUBFGB(I),I=1,8)

        WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FNAME(I),I=1,29),
     &                    (FTIMEB(I),I=1,10),(FSUBFGB(I),I=1,8 ) 
        READ(90,93,REC=1)FILENFGB
       ENDIF

!=== Open HDF/GRB daily output file     
       IF((LDAS%NUMOUTF.EQ.1.AND.LDAS%WHDF.EQ.1)) THEN

        WRITE(90,101,REC=1)'.FORgrid'
        READ(90,92,REC=1) (FSUBFG(I),I=1,8)

        WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FNAME(I),I=1,29),
     &                    (FTIME(I),I=1,8),(FSUBFG(I),I=1,8 ) 
        READ(90,93,REC=1)FILENFG

        IF(LDAS%WTIL.EQ.1)THEN
         WRITE(90,101,REC=1)'.FORtile'
         READ(90,92,REC=1) (FSUBFT(I),I=1,8)
         WRITE(90,92,REC=1)(FBASE(I),I=1,C), (FNAME(I),I=1,29),
     &                     (FTIME(I),I=1,8),(FSUBFT(I),I=1,8 ) 
         READ(90,93,REC=1)FILENFT
        ENDIF
 
        CLOSE(90)
!=== Forcing Data File Name Generation Complete

!=== Set up HDF Input Parameters =======================================
        PREC=0     !Data precision (0=32 bit, 1=64 bit) 
        TIMINC= 10000*INT(LDAS%WRITEINTF)
        DO C=1,2
         DO R=1,NVARS 
          VALID_RANGE(C,R)=LDAS%UDEF      !Do not use packing algorithm
          PACKING_RANGE(C,R)=LDAS%UDEF
         ENDDO
        ENDDO

        DO C=1,LDAS%NC
         DO R=1,LDAS%NR
          LAT(R)=GRID(C,R)%LAT
          LON(C)=GRID(C,R)%LON
         ENDDO
        ENDDO

        DO M=1,NVARS
         KMVARG(M)=0               !Number of levels for each variable in grid space; 0 for 2D variables
        ENDDO

        LEVSG=1.                   !Vertical level units for grid space

!=== Create HDF file for grid space output
        IF(LDAS%WTIL.EQ.1)THEN      !Set up tile space HDF input parameters    
         DO M=1,NVARS
          KMVART(M)=LDAS%MAXT       !Number of levels for each variable in tile space
         ENDDO
                       
         LEVST=0.
         DO M=1,LDAS%MAXT
          LEVST(M)=LEVST(M)+1       !Vertical level units for tile space    
         ENDDO   

!=== Create HDF file for tile space output   
        ENDIF

       ENDIF      !End daily HDF/Grib output file setup
!=== Write grid output in binary 
        IF(LDAS%WGRB.EQ.1)THEN
           if (LDAS%DOMAIN==1) then
              CALL GRIBOUT_FORCE (LDAS,GRID,FBASE,FYRMODIR)
              CALL GRIBOUT_FORCE_GRIB2(LDAS,GRID,FBASE,FYRMODIR)
           else 
              CALL GRIBOUT_GLBFORCE(LDAS,GRID,FBASE,FYRMODIR)
           end if
	ENDIF

       IF(LDAS%WBIN.EQ.1)THEN
        OPEN(56,FILE=FILENFGB,FORM='UNFORMATTED')
        WRITE(56)GRID%TOTALTMP
        WRITE(56)GRID%TOTALSPFH
        WRITE(56)GRID%TOTALPRES
        WRITE(56)GRID%TOTALUGRD
        WRITE(56)GRID%TOTALVGRD
        WRITE(56)GRID%TOTALDSWRF
        WRITE(56)GRID%TOTALDLWRF
        WRITE(56)GRID%TOTALAPCP
        WRITE(56)GRID%TOTALACPCP
        WRITE(56)GRID%TOTALPAR
        WRITE(56)GRID%TOTALODSWRF
        WRITE(56)GRID%TOTALSAPCP
        WRITE(56)GRID%TOTALOAPCP
	WRITE(56)GRID%TOTALOBRTTMP

        CLOSE(56)
       ENDIF

       IF(LDAS%WHDF.EQ.1)THEN  
!=== Write grid output in HDF format with the subroutine gfio_putvar
        KBEG=0                 !First level to write; if 2D grid kbeg=0
        KOUNT=1                !Number of levels to write


        IF((24-LDAS%GMT).LE.LDAS%WRITEINTF.or.
     &     LDAS%ENDTIME.EQ.1)THEN
         LDAS%NUMOUTF=0     !Reset output time counter
        ENDIF

!=== Write tile forcing data in HDF format with the subroutine gfio_putvar
        IF(LDAS%WTIL.EQ.1)THEN

         KBEG=1
         KOUNT=LDAS%MAXT
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(1),KBEG,KOUNT,T_TMP,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(2),KBEG,KOUNT,T_SPFH,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(7),KBEG,KOUNT,T_PRES,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(5),KBEG,KOUNT,T_UGRD,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(6),KBEG,KOUNT,T_VGRD,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(3),KBEG,KOUNT,T_DSWRF,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(4),KBEG,KOUNT,T_DLWRF,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(8),KBEG,KOUNT,T_APCP,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(9),KBEG,KOUNT,T_ACPCP,
     &    LDAS,GRID,TILE) 
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(10),KBEG,KOUNT,T_PAR,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(11),KBEG,KOUNT,T_ODSWRF,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(12),KBEG,KOUNT,T_SAPCP,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(13),KBEG,KOUNT,T_OAPCP,
     &    LDAS,GRID,TILE)
         CALL TSHDF(LDAS%FIDTF,LDAS%RCTF,VNAME(14),KBEG,KOUNT,T_BRTTMP,
     &    LDAS,GRID,TILE)



         IF((24-LDAS%GMT).LE.LDAS%WRITEINTF.OR.
     &      LDAS%ENDTIME.EQ.1)THEN


         ENDIF
        ENDIF

       ENDIF
!=== Reinitialize forcing total arrays to zero
       GRID%TOTALTMP=0.
       GRID%TOTALSPFH=0.
       GRID%TOTALDSWRF=0.
       GRID%TOTALDLWRF=0.
       GRID%TOTALUGRD=0.  
       GRID%TOTALVGRD=0.     
       GRID%TOTALPRES=0.
       GRID%TOTALAPCP=0.      
       GRID%TOTALACPCP=0.
       GRID%TOTALPAR=0.
       GRID%TOTALODSWRF=0.
       GRID%TOTALSAPCP=0.    
       GRID%TOTALOAPCP=0.     
       GRID%TOTALOBRTTMP=0.    

       GRID%COUNTFOR=0 

!=== Reinitialize tile space arrays to zero
       T_TMP=0.
       T_SPFH=0.
       T_DSWRF=0.
       T_DLWRF=0.
       T_UGRD=0.
       T_VGRD=0.
       T_PRES=0.
       T_APCP=0.
       T_ACPCP=0.
       T_PAR=0.
       T_ODSWRF=0.
       T_SAPCP=0.
       T_OAPCP=0.
       T_BRTTMP=0.


      ENDIF

      END SUBROUTINE FOR_OUT

                 
      SUBROUTINE TSHDF(FID,RC,VNAME,KBEG,KOUNT,VAR,LDAS,
     &  GRID,TILE)

! Declare modules and data structures
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use tile_module      ! LDAS non-model-specific tile variables
      use grid_module      ! LDAS non-model-specific grid variables

      implicit none
      type (ldasdec) ldas
      type (tiledec) tile(ldas%nch)
      type (griddec) grid(ldas%nc,ldas%nr)

!=== Local variables =====================================================
      INTEGER :: T,C,R,M,I
      REAL,DIMENSION(LDAS%NCH):: VAR
      REAL,DIMENSION(LDAS%NC,LDAS%NR,LDAS%MAXT) :: TEMPVARTS 

!=== HDF output variables
      INTEGER :: FID,RC,KBEG,KOUNT
      CHARACTER*80 :: VNAME

!=== End variable Definition ===============================================

!=== Prepare tile space output for the HDF subroutines
      DO M=1,LDAS%MAXT
       DO R=1,LDAS%NR
        DO C=1,LDAS%NC
         TEMPVARTS(C,R,M)=LDAS%UDEF
        ENDDO
       ENDDO
      ENDDO

      DO T=1,LDAS%NCH
       TEMPVARTS(TILE(T)%COL,TILE(T)%ROW,TILE(T)%PVEG)=VAR(T)
      ENDDO

!=== Write HDF tile space output

      END

