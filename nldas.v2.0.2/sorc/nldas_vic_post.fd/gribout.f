      SUBROUTINE GRIBOUT(NX,NY,pp_output,lats_table,LDASMASK,
     &        day,month,year,HOUR,  
     &        NSWRS,NLWRS,LHTFL,SHTFL,GFLUX,SNOHF,DSWRF,DLWRF,
     &        ASNOW,ARAIN,EVP,SSRUN,BGRUN,SNOM,
     &        SNOWT,AVSFT,RADT,ALBDO,WEASD,CWAT,
     &        soilt1,soilt2,soilt3,
     &        soilm,soilmr,soilm1m,
     &        soilm1,soilm2,soilm3,
     &        lsoil1,lsoil2,lsoil3,
     &        mstav,mstavr,
     &        evcw,trans,evbs,sbsno,acond,lai,
     &        snod,snoc,salbd)


C***********************************************************************
C  PROGRAM:  GRIBLDAS           WRITE LDAS OUTPUT IN GRIB FORMAT
C  PRGMMR:  MARSHALL & LOHMANN  ORG:  W/NP20
C
C  ABSTRACT:  TAKES HOURLY OUTPUT FILES FROM LDAS MODELS (IN BINARY
C             FORMAT - BIG ENDIAN) AND MAKES GRIB RECORDS FOR ALL FIELDS
C             IN THE FILE.  NOTE ORDER OF FIELDS IN EACH BINARY FILE
C             MUST CORRESPOND TO ORDER OF FIELDS IN KPDS.tbl.  TIME
C             STAMP INFORMATION MUST ALSO APPEAR (YYYYMMDDHH) IN THE
C             TITLE OF EACH BINARY HOURLY OUTPUT FILE.  THIS INFO IS
C             USED IN CONSTRUCTING THE TIME-RELATED PDS OCTETS WHICH
C             DO NOT APPEAR IN KPDS.tbl   
C
C  PROGRAM HISTORY LOG:
C  00-10-15 MARSHALL ORIGINAL TEMPLATE FOR NOAH MODEL
C  00-11-01 LOHMANN  FINALIZED AND CHANGED TO SUBROUTINE
C
C  Jul 2001 Sheffield - Removed hard coding of records between 13Z day 
C			before yesterday and 12Z yesterday.
C  15 October 2013 Xia - revised for GRIB2 output
C***********************************************************************

      IMPLICIT NONE

      include "postproc.h"

      INTEGER NX
      INTEGER NY
      character*256 lats_table, pp_output

      CHARACTER*8 YESTERDAY, TODAY, DAYBEFORE

      integer hour, day, month, year

      REAL    NSWRS(NX,NY),NLWRS(NX,NY),LHTFL(NX,NY),SHTFL(NX,NY) 
      REAL    GFLUX(NX,NY),SNOHF(NX,NY),DSWRF(NX,NY),DLWRF(NX,NY)
      REAL    ASNOW(NX,NY),ARAIN(NX,NY),EVP(NX,NY),SSRUN(NX,NY)
      REAL    BGRUN(NX,NY),SNOM(NX,NY)
      REAL    SNOWT(NX,NY),AVSFT(NX,NY),RADT(NX,NY),ALBDO(NX,NY)
      REAL    WEASD(NX,NY),CWAT(NX,NY)

      REAL    soilt1(NX,NY),soilt2(nx,ny),soilt3(nx,ny)
      real    soilm(nx,ny),soilmr(nx,ny),soilm1m(nx,ny)
      real    soilm1(nx,ny),soilm2(nx,ny),soilm3(nx,ny)
      real    lsoil1(nx,ny),lsoil2(nx,ny),lsoil3(nx,ny)
      real    mstav(nx,ny),mstavr(nx,ny)

      real    evcw(nx,ny),trans(nx,ny),evbs(nx,ny),sbsno(nx,ny)
      real    acond(nx,ny), lai(nx,ny)
     
      real    snod(nx,ny),snoc(nx,ny),salbd(nx,ny)
      REAL    DUMMY(NX,NY)

      INTEGER status, nwrote
C -----------------------13 values for IP Table ----------------------  
      INTEGER IPTABLE(13)

      INTEGER LUGB, I, J, K, N, IRET, JRET, NLDAS, IG, NFIELDS

      LOGICAL*1  LDASMASK(NX,NY)

      CHARACTER CENT*2, GRIBFILE*256
      CHARACTER CYEAR*4, CMONTH*2, CDAY*2, CHOUR*2
      CHARACTER file_ext*8

      CHARACTER*19 VANAM

C     the grib output file extension
      file_ext='.VIC.grb'
      nwrote = 0

C     read past header info in KPDS.tbl
      OPEN (UNIT = 30, FILE = lats_table, iostat=status, status='old') 
      if( status .gt. 0 ) then
        write(*,*) 'ERROR: Can''t open KPDS table: ',trim(lats_table)
        stop
      endif
      if(VERBOSE) then
        write(*,'(a,a)') 'Opened KPDS table: ',trim(lats_table)
      endif           
      DO K = 1, 42
         READ(30,*)
      END DO

C     MAKE OUTPUT FILENAME
      WRITE (CYEAR, '(i4.4)') YEAR
      WRITE (CMONTH, '(i2.2)') MONTH
      WRITE (CDAY, '(i2.2)') DAY
      WRITE (CHOUR, '(i2.2)') HOUR
      GRIBFILE = trim(pp_output)//
     &           CYEAR//CMONTH//CDAY//CHOUR//file_ext//char(0)

      if(VERBOSE) write(*,*) 'grib file: ',trim(GRIBFILE)
      WRITE(TODAY, '(i4.4,i2.2,i2.2)') YEAR, MONTH, DAY
      YESTERDAY = DAYBEFORE(TODAY)
      if(VERBOSE) write(*,*) 'today: ',TODAY," hour: ",HOUR,
     &                       ' yesterday: ',YESTERDAY

C     OPEN GRIB FILE

      LUGB = 9
      CALL BAOPEN (LUGB, GRIBFILE, IRET)

C     WRITE EACH FIELD (IN GRIB2) TO HOURLY OUTPUT FIELD 

      NLDAS = NX * NY

C NSWRS
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     NSWRS,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = NSWRS'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote NSWRS'
         nwrote = nwrote + 1
      END IF
     
C NLWRS
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     NLWRS,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = NLWRS'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote NLWRS'
         nwrote = nwrote + 1
      END IF

C LHTFL
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     LHTFL,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = LHTFL'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote LHTFL'
         nwrote = nwrote + 1
      END IF

C SHTFL
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     SHTFL,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SHTFL'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote SHTFL'
         nwrote = nwrote + 1
      END IF

C GFLUX
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     GFLUX,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = GFLUX'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote GFLUX'
         nwrote = nwrote + 1
      END IF

C SNOHF
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     SNOHF,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SNOHF'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote SNOHF'
         nwrote = nwrote + 1
      END IF

C DSWRF
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
       CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     DSWRF,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = DSWRF'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote DSWRF'
         nwrote = nwrote + 1
      END IF

C DLWRF
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
       CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     DLWRF,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = DLWRF'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote DLWRF'
         nwrote = nwrote + 1
      END IF

C ASNOW
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     ASNOW,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = ASNOW'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote ASNOW'
         nwrote = nwrote + 1
      END IF

C ARAIN
       READ (30,*) VANAM, (IPTABLE(I),I=1,13)
       CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     ARAIN,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = ARAIN'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote ARAIN'
         nwrote = nwrote + 1
      END IF

C EVP
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     EVP,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = EVP'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote EVP'
         nwrote = nwrote + 1
      END IF

C SSRUN
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     SSRUN,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SSRUN'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote SSRUN'
         nwrote = nwrote + 1
      END IF

C BGRUN
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     BGRUN,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = BGRUN'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote BGRUN'
         nwrote = nwrote + 1
      END IF

C SNOM
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     SNOM,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SNOM'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote SNOM'
         nwrote = nwrote + 1
      END IF

C SNOWT
       READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     SNOWT,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SNOWT'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote SNOWT'
         nwrote = nwrote + 1
      END IF

C AVSFT
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     AVSFT,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = AVSFT'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote AVSFT'
         nwrote = nwrote + 1
      END IF

C RADT
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     RADT,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = RADT'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote RADT'
         nwrote = nwrote + 1
      END IF

C ALBDO
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     ALBDO,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = ALBDO'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote ALBDO'
         nwrote = nwrote + 1
      END IF

C WEASD
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     WEASD,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = WEASD'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote WEASD'
         nwrote = nwrote + 1
      END IF

C CWAT
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +     CWAT,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = CWAT'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote CWAT'
         nwrote = nwrote + 1
      END IF

C soilt1
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilt1,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilt1'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilt1'
         nwrote = nwrote + 1
      END IF

C soilt2
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilt2,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilt2'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilt2'
         nwrote = nwrote + 1
      END IF

C soilt3
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilt3,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilt3'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilt3'
         nwrote = nwrote + 1
      END IF

C soilm
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilm,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilm'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilm'
         nwrote = nwrote + 1
      END IF

C soilmr
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilmr,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilmr'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilmr'
         nwrote = nwrote + 1
      END IF

C soilm1m
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilm1m,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilm1m'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilm1m'
         nwrote = nwrote + 1
      END IF

C soilm1
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilm1,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilm1'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilm1'
         nwrote = nwrote + 1
      END IF

C soilm2
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilm2,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilm2'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilm2'
         nwrote = nwrote + 1
      END IF

C soilm3
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        soilm3,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = soilm3'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote soilm3'
         nwrote = nwrote + 1
      END IF

C lsoil1
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        lsoil1,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = lsoil1'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote lsoil1'
         nwrote = nwrote + 1
      END IF

C lsoil2
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        lsoil2,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = lsoil2'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote lsoil2'
         nwrote = nwrote + 1
      END IF

C lsoil3
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        lsoil3,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = lsoil3'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote lsoil3'
         nwrote = nwrote + 1
      END IF

C mstav
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        mstav,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = mstav'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote mstav'
         nwrote = nwrote + 1
      END IF

C mstavr
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        mstavr,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = mstavr'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote mstavr'
         nwrote = nwrote + 1
      END IF

C evcw
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        evcw,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = evcw'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote evcw'
         nwrote = nwrote + 1
      END IF

C trans
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        trans,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = trans'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote trans'
         nwrote = nwrote + 1
      END IF

C evbs
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        evbs,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = evbs'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote evbs'
         nwrote = nwrote + 1
      END IF

C sbsno
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        sbsno,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = sbsno'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote sbsno'
         nwrote = nwrote + 1
      END IF

C acond
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        acond,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = acond'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote acond'
         nwrote = nwrote + 1
      END IF

C lai
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        lai,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = lai'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote lai'
         nwrote = nwrote + 1
      END IF

C snod
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        snod,IPTABLE,LUGB,IRET) 
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = snod'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote snod'
         nwrote = nwrote + 1
      END IF

C snoc
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        snoc,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = snoc'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote snoc'
         nwrote = nwrote + 1
      END IF

C salbd
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR,YESTERDAY,TODAY,
     +        salbd,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = salbd'
         STOP
      ELSE
         if(DEBUG) PRINT*, 'wrote salbd'
         nwrote = nwrote + 1
      END IF

      if(DEBUG) PRINT*, 'wrote ', nwrote, ' records'
      CALL BACLOSE (LUGB, JRET) 
      CLOSE (30)

      END
