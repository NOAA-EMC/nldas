      SUBROUTINE GRIBOUT(NX,NY,LDASMASK,NSOLD,HOUR,YESTERDAY,
     &        TODAY,NSWRS,NLWRS,LHTFL,SHTFL,GFLUX, 
     &        SNOHF, DSWRF,DLWRF,ASNOW,ARAIN,EVP,SSRUN,BGRUN,
     &        SNOM, AVSFT,ALBDO,WEASD,CWAT,SOILT,SOILM_TOTAL,
     &        SOILM_ROOT,SOILM_1M,SOILM,LSOIL,MSTAV,SOILW,
     &        EVCW,TRANS,EVBS,SBSNO,PEVPR,ACOND,SNOD,SNOC,
     &        RC,RCS,RCT,RCQ,RCSOIL,FX,RSMIN,LAI,GVF,DUMMY)

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
C  01-06-02 LOHMANN  ADDED OUTPUT VARIABLES FOR LDAS
C  13-09-25 YOULONG XIA REVISED for GRIB2 OUTPUT
C***********************************************************************

      IMPLICIT NONE

      INTEGER NX
      INTEGER NY
      INTEGER NSOLD

      CHARACTER*8 YESTERDAY, TODAY

      LOGICAL*1   LDASMASK(NX,NY)

      REAL    NSWRS(NX,NY) 
      REAL    NLWRS(NX,NY)
      REAL    LHTFL(NX,NY) 
      REAL    SHTFL(NX,NY) 
      REAL    GFLUX(NX,NY)
      REAL    SNOHF(NX,NY)
      REAL    DSWRF(NX,NY)
      REAL    DLWRF(NX,NY)
      REAL    ASNOW(NX,NY)
      REAL    ARAIN(NX,NY)

      REAL    EVP(NX,NY)
      REAL    SSRUN(NX,NY)
      REAL    BGRUN(NX,NY) 
      REAL    SNOM(NX,NY) 
      REAL    AVSFT(NX,NY)
      REAL    ALBDO(NX,NY) 
      REAL    WEASD(NX,NY)
      REAL    CWAT(NX,NY)
      REAL    SOILT(NSOLD,NX,NY)
      REAL    SOILM_TOTAL(NX,NY)

      REAL    SOILM_ROOT(NX,NY)
      REAL    SOILM_1M(NX,NY)
      REAL    SOILM(NSOLD,NX,NY)
      REAL    LSOIL(NSOLD,NX,NY)
      REAL    MSTAV(NX,NY)
      REAL    SOILW(NX,NY)
      REAL    EVCW(NX,NY)
      REAL    TRANS(NX,NY)
      REAL    EVBS(NX,NY)
      REAL    SBSNO(NX,NY)

      REAL    PEVPR(NX,NY)
      REAL    ACOND(NX,NY)
      REAL    RC(NX,NY)
      REAL    RCS(NX,NY)
      REAL    RCT(NX,NY)
      REAL    RCQ(NX,NY)
      REAL    RCSOIL(NX,NY)
      REAL    FX(NX,NY)
      REAL    RSMIN(NX,NY)
      REAL    LAI(NX,NY)
      REAL    GVF(NX,NY)
      REAL    SNOD(NX,NY)
      REAL    SNOC(NX,NY)

      REAL    DUMMY(NX,NY)

      INTEGER IPTABLE(13)
      INTEGER HOUR,HOUR1
      INTEGER LUGB, I, J, K, N, IRET, JRET, NLDAS, IG, NFIELDS

      CHARACTER CENT*2, CHOUR*2, GRIBFILE*256
      CHARACTER*19 VANAM 

C     read past header info in KPDS.tbl
      OPEN (UNIT = 30, FILE = 'KPDS.tbl') 
      DO K = 1, 42
         READ(30,*)
      END DO
      
C     MAKE OUTPUT FILENAME
      IF(HOUR.NE.24)THEN
      HOUR1=HOUR
      WRITE (CHOUR, '(i2.2)') HOUR
      GRIBFILE = YESTERDAY//CHOUR//'.NOAH.grb'//char(0)
      ELSE
      HOUR1=HOUR-24
      WRITE (CHOUR, '(i2.2)') HOUR1
      GRIBFILE = TODAY//CHOUR//'.NOAH.grb'//char(0)
      ENDIF

C     OPEN GRIB FILE

      LUGB = 40
      CALL BAOPEN(LUGB, GRIBFILE, IRET)    

C     WRITE EACH FIELD (IN GRIB) TO HOURLY OUTPUT FIELD

      NLDAS = NX * NY
       
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     NSWRS,IPTABLE,LUGB,IRET) 
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = NSWRS'
         STOP
      END IF
      
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     NLWRs,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = NLWRS'
         STOP
      ENDIF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     LHTFL,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = LHTFL'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SHTFL,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SHTFL'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     GFLUX,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = GFLUX'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SNOHF,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SNOHF'
         STOP
      END IF

       READ (30,*) VANAM, (IPTABLE(I),I=1,13)
       CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     DSWRF,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = DSWRF'
         STOP
      END IF

       READ (30,*) VANAM, (IPTABLE(I),I=1,13)
       CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     DLWRF,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = DLWRF'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     ASNOW,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = ASNOW'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     ARAIN,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = ARAIN'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     EVP,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = EVP'
         STOP
      END IF
      
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SSRUN,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SSRUN'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     BGRUN,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = BGRUN'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SNOM,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SNOM'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     AVSFT,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = AVSFT'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     ALBDO,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = ALBDO'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     WEASD,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = WEASD'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     CWAT,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = CWAT'
         STOP
      END IF
      
      DO N = 1, NSOLD
         DO J = 1,NY
            DO I = 1,NX
               DUMMY(I,J) = SOILT(N,I,J)
            END DO
         END DO
         READ (30,*) VANAM, (IPTABLE(I),I=1,13)
         CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +        DUMMY,IPTABLE,LUGB,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HOUR,',','FIELD = SOILT',N
            STOP
          END IF
      END DO
      
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SOILM_TOTAL,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD=SOILM_TOTAL'
         STOP
      END IF
      
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SOILM_ROOT,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD=SOILM_ROOT'
         STOP
      END IF
      
      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SOILM_1M,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SOILM_1M'
         STOP
      END IF
      
      DO N = 1, NSOLD
         DO J = 1,NY
            DO I = 1,NX
               DUMMY(I,J) = SOILM(N,I,J)
            END DO
         END DO
         READ (30,*) VANAM, (IPTABLE(I),I=1,13)
         CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     DUMMY,IPTABLE,LUGB,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HOUR,',','FIELD = SOILM',N
            STOP
         END IF
      END DO
      
      DO N = 1, NSOLD
         DO J = 1,NY
            DO I = 1,NX
               DUMMY(I,J) = LSOIL(N,I,J)
            END DO
         END DO
         READ (30,*) VANAM, (IPTABLE(I),I=1,13)
         CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     DUMMY,IPTABLE,LUGB,IRET)
         IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB FAILED FOR HOUR=',HOUR,',','FIELD = LSOIL',N
            STOP
         END IF
      END DO

       READ (30,*) VANAM, (IPTABLE(I),I=1,13)
       CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     MSTAV,IPTABLE,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = MSTAV'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SOILW,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SOILW'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     EVCW,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = EVCW'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     TRANS,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = TRANS'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     EVBS,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = EVBS'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SBSNO,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SBSNO'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     PEVPR,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = PEVPR'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     ACOND,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = ACOND'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SNOD,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SNOD'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     SNOC,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SNOC'
         STOP
      END IF


      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
        do i=1,nx
        do j=1,ny
        rc(i,j)=1./rc(i,j)
        enddo
        enddo
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     RC,iptable,LUGB,IRET)

      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = RC'
         STOP
      END IF


      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     RCS,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SNOC'
         STOP
      END IF


      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     RCT,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = RCT'
         STOP
      END IF


      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     RCQ,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = RCQ'
         STOP
      END IF


      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     RCSOIL,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = RCSOIL'
         STOP
      END IF


      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     RSMIN,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = RSMIN'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     LAI,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = LAI'
         STOP
      END IF

      READ (30,*) VANAM, (IPTABLE(I),I=1,13)
      CALL grib2_wrt_g2func(LDASMASK,HOUR1,YESTERDAY,TODAY,
     +     GVF,iptable,LUGB,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = GVF'
         STOP
      END IF
C     ------------------------ close grib file-----------------------
      CALL BACLOSE (LUGB, JRET)
      CLOSE (30)
      END
