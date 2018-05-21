      SUBROUTINE GRIB1OUT(NLDAS,LDASMASK,HOUR2,YESTERDAY,
     &        TODAY,TAIR,SPFH,PRES,UWIND,VWIND,DLWRF,
     &        FRAIN,CAPE,CONVPREC,PREC,DSWRF)

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
C  16-03-18 XIA      REVISED FOR NLDAS FORCING
C***********************************************************************

      IMPLICIT NONE

      CHARACTER*8 YESTERDAY, TODAY

      LOGICAL*1   LDASMASK(NLDAS)

      REAL    TAIR(NLDAS) 
      REAL    SPFH(NLDAS)
      REAL    PRES(NLDAS) 
      REAL    UWIND(NLDAS) 
      REAL    VWIND(NLDAS)
      REAL    DLWRF(NLDAS)
      REAL    FRAIN(NLDAS)
      REAL    CAPE(NLDAS)
      REAL    CONVPREC(NLDAS)
      REAL    PREC(NLDAS)
      REAL    DSWRF(NLDAS)

      INTEGER HOUR,HOUR1,HOUR2

      INTEGER LUGB, I, J, K, N, IRET, KPDS(25), KGDS(22),
     &        NLDAS, LENGDS, IG, JRET, NFIELDS

      CHARACTER CENT*2, CHOUR*2, GDS(400), GRIBFILE*256

C     NON-CHANGING PDS ELEMENTS

      KPDS(1)=7                 !ID FOR NCEP PRODUCTS
      KPDS(2)=84                !ID FOR NAMv4 MODEL
      KPDS(3)=110               !grid
      KPDS(4)=192               !BMS FLAG... don't worry about this.
      KPDS(12)=0                !assume output time minute always = 0
      KPDS(13)=1                !hourly output
      KPDS(17)=4                !number of time steps in hourly averaged variables
      KPDS(18)=0                !GRIB version -- left as 0 in NCEP products
      KPDS(19)=130                !version number of KPDS.tbl for LDAS.  
      KPDS(20)=0                !none missing from averages/accumulations (always4)
      KPDS(23)=12               !ID# Match to NASA frocing generation
      KPDS(24)=0                !does not apply to LDAS output
      KPDS(25)=0                !not used

C     DEFINE THE GDS FOR THE LDAS GRID (THIS ARRAY NEVER CHANGES FROM
C     FIELD TO FIELD. ONLY HAS TO BE DEFINED ONCE, UNLIKE THE PDS).

      IG = 110    !FOR LDAS GRID IN SUBROUTINE w3fi71.f in w3lib

      CALL MAKGDS(IG, KGDS, GDS, LENGDS, IRET)

C     read past header info in KPDS.force.ndas_nam.tbl
      DO K = 1, 42
         READ(60,*)
      END DO

C     MAKE OUTPUT FILENAME
      IF(HOUR2.LE.11)THEN
      HOUR1=HOUR2+12
      WRITE (CHOUR, '(i2.2)') HOUR1
      GRIBFILE = YESTERDAY//CHOUR//'.forcing.grb'//char(0)
      ELSE
      HOUR1=HOUR2-12
      WRITE (CHOUR, '(i2.2)') HOUR1
      GRIBFILE = TODAY//CHOUR//'.forcing.grb'//char(0)
      ENDIF

      HOUR=HOUR1

C     OPEN GRIB FILE

      LUGB = 40
      CALL BAOPEN(LUGB, GRIBFILE, IRET)

C     WRITE EACH FIELD (IN GRIB) TO HOURLY OUTPUT FIELD
 15   FORMAT (29x, 7I6) 

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)

C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,TAIR,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = TAIR'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,SPFH,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = SPFH'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,PRES,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = PRES'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,UWIND,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = UWIND'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,VWIND,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = VWIND'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DLWRF,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = DLWRF'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,FRAIN,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = FRAIN'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,CAPE,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = CAPE'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,CONVPREC,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = CONVPREC'
         STOP
      END IF

      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,PREC,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = PREC'
         STOP
      END IF
      READ (60, 15) KPDS(5), KPDS(6), KPDS(7), KPDS(14), 
     +              KPDS(15), KPDS(16), KPDS(22)
C     MAKE TIME DEPENDENT PDS PARAMETERS
      CALL MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      CALL PUTGB(LUGB,NLDAS,KPDS,KGDS,LDASMASK,DSWRF,IRET)
      IF (IRET .NE. 0) THEN
         PRINT*, 'PUTGB FAILED FOR HOUR = ',HOUR,',','FIELD = DSWRF'
         STOP
      END IF

      CALL BACLOSE (LUGB, JRET)
      CLOSE (60)

      END
