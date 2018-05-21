	PROGRAM MERGEFORCE
C
C***********************************************************************
C  PROGRAM:  MERGEFORCE	     MERGE LDAS FORCING FILES INTO SINGLE FILE
C  PRGMMR:  MARSHALL         ORG:  W/NP20
C
C  ABSTRACT:  MERGE DAILY FORCING FIELDS FOR LDAS INTO A SINGLE FILE.  
C  SEPARATE JOBS CREATE THE NDAS/NAM-BASED FORCING, GOES-BASED FORCING,
C  AND THE PRECIP FORCING, AND SEPARATE HOURLY FILES ARE WRITTEN OUT
C  FOR EACH TYPE.  THIS CODE TAKES THOSE THREE SEPARATE FILES, AND 
C  MERGES THEM INTO ONE FILE CONTAINING ALL THE GRIDS.  THE FILES 
C  THEN CONTAIN ALL THE LDAS FORCING FIELDS FOR A SINGLE HOUR.  THIS
C  CODE IS DESIGNED TO RUN DAILY AND PRODUCE 24 HOURLY FILES. 
C  
C  PROGRAM LOG
C  99-11-16 CURTIS M. UPDATED TO ADD STAGE IV TO FORCING FILES
C  16-10-06 XIA Y. UPDATED FOR NLDAS VERSION 3 TO ACHIEVE ACTUAL REALTIME RUN
C***********************************************************************
C
C  DECLARATIONS
C
	IMPLICIT NONE
C
	INTEGER I, IRET, N, LUGB, LUGI, J, NLDAS, JPDS(25),
     +          JGDS(22), KF, K, KPDS(25), KGDS(22), PARM(10), JRET
        INTEGER LVLTYPE(10),PARMLVL(10) 
        INTEGER L, HOUR, KK,KL,KN,KG 
C
	PARAMETER (NLDAS = 103936)
C
        INTEGER LAND_SEA(NLDAS)
	REAL F(NLDAS)

        REAL TAIR(NLDAS),PRES(NLDAS),SPFH(NLDAS)
        REAL DLWRF(NLDAS),UWIND(NLDAS),VWIND(NLDAS),CAPE(NLDAS)
        REAL CONVPREC(NLDAS),DSWRF(NLDAS)
        REAL NDAS_APCP(12,NLDAS),CONVPREC1(12,NLDAS)
        REAL FRAIN(NLDAS),PREC(NLDAS) 
C  	
        LOGICAL*1 LB(NLDAS) 	
        LOGICAL*1 LDASMASK(NLDAS)  
C
	CHARACTER OUTFILE*255, FILENAM*255, COUT*2, CLUGB*2
        CHARACTER*8  YESTERDAY, TODAY
        CHARACTER*11 ENVVAR, ENVVAR1
C  
C  CONSTANTS, DATA STATEMENTS,  AND PARAMETERS
	DATA PARM/11,1,204,51,205,33,34,157,154,63/
        DATA LVLTYPE /105,1,1,105,1,105,105,1,1,1/
        DATA PARMLVL /2,0,0,2,0,10,10,0,0,0/
C
C 11 - 2m air tempearture, 1-surface pressure, 204-DSWRF, 51-SPFH
C 205 -DLWRF, 33-UGRD, 34-VGRD, 157-CAPE,154-LSPA,63-ACPCP
C
C  BEGIN EXECUTABLE CODE
C READ NLDAS LAND-SEA MASK
        OPEN(99,form='unformatted')
        READ(99) LAND_SEA
        CLOSE(99)
C     DEFINE THE OUTPUT MASK FOR GRIB FILES.
        LDASMASK = .TRUE.
          WHERE (LAND_SEA .EQ. 0)
           LDASMASK = .FALSE.
        END WHERE
C READ BEFORE YESTESTERDAY AND YESTERDAY TO PRODUCE ACTUAL REALTIME 
C NLDAS-3 FORCING  
        OPEN(30, file = 'beforeyesterday')
        READ(30,'(A)') YESTERDAY
        CLOSE(30)

        OPEN(30, file = 'yesterday')
        READ(30,'(A)') TODAY
        CLOSE(30)
C
	DO I = 1, 24
C
C  FOR EACH HOUR, GET THE 16 FORCING GRIDS FROM THE 3 SEPARATE FILES
C  IN WHICH THEY ARE CONTAINED.  THESE THREE SEPARATE FILES  
C  ARE KNOWN BY THEIR UNIT NUMBERS (SEE DRIVING SCRIPT!).  THIS CODE
C  ASSUMES GRIDS 1-10 ARE IN THE NDAS FILES (25-48), GRIDS 11-14 COME
C  FROM THE GOES-BASED FILES (49-72), AND GRID 15 IS THE PRECIP FILE
C  (73-96).  
C
	  N = 1
   	  DO WHILE (N .LE. 10) 
	      LUGB = I+24
              JPDS(2)=84    ! NDAS/NAM DEFINITION
	    WRITE (CLUGB, '(i2.2)') LUGB
            ENVVAR1 = 'FORT'//CLUGB//char(0)
            CALL GET_ENVIRONMENT_VARIABLE(ENVVAR1,FILENAM)
C	    FILENAM = 'fort.'//CLUGB//char(0)
	    CALL BAOPEN (LUGB, FILENAM, IRET)
C ---------- CHECK IF NAM FILE IS OPENED SUCCESSFULY ---------------------
            IF(IRET.NE.0) THEN
            WRITE(*,*) 'CHECK NAM FILE OPEN ISSUE, STOP'
            STOP
            ENDIF
 
            J = 0
            LUGI = 0
            JPDS = -1
	    JPDS(5) = PARM(N)
            JPDS(6) = LVLTYPE(N)
            JPDS(7) = PARMLVL(N)
            IF(N.EQ.3.OR.N.EQ.5) THEN
            JPDS(16)=0       ! read simultaneous value
            ENDIF
	    CALL GETGB (LUGB, LUGI, NLDAS, J, JPDS, JGDS, KF, K, KPDS, 
     +                  KGDS, LB, F, IRET)

            IF(IRET.NE.0) THEN
            write(*,*) FILENAM
            WRITE(*,*) 'GETGB IRET=', IRET, 'N=',N
            WRITE(*,*) 'PLEASE CHECK NTH INDIVIDUAL FILE FOR MERGE PROCESS'  
            ENDIF
C  COPY DATA INTO ITS RESPECTIVE BINARY ARRAY. NOTE DATA ARE ASSIGNED 
C  TO PROPER ARRAY BASED ON LOOP INDEX (i).          

            IF (N.EQ.1) THEN
              DO L = 1, NLDAS
               TAIR(L) = F(L)
              END DO
            ELSE IF (N.EQ.2) THEN
              DO L = 1, NLDAS
                PRES(L) = F(L)
              END DO
            ELSE IF (N.EQ.3) THEN
            DO L=1,NLDAS
            DSWRF(L)=F(L)
            IF(DSWRF(L).LT.0.0) DSWRF(L)=0.0
            ENDDO
            ELSE IF (N.EQ.4) THEN
              DO L = 1, NLDAS
                SPFH(L) = F(L)
              END DO
            ELSE IF (N.EQ.5) THEN
            DO L=1,NLDAS
            DLWRF(L)=F(L)
            ENDDO
            ELSE IF (N.EQ.6) THEN
              DO L = 1, NLDAS
                UWIND(L) = F(L)
              END DO
            ELSE IF (N.EQ.7) THEN
              DO L = 1, NLDAS
                VWIND(L) = F(L)
              END DO
            ELSE IF (N.EQ.8) THEN
              DO L = 1, NLDAS
                CAPE(L) = F(L)
              END DO
            ELSE IF (N.EQ.9) THEN

C USING 12-HOUR BUCKET TO GET HOURLY PRECIPITATION
            KK=I-(INT(I/12)*12)
            IF(KK.EQ.0) KK=12
            DO L=1,NLDAS
            NDAS_APCP(KK,L)=F(L)
            ENDDO            
CI=1(13Z),12(01Z),AND 24(12Z)
            IF(KK.EQ.1) THEN
            DO L=1,NLDAS
            PREC(L)=NDAS_APCP(KK,L)
            IF(PREC(L).LT.0.0) PREC(L)=0.0
            ENDDO
            ELSE
            DO L=1,NLDAS
            PREC(L)=NDAS_APCP(KK,L)-NDAS_APCP(KK-1,L)
            IF(PREC(L).LT.0.0) PREC(L)=0.0
            ENDDO
            ENDIF

            ELSE IF (N.EQ.10) THEN

C USING 12-HOUR BUCKET TO GET HOURLY PRECIPITATION
            KG=I-(INT(I/12)*12)
            IF(KG.EQ.0) KG=12
            DO L=1,NLDAS
            CONVPREC1(KG,L)=F(L)
            ENDDO
CI=1(13Z),13(01Z),AND 24(12Z)
            IF(KG.EQ.1) THEN
            DO L=1,NLDAS
            CONVPREC(L)=CONVPREC1(KG,L)
            ENDDO
            ELSE
            DO L=1,NLDAS
            CONVPREC(L)=CONVPREC1(KG,L)-CONVPREC1(KG-1,L)
            ENDDO
            ENDIF

C ----------- CALCULATE FRACTION OF CONVECTION PRECIP --------
              DO L=1,NLDAS
               IF(PREC(L).EQ.-999.0.OR.CONVPREC(L).EQ.-999.0) THEN
               FRAIN(L)=-999.0
               ELSE
                 IF (PREC(L).NE.0.0) THEN
                 FRAIN(L)=CONVPREC(L)/PREC(L)
                ELSE
                 FRAIN(L)=0.0
                END IF 
               END IF
               IF(FRAIN(L).LT.0.0) FRAIN(L)=0.0
               IF(FRAIN(L).GT.1.0) FRAIN(L)=1.0 
               ENDDO
               END IF  
C CLOSE FILE AND READ NEXT FILE
	      N = N + 1
	      CALL BACLOSE (LUGB, JRET)
	
	  END DO        ! FIELD FILE LOOP
C
C  WRITE GRID TO ITS HOURLY OUTPUT FILE (I).

          HOUR=I
          CALL GRIB1OUT(NLDAS,LDASMASK,HOUR,YESTERDAY,
     &        TODAY,TAIR,SPFH,PRES,UWIND,VWIND,DLWRF,
     &        FRAIN,CAPE,CONVPREC,PREC,DSWRF)
          
          HOUR=I
          CALL GRIB2OUT(NLDAS,LDASMASK,HOUR,YESTERDAY,
     &        TODAY,TAIR,SPFH,PRES,UWIND,VWIND,DLWRF,FRAIN,CAPE,
     &        CONVPREC,PREC,DSWRF)

	END DO             ! HOUR LOOP
C
	STOP
	END

