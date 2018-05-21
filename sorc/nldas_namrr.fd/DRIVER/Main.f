	PROGRAM MERGEFORCE
C
C***********************************************************************
C  PROGRAM:  MERGEFORCE	     MERGE LDAS FORCING FILES INTO SINGLE FILE
C  PRGMMR:  MARSHALL         ORG:  W/NP20
C
C  ABSTRACT:  MERGE DAILY FORCING FIELDS FOR LDAS INTO A SINGLE FILE.  
C  SEPARATE JOBS CREATE THE NAMv4-BASED FORCING AND THE PRECIP FORCING, 
C  AND SEPARATE HOURLY FILES ARE WRITTEN OUT FOR EACH TYPE.  THIS CODE TAKES
C  THOSE THREE SEPARATE FILES, AND MERGES THEM INTO ONE FILE CONTAINING ALL THE C  GRIDS.  THE FILES THEN CONTAIN ALL THE LDAS FORCING FIELDS FOR A SINGLE HOUR.C  THIS CODE IS DESIGNED TO RUN DAILY AND PRODUCE 24 HOURLY FILES. 
C  
C  PROGRAM LOG
C  99-11-16 CURTIS M. UPDATED TO ADD STAGE IV TO FORCING FILES
C  16-10-06 XIA Y. UPDATED FOR NLDAS VERSION 3 TO ACHIEVE ACTUAL REALTIME RUN
C  17-05-15 XIA Y. DROP OFF GOES DOWNWARD SHORTWAVE RADIATION
C***********************************************************************
C
C  DECLARATIONS
C
	IMPLICIT NONE
C
	INTEGER I, IRET, N, LUGB, LUGI, J, NLDAS, JPDS(25),
     +          JGDS(22), KF, K, KPDS(25), KGDS(22), PARM(11), JRET
        INTEGER LVLTYPE(11),PARMLVL(11) 
        INTEGER L, HOUR 
C
	PARAMETER (NLDAS = 103936)
C
        INTEGER LAND_SEA(NLDAS)
	REAL F(NLDAS)

        REAL TAIR(NLDAS),PRES(NLDAS),NDAS_DSWRF(NLDAS),SPFH(NLDAS)
        REAL DLWRF(NLDAS),UWIND(NLDAS),VWIND(NLDAS),CAPE(NLDAS)
        REAL NDAS_APCP(NLDAS),CONVPREC(NLDAS),DSWRF(NLDAS)
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
	DATA PARM/11,1,204,51,205,33,34,157,154,63,61/
        DATA LVLTYPE /105,1,1,105,1,105,105,1,1,1,1/
        DATA PARMLVL /2,0,0,2,0,10,10,0,0,0,0/

C
C 11 - 2m air tempearture, 1-surface pressure, 204-DSWRF, 51-SPFH
C 205 -DLWRF, 33-UGRD, 34-VGRD, 157-CAPE,61-APCP,63-ACPCP,161-CPC GAUGE PRECIP
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
C  FOR EACH HOUR, GET THE 11 FORCING GRIDS FROM THE 2 SEPARATE FILES
C  IN WHICH THEY ARE CONTAINED.  THESE THREE SEPARATE FILES  
C  ARE KNOWN BY THEIR UNIT NUMBERS (SEE DRIVING SCRIPT!).  THIS CODE
C  ASSUMES GRIDS 1-10 ARE IN THE NAMV4 FILES (25-48),AND GRID 15 IS THE 
c  PRECIP FILE (73-96).  
C
	  N = 1
   	  DO WHILE (N .LE. 11)
            
            IF (N .LE. 10) THEN
	      LUGB = I+24
              JPDS(2)=84    ! NAMV4 DEFINITION
	    ELSE IF (N.EQ.11) THEN
	      LUGB = I+72
              JPDS(2)=155
	    END IF

	    WRITE (CLUGB, '(i2.2)') LUGB
            ENVVAR1 = 'FORT'//CLUGB//char(0)
            CALL GET_ENVIRONMENT_VARIABLE(ENVVAR1,FILENAM)
C	    FILENAM = 'fort.'//CLUGB//char(0)
	    CALL BAOPEN (LUGB, FILENAM, JRET)
            IF (N.LE.10.AND.JRET.NE.0) THEN
                    PRINT*, 'NDAS BAOPEN PROBLEM-STOP, JRET=', JRET
                    STOP  
            END IF
            
            J = 0
            LUGI = 0
            JPDS = -1
            JPDS(5) = PARM(N)
            JPDS(6) = LVLTYPE(N)
            JPDS(7) = PARMLVL(N)
            
C            write(*,*), 'N=', N, JPDS(5), JPDS(6), JPDS(7)
C USING instantaneous DOWNWARD RADIATION FROM NAMRR
            IF(N.EQ.3.OR.N.EQ.5) THEN
            JPDS(16)=0
            ENDIF

	    CALL GETGB (LUGB, LUGI, NLDAS, J, JPDS, JGDS, KF, K, KPDS, 
     +                  KGDS, LB, F, IRET)

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
              DO L = 1, NLDAS
                DSWRF(L) = F(L)
                IF(DSWRF(L).LT.0.0) DSWRF(L)=0.0 
              END DO
            ELSE IF (N.EQ.4) THEN
              DO L = 1, NLDAS
                SPFH(L) = F(L)
              END DO
            ELSE IF (N.EQ.5) THEN
              DO L = 1, NLDAS
                DLWRF(L) = F(L)
              END DO
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
              DO L = 1, NLDAS
                NDAS_APCP(L) = F(L)
              END DO
            ELSE IF (N.EQ.10) THEN
              DO L = 1, NLDAS
                CONVPREC(L) = F(L)
              END DO
C ----------- CALCULATE FRACTION OF CONVECTION PRECIP --------
              DO L=1,NLDAS
               IF(NDAS_APCP(L).EQ.-999.0.OR.CONVPREC(L).EQ.-999.0) THEN
               FRAIN(L)=-999.0
               ELSE
                 IF (NDAS_APCP(L).NE.0.0) THEN
                 FRAIN(L)=CONVPREC(L)/NDAS_APCP(L)
                ELSE
                 FRAIN(L)=0.0
                END IF 
               END IF
               IF(FRAIN(L).LT.0.0) FRAIN(L)=0.0
               IF(FRAIN(L).GT.1.0) FRAIN(L)=1.0 
               ENDDO
            ELSE IF (N.EQ.11) THEN
              IF (IRET.EQ.0) THEN
                DO L = 1, NLDAS
                 PREC(L)=F(L)
                END DO
              ELSE
                DO L = 1, NLDAS
                 PREC(L)=NDAS_APCP(L)
                END DO
              END IF
              END IF  
C CLOSE FILE AND READ NEXT FILE
	      N = N + 1
	      CALL BACLOSE (LUGB, JRET)
	
	  END DO        ! FIELD FILE LOOP
C
C  WRITE GRID TO ITS HOURLY OUTPUT FILE (I).
C
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

