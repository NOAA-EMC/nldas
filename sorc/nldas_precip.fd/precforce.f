	PROGRAM PRECFORCE
C
C***********************************************************************
C  PROGRAM:  PRECFORCE              MAKE PRECIP FORCING FIELDS FOR LDAS 
C  PRGMMR:  MARSHALL                ORG:  W/NP20
C
C  ABSTRACT:  MERGE STAGE IV HOURLY PRECIP DATA (MULTISENSOR ANALYSES) 
C  AND CPC DAILY TOTAL PRECIP (GAUGE ONLY) TO GENERATE HOURLY PRECIP
C  FORCING FIELDS TO BE USED FOR FORCING THE LDAS.  HOURLY STAGE IV
C  ESTIMATES ARE SUMMED OVER THE 24 HOUR PERIOD OF THE CPC ANALYSIS.  
C  WEIGHTS ARE DERIVED FOR EACH HOUR BY FINDING THE FRACTION OF 24 H
C  SUM OVER ALL HOURS. THESE HOURLY WEIGHTS, DERIVED ONLY FROM THE STAGE
C  IV, ARE THEN MULTIPLIED BY THE CPC DAILY TOTAL TO "SPLIT OUR" HOURLY
C  PRECIPITATION AMOUNTS.  NOTE THAT WHEN A STAGE IV DATA POINT IS    
C  MISSING, IT IS REPLACED WITH ITS NEAREST NEIGHBOR, OUT TO 240 KM.  
C  BEYOND THAT RANGE, NDAS HOURLY PRECIPS ARE USED.  PRACTICALLY, THIS
C  MEANS AREAS OUTSIDE THE EDGE OF THE STAGE IV DOMAIN, BUT INSIDE THE 
C  LDAS DOMAIN, ARE DOMINATED BY NDAS IN THE DERIVATION OF THE WEIGHTS.  
C
C  PROGRAM HISTORY LOG:
C  99-01-21  MARSHALL  ORIGINAL CODING
C  99-11-16  MARSHALL  ADDED STAGE II ON LDAS GRID TO OUTPUT FILES
C  02-02-01  LOHMANN   CHANGED TO PRISM BASED PRECIP
C  02-03-18  LOHMANN   CHANGED GRID on STAGE II
c  06-03-10  XIA       OUTPUT PROCESSED PRECIPITATION ONLY 
C***********************************************************************
C
C  DECLARATIONS
C
	IMPLICIT NONE
C
	INTEGER NLDAS, NX, NY, BOUND, JPDS(25), IP, IPOPT(20), I, N, 
     +          LUGB, IRET, LUGI, J, JGDS(22), KF, K, KPDS(25), KM,
     +          KGDSEDAS(22), KGDSS4(22), IBI, IBO, NO, DATE(25),
     +          KPDSOUT(24,25), X, Y, COUNT, DIST, XX, YY, NSTAGE4,
     +          KGDSLDAS(22), GDS(22), LENGDS, LDASGRD, JRET  
C
 	PARAMETER (NLDAS = 103936, NX = 464, NY = 224, BOUND=16,
     +             NSTAGE4 = 987601)
C
	REAL EDAS(NLDAS), S4IN(NSTAGE4), S4OUT(NLDAS), RLAT(NLDAS),
     +       RLON(NLDAS), S42DTEMP(NX,NY), EDAS2D(NX,NY), S42D(NX,NY),
     +       S4(24,NLDAS), SUM(NLDAS), WEIGHT(24,NLDAS), CPC(NLDAS),
     +       HOURPREC(NLDAS), HOURS4(NLDAS), WEIGHTSUM(NLDAS),
     +       CPCTMP(NX,NY), WTSUMTMP(NX,NY), WT2DTMP(24,NX,NY),
     +       WT2D(24,NX,NY), S4SAVE(24,NLDAS)
C
	LOGICAL*1 OCEANMASK(NLDAS), LB(NSTAGE4), LO(NLDAS), 
     +            MASKTEMP(NX,NY), FOUND, S4MASK(24,NLDAS),
     +            MASKS4(NLDAS)
	LOGICAL LEAP, S4AVAIL(24)
C
        CHARACTER  ENVVAR*11, ENVVAR1*11
	CHARACTER  FILENAM*255, CENT*2, CYR*2, CMON*2, CDAY*2, CLUGB*2,
     +             CHOUR*2
C
C  SOME CONSTANTS AND SUCH
C
	J = 0
	JPDS = -1
        IP = 3
	IPOPT(1) = 4
	IPOPT(2) = -1
        LDASGRD = 110
        IBI = 1
        KM = 1
C
C  MAKE A GDS FOR THE OUTPUT GRID
C
	CALL MAKGDS (LDASGRD,KGDSLDAS,GDS,LENGDS,IRET)
C
        DO N = 1, NLDAS
        SUM(N) = 0.0
        WEIGHTSUM(N)=0.0
        END DO

C  READ IN PRECIP GRID FROM EDAS-BASED FORCING FILE    
C
	DO I = 1, 24
        JPDS(5) = 154  ! NAMRR LSPA
C          print*, I
          LUGB = 10+I
          WRITE (CLUGB, '(i2.2)') LUGB
          ENVVAR = 'FORT'//CLUGB//char(0)
          CALL GET_ENVIRONMENT_VARIABLE(ENVVAR,FILENAM)
C          FILENAM = 'fort.'//CLUGB//char(0)

          CALL BAOPEN (LUGB, FILENAM, IRET)
          IF(IRET.NE.0) THEN
          WRITE(*,*) 'NDAS FILE OPEN ISSUE, STOP TO CHECK'
          STOP
          ENDIF
  
          CALL GETGB (LUGB,LUGI,NLDAS,J,JPDS,JGDS,KF,K,KPDS,KGDSEDAS,
     +                OCEANMASK,EDAS,IRET)
          IF(IRET.NE.0) THEN
           WRITE(*,*) 'NDAS FILE IRET=', IRET, 'check ndas precipitation'
          ENDIF

	  CALL BACLOSE (LUGB, JRET)
C          print*, 'got ndas'
          do n=1,nldas
           if(edas(n).lt.0)then
           edas(n)=0
           endif
          enddo
C
C  SAVE PDS INFO FOR USE IN OUTPUT FILES.  
C
	  DO N = 1, 25
	    KPDSOUT(I,N) = KPDS(N)
          END DO
C
C   READ IN THE STAGE IV GRID AND INTERPOLATE TO LDAS GRID.  IF ENTIRE
C   STAGE IV GRID IS MISSING (I.E., THE HOURLY FILE IS MISSING) THEN
C   REPLACE GRID WITH EDAS GRID.  ELSE, PROCEED AND FILL HOLES IN 
C   STAGE IV GRID WITH EDAS VALUES BELOW.
C
          JPDS(5) = 61  ! STAGE II precipitation
	  LUGB = 34+I
          WRITE (CLUGB, '(i2.2)') LUGB
C          FILENAM = 'fort.'//CLUGB//char(0)
          ENVVAR1 = 'FORT'//CLUGB//char(0)
          CALL GET_ENVIRONMENT_VARIABLE(ENVVAR1,FILENAM)
          CALL BAOPEN (LUGB, FILENAM, IRET)
          CALL GETGB (LUGB,LUGI,NSTAGE4,J,JPDS,JGDS,KF,K,KPDS,KGDSS4,
     +                LB,S4IN,IRET)
          IF (IRET.NE.0) THEN
          S4AVAIL(I) = .FALSE.
          print *, 'no stage4, IRET=', IRET
          DO N = 1, NLDAS
            S4(I,N)=EDAS(N)
          END DO
          ELSE 
          S4AVAIL(I) = .TRUE.
C	  print*, 'got stage4'
          CALL IPOLATES (IP,IPOPT,KGDSS4,KGDSLDAS,NSTAGE4,NLDAS,KM,IBI,
     +                   LB,S4IN,NO,RLAT,RLON,IBO,LO,S4OUT,IRET)
          END IF
          CALL BACLOSE (LUGB,JRET)
C
C  SAVE S4OUT AND MASK FOR LATER OUTPUT IF STAGE IV AVAILABLE.
C
          IF (S4AVAIL(I)) THEN
          DO N = 1, NLDAS
            S4SAVE(I,N)=S4OUT(N)
            S4MASK(I,N)=LO(N) 
          END DO
C
C  PUT 1D ARRAYS OF DATA ONTO 2D LDAS GRID (FOR NEIGHBOR SEARCHING 
C  ALGORITHM BELOW).
C
          COUNT = 0
          DO Y = 1, NY
            DO X = 1, NX
              S42DTEMP(X,Y) = S4OUT(X+COUNT)
              EDAS2D(X,Y) = EDAS(X+COUNT)
              MASKTEMP(X,Y) = LO(X+COUNT)
            END DO
            COUNT = COUNT + NX
          END DO
C
C  REPLACE MISSING DATA POINTS WITH NEAREST NEIGHBOR IF NEAREST NEIGHBOR
C  WITHIN ABOUT 240 KM ("BOUND").  ELSE, REPLACE WITH EDAS VALUE.
C 
	  DO Y = 1, NY
            DO X = 1, NX
              FOUND = .TRUE.
	      IF (.NOT. MASKTEMP(X,Y)) THEN
                FOUND = .FALSE.
                IF ((X.LE.(NX-BOUND)).AND.(X.GT.BOUND).AND.   
     +              (Y.LE.(NY-BOUND)).AND.(Y.GT.BOUND)) THEN
                  DIST = 1
                  DO WHILE ((DIST .LE. BOUND).AND.(.NOT. FOUND))    
                    YY = -1*DIST
                    DO WHILE ((YY .LE. DIST).AND.(.NOT. FOUND))
                      XX = -1*DIST
                      DO WHILE ((XX .LE. DIST).AND.(.NOT. FOUND))
	                IF (MASKTEMP(X+XX,Y+YY)) THEN 
	                  S42D(X,Y) = S42DTEMP(X+XX,Y+YY)
	                  FOUND = .TRUE.
	                END IF
	                XX = XX+1
	              END DO
                      YY = YY+1
	            END DO
                  DIST = DIST+1
	          END DO
                ELSE
                  S42D(X,Y)=EDAS2D(X,Y)        
                  FOUND = .TRUE.
                END IF
              ELSE
	        S42D(X,Y)=S42DTEMP(X,Y)	
	      END IF
              IF (.NOT. FOUND) THEN
                S42D(X,Y)=EDAS2D(X,Y)
              END IF
            END DO
          END DO
C
C  WRITE THE HOLE-FILLED GRID BACK OUT INTO A 1D ARRAY
C
	  COUNT = 0
	  DO Y = 1, NY
	    DO X = 1, NX
	      S4(I,X+COUNT) = S42D(X,Y)
            END DO
	    COUNT = COUNT+NX
	  END DO 
          END IF
C
C  SUM HOURLY STAGE IV PRECIPS (WITH NO MISSING DATA) OVER 24H PERIOD  
C
	  DO N = 1, NLDAS
	    SUM(N) = SUM(N)+S4(I,N)
	  END DO
	END DO
        print*, 'sums complete'

C  This call to convert global product to NLDAS domain

	CALL GRIDCPC (CPC,NLDAS) 

	DO I = 1, NLDAS
	   CPC(I) = CPC(I) * 25.4
C If CPC missing using STAGE IV or EDAS data 
 	   IF (CPC(I) .LT. 0.0) CPC(I) = SUM(I) 
	END DO
C
C  COMPUTE WEIGHT FOR EACH HOUR.  WATCH OUT FOR SINGULARITY!  FIND SUM
C  OF HOURLY WEIGHTS.  SHOULD ALWAYS BE ZERO OR ONE!  
C
	DO I = 1, 24
          DO N = 1, NLDAS
            IF (SUM(N).GT.0) THEN
	      WEIGHT(I,N) = S4(I,N)/SUM(N)
            ELSE
              WEIGHT(I,N) = 0
            END IF
            WEIGHTSUM(N)=WEIGHTSUM(N)+WEIGHT(I,N)
          END DO
 	END DO
        print*, 'got weights'
C
C  PUT WEIGHTS AND WEIGHTSUM ON 2D ARRAYS FOR SEARCHING ALGORITHM BELOW.
C
        COUNT = 0
        DO Y = 1, NY
          DO X = 1, NX
            CPCTMP(X,Y) = CPC(X+COUNT)
            WTSUMTMP(X,Y) = WEIGHTSUM(X+COUNT)
            DO I = 1, 24
              WT2DTMP(I,X,Y)=WEIGHT(I,X+COUNT)
            END DO
          END DO
          COUNT = COUNT + NX
        END DO
C
C  SEARCH FOR OCCURRENCES WHERE WEIGHTSUM IS ZERO, BUT CPC PRECIP IS
C  NON-ZERO.  IN THESE INSTANCES, TAKE HOURLY WEIGHTS FROM NEAREST NEIGHBOR.
C  IN THIS MANNER, THE SUM OF THE 24 HOURLY PRECIPS TO BE DERIVED BELOW
C  WILL ALWAYS EQUAL THE CPC 24H ANALYSIS AT ALL GRID POINTS.  
C
	DO Y = 1, NY
          DO X = 1, NX
            FOUND = .TRUE.
            IF ((WTSUMTMP(X,Y).EQ.0).AND.(CPCTMP(X,Y).NE.0)) THEN
              FOUND = .FALSE.
              IF ((X.LE.(NX-BOUND)).AND.(X.GT.BOUND).AND.   
     +            (Y.LE.(NY-BOUND)).AND.(Y.GT.BOUND)) THEN
                DIST = 1
                DO WHILE ((DIST .LE. BOUND).AND.(.NOT. FOUND))    
                  YY = -1*DIST
                  DO WHILE ((YY .LE. DIST).AND.(.NOT. FOUND))
                    XX = -1*DIST
                    DO WHILE ((XX .LE. DIST).AND.(.NOT. FOUND))
	              IF (WTSUMTMP(X+XX,Y+YY).NE.0) THEN 
                        DO I = 1, 24
	                  WT2D(I,X,Y) = WT2DTMP(I,X+XX,Y+YY)
	                END DO
	                FOUND = .TRUE.
	              END IF
	              XX = XX+1
	            END DO
                    YY = YY+1
	          END DO
                  DIST = DIST+1
	        END DO
              ELSE
                DO I = 1, 24
                  WT2D(I,X,Y)=1/24        
                END DO
                FOUND = .TRUE.
              END IF
            ELSE
              DO I = 1, 24
                WT2D(I,X,Y)=WT2DTMP(I,X,Y)
              END DO
	    END IF
            IF (.NOT. FOUND) THEN
              DO I = 1, 24
                WT2D(I,X,Y) = 1/24
              END DO 
            END IF
          END DO
        END DO
C
C  NOW WRITE WEIGHTS BACK OUT INTO 1D ARRAY
C
        COUNT = 0
        DO Y = 1, NY
          DO X = 1, NX
            DO I = 1, 24
              WEIGHT(I,X+COUNT)=WT2D(I,X,Y)
            END DO
          END DO
          COUNT = COUNT + NX
        END DO
        print*, 'done checking weights'
C
C  COMPUTE HOURLY PRECIPS FROM WEIGHTS AND CPC TOTAL
C
	DO I = 1, 24
          DO N = 1, NLDAS
	    HOURPREC(N)=WEIGHT(I,N)*CPC(N)
            HOURS4(N)=S4SAVE(I,N)
            MASKS4(N)=S4MASK(I,N) 
          END DO 
C 
C  GET PDS INFO, OPEN FILE AND WRITE DATA TO HOURLY GRIB FILES.  
C 
          LUGB = 59+I 
          DO N = 1, 25 
            KPDS(N) = KPDSOUT(I,N) 
            if (KPDS(8).eq.0) then
              KPDS(8) = 100
            end if
	    DATE(N) = KPDS(N)
          END DO
	  LEAP = (MOD(DATE(8),4).EQ.0)
	  DATE(11) = DATE(11)+1
	  CALL UPDATE(DATE, LEAP)
          if (DATE(8).eq.100)then
            WRITE (CENT, '(i2.2)') DATE(21)
            WRITE (CYR, '(i2.2)') 0 
          else
            WRITE (CENT, '(i2.2)') DATE(21)-1
            WRITE (CYR, '(i2.2)') DATE(8)
          end if
          WRITE (CMON, '(i2.2)') DATE(9)
          WRITE (CDAY, '(i2.2)') DATE(10)
          WRITE (CHOUR, '(i2.2)') DATE(11)
          FILENAM  = CENT//CYR//CMON//CDAY//CHOUR//'.PREC'
C          print*, FILENAM,S4AVAIL(I)
          KPDS(2) = 155
          KPDS(5) = 61
          KPDS(22) = 2
          CALL BAOPEN (LUGB, FILENAM, IRET)
          CALL PUTGB (LUGB,NLDAS,KPDS,KGDSLDAS,OCEANMASK,HOURPREC,IRET)
          IF (IRET .NE. 0) THEN
            PRINT*, 'PUTGB LDAS PREC FAILED, IRET=', IRET
          END IF
C
          CALL BACLOSE (LUGB, JRET)
	END DO
C
	STOP
	END

	  
