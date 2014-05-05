      SUBROUTINE CHECK_INITIAL(NX,NY,NSOLD,LAND_SEA,SMC,
     &     SOILTYP,VEGTYP,SLOPETYP,ALBEDO,SHDFAC,
     &     NSOIL,NMONTH,SOILDEPTH)

      IMPLICIT NONE

C     READS INITIAL CONDITIONS

      INTEGER NX
      INTEGER NY
      INTEGER NSOLD
      INTEGER NMONTH
      INTEGER LAND_SEA(NX,NY)
      INTEGER SOILTYP(NX,NY)
      INTEGER VEGTYP(NX,NY)
      INTEGER SLOPETYP(NX,NY)
      INTEGER NSOIL(NX,NY)

      REAL    SOILDEPTH(NSOLD,NX,NY)
      REAL    SMC(NSOLD,NX,NY)

      REAL    ALBEDO(NMONTH,NX,NY)
      REAL    SHDFAC(NMONTH,NX,NY)

      REAL    MAXSMC(9)
      REAL    WLTSMC(9)
      REAL    MAX_M, MIN_M

      INTEGER I, J, K

      DATA MAXSMC/0.421, 0.464, 0.468, 0.434, 0.406, 0.465, 0.404,
     &     0.439, 0.421/
      DATA WLTSMC/0.029, 0.119, 0.139, 0.047, 0.020, 0.103, 0.069,
     &     0.066, 0.029/

      DO J = 1,NY
         DO I = 1,NX
            IF (LAND_SEA(I,J) .EQ. 1) THEN
               
               IF (SOILTYP(I,J) .LE. 0) THEN
                  WRITE(*,*) 'SOILTYP is smaller than 1, ', I, J
                  SOILTYP(I,J) = 2
               END IF
               IF (SOILTYP(I,J) .GE. 10) THEN
                  WRITE(*,*) 'SOILTYP is larger than 10, ', I, J
                  SOILTYP(I,J) = 2
               END IF
               MAX_M = MAXSMC(SOILTYP(I,J))
               MIN_M = WLTSMC(SOILTYP(I,J))
               
               IF (VEGTYP(I,J) .LE. 0) THEN
                  WRITE(*,*) 'VEGTYP is smaller than 1, ', I, J
                  VEGTYP(I,J) = 7
               END IF
               IF (VEGTYP(I,J) .GE. 14) THEN
                  WRITE(*,*) 'VEGTYP is larger than 13, ', I, J
                  VEGTYP(I,J) = 7
               END IF
               
               IF (SLOPETYP(I,J) .LE. 0) THEN
                  WRITE(*,*) 'SLOPETYP is smaller than 1, ', I, J
                  SLOPETYP(I,J) = 1
               END IF
               IF (SLOPETYP(I,J) .GE. 10) THEN
                  WRITE(*,*) 'SLOPETYP is larger than 9, ', I, J
                  SLOPETYP(I,J) = 1
               END IF
               
               DO K = 1, NMONTH
                  IF (ALBEDO(K,I,J) .LT. 0.0) THEN
                     ALBEDO(K,I,J) = 0.0
                     WRITE(*,*) 'Warning ALBEDO < 0.0 at ', I, J 
                  END IF
                  IF (ALBEDO(K,I,J) .GT. 1.0) THEN
                     ALBEDO(K,I,J) = 1.0
                     WRITE(*,*) 'Warning ALBEDO > 1.0 at ', I, J 
                  END IF
                  
                  IF (SHDFAC(K,I,J) .LT. 0.0) THEN
                     SHDFAC(K,I,J) = 0.0
                     WRITE(*,*) 'Warning SHDFAC < 0.0 at ', I, J 
                  END IF
                  IF (SHDFAC(K,I,J) .GT. 1.0) THEN
                     SHDFAC(K,I,J) = 1.0
                     WRITE(*,*) 'Warning SHDFAC > 1.0 at ', I, J 
                  END IF
               END DO

c               IF (NSOIL(I,J) .LT. 1) THEN
c                  NSOIL(I,J) = 4
c                  WRITE(*,*) 'Warning NSOIL < 1 at ', I, J 
c               END IF
c               IF (NSOIL(I,J) .GT. NSOLD) THEN
c                  NSOIL(I,J) = NSOLD
c                  WRITE(*,*) 'Warning NSOIL >',NSOLD,' at ', I, J 
c               END IF

C               DO K = 1, NSOLD
C                  IF (SMC(K,I,J) .LT. MIN_M) THEN
C                     SMC(K,I,J) = MIN_M
C                     WRITE(*,*) 'Warning SMC(',K,') <',MIN_M,' at ',I,J 
C                  END IF
C                  IF (SMC(K,I,J) .GT. MAX_M) THEN
C                     SMC(K,I,J) = MAX_M
C                     WRITE(*,*) 'Warning SMC(',K,') >',MAX_M,' at ',I,J 
C                  END IF
C               END DO

c               DO K = 1, NSOLD
c                  IF (SOILDEPTH(K,I,J) .LT. 0.0) THEN
c                     WRITE(*,*) 'SOILDEPTH(',K,') < 0.0)'
c                     PAUSE
c                  END IF
c               END DO

            END IF
         END DO
      END DO
         
      END
