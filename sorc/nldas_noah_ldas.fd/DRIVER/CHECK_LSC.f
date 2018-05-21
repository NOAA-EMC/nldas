      SUBROUTINE CHECK_LSC(NX,NY,NSOLD,NMONTH,NSOIL,SOILTYP,
     &     LAND_SEA,VEGTYP,SLOPETYP,SOILDEPTH,ALBEDO,SHDFAC,
     &     TBOT,MAXSNOWALB,PRCP_MASK)

      IMPLICIT NONE

C     CHECKS LAND SURFACE CHARACTERISTICS

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
      REAL    TBOT(NX,NY)
      REAL    MAXSNOWALB(NX,NY)
      REAL    ALBEDO(NMONTH,NX,NY)
      REAL    SHDFAC(NMONTH,NX,NY)
      REAL    PRCP_MASK(NX,NY)

      INTEGER I, J, K

      DO J = 1,NY
         DO I = 1,NX
            IF (LAND_SEA(I,J) .EQ. 1) THEN
               IF (SOILTYP(I,J) .LE. 0) THEN
                  WRITE(*,*) 'SOILTYP is smaller than 1, ', I, J
                  SOILTYP(I,J) = 2
               END IF
               IF (SOILTYP(I,J) .GE. 19) THEN
                  WRITE(*,*) 'SOILTYP is larger than 19, ', I, J
                  SOILTYP(I,J) = 2
               END IF
               IF (VEGTYP(I,J) .LE. 0) THEN
                  WRITE(*,*) 'VEGTYP is smaller than 1, ', I, J, VEGTYP(i,j)
                  VEGTYP(I,J) = 12
               END IF
               IF (VEGTYP(I,J) .GE. 14) THEN
                  WRITE(*,*) 'VEGTYP is larger than 14, ', I, J,VEGTYP(i,j)
                  VEGTYP(I,J) = 12
               END IF
               
               IF (SLOPETYP(I,J) .LE. 0) THEN
                  WRITE(*,*) 'SLOPETYP is smaller than 1, ', I, J
                  SLOPETYP(I,J) = 1
               END IF
               IF (SLOPETYP(I,J) .GE. 10) THEN
                  WRITE(*,*) 'SLOPETYP is larger than 9, ', I, J
                  SLOPETYP(I,J) = 1
               END IF
               
               IF (TBOT(I,J) .LT. 230.0) THEN
                  TBOT(I,J) = 230.0
                  WRITE(*,*) 'Warning TBOT < 230.0K at ', I, J 
               END IF
               IF (TBOT(I,J) .GT. 350.0) THEN
                  TBOT(I,J) = 350.0
                  WRITE(*,*) 'Warning TBOT > 350.0K at ', I, J 
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

               IF (NSOIL(I,J) .LT. 1) THEN
                  NSOIL(I,J) = 4
                  WRITE(*,*) 'Warning NSOIL < 1 at ', I, J 
               END IF
               IF (NSOIL(I,J) .GT. NSOLD) THEN
                  NSOIL(I,J) = NSOLD
                  WRITE(*,*) 'Warning NSOIL >',NSOLD,' at ', I, J 
               END IF

               DO K = 1, NSOLD
                  IF (SOILDEPTH(K,I,J) .LT. 0.0) THEN
                     WRITE(*,*) 'SOILDEPTH(',K,') < 0.0)'
                     PAUSE
                  END IF
               END DO

               IF (MAXSNOWALB(I,J) .LT. 0.0) THEN
                  MAXSNOWALB(I,J) = 0.0
                  WRITE(*,*) 'Warning MAXSNOWALB < 0.0 at ', I, J 
               END IF
               IF (MAXSNOWALB(I,J) .GT. 1.0) THEN
                  MAXSNOWALB(I,J) = 1.0
                  WRITE(*,*) 'Warning MAXSNOWALB > 1.0 at ', I, J 
               END IF               

               IF (PRCP_MASK(I,J) .LT. 0.0) THEN
                  PRCP_MASK(I,J) = 0.0
                  WRITE(*,*) 'Warning PRCP_MASK < 0.0 at ', I, J 
               END IF
               IF (PRCP_MASK(I,J) .GT. 1.0) THEN
                  PRCP_MASK(I,J) = 1.0
                  WRITE(*,*) 'Warning PRCP_MASK > 1.0 at ', I, J 
               END IF               

            END IF
         END DO
      END DO

      END


