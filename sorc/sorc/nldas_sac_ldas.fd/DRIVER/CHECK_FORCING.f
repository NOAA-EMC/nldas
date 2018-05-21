      SUBROUTINE CHECK_FORCING(PRCP,NX,NY,NHOUR,LAND_SEA)

      IMPLICIT NONE

      INTEGER NX
      INTEGER NY
      INTEGER NHOUR
      REAL    PRCP(NX,NY,NHOUR)
      INTEGER LAND_SEA(NX,NY)

      INTEGER I, J, K

      DO J = 1,NY
         DO I = 1,NX
            IF (LAND_SEA(I,J) .EQ. 1) THEN
               DO K = 1, NHOUR
                  IF (PRCP(I,J,K) .LT. 0.0) THEN
                     PRCP(I,J,K) = 0.0
                     WRITE(*,*) 'Warning PRCP(',K,') < 0.0mm at ', I, J 
                  END IF
                  IF (PRCP(I,J,K) .GT. 50.0) THEN
                     PRCP(I,J,K) = 0.0
                     WRITE(*,*) 'Warning PRCP(',K,') > 50.0mm at ',I,J 
                  END IF
               END DO
            END IF
         END DO
      END DO
         
      END
