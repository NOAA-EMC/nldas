      SUBROUTINE MAKEPDS(TODAY, KPDS, HOUR)
      
      IMPLICIT NONE

      CHARACTER*8 TODAY
      INTEGER     KPDS(25), HOUR

C     SET TIME-RELATED KPDS OCTETS
            READ (TODAY,'(4(I2))') KPDS(21), KPDS(8), 
     &           KPDS(9), KPDS(10)
            KPDS(11) = HOUR-1
      
      IF (KPDS(8) .EQ. 0) THEN 
         KPDS(8) = 100
      ELSE
         KPDS(21) = KPDS(21) + 1
      END IF
      END
