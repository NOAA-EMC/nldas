      SUBROUTINE MAKEPDS(YESTERDAY, TODAY, KPDS, HOUR)
      
      IMPLICIT NONE

      CHARACTER*8 YESTERDAY, TODAY
      INTEGER     KPDS(25), HOUR

C     SET TIME-RELATED KPDS OCTETS
      IF (KPDS(16) .NE. 0) THEN

         IF(HOUR.GE.13.AND.HOUR.LE.23) THEN 
            READ (YESTERDAY,'(4(I2))') KPDS(21), KPDS(8), 
     &           KPDS(9), KPDS(10)
            KPDS(11) = HOUR-1
          ENDIF

          IF(HOUR.EQ.0) THEN
            READ (YESTERDAY,'(4(I2))') KPDS(21), KPDS(8),
     &           KPDS(9), KPDS(10)
            KPDS(11) = 23
           ENDIF

           IF(HOUR.GE.1.AND.HOUR.LE.12) THEN
            READ (TODAY,'(4(I2))') KPDS(21), KPDS(8),
     &           KPDS(9), KPDS(10)
            KPDS(11) = HOUR-1
            ENDIF
      ELSE
         IF (HOUR.GE.0.AND.HOUR.LE.12) THEN 
            READ (TODAY,'(4(I2))') KPDS(21), KPDS(8),
     &           KPDS(9), KPDS(10)
            KPDS(11) = HOUR
           ENDIF 
         IF(HOUR.GE.13.AND.HOUR.LE.23) THEN 
         READ (YESTERDAY,'(4(I2))') KPDS(21), KPDS(8), 
     &        KPDS(9), KPDS(10)
         KPDS(11) = HOUR
         ENDIF

      END IF
      
      IF (KPDS(8) .EQ. 0) THEN 
         KPDS(8) = 100
      ELSE
         KPDS(21) = KPDS(21) + 1
      END IF

      END
