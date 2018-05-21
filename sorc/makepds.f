      SUBROUTINE MAKEPDS(TODAY, YESTERDAY, HOUR, KPDS)
      
      IMPLICIT NONE

      include 'postproc.h'

      CHARACTER*8 YESTERDAY, TODAY
      INTEGER     KPDS(25), HOUR

C     SET TIME-RELATED KPDS OCTETS
      if(DEBUG) then
        write(*,*) 'kpds(16) = ', kpds(16)
        write(*,*) 'today = ', TODAY
        write(*,*) 'yester = ', YESTERDAY
        write(*,*) 'hour = ', HOUR
      endif
C     it's a time-step accumulated variable
      IF (KPDS(16) .NE. 0) THEN
        READ (TODAY,'(4(I2))') KPDS(21), KPDS(8), 
     &                                KPDS(9), KPDS(10)
        KPDS(11) = HOUR - 1
        IF (HOUR .LE. 0) THEN
          READ (YESTERDAY,'(4(I2))') KPDS(21), KPDS(8), 
     &                                    KPDS(9), KPDS(10)
          KPDS(11) = HOUR + 23
        END IF
C     it's a normal variable
      ELSE
        READ (TODAY,'(4(I2))') KPDS(21), KPDS(8), 
     &                         KPDS(9), KPDS(10)
        KPDS(11) = HOUR
      END IF
      
      IF (KPDS(8) .EQ. 0) THEN 
         KPDS(8) = 100
      ELSE
         KPDS(21) = KPDS(21) + 1
      END IF

      if(DEBUG) then
        write(*,*) 'kpds = ', kpds(21),KPDS(8),kpds(9),kpds(10),kpds(11)
      endif

      END
