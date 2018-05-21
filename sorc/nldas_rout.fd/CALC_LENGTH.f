      SUBROUTINE CALC_LENGTH(line,length)

      IMPLICIT NONE

      CHARACTER*100 line
      INTEGER       length

      length = 1
 100  CONTINUE
      IF ((line(length:length)) .NE. ' ') THEN
         length = length + 1
      ELSE
         GOTO 200
      END IF
      GOTO 100
 200  CONTINUE

      length = length - 1

      END
