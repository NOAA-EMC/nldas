      SUBROUTINE SKIP_COMMENTS(fn1,one_line)

      IMPLICIT NONE

      INTEGER fn1
      CHARACTER*100 one_line
      
 100  CONTINUE
      READ(fn1,'(A100)',end=200) one_line
      IF (one_line(1:1) .NE. '#') THEN
         GOTO 200
      END IF
      GOTO 100
 200  CONTINUE

      END
