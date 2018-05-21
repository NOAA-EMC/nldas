      SUBROUTINE SKIP_COMMENTS(fn1,one_line)

      IMPLICIT NONE

      INTEGER fn1
      CHARACTER*100 one_line

      one_line(1:1)='#'
      DO WHILE (one_line(1:1) .EQ. '#')
      READ(fn1,'(A100)') one_line
      ENDDO
      
      END
