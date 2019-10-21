      SUBROUTINE READCNTL(CNTRFL,ICE,DT,Z)

      IMPLICIT NONE

C     READS A CONTROL FILE

      INTEGER ICE

      REAL  DT
      REAL  Z

      CHARACTER*11  CNTRFL
      CHARACTER*100 one_line
      
      OPEN (UNIT=99, FILE=CNTRFL, STATUS='OLD',
     &     FORM='FORMATTED')
      
      CALL SKIP_COMMENTS(99,one_line)
      READ(one_line,*) DT
      CALL SKIP_COMMENTS(99,one_line)
      READ(one_line,*) Z
      CALL SKIP_COMMENTS(99,one_line)
      READ(one_line,*) ICE

      CLOSE(99)

      END
