      SUBROUTINE READCNTL(CNTRFL,ICE,DT,Z,SNOALB,INOAHETP)

      IMPLICIT NONE

C     READS A CONTROL FILE

      INTEGER ICE
      INTEGER INOAHETP

      REAL    DT
      REAL    Z
      REAL    SNOALB

      CHARACTER*11  CNTRFL
      CHARACTER*100 one_line
      
      OPEN (UNIT=99, FILE=CNTRFL, STATUS='OLD',
     &     FORM='FORMATTED')
      
      CALL SKIP_COMMENTS(99,one_line)
      READ(one_line,*) DT
      CALL SKIP_COMMENTS(99,one_line)
      READ(one_line,*) Z
      CALL SKIP_COMMENTS(99,one_line)
      READ(one_line,*) SNOALB
      CALL SKIP_COMMENTS(99,one_line)
      READ(one_line,*) ICE
      CALL SKIP_COMMENTS(99,one_line)
      READ(one_line,*) INOAHETP

      CLOSE(99)

      END
