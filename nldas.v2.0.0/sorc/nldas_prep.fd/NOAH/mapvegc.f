      SUBROUTINE MAPVEGC(N,VEGT) 

!  This subroutine converts the UMD classes to the SIB classes 
!  used by NOAH LSM (v 2.5).
!  (Originally from Dag Lohmann at NCEP)
!
!  28 Apr 2002, K Arsenault:  Added NOAH LSM to LDAS.

      IMPLICIT NONE

!  Local Variables

      INTEGER :: N
      INTEGER :: VEGT, SIBVEG

!  Convert UMD Classes to SIB Classes.

        IF (VEGT .EQ. 1)  SIBVEG = 4
        IF (VEGT .EQ. 2)  SIBVEG = 1
        IF (VEGT .EQ. 3)  SIBVEG = 5
        IF (VEGT .EQ. 4)  SIBVEG = 2
        IF (VEGT .EQ. 5)  SIBVEG = 3
        IF (VEGT .EQ. 6)  SIBVEG = 3
        IF (VEGT .EQ. 7)  SIBVEG = 6
        IF (VEGT .EQ. 8)  SIBVEG = 8
        IF (VEGT .EQ. 9)  SIBVEG = 9
        IF (VEGT .EQ. 10) SIBVEG = 7
        IF (VEGT .EQ. 11) SIBVEG = 12
        IF (VEGT .EQ. 12) SIBVEG = 11
        IF (VEGT .EQ. 13) SIBVEG = 11
        IF (VEGT .GT. 13) THEN
            SIBVEG = 7
        END IF

      VEGT=SIBVEG

      RETURN
      END 

