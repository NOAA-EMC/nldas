
      SUBROUTINE CALC_AREA(NX, NY, AREA, DT)

      IMPLICIT NONE

      INTEGER NX, NY
      REAL    AREA(NX,NY)
      REAL    PI, DPHI, RERD, E, SOUTH, LAT
      INTEGER DT

      PARAMETER (RERD  = 6378.136)
      PARAMETER (DPHI  = 0.125)
      PARAMETER (E     = 0.00669447)
      PARAMETER (SOUTH = 25.0625)
      INTEGER I, J

      PI = ATAN(1.0) * 4.0

      DO J = 1, NY
         DO I = 1, NX
            LAT = SOUTH + (J-1) * DPHI
            AREA(I,J) = (2.0*PI*RERD*DPHI/360.0)**2.0 * 
     &           COS((LAT)*PI/180.0)

         END DO
      END DO

      END
