
      PROGRAM CALC_AREA

      INTEGER NX, NY
      REAL    AREA, AREA2
      REAL    PI, DPHI, RERD, E, SOUTH, LAT

      PARAMETER (RERD  = 6378136.0)
      PARAMETER (DPHI  = 0.125)
      PARAMETER (E     = 0.00669447)
      PARAMETER (SOUTH = -38.0)
      PARAMETER (NX = 1)
      PARAMETER (NY = 1)

      INTEGER I, J

      PI = ATAN(1.0) * 4.0

      DO J = 1, NY
         DO I = 1, NX
            LAT = SOUTH + (J-1) * DPHI
            AREA   = (2.0*PI*RERD*DPHI/360.0)**2.0 * 
     &           COS((LAT)*PI/180.0) / (1000.0*1000.0)

            AREA2 = (2.0*PI*RERD*DPHI/360.0)**2.0 * 
     &           COS((LAT)*PI/180.0) / 1000**2.0

            WRITE(*,*) LAT, AREA, AREA2
         END DO
      END DO

      END
