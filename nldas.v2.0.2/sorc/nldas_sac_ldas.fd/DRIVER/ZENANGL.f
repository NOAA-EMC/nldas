      SUBROUTINE ZENANGL(NT,IT,JULDAY,DT,XLAT,ZEN,ZENAVG,DAYLIGHT)
C
C     NT - Number of model time steps in a day
C     IT - Time step index from local 0 hour.
C     JULDAY - Julian day index
C     DT - Time step in hours
C     ZEN - Zenith angle
C     ZENAVG - Average Zenith angle
C     DAYLIGHT - Daylight in hours
C
      IMPLICIT NONE
      INTEGER NT,IT,JULDAY,K
      REAL COSZEN(NT),DT,XLAT,ZEN,ZENAVG,DAYLIGHT,PI,PID12,PID180
      REAL SLAT,CLAT,SDEC,SDEL,CDEL,TIME,HRANGL,ZENSUM
      DATA PI     /3.14159265/
      DATA PID12  /0.2617993878/
      DATA PID180 /0.0174532925/
C
      SLAT = SIN (XLAT * PID180)
      CLAT = COS (XLAT * PID180)
      SDEC = 23.5 * PID180 * SIN(2.0 * PI * (JULDAY-80.0)/365.0)
      SDEL = SIN (SDEC)
      CDEL = COS (SDEC)
      DAYLIGHT = 0.0
      ZENSUM = 0.0
      DO K = 1, NT
        TIME = (K-1) * DT
        HRANGL = (12.0-TIME) * PID12
        COSZEN(K) = (SLAT * SDEL) + ( CLAT * CDEL * COS(HRANGL) )
        IF (COSZEN(K) .GE. 0.01) THEN
          DAYLIGHT = DAYLIGHT + DT
          ZENSUM = ZENSUM + DT * COSZEN(K)
        ELSE 
          COSZEN(K) = 0.0
        END IF
      END DO
      ZENAVG = ZENSUM / DAYLIGHT
      ZEN = COSZEN(IT)
C
      RETURN
      END
