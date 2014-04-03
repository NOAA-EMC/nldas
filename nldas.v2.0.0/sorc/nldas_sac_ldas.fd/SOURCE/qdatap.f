      SUBROUTINE QDATAP (T,P,RH,QD,QS,ES) 

      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  PURPOSE:  OBTAIN SPECIFIC HUMIDITY (q) FROM RELATIVE HUMIDITY 
CC            AND GIVEN PRESSURE AND TEMPERATURE.
CC            
CC
CC            FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989: 'A 
CC            SHORT COURSE IN CLOUD PHYSICS', PERGAMON PRESS, 3rd ED.
CC
CC                                   Pablo J. Grunmann, 3/6/98.
CC                Updated to eliminate subroutine SVP, 6/24/98.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C----------------------------------------
C In:
C        T    Temperature (K)
C        P    Pressure (Pa)
C        RH   Relative humidity (%)
C-----------------------------------------
C Out:
C        QD   Specific humidity (Kg/Kg)
C        QS   Saturation Specific humidity (Kg/Kg)
C        ES   Saturation vapor pressure for water (Pa)
C----------------------------------------
      REAL T
      REAL P
      REAL RH
      REAL RHF
      REAL QD
      REAL QS
      REAL ES
      REAL EP
      REAL EPS
      REAL E

      PARAMETER (eps=0.622 )
C 
C     ABOUT THE PARAMETER:
C      
C     eps ---------- (Water)/(dry air) molecular mass ratio, epsilon
C _____________________________________________________________________
C
C    function E(T) = Sat. vapor pressure (in Pascal) at 
C                    temperature T (uses Clausius-Clapeyron).
          Es = E(T)
C  CONVERT REL. HUMIDITY (%) TO THE FRACTIONAL VALUE
          RHF = RH/100.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC      CALCULATE SATURATION MIXING RATIO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      QS = 0.622 * ES /P   was substituted by a more precise
C formula:                              -PABLO J. GRUNMANN, 05/28/98.
        QS = 0.622 * ES /(P - (1.-0.622)*ES)

C
C  CONVERSION FROM REL. HUMIDITY:
C     (Rogers, pg. 17)
C
        EP = (P*Es*RHF)/(P - Es*(1. - RHF))
        QD = eps*EP/(P - (1. - eps)*EP)
C     
          RETURN
          END
