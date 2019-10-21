      SUBROUTINE QDATAP (T,P,QS) 

      IMPLICIT NONE

C----------------------------------------
C In:
C	T    Temperature (K)
C	P    Pressure (Pa)
C	QD   Specific humidity (Kg/Kg)
C Out:
C	QS   Saturation Specific humidity (Kg/Kg)
C----------------------------------------
      REAL  T
      REAL  P
      REAL  QS
      REAL  ES
      REAL  E

C     ABOUT THE PARAMETER:
C      
C     eps ---------- (Water)/(dry air) molecular mass ratio, epsilon
C _____________________________________________________________________
C
C    function E(T) = Sat. vapor pressure (in Pascal) at 
C                    temperature T (uses Clausius-Clapeyron).

      ES = E(T)
      QS = 0.622 * ES /(P - (1.-0.622)*ES)

      END
