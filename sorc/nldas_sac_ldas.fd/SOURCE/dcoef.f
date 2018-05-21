      SUBROUTINE DCOEF ( Z, Z0, T1V, TH2V, SFCSPD, CM, CH )

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    NAME:  DETERMINE COEFFICIENTS (DCOEF)       VERSION: N/A
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL A
      REAL AH
      REAL AM
      REAL BETA
      REAL B1
      REAL B2

      DATA B1     / 9.4 /
      DATA B2     / 15.0 /
      DATA CUS    / 7.4 /
      DATA EXMCH  / -1. /
      DATA G      / 9.806 /
C
C      DATA PR     / .74 /
C   Set PR=1 the same as in Ek's version for PILPS 1994
C
      DATA PR     / 1.0 /
      DATA VK     / .4 /

C
C   Set Z0H as function of Z0 as in Beljaars and Betts (1992?)
C
      Z0H = Z0/10.0
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CALC A FRICTION VELOCITY FOR USE IN CALCULATING THE DRAG
C     COEF FOR MOMENTUM AND ONE FOR USE IN CALCULATING THE DRAG
C     COEF FOR HEAT.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      A = VK / ALOG( Z / Z0 )
      AM = A * A
      AH = A * VK / ALOG( Z / Z0H )

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CALC A BULK RICHARDSON NUMBER.  CONSTRAIN ITS VALUE IN THE
C     STABLE CASE TO A MAXIMUM OF 1.0 TO AVOID CREATING EXCHANGE
C     COEFFICIENTS THAT APPROACH ZERO.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RIB = G * Z * ( TH2V - T1V ) / ( TH2V * SFCSPD * SFCSPD )
      RIB = AMIN1( RIB, 1.0 )

      IF ( RIB .GE. 0. ) THEN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       IF THE RICHARDSON NUMBER IS .GE. ZERO, THE AIR IS STABLY
C       STRATIFIED.  CALC THE DRAG COEFFICIENTS USING A METHOD
C       DEVELOPED BY MAHRT (MONTHLY WEATHER REVIEW, 1987).  THE
C       SMALLEST ALLOWABLE CM VALUE IS 1.0E-6.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        CM = AM * SFCSPD * EXP( EXMCH * RIB )
        IF ( CM .LT. 1.0E-6) THEN
          CM = 1.0E-6
          CH = 1.0E-6
        ELSE
          CH = ( CM * AH / AM ) / PR
        ENDIF
      ELSE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       OTHERWISE, THE AIR IS UNSTABLY STRATIFIED AND THE DRAG
C       COEFFICIENTS WILL BE CALCULATED AS FOLLOWS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        CUP = CUS * AM * B1 * SQRT( -RIB * Z / Z0 )
        CTP = CUS * AH * B1 * SQRT( -RIB * Z / Z0H )
        CM = ( 1.0 - (B1 * RIB)/(1.0 + CUP) ) * AM * SFCSPD
        CH = ( 1.0 - (B2 * RIB)/(1.0 + CTP) ) * AH * SFCSPD / PR
      END IF

      RETURN
      END
