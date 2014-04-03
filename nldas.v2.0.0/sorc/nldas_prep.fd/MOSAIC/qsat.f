C**** ------------------------------------------------------------------
C**** //////////////////////////////////////////////////////////////////
C**** ------------------------------------------------------------------
C****
      REAL FUNCTION QSAT(T,PR,ALHX)
C**** 
      IMPLICIT NONE

      REAL  T, PR, ALHX, C1, C2, C3, ALHE, ALHS, ALHM, C1LOG
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      DATA C1/3.797915/, C2/7.93252E-6/, C3/2.166847E-3/,
     &           C1LOG/1.33445/
C****
c      QSAT = C1*EXP(ALHX*(C2-C3/T))/PR
      QSAT = C1*EXP((ALHX/ALHE)*(21.18123-C1LOG-5418./T))/PR
C****
      RETURN
      END

C**** ------------------------------------------------------------------
C**** //////////////////////////////////////////////////////////////////
C**** ------------------------------------------------------------------
C****
      REAL FUNCTION DQSAT(T,PR,ALHX)
C**** 
      IMPLICIT NONE

      REAL  T, PR, ALHX, C1, C2, C3, QS, ALHE, ALHS, ALHM
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      DATA C1/3.797915/, C2/7.93252E-6/, C3/2.166847E-3/
C****
C      QS = C1*EXP(ALHX*(C2-C3/T))/PR
      QS = C1*EXP((ALHX/ALHE)*(21.18123-ALOG(C1)-5418./T))/PR
C      DQSAT = QS*C3*ALHX/(T*T)
      DQSAT = QS * 5418. / ( T * T )
      
C****
      RETURN
      END
