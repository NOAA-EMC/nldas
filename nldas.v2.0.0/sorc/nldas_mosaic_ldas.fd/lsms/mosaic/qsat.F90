!-------------------------------------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
! Released October 2005
!
! See SOFTWARE DISTRIBUTION POLICY for software distribution policies
!
! The LIS source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See COPYRIGHT.TXT for copyright details.
!
!-------------------------------------------------------------------------
!**** ------------------------------------------------------------------
!**** //////////////////////////////////////////////////////////////////
!**** ------------------------------------------------------------------
!****
      REAL FUNCTION QSAT(T,PR,ALHX)
!**** 
      IMPLICIT NONE

      REAL  T, PR, ALHX, C1, C2, C3, ALHE, ALHS, ALHM, C1LOG
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      DATA C1/3.797915/, C2/7.93252E-6/, C3/2.166847E-3/, &
                C1LOG/1.33445/
!****
!      QSAT = C1*EXP(ALHX*(C2-C3/T))/PR
      QSAT = C1*EXP((ALHX/ALHE)*(21.18123-C1LOG-5418./T))/PR
!****
      RETURN
      END

!**** ------------------------------------------------------------------
!**** //////////////////////////////////////////////////////////////////
!**** ------------------------------------------------------------------
!****
      REAL FUNCTION DQSAT(T,PR,ALHX)
!**** 
      IMPLICIT NONE

      REAL  T, PR, ALHX, C1, C2, C3, QS, ALHE, ALHS, ALHM
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      DATA C1/3.797915/, C2/7.93252E-6/, C3/2.166847E-3/
!****
!      QS = C1*EXP(ALHX*(C2-C3/T))/PR
      QS = C1*EXP((ALHX/ALHE)*(21.18123-ALOG(C1)-5418./T))/PR
!      DQSAT = QS*C3*ALHX/(T*T)
      DQSAT = QS * 5418. / ( T * T )
      
!****
      RETURN
      END
