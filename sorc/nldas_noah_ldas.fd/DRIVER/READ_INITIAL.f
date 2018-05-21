      SUBROUTINE READ_INITIAL(NX,NY,NSOLD,T1,STC,SMC,SH2O,
     &     CMC,SNOWH,SNEQV,CH,CM,LSTSNW1)

      IMPLICIT NONE

C     READS INITIAL CONDITIONS

      INTEGER NX
      INTEGER NY
      INTEGER NSOLD
      integer I,J,K

      REAL    T1(NX,NY)
      REAL    STC(NSOLD,NX,NY)
      REAL    SMC(NSOLD,NX,NY)
      REAL    SH2O(NSOLD,NX,NY)
      REAL    CMC(NX,NY)
      REAL    SNOWH(NX,NY)
      REAL    SNEQV(NX,NY)
      REAL    CH(NX,NY)
      REAL    CM(NX,NY)
      REAL    WA(NX,NY)
      REAL    WT(NX,NY)
      REAL    ZWT(NX,NY)
      REAL    LSTSNW1(NX,NY)

      READ(69) ((LSTSNW1(J,K),J=1,NX),K=1,NY)
      READ(70) ((T1(J,K),J=1,NX),K=1,NY)
      READ(71) (((STC(I,J,K),I=1,NSOLD),J=1,NX),K=1,NY)
      READ(72) (((SMC(I,J,K),I=1,NSOLD),J=1,NX),K=1,NY)
      READ(73) (((SH2O(I,J,K),I=1,NSOLD),J=1,NX),K=1,NY)
      READ(74) ((CMC(J,K),J=1,NX),K=1,NY)
      READ(75) ((SNOWH(J,K),J=1,NX),K=1,NY)
      READ(76) ((SNEQV(J,K),J=1,NX),K=1,NY)
      READ(77) ((CH(J,K),J=1,NX),K=1,NY)
      READ(78) ((CM(J,K),J=1,NX),K=1,NY)
      END
