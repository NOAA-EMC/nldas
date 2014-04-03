      SUBROUTINE READ_LSC(NX,NY,NSOLD,NMONTH,NSOIL,SOILTYP,
     &     LAND_SEA,VEGTYP,SLOPETYP,SOILDEPTH,ALBEDO,SHDFAC,
     &     TBOT,MAXSNOWALB,PRCP_MASK)

      IMPLICIT NONE
      
C     READS SPATIAL INPUT PARAMETERS FOR NOAH - LSM
      
      INTEGER NX
      INTEGER NY
      INTEGER NSOLD
      INTEGER NMONTH
      
      INTEGER LAND_SEA(NX,NY)
      INTEGER NSOIL(NX,NY)
      INTEGER SOILTYP(NX,NY)
      INTEGER VEGTYP(NX,NY)
      INTEGER SLOPETYP(NX,NY)

      REAL    SOILDEPTH(NSOLD,NX,NY)
      REAL    ALBEDO(NMONTH,NX,NY)
      REAL    SHDFAC(NMONTH,NX,NY)
      REAL    TBOT(NX,NY)
      REAL    MAXSNOWALB(NX,NY)
      REAL    PRCP_MASK(NX,NY)

      INTEGER NREAD, I, J, K

      READ(50) ((LAND_SEA(J,K),J=1,NX),K=1,NY)
      READ(51) ((NSOIL(J,K),J=1,NX),K=1,NY)
      READ(53) ((SOILTYP(J,K),J=1,NX),K=1,NY)
      READ(54) ((VEGTYP(J,K),J=1,NX),K=1,NY)
      READ(55) ((SLOPETYP(J,K),J=1,NX),K=1,NY)
      READ(56) (((SOILDEPTH(I,J,K),I=1,NSOLD),J=1,NX),K=1,NY)
      READ(57) ((TBOT(J,K),J=1,NX),K=1,NY)
      READ(58) (((ALBEDO(I,J,K),I=1,NMONTH),J=1,NX),K=1,NY)
      READ(59) (((SHDFAC(I,J,K),I=1,NMONTH),J=1,NX),K=1,NY)
      READ(60) ((MAXSNOWALB(J,K),J=1,NX),K=1,NY)
      READ(61) ((PRCP_MASK(J,K),J=1,NX),K=1,NY)

      END
