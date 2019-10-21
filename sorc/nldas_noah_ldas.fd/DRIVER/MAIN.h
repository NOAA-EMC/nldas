      INTEGER NX
      INTEGER NY
      INTEGER NSOLD
      INTEGER NMONTH
      INTEGER NHOUR
      PARAMETER (NX = 464)
      PARAMETER (NY = 224)
      PARAMETER (NSOLD = 4)
      PARAMETER (NMONTH = 12)
      PARAMETER (NHOUR  = 25)

      REAL    TFREEZ
      PARAMETER (TFREEZ = 273.15)

      INTEGER NSOIL(NX,NY)
      INTEGER SOILTYP(NX,NY)
      INTEGER VEGTYP(NX,NY)
      INTEGER SLOPETYP(NX,NY)
      INTEGER LAND_SEA(NX,NY)

      LOGICAL*1 LDASMASK(NX,NY)

      REAL    MAXSNOWALB(NX,NY)
      REAL    SOILDEPTH(NSOLD,NX,NY)
      REAL    STC(NSOLD,NX,NY)
      REAL    SMC(NSOLD,NX,NY)
      REAL    SH2O(NSOLD,NX,NY)
      REAL    SHDFAC(NMONTH,NX,NY)
      REAL    SHDMAX(NX,NY)
      REAL    SHDMIN(NX,NY)
      REAL    TBOT(NX,NY)
      REAL    T1(NX,NY)
      REAL    CMC(NX,NY)
      REAL    SNOWH(NX,NY)
      REAL    SNEQV(NX,NY)
      REAL    DUMMY(NX,NY)

      REAL    EC(NX,NY)
      REAL    EDIR(NX,NY)
      REAL    ET(NSOLD,NX,NY)
      REAL    ETT(NX,NY)
      REAL    ESNOW(NX,NY)
      REAL    DRIP(NX,NY)
      REAL    DEW(NX,NY)

      REAL    FLX1(NX,NY)
      REAL    FLX2(NX,NY)
      REAL    FLX3(NX,NY)

      REAL    SMCWLT(NX,NY)
      REAL    SMCDRY(NX,NY)
      REAL    SMCREF(NX,NY)
      REAL    SMCMAX(NX,NY)

      CHARACTER*100 CNTRFL, DIREC
      CHARACTER*150 FILENAME
      CHARACTER*2   HOUR_CH
      CHARACTER*8   YESTERDAY, TODAY
      INTEGER       LENGTH
      INTEGER       LUGB
      INTEGER       JULDAY
      INTEGER       SUB_DT

      REAL    W1
      REAL    W2
      REAL    W_PRECIP(4)
      DATA    W_PRECIP  /0.25,0.25,0.25,0.25/

      REAL    W_now(4)
      REAL    W_minus(4)
      REAL    W_plus(4)
      DATA    W_now   /0.625,0.875,0.875,0.625/
      DATA    W_minus /0.375,0.125,0.000,0.000/
      DATA    W_plus  /0.000,0.000,0.125,0.375/


      INTEGER M1
      INTEGER M2
      INTEGER ICE(NX,NY)
      INTEGER NT
      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER NROOT(NX,NY)

      REAL    DT(NX,NY)
      REAL    ZLVL(NX,NY)
      REAL    SNODEP
      REAL    SOLNET(NX,NY)

      REAL    ICE1
      REAL    DT1
      REAL    ZLVL1

      REAL    SFCSPD(NX,NY)
      REAL    PTU(NX,NY)
      REAL    CH(NX,NY)
      REAL    CM(NX,NY)
      REAL    ALBEDO(NMONTH,NX,NY)
      REAL    ALBEDO_D(NX,NY)
      REAL    SHDFAC_D(NX,NY)

      INTEGER   NLDAS
      PARAMETER (NLDAS = NX*NY)
      REAL      F(NX,NY)
      REAL      TAIR(NX,NY,NHOUR)
      REAL      SPFH(NX,NY,NHOUR)
      REAL      PSFC(NX,NY,NHOUR)
      REAL      UWIND(NX,NY,NHOUR)
      REAL      VWIND(NX,NY,NHOUR)
      REAL      EDASSWRD(NX,NY)
      REAL      LWDN(NX,NY,NHOUR)
      REAL      EDASPREC(NX,NY)
      REAL      CONVPREC(NX,NY)
      REAL      CAPE(NX,NY)
      REAL      SOLDN(NX,NY,NHOUR)
      REAL      SKIN(NX,NY)
      REAL      PRCP(NX,NY,NHOUR)
      REAL      PRCP_MASK(NX,NY)
      REAL      SNOTEL_MASK(NX,NY)
      REAL      SCALE_PREC(3,NX,NY)
      REAL      PAR(NX,NY)
      LOGICAL*1 LB(NX,NY)

      REAL      DT_TAIR(NX,NY)
      REAL      DT_SPFH(NX,NY)
      REAL      DT_PSFC(NX,NY)
      REAL      DT_UWIND(NX,NY)
      REAL      DT_VWIND(NX,NY)
      REAL      DT_LWDN(NX,NY)
      REAL      DT_SOLDN(NX,NY)
      REAL      DT_SOLDNB(NX,NY)
      REAL      DT_PRCP(NX,NY)
      REAL      DT_SPFH_SAT(NX,NY)
      REAL      DQSDT2(NX,NY)
      REAL      DQSDT
      REAL      TH2(NX,NY)
      REAL      RC(NX,NY)
      REAL      PC(NX,NY)
      REAL      XLAI(NX,NY)
      REAL      RSMIN(NX,NY)
      REAL      RGL(NX,NY)
      REAL      HS(NX,NY)
      REAL      RCS(NX,NY)
      REAL      RCT(NX,NY)
      REAL      RCQ(NX,NY)
      REAL      RCSOIL(NX,NY)
      REAL      FX(NX,NY)
      REAL      BETA(NX,NY)

      REAL      ETP(NX,NY)
      REAL      ETA(NX,NY)
      REAL      H(NX,NY)
      REAL      S(NX,NY)
      REAL      RUNOFF1(NX,NY)
      REAL      RUNOFF2(NX,NY)
      REAL      RUNOFF3(NX,NY)
      REAL      Q1(NX,NY)
      REAL      SNOMLT(NX,NY)
      REAL      SNCOVR(NX,NY)
      REAL      SOILM(NX,NY)
      REAL      SOILW(NX,NY)
      REAL      ALBSNO(NX,NY)
      REAL      SOILM_OUT(NSOLD,NX,NY)
      REAL      SOILL_OUT(NSOLD,NX,NY)

      REAL      HOUR_DEW(NX,NY)
      REAL      HOUR_NSWRS(NX,NY)
      REAL      HOUR_LWDN_NET(NX,NY)
      REAL      HOUR_ETA(NX,NY)
      REAL      HOUR_H(NX,NY)
      REAL      HOUR_S(NX,NY)
      REAL      DelColdCont(NX,NY)
      REAL      HOUR_DSWRF(NX,NY)
      REAL      HOUR_LWDN(NX,NY)
      REAL      HOUR_ASNOW(NX,NY)
      REAL      HOUR_ARAIN(NX,NY)
      REAL      HOUR_RUNOFF1(NX,NY)
      REAL      HOUR_RUNOFF2(NX,NY)
      REAL      HOUR_SNOMLT(NX,NY)
      REAL      HOUR_CMC(NX,NY)
      REAL      HOUR_EVP(NX,NY)
      REAL      HOUR_ETP(NX,NY)
      REAL      HOUR_EDIR(NX,NY)
      REAL      HOUR_EC(NX,NY)
      REAL      HOUR_ETT(NX,NY)
      REAL      HOUR_ESNOW(NX,NY)
      REAL      HOUR_ACOND(NX,NY)
      REAL      HOUR_RC(NX,NY)
      REAL      HOUR_RCS(NX,NY)
      REAL      HOUR_RCT(NX,NY)
      REAL      HOUR_RCQ(NX,NY)
      REAL      HOUR_RCSOIL(NX,NY)
      REAL      HOUR_FX(NX,NY)
      REAL      HOUR_RSMIN(NX,NY)
      REAL      HOUR_LAI(NX,NY)
      REAL      HOUR_GVF(NX,NY)
      REAL      SOILROOT(NX,NY)
      REAL      SOIL1M(NX,NY)
      REAL      MSTAV(NX,NY)
      REAL      HOUR_ALB(NX,NY)

      REAL      LAT(NX,NY)
      REAL      LON(NX,NY)
