      SUBROUTINE energy1 (
     I                 NCH, DTSTEP, ITYP, UM, RCIN,
     I                 ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,
     I                 QM,     RA,   SWNET,  HLWDWN, PSUR,
     I                 RDC,    U2FAC, HFTDS, DHFTDS,
     I                 QSATTC, DQSDTC, ALWRAD, BLWRAD,
     I                 EMAXRT,CSOIL,SWSRF,POTFRC,BUG,
     U                 TC, QA,
     O                 EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX
     &                 )

C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MXCHP, MemFac
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (NTYPS = 10, FRSTCH = 1, MXCHP = 3500, MemFac = 5)
      PARAMETER (NLAY = 3)
      PARAMETER (SFCLY = 1, ROOTLY = SFCLY + 1, RECHLY = ROOTLY + 1)

      PARAMETER (ZERO = 0., ONE = 1., PIE = 3.14159265)
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      PARAMETER (TF = 273.16)
      PARAMETER (STEFAN = 5.669E-8)
      PARAMETER (RGAS = .286*1003.5)
      PARAMETER (SHW = 4200., SHI = 2060.)
      PARAMETER (RHOW = 1000.)
      PARAMETER (GRAV = 9.81)
      PARAMETER (EPSILON = 18.01/28.97)
      PARAMETER (NOSNOW = 0.)

C****
      INTEGER NCH
      INTEGER ITYP(NCH)
      REAL  DTSTEP,  UM(NCH),
     &       ETURB(NCH),  DEDQA(NCH), HSTURB(NCH), DHSDTC(NCH),
     &      DHSDQA(NCH),  QM (NCH),  SWNET(NCH),
     &      HLWDWN(NCH),   PSUR(NCH), DEDTC(NCH), RCIN(NCH)
      REAL  RDC(NCH),  U2FAC(NCH),
     &      QSATTC(NCH), DQSDTC(NCH), ALWRAD(NCH), BLWRAD(NCH),
     &          TC(NCH),     QA(NCH), POTFRC(NCH), SNWFRC(NCH),
     &        EVAP(NCH), SHFLUX(NCH)
      REAL   HLWUP(NCH),    RX1(NCH), RX2(NCH), RA(NCH), GHFLUX(NCH)
      REAL EMAXRT(NCH), CSOIL(NCH),    SWSRF(NCH)
      LOGICAL*1 BUG
C****
      INTEGER ChNo
      REAL   DELTC,  DELEA
C****
      REAL VPDSTR(MXCHP),ESATTX(MXCHP),VPDSTX(MXCHP)
      REAL FTEMP(MXCHP), RC(MXCHP),  EAX(MXCHP),   TX(MXCHP),   
     &     RCX(MXCHP),  DRCDTC(MXCHP),   DUMMY(MXCHP)
      REAL FTEMPX(MXCHP),  DRCDEA(MXCHP)
      REAL  DEDEA(MXCHP),  DHSDEA(MXCHP),  EM(MXCHP), ESATTC(MXCHP),
     &     DESDTC(MXCHP),      EA(MXCHP),
     &      EPFRC, RTBS
      REAL RC1(MXCHP), RC2(MXCHP), RC3(MXCHP), RC4(MXCHP),
     &     RC1X(MXCHP),RC2X(MXCHP),RC3X(MXCHP),RC4X(MXCHP)
      REAL ADEN1,ADEN2,ADEN3,ADEN4
      REAL HFTDS(MXCHP),DHFTDS(MXCHP),TCOLD(MXCHP)


      real tstart,tfinal,tenergy,tsnow,tfluxes,tground
      real tstartx,tfinalx,t1en,t2en,t3en
      common/timing/tcatch,tenergy,tsnow,tfluxes,tground
      common/cenergy/t1en,t2en,t3en

C**** - - - - - - - -
C****
      DATA DELTC /0.01/, DELEA /0.001/
C****
C**** --------------------------------------------------------------------- 
C****
C**** Expand data as specified by ITYP into arrays of size NCH.
C****

C****
C**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)

 100  CONTINUE

C****
C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
C****          ( VPDFAC  computes the vapor pressure deficit stress term,
C****           TMPFAC  computes the temperature stress term,
C****           RSURFP  computes rc given a parallel resist. from the
C****                   surface,
C****           RCANOP  computes rc corrected for snow and interception.)
C****

      CALL VPDFAC (
     I             NCH,  ITYP,  ESATTC, EA,   
     O             VPDSTR
     &            )


      CALL TMPFAC (
     I             NCH,  ITYP, TC,
     O             FTEMP
     &            )


      DO N=1,NCH
        RC(N)=RCIN(N)/(VPDSTR(N)*FTEMP(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (
     I             NCH, UM, U2FAC, RDC, SWSRF,
     I             ESATTC, EA,
     U             RC,
     O             RX1, RX2
     &            )


      CALL RCANOP (
     I             NCH, RA, ETURB, POTFRC,
     U             RC
     &            )

C****
C**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
C****
C**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
C****

      DO ChNo = 1, NCH
        TX(ChNo) = TC(ChNo) + DELTC
        ESATTX(ChNo) = ESATTC(ChNo) + DESDTC(CHNO) * DELTC
        EAX(ChNo) = EA(ChNo) + DELEA
        ENDDO

C****
C**** temperature:
      CALL VPDFAC (NCH, ITYP, ESATTX, EA, VPDSTX)
      CALL TMPFAC (NCH, ITYP, TX, FTEMPX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMPX(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (NCH, UM, U2FAC, RDC, SWSRF,
     I             ESATTX, EA,
     I             RCX,
     O             DUMMY,DUMMY
     &            )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
C****
      DO  ChNo = 1, NCH
        DRCDTC(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELTC
        ENDDO
C****

C**** vapor pressure:
      CALL VPDFAC (NCH, ITYP, ESATTC, EAX, VPDSTX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMP(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (NCH, UM, U2FAC, RDC, SWSRF,
     I             ESATTC, EAX,
     I             RCX,
     O             DUMMY,DUMMY
     &            )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
C****
      DO ChNo = 1, NCH
         DRCDEA(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELEA
         ENDDO
C****
C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****       STEP 2: Solve the energy balance at the surface.
C****
      CALL FLUXES (
     I                NCH,   ITYP, DTSTEP, ESATTC, DESDTC,
     I              ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC, 
     I                 RC, DRCDEA, DRCDTC, 
     I              SWNET, HLWDWN, ALWRAD, BLWRAD,
     I                 EM,  CSOIL,   PSUR, EMAXRT, 
     I              HFTDS, DHFTDS,
     U                 TC,     EA,
     O             EVAP, SHFLUX,  HLWUP, GHFLUX
     &            )

c****
C****
C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****
C****
C**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE

C****
      RETURN
      END

C****
C**** [ END CHIP ]
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
      SUBROUTINE energy2 (
     I                 NCH, DTSTEP, ITYP, UM, RCIN,
     I                 ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,
     I                 QM,     RA,   SWNET,  HLWDWN, PSUR,
     I                 RDC,    U2FAC, HFTDS, DHFTDS,
     I                 QSATTC, DQSDTC, ALWRAD, BLWRAD,
     I                 EMAXRT,CSOIL,SWSRF,POTFRC,BUG,
     U                 TC, QA,
     O                 EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX
     &                 )

C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MXCHP, MemFac
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (NTYPS = 10, FRSTCH = 1, MXCHP = 3500, MemFac = 5)
      PARAMETER (NLAY = 3)
      PARAMETER (SFCLY = 1, ROOTLY = SFCLY + 1, RECHLY = ROOTLY + 1)

      PARAMETER (ZERO = 0., ONE = 1., PIE = 3.14159265)
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      PARAMETER (TF = 273.16)
      PARAMETER (STEFAN = 5.669E-8)
      PARAMETER (RGAS = .286*1003.5)
      PARAMETER (SHW = 4200., SHI = 2060.)
      PARAMETER (RHOW = 1000.)
      PARAMETER (GRAV = 9.81)
      PARAMETER (EPSILON = 18.01/28.97)
      PARAMETER (NOSNOW = 0.)

C****
      INTEGER NCH
      INTEGER ITYP(NCH)
      REAL  DTSTEP,  UM(NCH),
     &       ETURB(NCH),  DEDQA(NCH), HSTURB(NCH), DHSDTC(NCH),
     &      DHSDQA(NCH),  QM (NCH),  SWNET(NCH),
     &      HLWDWN(NCH),   PSUR(NCH), DEDTC(NCH), RCIN(NCH)
      REAL  RDC(NCH),  U2FAC(NCH),
     &      QSATTC(NCH), DQSDTC(NCH), ALWRAD(NCH), BLWRAD(NCH),
     &          TC(NCH),     QA(NCH), POTFRC(NCH), SNWFRC(NCH),
     &        EVAP(NCH), SHFLUX(NCH)
      REAL   HLWUP(NCH),    RX1(NCH), RX2(NCH), RA(NCH), GHFLUX(NCH)
      REAL EMAXRT(NCH), CSOIL(NCH),    SWSRF(NCH)
      LOGICAL*1 BUG
C****
      INTEGER ChNo
      REAL   DELTC,  DELEA
C****
      REAL VPDSTR(MXCHP),ESATTX(MXCHP),VPDSTX(MXCHP)
      REAL FTEMP(MXCHP), RC(MXCHP),  EAX(MXCHP),   TX(MXCHP),   
     &     RCX(MXCHP),  DRCDTC(MXCHP),   DUMMY(MXCHP)
      REAL FTEMPX(MXCHP),  DRCDEA(MXCHP)
      REAL  DEDEA(MXCHP),  DHSDEA(MXCHP),  EM(MXCHP), ESATTC(MXCHP),
     &     DESDTC(MXCHP),      EA(MXCHP),
     &      EPFRC, RTBS
      REAL RC1(MXCHP), RC2(MXCHP), RC3(MXCHP), RC4(MXCHP),
     &     RC1X(MXCHP),RC2X(MXCHP),RC3X(MXCHP),RC4X(MXCHP)
      REAL ADEN1,ADEN2,ADEN3,ADEN4
      REAL HFTDS(MXCHP),DHFTDS(MXCHP),TCOLD(MXCHP)


      real tstart,tfinal,tenergy,tsnow,tfluxes,tground
      real tstartx,tfinalx,t1en,t2en,t3en
      common/timing/tcatch,tenergy,tsnow,tfluxes,tground
      common/cenergy/t1en,t2en,t3en

C**** - - - - - - - -
C****
      DATA DELTC /0.01/, DELEA /0.001/
C****
C**** --------------------------------------------------------------------- 
C****
C**** Expand data as specified by ITYP into arrays of size NCH.
C****

C****
C**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)

 100  CONTINUE

C****
C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
C****          ( VPDFAC  computes the vapor pressure deficit stress term,
C****           TMPFAC  computes the temperature stress term,
C****           RSURFP  computes rc given a parallel resist. from the
C****                   surface,
C****           RCANOP  computes rc corrected for snow and interception.)
C****

      CALL VPDFAC (
     I             NCH,  ITYP,  ESATTC, EA,   
     O             VPDSTR
     &            )


      CALL TMPFAC (
     I             NCH,  ITYP, TC,
     O             FTEMP
     &            )


      DO N=1,NCH
        RC(N)=RCIN(N)/(VPDSTR(N)*FTEMP(N)+1.E-20)
        ENDDO

      CALL RSURFP2 (
     I             NCH, UM, U2FAC, RDC, SWSRF,
     I             ESATTC, EA,
     U             RC,
     O             RX1, RX2
     &            )


      CALL RCANOP (
     I             NCH, RA, ETURB, POTFRC,
     U             RC
     &            )

C****
C**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
C****
C**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
C****

      DO ChNo = 1, NCH
        TX(ChNo) = TC(ChNo) + DELTC
        ESATTX(ChNo) = ESATTC(ChNo) + DESDTC(CHNO) * DELTC
        EAX(ChNo) = EA(ChNo) + DELEA
        ENDDO

C****
C**** temperature:
      CALL VPDFAC (NCH, ITYP, ESATTX, EA, VPDSTX)
      CALL TMPFAC (NCH, ITYP, TX, FTEMPX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMPX(N)+1.E-20)
        ENDDO

      CALL RSURFP2 (NCH, UM, U2FAC, RDC, SWSRF,
     I             ESATTX, EA,
     I             RCX,
     O             DUMMY,DUMMY
     &            )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
C****
      DO  ChNo = 1, NCH
        DRCDTC(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELTC
        ENDDO
C****

C**** vapor pressure:
      CALL VPDFAC (NCH, ITYP, ESATTC, EAX, VPDSTX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMP(N)+1.E-20)
        ENDDO

      CALL RSURFP2 (NCH, UM, U2FAC, RDC, SWSRF,
     I             ESATTC, EAX,
     I             RCX,
     O             DUMMY,DUMMY
     &            )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
C****
      DO ChNo = 1, NCH
         DRCDEA(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELEA
         ENDDO
C****
C****
C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****       STEP 2: Solve the energy balance at the surface.
C****
      CALL FLUXES (
     I                NCH,   ITYP, DTSTEP, ESATTC, DESDTC,
     I              ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC, 
     I                 RC, DRCDEA, DRCDTC, 
     I              SWNET, HLWDWN, ALWRAD, BLWRAD,
     I                 EM,  CSOIL,   PSUR, EMAXRT, 
     I              HFTDS, DHFTDS,
     U                 TC,     EA,
     O             EVAP, SHFLUX,  HLWUP, GHFLUX
     &            )

c****
C****
C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****
C****
C**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE
C****

      RETURN
      END

C****
C**** [ END CHIP ]
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
      SUBROUTINE energy4 (
     I                 NCH, DTSTEP, ITYP, UM, RCIN,
     I                 ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,
     I                 QM,     RA,   SWNET,  HLWDWN, PSUR,
     I                 RDC,    U2FAC, HFTDS, DHFTDS,
     I                 QSATTC, DQSDTC, ALWRAD, BLWRAD,
     I                 EMAXRT,CSOIL,SWSRF,POTFRC,BUG,
     U                 TC, QA,
     O                 EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX
     &                 )

C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MXCHP, MemFac
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (NTYPS = 10, FRSTCH = 1, MXCHP = 3500, MemFac = 5)
      PARAMETER (NLAY = 3)
      PARAMETER (SFCLY = 1, ROOTLY = SFCLY + 1, RECHLY = ROOTLY + 1)

      PARAMETER (ZERO = 0., ONE = 1., PIE = 3.14159265)
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      PARAMETER (TF = 273.16)
      PARAMETER (STEFAN = 5.669E-8)
      PARAMETER (RGAS = .286*1003.5)
      PARAMETER (SHW = 4200., SHI = 2060.)
      PARAMETER (RHOW = 1000.)
      PARAMETER (GRAV = 9.81)
      PARAMETER (EPSILON = 18.01/28.97)
      PARAMETER (NOSNOW = 0.)

C****
      INTEGER NCH
      INTEGER ITYP(NCH)
      REAL  DTSTEP,  UM(NCH),
     &       ETURB(NCH),  DEDQA(NCH), HSTURB(NCH), DHSDTC(NCH),
     &      DHSDQA(NCH),  QM (NCH),  SWNET(NCH),
     &      HLWDWN(NCH),   PSUR(NCH), DEDTC(NCH), RCIN(NCH)
      REAL  RDC(NCH),  U2FAC(NCH),
     &      QSATTC(NCH), DQSDTC(NCH), ALWRAD(NCH), BLWRAD(NCH),
     &          TC(NCH),     QA(NCH), POTFRC(NCH), SNWFRC(NCH),
     &        EVAP(NCH), SHFLUX(NCH)
      REAL   HLWUP(NCH),    RX1(NCH), RX2(NCH), RA(NCH), GHFLUX(NCH)
      REAL EMAXRT(NCH), CSOIL(NCH),    SWSRF(NCH)
      LOGICAL*1 BUG
C****
      INTEGER ChNo
      REAL   DELTC,  DELEA
C****
      REAL VPDSTR(MXCHP),ESATTX(MXCHP),VPDSTX(MXCHP)
      REAL FTEMP(MXCHP), RC(MXCHP),  EAX(MXCHP),   TX(MXCHP),   
     &     RCX(MXCHP),  DRCDTC(MXCHP),   DUMMY(MXCHP)
      REAL FTEMPX(MXCHP),  DRCDEA(MXCHP)
      REAL  DEDEA(MXCHP),  DHSDEA(MXCHP),  EM(MXCHP), ESATTC(MXCHP),
     &     DESDTC(MXCHP),      EA(MXCHP),
     &      EPFRC, RTBS
      REAL RC1(MXCHP), RC2(MXCHP), RC3(MXCHP), RC4(MXCHP),
     &     RC1X(MXCHP),RC2X(MXCHP),RC3X(MXCHP),RC4X(MXCHP)
      REAL ADEN1,ADEN2,ADEN3,ADEN4
      REAL HFTDS(MXCHP),DHFTDS(MXCHP),TCOLD(MXCHP)


      real tstart,tfinal,tenergy,tsnow,tfluxes,tground
      real tstartx,tfinalx,t1en,t2en,t3en
      common/timing/tcatch,tenergy,tsnow,tfluxes,tground
      common/cenergy/t1en,t2en,t3en

C**** - - - - - - - -
C****
      DATA DELTC /0.01/, DELEA /0.001/
C****
C**** --------------------------------------------------------------------- 
C****
C**** Expand data as specified by ITYP into arrays of size NCH.
C****

C****
C**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)

 100  CONTINUE

C****
C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
C****          ( VPDFAC  computes the vapor pressure deficit stress term,
C****           TMPFAC  computes the temperature stress term,
C****           RSURFP  computes rc given a parallel resist. from the
C****                   surface,
C****           RCANOP  computes rc corrected for snow and interception.)
C****

      DO N=1,NCH
        RC(N)=RCIN(N)
        ENDDO

      CALL RSURFP2 (
     I             NCH, UM, U2FAC, RDC, SWSRF,
     I             ESATTC, EA,
     U             RC,
     O             RX1, RX2
     &            )

      CALL RCANOP (
     I             NCH, RA, ETURB, POTFRC,
     U             RC
     &            )

C****
C**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
C****
C**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
C****

      DO ChNo = 1, NCH
        DRCDTC(CHNO)=0.
        DRCDEA(CHNO)=0.
        ENDDO

C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****       STEP 2: Solve the energy balance at the surface.
C****
      CALL FLUXES (
     I                NCH,   ITYP, DTSTEP, ESATTC, DESDTC,
     I              ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC, 
     I                 RC, DRCDEA, DRCDTC, 
     I              SWNET, HLWDWN, ALWRAD, BLWRAD,
     I                 EM,  CSOIL,   PSUR, EMAXRT, 
     I              HFTDS, DHFTDS,
     U                 TC,     EA,
     O             EVAP, SHFLUX,  HLWUP, GHFLUX
     &            )

c****
C****
C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****
C****
C**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE
C****

      RETURN
      END

C****
C**** [ END CHIP ]
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C**** [ BEGIN FLUXES ]
C****
      SUBROUTINE FLUXES (
     I                      NCH,   ITYP, DTSTEP, ESATTC, DESDTC,
     I                    ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,
     I                       RC, DRCDEA, DRCDTC,
     I                    SWNET, HLWDWN, ALWRAD, BLWRAD,
     I                       EM,  CSOIL,   PSUR, EMAXRT,
     I                    HFTDS, DHFTDS,
     U                       TC,     EA,
     O                     EVAP, SHFLUX,  HLWUP, GHFLUX
     &                  )
C****
C**** This subroutine computes the fluxes of latent and sensible heat
C**** from the surface through an energy balance calculation.
C****
      IMPLICIT NONE
C****
C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MemFac
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (NTYPS = 10, FRSTCH = 1, MemFac = 5)
      PARAMETER (NLAY = 3)
      PARAMETER (SFCLY = 1, ROOTLY = SFCLY + 1, RECHLY = ROOTLY + 1)

      PARAMETER (ZERO = 0., ONE = 1., PIE = 3.14159265)
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      PARAMETER (TF = 273.16)
      PARAMETER (STEFAN = 5.669E-8)
      PARAMETER (RGAS = .286*1003.5)
      PARAMETER (SHW = 4200., SHI = 2060.)
      PARAMETER (RHOW = 1000.)
      PARAMETER (GRAV = 9.81)
      PARAMETER (EPSILON = 18.01/28.97)
      PARAMETER (NOSNOW = 0.)

C****
      INTEGER NCH
      INTEGER ITYP(NCH), ChNo

      REAL       DTSTEP, ESATTC(NCH), DESDTC(NCH),
     &       ETURB(NCH),  DEDEA(NCH),  DEDTC(NCH),
     &      HSTURB(NCH), DHSDEA(NCH), DHSDTC(NCH),
     &          RC(NCH), DRCDEA(NCH), DRCDTC(NCH),
     &       SWNET(NCH), HLWDWN(NCH), ALWRAD(NCH), BLWRAD(NCH),
     &          EM(NCH),  CSOIL(NCH),   PSUR(NCH),
     &      EMAXRT(NCH),  HFTDS(NCH), DHFTDS(NCH),
     &          TC(NCH),     EA(NCH),
     &        EVAP(NCH), SHFLUX(NCH), HLWUP(NCH), GHFLUX(NCH)

      REAL     HLWTC,  CDEEPS,      Q0,  RHOAIR,   CONST,  DHLWTC,
     &        EPLANT,     A11,     A12,     A21,     A22,      F0,
     &           DEA,     DTC,
     &         EANEW,  ESATNW,  EHARMN, DETERM, DENOM, EDIF

      real tstart,tfinal,tcatch,tenergy,tsnow,tfluxes,tground
      common/timing/tcatch,tenergy,tsnow,tfluxes,tground

      LOGICAL*1 DEBUG, CHOKE
      DATA DEBUG /.FALSE./

C****
C**** -------------------------------------------------------------------

      DO 200 ChNo = 1, NCH
C****
      HLWTC = ALWRAD(CHNO) + BLWRAD(CHNO) * TC(CHNO)
      RHOAIR = PSUR(ChNo) * 100. / (RGAS * TC(ChNo))
      CONST = RHOAIR * EPSILON / PSUR(ChNo)
      DHLWTC = BLWRAD(CHNO)
C****
C**** Compute matrix elements A11, A22, AND Q0 (energy balance equation).
C****
      A11 = CSOIL(ChNo)/DTSTEP +
     &        DHLWTC +
     &        DHSDTC(ChNo) +
     &        ALHE*DEDTC(ChNo) +
     &        DHFTDS(CHNO)
      A12 = DHSDEA(ChNo) + ALHE * DEDEA(ChNo)
      Q0 =  SWNET(ChNo) +
     &        HLWDWN(ChNo) -
     &        HLWTC -
     &        HSTURB(ChNo) -
     &        ALHE * ETURB(ChNo) -
     &        HFTDS(CHNO)
C****
C**** Compute matrix elements A21, A22, and F0 (canopy water budget  
C**** equation) and solve for fluxes.  Three cases are considered:
C****
C**** 1. Standard case: RC>0.
C**** 2. RC = 0.  Can only occur if CIR is full or ETURB is negative.
C****
      CHOKE = .TRUE.

      IF( RC(CHNO) .GT. 0.) THEN 
          EPLANT = CONST * (ESATTC(ChNo) - EA(ChNo)) / RC(ChNo)
          IF(EPLANT*ETURB(ChNo).GT.0.) THEN
              EHARMN = 2.*EPLANT*ETURB(CHNO) / (EPLANT + ETURB(ChNo))
            ELSE
              EHARMN=0.
            ENDIF
C****
C****            Some limitations to A21 and A22 are applied:
C****            we assume that the increase in plant evaporation
C****            due to an increase in either TC or EA balances 
C****            or outweighs any decrease due to RC changes.
C****

          A21 =  -DEDTC(ChNo)*RC(ChNo) +
     &      amax1(0., CONST*DESDTC(ChNo) - EHARMN*DRCDTC(ChNo) )
          A22 = -( RC(ChNo)*DEDEA(ChNo) +
     &               amax1( 0., CONST + EHARMN*DRCDEA(ChNo) )   )

          F0 = RC(ChNo) * (ETURB(ChNo) - EPLANT)
          DETERM = AMIN1( A12*A21/(A11*A22) - 1., -0.1 )
          DEA = ( Q0*A21 - A11*F0 ) / ( DETERM * A11*A22 )
          DTC = ( Q0 - A12*DEA ) / A11
          EVAP(ChNo) = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
          SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA +
     &                                            DHSDTC(ChNo)*DTC
          DENOM = DETERM * A11*A22
        ELSE
          CHOKE = .FALSE.
          A21 = -DESDTC(ChNo)
          A22 = 1.
          F0 = ESATTC(ChNo) - EA(ChNo)
          DEA = ( Q0*A21 - A11*F0 ) / ( A12*A21 - A11*A22 )
          DTC = ( Q0 - A12*DEA ) / A11
          EVAP(ChNo) = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
          SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA +
     &                                            DHSDTC(ChNo)*DTC
          DENOM = A12 * A21 - A11*A22
        ENDIF

C**** - - - - - - - - - - - - - - - - - - - - - - -
C**** Adjustments

C**** 1. Adjust deltas and fluxes if all available water evaporates
C****    during time step:
C****
      IF( EVAP(CHNO) .GT. EMAXRT(CHNO) ) THEN
        CHOKE = .FALSE.
        DEA = EM(CHNO) - EA(CHNO)
        DTC = 
     &   (Q0 + ALHE*(ETURB(ChNo)-EMAXRT(CHNO)) - DHSDEA(CHNO)*DEA) 
     &          /  ( A11 - ALHE*DEDTC(ChNo) )
        EVAP(CHNO) = EMAXRT(CHNO)
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA +
     &                                            DHSDTC(ChNo)*DTC
        ENDIF
C****
C**** Adjust DEA and DTC if solutions were pathological:
C****
      ESATNW = ESATTC(CHNO)+DESDTC(CHNO)*DTC
      EANEW = EA(CHNO) + DEA



C**** 2. Pathological cases. 

C**** Case 1: EVAP is positive in presence of negative gradient.
C**** Case 2: EVAP and ETURB have opposite sign, implying that 
C**** "virtual effects" derivatives are meaningless and thus that we
C**** don't know the proper tendency terms.
C**** In both cases, assume zero evaporation for the time step.

C      IF( ( EVAP(CHNO) .LT. 0. .AND. EM(CHNO).LT.ESATNW )
C     &    .OR. ( EVAP(CHNO)*ETURB(CHNO) .LT. 0. ) )  THEN 
      IF( EVAP(CHNO) .LT. 0. .AND. EM(CHNO).LT.ESATNW )  THEN 
        CHOKE = .FALSE.
        DEA = EM(CHNO) - EA(CHNO)
        DTC = ( Q0 + ALHE*ETURB(ChNo) - DHSDEA(CHNO)*DEA ) / 
     &            ( A11 - ALHE*DEDTC(ChNo) )
        EVAP(CHNO) = 0.
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA +
     &                                            DHSDTC(ChNo)*DTC
        ENDIF



C**** 3. Excessive dea change: apply "choke".
c -- (03.09.98) : changed to correct the conservation.

      IF( CHOKE .AND. ABS(DEA) .GT. 0.5*EA(CHNO) ) THEN
        DEA = SIGN(.5*EA(CHNO),DEA)
        DTC = ( Q0 - A12*DEA ) / A11
        EVAP(ChNo) = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA +
     &                                            DHSDTC(ChNo)*DTC

        IF(EVAP(CHNO) .GT. EMAXRT(CHNO)) THEN
          EDIF=EVAP(CHNO)-EMAXRT(CHNO)
          EVAP(CHNO) = EMAXRT(CHNO)
          SHFLUX(ChNo) = SHFLUX(CHNO)+EDIF*ALHE
          ENDIF

        ENDIF

C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      TC(ChNo) = TC(ChNo) + DTC
      EA(ChNo) = EA(ChNo) + DEA
      HLWUP(CHNO) = HLWTC + DHLWTC*DTC

      HLWTC = ALWRAD(CHNO) + BLWRAD(CHNO) * TC(CHNO)
c warning: this ghflux is the real ground heat flux, and does not include
c the temperature variation
      GHFLUX(CHNO)=HFTDS(CHNO)+DHFTDS(CHNO)*DTC

C**** Make sure EA remains positive

      EA(CHNO) = AMAX1(EA(CHNO), 0.0)

  200 CONTINUE

      RETURN
      END
C****
C**** [ END FLUXES ]
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
C**** [ BEGIN VPDFAC ]
C****
      SUBROUTINE VPDFAC (
     I                   NCH, ITYP, ESATTC, EA,
     O                   VPDSTR
     &                  )
C****
C**** This subroutine computes the vapor pressure deficit stress.
C****
      IMPLICIT NONE
C****
C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MemFac
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (NTYPS = 10, FRSTCH = 1, MemFac = 5)
      PARAMETER (NLAY = 3)
      PARAMETER (SFCLY = 1, ROOTLY = SFCLY + 1, RECHLY = ROOTLY + 1)

      PARAMETER (ZERO = 0., ONE = 1., PIE = 3.14159265)
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      PARAMETER (TF = 273.16)
      PARAMETER (STEFAN = 5.669E-8)
      PARAMETER (RGAS = .286*1003.5)
      PARAMETER (SHW = 4200., SHI = 2060.)
      PARAMETER (RHOW = 1000.)
      PARAMETER (GRAV = 9.81)
      PARAMETER (EPSILON = 18.01/28.97)
      PARAMETER (NOSNOW = 0.)

C****
      INTEGER NCH
      INTEGER ITYP(NCH), ChNo
      REAL    ESATTC(NCH), EA(NCH),  VPDSTR(NCH)
      REAL  VGDFAC(NTYPS)
C****
      DATA VGDFAC /   .0273,    .0357,    .0310,    .0238,
     5                .0275,    .0275,       0.,       0.,
     9                   0.,       0. /
C****
C**** -----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
C****
c      VPDSTR(ChNo) = 1. - (ESATTC(ChNo)-EA(ChNo)) * VGDFAC(ITYP(ChNo))
c      VPDSTR (ChNo) = AMIN1( 1., AMAX1( VPDSTR(ChNo), 1.E-10 ) )
      VPDSTR(CHNO) = 1.
C****
 100  CONTINUE
C****
      RETURN
      END
C****
C**** [ END VPDFAC ]
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
C**** [ BEGIN TMPFAC ]
C****
      SUBROUTINE TMPFAC (
     I                   NCH,  ITYP, TC,
     O                   FTEMP
     &                  )
C****
C**** Compute temperature stress factor.
C****
      IMPLICIT NONE
C****
C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MemFac
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (NTYPS = 10, FRSTCH = 1, MemFac = 5)
      PARAMETER (NLAY = 3)
      PARAMETER (SFCLY = 1, ROOTLY = SFCLY + 1, RECHLY = ROOTLY + 1)

      PARAMETER (ZERO = 0., ONE = 1., PIE = 3.14159265)
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      PARAMETER (TF = 273.16)
      PARAMETER (STEFAN = 5.669E-8)
      PARAMETER (RGAS = .286*1003.5)
      PARAMETER (SHW = 4200., SHI = 2060.)
      PARAMETER (RHOW = 1000.)
      PARAMETER (GRAV = 9.81)
      PARAMETER (EPSILON = 18.01/28.97)
      PARAMETER (NOSNOW = 0.)

C****
      INTEGER NCH
      INTEGER ITYP(NCH), ChNo, TypPtr
      REAL TC(NCH), FTEMP(NCH)
      REAL    VGTLL(MemFac*NTYPS),    VGTU(MemFac*NTYPS),
     &       VGTCF1(MemFac*NTYPS),  VGTCF2(MemFac*NTYPS),
     &       VGTCF3(MemFac*NTYPS)
C****
      DATA VGTLL /MemFac*273., MemFac*273., MemFac*268., MemFac*283.,
     5            MemFac*283., MemFac*273., MemFac*  0., MemFac*  0., 
     9            MemFac*  0., MemFac*  0. /
      DATA VGTU /MemFac*318., MemFac*318., MemFac*313., MemFac*328.,
     5           MemFac*323., MemFac*323., MemFac*  0., MemFac*  0.,
     9           MemFac*  0., MemFac*  0. /
      DATA VGTCF1 / MemFac*-1.43549E-06,  MemFac*-6.83584E-07,
     3              MemFac* 1.67699E-07,  MemFac*-1.43465E-06,
     5              MemFac*-2.76097E-06,  MemFac*-1.58094E-07,
     7                       MemFac* 0.,           MemFac* 0.,
     9                       MemFac* 0.,           MemFac* 0. /
      DATA VGTCF2 / MemFac* 7.95859E-04,  MemFac* 3.72064E-04,
     3              MemFac*-7.65944E-05,  MemFac* 8.24060E-04,
     5              MemFac* 1.57617E-03,  MemFac* 8.44847E-05, 
     7                       MemFac* 0.,           MemFac* 0.,
     9                       MemFac* 0. ,          MemFac* 0./
      DATA VGTCF3 / MemFac*-1.11575E-01,  MemFac*-5.21533E-02,
     3              MemFac* 6.14960E-03,  MemFac*-1.19602E-01,
     5              MemFac*-2.26109E-01,  MemFac*-1.27272E-02,
     7                       MemFac* 0.,           MemFac* 0.,
     9                       MemFac* 0.,           MemFac* 0. /
C****
C**** ----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
C****
      TypPtr = MOD(ChNo,MemFac) + (ITYP(ChNo)-1)*MemFac + 1
      FTEMP(ChNo) = (TC(ChNo) - VGTLL(TypPtr)) *
     &              (TC(ChNo) - VGTU(TypPtr)) *
     &                    ( VGTCF1(TypPtr)*TC(ChNo)*TC(ChNo) +
     &                      VGTCF2(TypPtr)*TC(ChNo) +
     &                      VGTCF3(TypPtr) )
      IF ( TC(ChNo) .LE. VGTLL(TypPtr) .OR. TC(ChNo) .GE. VGTU(TypPtr) )
     &      FTEMP (ChNo) = 1.E-10
      FTEMP(CHNO) = AMIN1( 1., AMAX1( FTEMP(ChNo), 1.E-10 ) )
C****
 100  CONTINUE
C****
      RETURN
      END
C****
C**** [ END TMPFAC ]
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
C**** [ BEGIN RSURFP ]
C****
      SUBROUTINE RSURFP1 (
     I                   NCH, UM, U2FAC, RDC, WET, ESATTC, EA,
     U                   RC,
     O                   RX1, RX2
     &                  )
C****
      IMPLICIT NONE
      INTEGER NCH
      INTEGER ChNo
      REAL     UM(NCH),   U2FAC(NCH),    RDC(NCH), 
     &        WET(NCH),  ESATTC(NCH),     EA(NCH),
     &         RC(NCH),     RX1(NCH),    RX2(NCH)
      REAL  U2, RSURF, HESAT
C****
C**** -----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
C****
      U2 = UM(ChNo) * U2FAC(ChNo)
C      RSURF = RDC(ChNo) / U2 + 26. + 6. / (1.E-10 + WET(ChNo))**2
      RSURF = RDC(ChNo) / U2

      RX1(CHNO)=RC(CHNO)
      RX2(CHNO)=RSURF

      RC(ChNo) = RC(CHNO) * RSURF / ( RC(ChNo) + RSURF )
C****
 100  CONTINUE
C****
      RETURN
      END
C****
C**** [ END RSURFP ]
C****
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
C**** [ BEGIN RSURFP ]
C****
      SUBROUTINE RSURFP2 (
     I                   NCH, UM, U2FAC, RDC, WET, ESATTC, EA,
     U                   RC,
     O                   RX1, RX2
     &                  )
C****
      IMPLICIT NONE
      INTEGER NCH
      INTEGER ChNo
      REAL     UM(NCH),   U2FAC(NCH),    RDC(NCH), 
     &        WET(NCH),  ESATTC(NCH),     EA(NCH),
     &         RC(NCH),     RX1(NCH),    RX2(NCH)
      REAL  U2, RSURF, HESAT
C****
C**** -----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
C****
      U2 = UM(ChNo) * U2FAC(ChNo)
c      RSURF = RDC(ChNo) / U2 + 30. / (1.E-20 + WET(ChNo))
c      RSURF = RDC(ChNo) / U2 + 26. + 6. / (1.E-20 + WET(ChNo))**2
c the second change was previously done in MOSAIC
      RSURF = RDC(ChNo) / U2 + 26. + 6. / (1.E-10 + WET(ChNo))**2

C**** Account for subsaturated humidity at soil surface:
C****
c      HESAT = ESATTC(CHNO) * MIN( 1., WET(CHNO)*2. )
c      IF( EA(CHNO) .LT. HESAT ) THEN
c          RSURF=RSURF*( 1. + (ESATTC(CHNO)-HESAT)/(HESAT-EA(CHNO)) )
c        ELSE
c          RSURF=1.E10
c        ENDIF


      RX1(CHNO)=RC(CHNO)
      RX2(CHNO)=RSURF

      RC(ChNo) = RC(CHNO) * RSURF / ( RC(ChNo) + RSURF )
C****
 100  CONTINUE
C****
      RETURN
      END
C****
C**** [ END RSURFP ]
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
C**** [ BEGIN RCANOP ]
C****
      SUBROUTINE RCANOP (
     I                   NCH, RA, ETURB, POTFRC,
     U                   RC
     &                  )
C****
C**** The effective latent heat resistance RC depends on the quantity 
C**** of interception reservoir water and the snow cover.  POTFRC
C**** is the fraction of the tile from which potential evaporation
C**** occurs.
C****
      IMPLICIT NONE
      INTEGER NCH,N

      REAL RA(NCH), ETURB(NCH), RC(NCH), POTFRC(NCH)
      REAL ETCRIT,RAMPFC

C**** (Note: ETCRIT arbitrarily set to ~-5 W/m2, or -2.e-6 mm/sec.)
      DATA ETCRIT/ -2.E-6 /
C****
C**** -----------------------------------------------------------------

      DO N = 1, NCH

        RC(N)=RC(N)*(1.-POTFRC(N))/
     &                ( 1.+POTFRC(N)*RC(N)/RA(N) )

C**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C****   Assume RC=0 for condensation (dew).
C****   RAMPFC is used to ensure continuity in RC.

        RAMPFC=ETURB(N)/ETCRIT
        IF ( RAMPFC .GE. 0. ) RC(N) = RC(N)*(1.-RAMPFC)
        IF ( RAMPFC .GT. 1. ) RC(N) = 0.
C****
        ENDDO

      RETURN
      END
C****
C**** [ END RCANOP ]
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------

