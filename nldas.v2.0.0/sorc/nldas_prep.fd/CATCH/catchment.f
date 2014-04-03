      SUBROUTINE CATCHMENT (
     I               NCH, DTSTEP, SFRAC, ITYP,TRAINC,TRAINL, TSNOW,  UM,
     I               ETURB1, DEDQA1, DEDTC1, HSTURB1,DHSDQA1, DHSDTC1,
     I               ETURB2, DEDQA2, DEDTC2, HSTURB2,DHSDQA2, DHSDTC2,
     I               ETURB4, DEDQA4, DEDTC4, HSTURB4,DHSDQA4, DHSDTC4,
     I               ETURBS, DEDQAS, DEDTCS, HSTURBS,DHSDQAS, DHSDTCS,
     I               TM, QM, CD1, CD2, CD4, CDS, SUNANG, PARDIR, PARDIF,
     I               SWNETF,SWNETS,  HLWDWN, PSUR,  ZLAI,   GREEN,  Z2,
     I               SQSCAT, RSOIL1, RSOIL2,   RDC,    U2FAC,
     I               QSAT1, DQS1, ALW1, BLW1,
     I               QSAT2, DQS2, ALW2, BLW2,
     I               QSAT4, DQS4, ALW4, BLW4,
     I               QSATS, DQSS, ALWS, BLWS,
     I               BF1, BF2, BF3,VGWMAX,
     I               CDCR1,CDCR2, psis, bee, poros, wpwet, cond, gnu,
     I               ARS1,ARS2,ARS3,ARA1,ARA2,ARA3,ARA4,
     I               ARW1,ARW2,ARW3,ARW4,
     I               tsa1,tsa2,tsb1,tsb2,BUG,assim,
     U               TC1, TC2, TC4, QA1, QA2, QA4, CAPAC,
     U               CATDEF, RZEXC, srfexc, GHTCNT, TSURF,
     U               WESNN, HTSNNN, SNDZN,
     O               EVAP, SHFLUX, RUNOFF,
     O               EINT,   ESOI,   EVEG,   ESNO,
     O               BFLOW,RUNSRF,SMELT, 
     O               HLWUP,HLATN,QINFIL,AR1, AR2, RZEQ, 
     O               GHFLUX, TPSN1, ASNOW0, TP1, TP2,
     O               TP3, TP4, TP5, TP6, 
     o               a11, a12, a21, a22, a23, a32, a33, q11, q22, q33,
     o               es,ev,et,qin,srflw,rzflw,bflw,adj,cor,exc1,exc2
     &                     )

      IMPLICIT NONE
      INTEGER NCH,MXCHP,K,LAYER

C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MEMFAC
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY
      INTEGER N

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (NTYPS = 10, FRSTCH = 1, MXCHP = 3500, MEMFAC = 5)
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
      INTEGER ITYP(NCH)
      REAL  DTSTEP,SFRAC, TRAINL(NCH), TRAINC(NCH), TSNOW(NCH), UM(NCH),
     &     TM (NCH), CD1(NCH),CD2(NCH),CD4(NCH), CDS(NCH), SUNANG(NCH),
     &       QM(NCH), PARDIR(NCH), PARDIF(NCH), SWNETF(NCH),SWNETS(NCH),
     &      HLWDWN(NCH),   PSUR(NCH),   ZLAI(NCH),  GREEN(NCH),
     &          Z2(NCH), SQSCAT(NCH)
      REAL  RSOIL1(NCH), RSOIL2(NCH),    RDC(NCH),  U2FAC(NCH),
     &     QSATTC(NCH), DQSDTC(NCH),  ALWRAD(NCH), BLWRAD(NCH),
     &    TC1(NCH), TC2(NCH), TC4(NCH), QA1(NCH), QA2(NCH), QA4(NCH), 
     &        CAPAC(NCH), VGWMAX(NCH), CDCR1(NCH), CDCR2(NCH),
     &         EVAP(NCH), SHFLUX(NCH),RUNOFF(NCH), TSURF(NCH)
      REAL   EINT(NCH),    ESOI(NCH),   EVEG(NCH),   ESNO(NCH),
     &      SMELT(NCH),   HLATN(NCH),  HLWUP(NCH),
     &     RUNSRF(NCH),  QINFIL(NCH),    AR1(NCH),  BFLOW(NCH),
     &       ARS1(NCH),    ARS2(NCH),   ARS3(NCH),
     &       ARA1(NCH),    ARA2(NCH),   ARA3(NCH),   ARA4(NCH),
     &       ARW1(NCH),    ARW2(NCH),   ARW3(NCH),   ARW4(NCH),
     &       tsa1(NCH),    tsa2(NCH),   tsb1(NCH),   tsb2(NCH),   
     &        BEE(NCH),   PSIS(NCH),
     &      poros(nch),    WPWET(NCH),  cond(nch),   gnu(nch)
      REAL RZEQ(NCH), GHFLUX(NCH),  ASNOW(NCH)
      REAL BF1(NCH),    BF2(NCH),    BF3(NCH)
      REAL CATDEF(NCH),  RZEXC(NCH), SRFEXC(NCH), TP1(nch),TP2(nch),
     &  TP3(nch), TP4(nch), TP5(nch), TP6(nch), GHTCNT(6,NCH),
     &     TPSN1(NCH), TPSN2(NCH), TPSN3(NCH), WESNN(3,NCH),
     &     HTSNNN(3,NCH),SNDZN(3,NCH)

      REAL ETURB1(NCH),DEDQA1(NCH),HSTURB1(NCH),DHSDTC1(NCH),
     &     DHSDQA1(NCH),DEDTC1(NCH),QSAT1(NCH),DQS1(NCH),
     &     ALW1(NCH),BLW1(NCH)
      REAL ETURB2(NCH),DEDQA2(NCH),HSTURB2(NCH),DHSDTC2(NCH),
     &     DHSDQA2(NCH),DEDTC2(NCH),QSAT2(NCH),DQS2(NCH),
     &     ALW2(NCH),BLW2(NCH)
      REAL ETURB4(NCH),DEDQA4(NCH),HSTURB4(NCH),DHSDTC4(NCH),
     &     DHSDQA4(NCH),DEDTC4(NCH),QSAT4(NCH),DQS4(NCH),
     &     ALW4(NCH),BLW4(NCH)
      REAL ETURBS(NCH),DEDQAS(NCH),HSTURBS(NCH),DHSDTCS(NCH),
     &     DHSDQAS(NCH),DEDTCS(NCH),QSATS(NCH),DQSS(NCH),
     &     ALWS(NCH),BLWS(NCH)
C
C****
      REAL SNWMID(NTYPS)
C****

      REAL CSOIL(MXCHP)
      REAL EMAXRT(MXCHP),     RC(MXCHP), SATCAP(MXCHP)
      REAL RA1(MXCHP),RA2(MXCHP),RA4(MXCHP),RAS(MXCHP)  
      REAL     RX1(MXCHP),    RX2(MXCHP),
     &     SNWFRC(MXCHP), POTFRC(MXCHP), ESNFRC(MXCHP), EIRFRC(MXCHP),
     &       FCAN(MXCHP),   THRU(MXCHP), RZEQOL(MXCHP)
      REAL EVSNOW(MXCHP), SHFLUXS(MXCHP), HLWUPS(MXCHP)
      REAL HFTDS1(MXCHP),HFTDS2(MXCHP),HFTDS4(MXCHP),DHFT1(MXCHP),
     &    DHFT2(MXCHP),DHFT4(MXCHP)
      REAL SWSRF1(MXCHP),SWSRF2(MXCHP),SWSRF3(MXCHP),SWSRF4(MXCHP),
     &        AR2(MXCHP),   AR3(MXCHP),   AR4(MXCHP)
      REAL RX11(MXCHP),RX21(MXCHP),RX12(MXCHP),RX22(MXCHP),RX13(MXCHP),
     &     RX23(MXCHP),RX14(MXCHP),RX24(MXCHP)
      REAL TSC1(MXCHP),frice(mxchp),srfmx(mxchp),srfmn(mxchp),
     &     RCST(MXCHP),EVAPFR(MXCHP)
      REAL RCUN(MXCHP),PAR(MXCHP),PDIR(MXCHP)
      REAL RTBS1,RTBS2,RTBS4
      REAL EVAP1(MXCHP),EVAP2(MXCHP),EVAP4(MXCHP),SHFLUX1(MXCHP),
     &  SHFLUX2(MXCHP),SHFLUX4(MXCHP),HLWUP1(MXCHP),HLWUP2(MXCHP),
     &  HLWUP4(MXCHP),GHFLUX1(MXCHP),GHFLUX2(MXCHP),GHFLUX4(MXCHP),
     &  GHFLUXS(MXCHP),asnow0(mxchp),EPFRC1,EPFRC2,EPFRC4,SUMEP,SUME
      REAL TC1SF(MXCHP),TC2SF(MXCHP),TC4SF(MXCHP),TC1SN,TC2SN,TC4SN
      REAL T1(3),AREA(3),TPSN(3),ZBAR,THETAF,HT(7),TP(7),XFICE,
     &  f2ice(mxchp),fice(7),
     &  FH21,FH21W,FH21I,FH21D,DFH21W,DFH21I,DFH21D,
     &  EVSN,SHFLS,HUPS
      REAL SWNET0,HLWDWN0,TMPSNW,HLWTC,DHLWTC,HSTURB,DHSDEA,DHSDTC,
     &  ESATTC,ETURB,DEDEA,DEDTC,SNOWF,TS,fh31w,fh31i,fh31d,pr,ea,
     &  desdtc,areasc,pre
      real xklh0(7,mxchp),xklh(7),zc0(7,mxchp),zc(7),xklhw0(mxchp),xklhw
      REAL WESN(3),HTSNN(3),SNDZ(3),dummy1,dummy2,dummy3,areasc0,EDIF,
     &     EINTX

c for procedures timing
      real tstart,tfinal,tcatch,tenergy,tsnowr,tfluxes,tground
      real tstartx,tfinalx,tstart1,tfinal1,tstarts,tfinals
      real tpart1,tpart2,tpart3,tpart4,tpart5
      common/timing/tcatch,tenergy,tsnowr,tfluxes,tground
      common/timing2/tpart1,tpart2,tpart3,tpart4,tpart5

c for the correction of the unbalanced budgets (12-17-97)
      REAL ADJUST
c
cc jpw 8/4/99 ----------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
      real a11(nch),a12(nch),a21(nch),a22(nch),a23(nch),a32(nch)
      real a33(nch),q11(nch),q22(nch),q33(nch)
      logical*1 assim
c
      real es(nch),ev(nch),et(nch),qin(nch),srflw(nch),rzflw(nch)
      real bflw(nch),adj(nch),cor(nch),exc1(nch),exc2(nch)
c
cc ---------------------------------------------------------------------
c

C****

c LIR is a flag for vegetation as opposed to non vegetation cover types
c in MOSAIC, vegetation is assessed with a test on SATCAP, compared to .001
      LOGICAL*1 LIR(NTYPS),LIRX(MXCHP),BUG

C****

c     DATA TSC1 /MXCHP*432000./      ! 5 DAY TIME SCALE 
c     DATA TSC1 /MXCHP*86400./       ! 1 DAY TIME SCALE
      DATA TSC1 /MXCHP*864000./      ! 10 DAY TIME SCALE

      DATA SNWMID /50.,50.,50.,2.,50.,2.,2.,2.,2.,2./

      DATA LIR / .TRUE.,  .TRUE.,  .TRUE.,  .TRUE., 
     &           .TRUE.,  .TRUE., .FALSE., .FALSE.,
     &          .FALSE., .FALSE. /
c      DATA VGPH1 /     -100.,     -190.,     -200.,     -120.,
c     5                 -200.,     -200.,     -200.,      -10.,
c     9                  -10.,      -10. /
c      DATA VGPH2 /     -500.,     -250.,     -250.,     -230.,
c     5                 -400.,     -400.,     -250.,     -100.,
c     9                 -100.,     -100. /

      DATA RCST /MXCHP*1.E10/

C**** ---------------------------------------------------
C**** PRE-PROCESS DATA AS NECESSARY:
C****

      DO N=1,NCH
c        SATCAP(N) = 0.1 * ZLAI(N)
c change for convergence towards MOSAIC
         SATCAP(N) = 0.2 * ZLAI(N)

         CSOIL(N)  = 70000.
         FCAN(N) = AMIN1( 1., AMAX1(0.,CAPAC(N)/SATCAP(N)) )
         POTFRC(N)=FCAN(N)
c (02.02.98) EMAXRT(N) = (VGWMAX(N)+SATCAP(N)) / DTSTEP
         RA1(N)     = ONE / ( CD1(N) * UM(N) )
         RA2(N)     = ONE / ( CD2(N) * UM(N) )
         RA4(N)     = ONE / ( CD4(N) * UM(N) )
         RAS(N)     = ONE / ( CDS(N) * UM(N) ) 
         LIRX(N)=LIR(ITYP(N))
         RUNSRF(N)=0.
      ENDDO

C**** ---------------------------------------------------
C**** DETERMINE INITIAL VALUE OF RZEQ:

      CALL RZEQUIL (
     I              NCH, ITYP, CATDEF, VGWMAX,CDCR1,CDCR2,WPWET,
     I              ars1,ars2,ars3,ara1,ara2,ara3,ara4,
     I              arw1,arw2,arw3,arw4,
     O              RZEQOL
     &             )

      IF (BUG) THEN
         WRITE(*,*) 'RZEQUIL OK'
      ENDIF

c (02.02.08)
      DO N=1,NCH
c rk 10.26.99
         EMAXRT(N)=(CAPAC(N)+(RZEQOL(N)+RZEXC(N)-WPWET(N)*VGWMAX(N))+
     &             SRFEXC(N))/DTSTEP
      ENDDO


C**** PARTITION CATCHMENT INTO EVAPORATION SUBREGIONS:
C****
      CALL PARTITION (
     I                NCH,ITYP,RZEXC,
     I                RZEQOL,VGWMAX,CDCR1,CDCR2,
     I                PSIS,BEE,poros,WPWET,
     I                ars1,ars2,ars3,ara1,ara2,ara3,ara4,
     I                arw1,arw2,arw3,arw4,BUG,
     u                SRFEXC,CATDEF,
     O                AR1, AR2, AR4,srfmx,srfmn, 
     O                SWSRF1,SWSRF2,SWSRF4,
     o                cor
     &               )

      IF (BUG) THEN
         WRITE(*,*) 'PARTITION OK'
      ENDIF

C**** ========================================================
C**** ENERGY BALANCES.

C**** COMPUTE "INITIAL ESTIMATE" OF HEAT FLUX TO DEEP SOIL (HFTDS)
C**** AND ITS DERIVATIVE WITH RESPECT TO TEMPERATURE (DHFTDS):

      DO N=1,NCH
         AREA(1)=AR1(N)
         AREA(2)=AR2(N)+AR3(N)
         AREA(3)=AR4(N)
         T1(1)=TC1(N)
         T1(2)=TC2(N)
         T1(3)=TC4(N)
         ZBAR=-CATDEF(N)/100.
         THETAF=.5
         DO LAYER=2,7
            HT(LAYER)=GHTCNT(LAYER-1,N)
         ENDDO

         CALL GNDTP0(
     I              T1,ZBAR,THETAF,
     U              HT,
     O              fh21w,fH21i,fh21d,dfh21w,dfh21i,dfh21D,
     O              tp, xklh,zc,xklhw
     &             )

         HFTDS1(N)=-FH21W
         HFTDS2(N)=-FH21I
         HFTDS4(N)=-FH21D
         DHFT1(N)=-DFH21W
         DHFT2(N)=-DFH21I
         DHFT4(N)=-DFH21D
         DO LAYER=1,6
            GHTCNT(LAYER,N)=HT(LAYER+1)
         ENDDO
         tp1(n)=tp(2)
         tp2(n)=tp(3)
         tp3(n)=tp(4)
         tp4(n)=tp(5)
         tp5(n)=tp(6)
         tp6(n)=tp(7)
         frice(n)=xfice

         do k=1,7
            xklh0(k,n)=xklh(k)
            zc0(k,n)=zc(k)
         enddo
         xklhw0(n)=xklhw
      ENDDO 

      IF (BUG) THEN
         WRITE(*,*) 'HEAT FLUX INITIAL ESTIMATE OK'
      ENDIF

C**** -------------------------------------------------------------
C**** A. SNOW-FREE FRACTION.
C**** DETERMINE EVAPORATION, SENSIBLE HEAT FLUXES; UPDATE TEMPS:

      DO N=1,NCH
        PAR(N)    = PARDIR(N) + PARDIF(N) + 1.E-20
        PDIR(N)   = PARDIR(N) / PAR(N)
        TC1SF(N)  = TC1(N)
        TC2SF(N)  = TC2(N)
        TC4SF(N)  = TC4(N)
        ENDDO

      CALL RCUNST (
     I             NCH, ITYP, SUNANG, SQSCAT, PDIR,
     I             PAR, ZLAI, GREEN, BUG,
     O             RCUN
     &            )

      IF (BUG) THEN
         WRITE(*,*) 'RCUNST OK'
      ENDIF

C**** 1. SATURATED FRACTION

      CALL ENERGY1 (
     I             NCH, DTSTEP, ITYP, UM, RCUN,
     I             ETURB1, DEDQA1, DEDTC1, HSTURB1, DHSDQA1, DHSDTC1,
     I             QM,     RA1,   SWNETF,  HLWDWN, PSUR,
     I             RDC,    U2FAC, HFTDS1, DHFT1,
     I             QSAT1, DQS1, ALW1, BLW1,
     I             EMAXRT,CSOIL,SWSRF1,POTFRC,.false.,
     U             TC1SF, QA1,
     O             EVAP1, SHFLUX1, HLWUP1, RX11, RX21, GHFLUX1
     &            )

      IF (BUG) THEN
         WRITE(*,*) 'ENERGY1 OK'
      ENDIF

C**** 2. SUBSATURATED BUT UNSTRESSED FRACTION

      CALL ENERGY2 (
     I             NCH, DTSTEP, ITYP, UM, RCUN,
     I             ETURB2, DEDQA2, DEDTC2, HSTURB2, DHSDQA2, DHSDTC2,
     I             QM,     RA2,   SWNETF,  HLWDWN, PSUR,
     I             RDC,    U2FAC, HFTDS2, DHFT2,
     I             QSAT2, DQS2, ALW2, BLW2,
     I             EMAXRT,CSOIL,SWSRF2,POTFRC,.false.,
     U             TC2SF, QA2,
     O             EVAP2, SHFLUX2, HLWUP2, RX12, RX22, GHFLUX2
     &            )

      IF (BUG) THEN
         WRITE(*,*) 'ENERGY2 OK'
      ENDIF

C**** 3. WILTING FRACTION
      CALL ENERGY4 (
     I             NCH, DTSTEP, ITYP, UM, RCST,
     I             ETURB4, DEDQA4, DEDTC4, HSTURB4, DHSDQA4, DHSDTC4,
     I             QM,     RA4,   SWNETF,  HLWDWN, PSUR,
     I             RDC,    U2FAC, HFTDS4, DHFT4,
     I             QSAT4, DQS4, ALW4, BLW4,
     I             EMAXRT,CSOIL,SWSRF4,POTFRC,.false.,
     U             TC4SF, QA4,
     O             EVAP4, SHFLUX4, HLWUP4, RX14, RX24, GHFLUX4
     &            )

      IF (BUG) THEN
         WRITE(*,*) 'ENERGY3 OK'
      ENDIF

C**** COMPUTE EIRFRC
      DO N=1,NCH
         
         RTBS1=RX11(N)*RX21(N)/(RX11(N)+RX21(N)+1.E-20)
         EPFRC1=POTFRC(N) * ( RA1(N) + RTBS1 ) /
     &        ( RA1(N) + POTFRC(N)*RTBS1 )
         
         RTBS2=RX12(N)*RX22(N)/(RX12(N)+RX22(N)+1.E-20)
         EPFRC2=POTFRC(N) * ( RA2(N) + RTBS2 ) /
     &        ( RA2(N) + POTFRC(N)*RTBS2 )
         
         RTBS4=RX14(N)*RX24(N)/(RX14(N)+RX24(N)+1.E-20)
         EPFRC4=POTFRC(N) * ( RA4(N) + RTBS4 ) /
     &        ( RA4(N) + POTFRC(N)*RTBS4 )
         
         SUMEP=EPFRC1*EVAP1(N)*AR1(N)+EPFRC2*EVAP2(N)*AR2(N)+
     &        EPFRC4*EVAP4(N)*AR4(N)
         SUME=EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N)
         EIRFRC(N)=SUMEP/(SUME+1.E-20)
         
      ENDDO

      IF (BUG) THEN
         WRITE(*,*) 'EIRFRC OK'
      ENDIF

C**** --------------------------------------------------------
C**** B. SNOW-COVERED FRACTION.

      DO N=1,NCH

         if (bug) write(*,*) 'snow: n=',n

         ts=tm(n)
         T1(1)=TC1(N)-273.16
         T1(2)=TC2(N)-273.16
         T1(3)=TC4(N)-273.16
         AREA(1)=AR1(N)
         AREA(2)=AR2(N)+AR3(N)
         AREA(3)=AR4(N)
         pr=.001*(trainc(n)+trainl(n)+tsnow(n)) !  in m/s
         snowf=.001*(tsnow(n))  !  in m/s
         eturb=eturbs(n)
         dedtc=dedtcs(n)
         dedea=dedqas(n)*epsilon/psur(n)
         dhsdtc=dhsdtcs(n)
         dhsdea=dhsdqas(n)*epsilon/psur(n)
         ea=qm(n)*psur(n)/epsilon

c**** ask Marc: use (dedqas,dhsdqas) or (dedeas,dhsdeas)?
         hsturb=hsturbs(n)
c**** ask Marc: use (qsats, dqss) or (esats, dess)?
         esattc=qsats(n)*psur(n)/epsilon
         desdtc=dqss(n)*psur(n)/epsilon
         hlwtc=ALWS(N) + BLWS(N)*TPSN1(N)
         dhlwtc=BLWS(N)
         swnet0=swnets(n)
         hlwdwn0=hlwdwn(n)
         do k=1,7
            xklh(k)=xklh0(k,n)
            zc(k)=zc0(k,n)
         enddo
         xklhw=xklhw0(n)
          
         wesn(1)=wesnn(1,n)
         wesn(2)=wesnn(2,n)
         wesn(3)=wesnn(3,n)
         htsnn(1)=htsnnn(1,n)
         htsnn(2)=htsnnn(2,n)
         htsnn(3)=htsnnn(3,n)
         sndz(1)=sndzn(1,n)
         sndz(2)=sndzn(2,n)
         sndz(3)=sndzn(3,n)

C**** 1. RUN SNOW MODEL:

        CALL SNOWRT(
     I             ts,t1,area,pr,snowf,DTSTEP,eturb,dedtc,dedea,
     I             dhsdtc,dhsdea,esattc,desdtc,dhlwtc,ea,
     I             swnet0,hlwdwn0,hlwtc,hsturb,
     I             xklhw,xklh,zc,bug,
     U             wesn,htsnn,sndz,
     O             tpsn,areasc,pre,fh31w,fh31i,fh31d,HUPS,
     O             EVSN,SHFLS,areasc0
     &            )         

        TPSN1(N)=TPSN(1)+TF
        TPSN2(N)=TPSN(2)+TF
        TPSN3(N)=TPSN(3)+TF
        asnow(n)=areasc
        asnow0(n)=areasc0
        areasc0=areasc

        wesnn(1,n)=wesn(1)
        wesnn(2,n)=wesn(2)
        wesnn(3,n)=wesn(3)
        htsnnn(1,n)=htsnn(1)
        htsnnn(2,n)=htsnn(2)
        htsnnn(3,n)=htsnn(3)
        sndzn(1,n)=sndz(1)
        sndzn(2,n)=sndz(2)
        sndzn(3,n)=sndz(3)
        SMELT(N)=PRE*1000.

        trainc(n)=trainc(n)*(1.-areasc)
        trainl(n)=trainl(n)*(1.-areasc)

C
C**** 2. UPDATE SURFACE TEMPERATURE
C
        TC1SN=TC1(N)+
     &        (-(FH31W/(area(1)+1.e-20))-HFTDS1(N))*DTSTEP/CSOIL(N)
        TC2SN=TC2(N)+
     &        (-(FH31I/(area(2)+1.e-20))-HFTDS2(N))*DTSTEP/CSOIL(N)
        TC4SN=TC4(N)+
     &        (-(FH31D/(area(3)+1.e-20))-HFTDS4(N))*DTSTEP/CSOIL(N)
        TC1(N)=TC1SF(N)*(1-AREASC0)+TC1SN*AREASC0
        TC2(N)=TC2SF(N)*(1-AREASC0)+TC2SN*AREASC0
        TC4(N)=TC4SF(N)*(1-AREASC0)+TC4SN*AREASC0
        
        EVSNOW(N)=EVSN/alhs
        esno(n)=evsnow(n)*asnow(n)*DTSTEP ! to have esno in mm/20min (03-17-99)
        SHFLUXS(N)=SHFLS
        HLWUPS(N)=HUPS
        GHFLUXS(N)=(AREA(1)*HFTDS1(N)+AREA(2)*HFTDS2(N)+
     &            AREA(3)*HFTDS4(N))
        ENDDO

      IF (BUG) THEN
         WRITE(*,*) 'SNOW FRACTION OK'
      ENDIF

      DO N=1,NCH
         HLATN(N)=(1.-ASNOW(N))*
     &        (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))*ALHE
     &        +ASNOW(N)*EVSNOW(N)*ALHS
         EVAP(N)=(1.-ASNOW(N))*
     &        (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))
     &        +ASNOW(N)*EVSNOW(N)
         EVAPFR(N)=(1.-ASNOW(N))*
     &        (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))
         SHFLUX(N)=(1.-ASNOW(N))*
     &        (SHFLUX1(N)*AR1(N)+SHFLUX2(N)*AR2(N)+SHFLUX4(N)*AR4(N))
     &        +ASNOW(N)*SHFLUXS(N)
         HLWUP(N)=(1.-ASNOW(N))*
     &        (HLWUP1(N)*AR1(N)+HLWUP2(N)*AR2(N)+HLWUP4(N)*AR4(N))
     &        +ASNOW(N)*HLWUPS(N)
         GHFLUX(N)=(1.-ASNOW(N))*
     &        (GHFLUX1(N)*AR1(N)+GHFLUX2(N)*AR2(N)+GHFLUX4(N)*AR4(N))
     &        +ASNOW(N)*GHFLUXS(N)

      ENDDO

      IF (BUG) THEN
         WRITE(*,*) 'ENERGY FLUXES OK'
      ENDIF

C****
C**** NOW ALLOW DEEPER SOIL TEMPERATURES TO BE UPDATED:

      DO N=1,NCH
         ZBAR=-CATDEF(N)/100.
         THETAF=.5
         DO LAYER=2,7
            HT(LAYER)=GHTCNT(LAYER-1,N)
         ENDDO
         FH21=-GHFLUX(N)

         CALL GNDTMP(
     I        dtstep,zbar,thetaf,fh21,
     U        ht,
     O        xfice,tp)

         DO LAYER=1,6
            GHTCNT(LAYER,N)=HT(LAYER+1)
         ENDDO
         tp1(n)=tp(2)
         tp2(n)=tp(3)
         tp3(n)=tp(4)
         tp4(n)=tp(5)
         tp5(n)=tp(6)
         tp6(n)=tp(7)
         frice(n)=xfice
         
      ENDDO

      IF (BUG) THEN
         WRITE(*,*) 'DEEPER SOIL TEMPERATURES UPDATE OK'
      ENDIF

C**** ========================================================

C**** REMOVE EVAPORATED WATER FROM SURFACE RESERVOIRS:
C****
C**** (FIRST CORRECT FOR EXCESSIVE INTERCEPTION LOSS)

      DO N=1,NCH
        EINTX=EIRFRC(N)*EVAPFR(N)*DTSTEP
        IF(EINTX .GT. CAPAC(N)) THEN
          EDIF=(EINTX-CAPAC(N))/DTSTEP
          EVAPFR(N)=EVAPFR(N)-EDIF
          EVAP(N)=EVAP(N)-EDIF
          HLATN(N)=HLATN(N)-EDIF*ALHE
          SHFLUX(N)=SHFLUX(N)+EDIF*ALHE
          EIRFRC(N)=CAPAC(N)/(EVAPFR(N)*DTSTEP)
          ENDIF
        ENDDO
c
cc jpw 9/15/99 ---------------------------------------------------------
c  added variables assim, q11, q22 and q33
c ----------------------------------------------------------------------
c
      CALL WUPDAT (
     I               NCH,   ITYP, DTSTEP,
     I               EVAPFR, SATCAP, TC1, RA1, RC, 
     I               RX11,RX21,RX12,RX22,RX13,RX23,RX14,RX24,
     I               AR1,AR2,AR3,AR4,CDCR1,
     I               RX1, RX2,EIRFRC,RZEQOL,srfmn,WPWET,VGWMAX,assim,
     U               CAPAC, RZEXC, CATDEF, SRFEXC,
     O               EINT, ESOI, EVEG,
     o               q11, q22, q33,
     o               es, ev, et
     &              )
c
cc ---------------------------------------------------------------------
c
      IF (BUG) THEN
         WRITE(*,*) 'WUPDAT OK'
      ENDIF

C**** REDISTRIBUTE MOISTURE BETWEEN RESERVOIRS:

c
cc jpw 8/4/99 ----------------------------------------------------------
c  added variables assim, a11, a12, a21, a22, a23, a32 and a33
c ----------------------------------------------------------------------
c
      CALL RZDRAIN (
     I              NCH,DTSTEP,TSC1,VGWMAX,SATCAP,RZEQOL,AR1,
     I              tsa1,tsa2,tsb1,tsb2, CDCR2,BUG,assim,
     U              CAPAC,RZEXC,SRFEXC,CATDEF,RUNSRF,
     O              a11,a12,a21,a22,a23,a32,a33,
     o              srflw,rzflw,exc1,exc2
     &              )
c
cc ---------------------------------------------------------------------
c
      IF (BUG) THEN
         WRITE(*,*) 'RZDRAIN OK'
      ENDIF

C**** COMPUTE BASEFLOW FROM TOPMODEL EQUATIONS

c
cc jpw 8/4/99 ----------------------------------------------------------
c  added variables assim and a33
c ----------------------------------------------------------------------
c
      CALL BASE (
     I           NCH, DTSTEP,BF1, BF2, BF3, CDCR1, FRICE, COND, GNU,
     I           CDCR2,
     i           assim,
     U           CATDEF,
     u           a33,
     O           BFLOW,bflw
     &          )
c
cc ---------------------------------------------------------------------
c
      IF (BUG) THEN
         WRITE(*,*) 'BASE OK'
      ENDIF

C**** UPDATE CANOPY INTERCEPTION; DETERMINE THROUGHFALL RATES.

      CALL INTERC (
     I             NCH, ITYP, DTSTEP, TRAINL, TRAINC, SMELT, 
     I             SATCAP, CSOIL, SFRAC,BUG,
     U             TC1, CAPAC, GHFLUX,
     O             THRU
     &            )

      IF (BUG) THEN
         WRITE(*,*) 'INTERC OK'
      ENDIF

C**** DETERMINE SURFACE RUNOFF AND INFILTRATION RATES:

c
cc jpw 9/15/99 ---------------------------------------------------------
c  added variables assim and q11
c ----------------------------------------------------------------------
c
      CALL SRUNOFF (
     I              NCH,DTSTEP,AR1,ar2,ar4,
     I              THRU,frice,tp1,srfmx,BUG,assim,
     U              RZEXC,SRFEXC,RUNSRF,
     u              q11,
     O              QINFIL,qin
     &             )
c     write(15,'(10(f10.5,1x))') a11(1),a12(1),a21(1),a22(1),a23(1),
c    &     a32(1),a33(1),q11(1),q22(1),q33(1)
c
cc ---------------------------------------------------------------------
c
      IF (BUG) THEN
         WRITE(*,*) 'SRUNOFF'
      ENDIF

C**** (ADD CHECK TO ENSURE RZEXC KEPT WITHIN BOUNDS IN SRUNOFF)
      
C**** RECOMPUTE RZEXC:

      CALL RZEQUIL (
     I              NCH, ITYP, CATDEF, VGWMAX,CDCR1,CDCR2,WPWET,
     I              ars1,ars2,ars3,ara1,ara2,ara3,ara4,
     I              arw1,arw2,arw3,arw4,
     O              RZEQ
     &             )

      IF (BUG) THEN
         WRITE(*,*) 'RZEQUIL'
      ENDIF

      DO N=1,NCH
         ADJ(N)=0.5*(RZEQOL(N)-RZEQ(N))
         RZEXC(N)=RZEXC(N)+ADJ(N)
         CATDEF(N)=CATDEF(N)+ADJ(N)
      ENDDO

C**** ---------------------------------------------------
C**** PROCESS DATA AS NECESSARY PRIOR TO RETURN:
C****
C**** ---------------------------------------------------

      DO N=1,NCH

         IF(CAPAC(N).LT.1.E-10) CAPAC(N) = 0.0
c        IF(SNOW (N).LT.1.E-10) SNOW (N) = 0.0
         RUNOFF(N) = RUNSRF(N)+BFLOW(N)
         EINT(N) = EINT(N) * ALHE / DTSTEP
         ESOI(N) = ESOI(N) * ALHE / DTSTEP
         EVEG(N) = EVEG(N) * ALHE / DTSTEP
         ESNO(N) = ESNO(N) * ALHS / DTSTEP
         
         TSURF(N)=AR1(N)*TC1(N)+AR2(N)*TC2(N)+AR4(N)*TC4(N)
         TSURF(N)=(1.-ASNOW(N))*TSURF(N)+ASNOW(N)*TPSN1(N)

         if(asnow(n) .eq. 0) then
            tpsn1(n)=TSURF(N)
            tpsn2(n)=TSURF(N)
            tpsn3(n)=TSURF(N)
         endif

      ENDDO     
      
      RETURN
      END

C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
C**** [ BEGIN RCUNST ]
C****
      SUBROUTINE RCUNST (
     I                   NCH, ITYP, SUNANG, SQSCAT, PDIR,
     I                   PAR, ZLAI, GREEN,BUG,
     O                   RCUN
     &                  )
C****
C****     This subroutine calculates the unstressed canopy resistance.
C**** (p. 1353, Sellers 1985.)  Extinction coefficients are computed first.
C****
      IMPLICIT NONE
C****
C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MemFac
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON

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

C****
      LOGICAL*1 BUG
      INTEGER NCH
      INTEGER ITYP(NCH), ChNo

      REAL  SUNANG(NCH),   PDIR(NCH),   PAR(NCH),   ZLAI(NCH),
     &      SQSCAT(NCH),  GREEN(NCH),  RCUN(NCH)

      REAL   VGCHIL(NTYPS), VGZMEW(NTYPS),         
     &       VGRST1(NTYPS), VGRST2(NTYPS), VGRST3(NTYPS)

      REAL            RHO4,         EXTK1,         EXTK2,
     &               RCINV,         GAMMA,          EKAT,    DUM1,
     &                DUM2,          DUM3,            AA,      BB,
     &                  ZK,            CC


      DATA VGCHIL /        0.1,        0.25,        0.01,        -0.3,
     5                    0.01,        0.20,         0.0,         0.0,
     9                     0.0,         0.0 /

      DATA VGZMEW/      0.9809,      0.9638,      0.9980,      1.0773,
     5                  0.9980,      0.9676,       1.000,       1.000,
     9                   1.000,       1.000 /

      DATA VGRST1 /     2335.9,      9802.2,      2869.7,      2582.0,
     5                 93989.4,      9802.2,         0.0,         0.0,
     9                     0.0,         0.0 /

      DATA VGRST2 /        0.0,        10.6,         3.7,         1.1,
     5                    0.01,        10.6,         0.0,         0.0,
     9                     0.0,         0.0 /

      DATA VGRST3 /      153.5,       180.0,       233.0,       110.0,
     5                   855.0,       180.0,         1.0,         1.0,
     9                     1.0,         1.0 /




      DO 100 ChNo = 1, NCH

C**** First compute optical parameters.
C**** (Note: CHIL is constrained to be >= 0.01, as in SiB calcs.)

      AA = 0.5 - (0.633 + 0.330*VGCHIL(ITYP(ChNo)))*VGCHIL(ITYP(ChNo))
      BB = 0.877 * ( ONE - 2.*AA )
      CC =  ( AA + BB*SUNANG(ChNo) ) / SUNANG(ChNo)

      EXTK1 =  CC * SQSCAT(ChNo)
      EXTK2 = (ONE / VGZMEW(ITYP(ChNo))) * SQSCAT(ChNo)

      DUM1 =      PDIR(ChNo)  *   CC
      DUM2 = (ONE-PDIR(ChNo)) * ( BB*(ONE/3.+PIE/4.) + AA*1.5 )

C**** Bound extinction coefficient by 50./ZLAI:

      ZK =     PDIR(ChNo) *AMIN1( EXTK1, 50./ZLAI(ChNo) ) +
     &    (ONE-PDIR(ChNo))*AMIN1( EXTK2, 50./ZLAI(ChNo) )

C**** Now compute unstressed canopy resistance:

      GAMMA = VGRST1(ITYP(ChNo)) / VGRST3(ITYP(ChNo)) +
     &        VGRST2(ITYP(ChNo))

      EKAT = EXP( ZK*ZLAI(ChNo) )
      RHO4 = GAMMA / (PAR(ChNo) * (DUM1 + DUM2))

      DUM1 = (VGRST2(ITYP(ChNo)) - GAMMA) / (GAMMA + 1.E-20)
      DUM2 = (RHO4 * EKAT + ONE) / (RHO4 + ONE)
      DUM3 = ZK * VGRST3(ITYP(ChNo))

      RCINV = ( DUM1*ALOG(DUM2) + ZK*ZLAI(ChNo) ) / DUM3         

      RCUN(ChNo) = ONE / (RCINV * GREEN(ChNo) + 1.E-10)

 100  CONTINUE


      RETURN
      END
C****
C**** [ END RCUNST ]

C**** ===================================================
C**** ///////////////////////////////////////////////////
C**** ===================================================

      SUBROUTINE SRUNOFF (
     I                    NCH,DTSTEP,AR1,ar2,ar4,
     I                    THRU,frice,tp1,srfmx,BUG,
     i                    assim,
     U                    RZEXC,SRFEXC,RUNSRF,
     u                    q11,
     O                    QINFIL,QIN
     &                   )

      IMPLICIT NONE
      INTEGER NCH,N
      REAL DTSTEP,AR1(NCH),ar2(nch),ar4(nch),
     &     THRU(NCH),frice(nch),tp1(nch),
     &     RZEXC(NCH),SRFEXC(NCH),RUNSRF(NCH),QINFIL(NCH),srfmx(nch)
      REAL PTOTAL,srun0,frun
      LOGICAL*1 BUG
c
cc jpw 9/15/99 ---------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
      real q11(nch)
      logical*1 assim
c
      real qin(nch)
c
cc ---------------------------------------------------------------------
c

C**** - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO N=1,NCH

         PTOTAL=THRU(N)
         frun=AR1(N)
         if(srfexc(n) .gt. 0.) then
            frun=frun+ar2(n)*(srfexc(n)/(srfmx(n)+1.e-20))
c           frun=frun+ar4(n)*(srfexc(n)/(srfmx(n)+1.e-20))**2
            endif
         frun=frun+(1-frun)*frice(n)
         srun0=PTOTAL*frun

c**** Comment out this line in order to allow moisture
c**** to infiltrate soil (until snow model is fixed):
         if(tp1(n) .lt. 0.) srun0=ptotal

         if(ptotal-srun0 .gt. srfmx(n)-srfexc(n)) 
     &                srun0=ptotal-(srfmx(n)-srfexc(n)) 

         if (srun0 .gt. ptotal) then
            write(*,*) 'srun0 > ptotal: N=',N
            write(*,*) ' frice=',frice(n),' ar1=',ar1(n),' ptotal=',
     &           ptotal,' tp1=',tp1(n)
            write(*,*) ' ar2=',ar2(n),' ar4=',ar4(n),' srfexc=',
     &           srfexc(n),' srfmx=',srfmx(n),' thru=',thru(n),
     &           ' rzexc=',rzexc(n)
            write(*,*) '=====> CORRECTION'
            srun0=ptotal
         endif

         RUNSRF(N)=RUNSRF(N)+srun0
         QIN(N)=PTOTAL - srun0

         if (qin(n) .lt. 0.) then
            write(*,*) 'Bomb: qin lt 0.=',qin(n),' frice=',
     &           frice(n),' ar1=',ar1(n),' ptotal=',ptotal, ' NCH=',NCH,
     &           ' N=',N,' tp1=',tp1(n)
            stop
         endif

c        if (bug .eq. .true.) then
c            write(*,*) 'Bomb: qin ge 0.=',qin(n),' frice=',
c     &           frice(n),' ar1=',ar1(n),' ptotal=',ptotal, ' NCH=',NCH
c            write(*,*) 'It went through!!!'
c         endif

         SRFEXC(N)=SRFEXC(N)+QIN(N)
         RUNSRF(N)=RUNSRF(N)/DTSTEP
         QINFIL(N)=QIN(N)/DTSTEP
c
cc jpw 9/15/99 ---------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
         if (assim) then
           q11(n)=sqrt(q11(n)**2.+qinfil(n)**2.)
         end if
c
cc ---------------------------------------------------------------------
c
         
      ENDDO
      
      RETURN
      END

C**** ===================================================
C**** ///////////////////////////////////////////////////
C**** ===================================================

      SUBROUTINE RZDRAIN (
     I                    NCH,DTSTEP,TSC1,VGWMAX,SATCAP,RZEQ,AR1,
     I                    tsa1,tsa2,tsb1,tsb2,CDCR2,BUG,assim,
     U                    CAPAC,RZEXC,SRFEXC,CATDEF,RUNSRF,
     O                    a11,a12,a21,a22,a23,a32,a33,
     o                    srflw,rzflw,exc1,exc2
     &                   )

c-----------------------------------------------------------------
c        defines drainage timescales:
c             - tsc0, between srfex and rzex
c             - tsc2, between rzex and catdef
c        then defines correponding drainages
c        and updates the water contents
c-----------------------------------------------------------------

      IMPLICIT NONE
      INTEGER NCH,N
      REAL DTSTEP,TSC1(NCH),RZEXC(NCH),SRFEXC(NCH),CATDEF(NCH),
     &    VGWMAX(NCH),CAPAC(NCH),SATCAP(NCH),RUNSRF(NCH),RZEQ(NCH)
      REAL tsa1(NCH),tsa2(NCH),tsb1(NCH),tsb2(NCH),CDCR2(NCH),AR1(NCH)
      REAL FLOW,EXCESS,TSC0,tsc2,rzave,rz0,wanom,alog2,rztot,
     &      a1,a2,b1,b2,rzx,ax,bx

      logical*1 BUG
      data alog2/0.693147/
c
cc jpw 8/4/99 ----------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
      integer flag
      real a11(nch),a12(nch),a21(nch),a22(nch),a23(nch),a32(nch)
      real a33(nch)
      logical*1 assim
c
      real aa,bb,cc,x1
c
      real srflw(nch),rzflw(nch),exc1(nch),exc2(nch)
c
cc ---------------------------------------------------------------------
c

C**** - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO N=1,NCH

c        TSC0=86400./24. !! this is the tuned newp2
c        if(srfexc(n) .lt. 0.) tsc0=86400.

c****   Compute equivalent of root zone excess in non-saturated area:
        rztot=rzeq(n)+rzexc(n)
        if(ar1(n).ne.1.) then
            rzave=(rztot-ar1(n)*vgwmax(n))/(1.-ar1(n))
            rzave=rzave*0.4/vgwmax(n)
          else
            rzave=0.4
          endif

c**** tsc0, between srfex and rzex
        if(srfexc(n).ge.0.) then
          rz0=amax1(.001,rzave+srfexc(n)/(1000.*0.0329821))
          tsc0=0.0529276/(rz0**2.5)
          endif
        if(srfexc(n).lt.0.) then
          rz0=amax1(.001,rzave+srfexc(n)/(1000.*0.107143))
          tsc0=0.684850/(rz0**2)
          endif

        tsc0=tsc0*3600.
        if(tsc0.lt.dtstep) tsc0=dtstep
c
cc jpw 9/7/99 ----------------------------------------------------------
c evaluate derivative of flow0 with respect to srfexc (a21) and rzexc 
c (a12)
c ----------------------------------------------------------------------
c
        if (assim) then
          if (tsc0 .lt. dtstep) then
            a12(n)=0.
            a21(n)=1.
          else if (rz0 .eq. .001) then
            a12(n)=0.
            a21(n)=dtstep/tsc0
          else
            if (srfexc(n) .ge. 0.) then
              aa=0.0529276
              bb=32.9821
              cc=2.5
            else
              aa=0.6848500
              bb=107.1430
              cc=2.
            end if
            if (ar1(n) .ne. 1.) then
              x1=srfexc(n)/bb+0.4*(rzeq(n)+rzexc(n)-ar1(n)*vgwmax(n))/
     &               ((1.-ar1(n))*vgwmax(n))
              a12(n)=(0.4*cc*dtstep*srfexc(n)*x1**(cc-1.))/
     &               (3600.*aa*(1.-ar1(n))*vgwmax(n))
              a21(n)=(cc*dtstep*srfexc(n)*x1**(cc-1.))/(3600.*aa*bb)+
     &               (dtstep*x1**cc)/(3600.*aa)
            else
              a12(n)=0.
              a21(n)=(cc*dtstep*srfexc(n)*(0.4+srfexc(n)/bb)**(cc-1.))/
     &               (3600.*aa*bb)+(dtstep*(0.4+srfexc(n)/bb)**cc)/
     &               (3600.*aa)
            end if
          end if
        end if
c
cc ---------------------------------------------------------------------
c
        SRFLW(N)=SRFEXC(N)*DTSTEP/TSC0
        RZEXC(N)=RZEXC(N)+SRFLW(N)
        SRFEXC(N)=SRFEXC(N)-SRFLW(N)

c**** Topography-dependent tsc2, between rzex and catdef

        rzx=rzexc(n)/vgwmax(n)

        if(rzx .gt. .01) then
            ax=tsa1(n)
            bx=tsb1(n)
          elseif(rzx .lt. -.01) then
            ax=tsa2(n)
            bx=tsb2(n)
          else
            ax=tsa2(n)+(rzx+.01)*(tsa1(n)-tsa2(n))/.02
            bx=tsb2(n)+(rzx+.01)*(tsb1(n)-tsb2(n))/.02
          endif

        tsc2=exp(ax+bx*catdef(n))

c        tsc2=(drpar1(n)+drpar2(n)*100.*(rzexc(n)/vgwmax(n)))
c     &        *1000./catdef(n)
c        tsc2=amax1(0., amin1( 1., tsc2 ))
c this is a test !!!
c        tsc2=tsc2*exp(-catdef(n)/1000.)
c        tsc2=tsc2*exp(-catdef(n)/500.)

        rzflw(n)=rzexc(n)*tsc2*dtstep/3600.
c sensitivity test (08-11-98)
c        rzflw(n)=AMAX1(rzflw(n),0.) 

C        if(rzflw(n).lt.0.) rzflw(n)=rzflw(n)/100.

c 05.12.98: first attempt to include bedrock
c
cc jpw 8/4/99 ----------------------------------------------------------
c added flag
c ----------------------------------------------------------------------
c
        IF (CATDEF(N)-RZFLW(N) .GT. CDCR2(N)) then
          RZFLW(N)=CATDEF(N)-CDCR2(N)
          flag=0
        else
          flag=1
        end if
c
cc jpw 9/7/99 ----------------------------------------------------------
c evaluate derivative of flow2 with respect to rzexc (a32) and catdef 
c (a23)
c ----------------------------------------------------------------------
c
        if (assim) then
          if (flag .eq. 0) then
            a23(n)=0.
            a32(n)=tsc2*dtstep/3600.
          else
            if (rzx.ge.-0.01 .and. rzx.le.0.01) then
              x1=exp(ax+bx*catdef(n))
              a23(n)=bx*dtstep*x1*rzexc(n)/3600.
              a32(n)=(dtstep/3600.)*(x1+x1*rzexc(n)*50.*
     &               ((tsa1(n)-tsa2(n))+catdef(n)*(tsb1(n)-tsb2(n)))/
     &               vgwmax(n))
            else
              x1=exp(ax+bx*catdef(n))
              a23(n)=bx*dtstep*x1*rzexc(n)/3600.
              a32(n)=dtstep*x1/3600.
            end if
          end if
        end if
c
cc ---------------------------------------------------------------------
c
        CATDEF(N)=CATDEF(N)-RZFLW(N)
        RZEXC(N)=RZEXC(N)-RZFLW(N)

C****   REMOVE ANY EXCESS FROM MOISTURE RESERVOIRS:

        IF(CAPAC(N) .GT. SATCAP(N)) THEN
          EXC1(N)=CAPAC(N)-SATCAP(N)
          CAPAC(N)=SATCAP(N)
          RZEXC(N)=RZEXC(N)+EXC1(N)
        ELSE 
          EXC1(N)=0.
          ENDIF

        IF(RZEQ(N) + RZEXC(N) .GT. VGWMAX(N)) THEN
          EXCESS=RZEQ(N)+RZEXC(N)-VGWMAX(N)
          RZFLW(N)=RZFLW(N)-EXCESS
          RZEXC(N)=VGWMAX(N)-RZEQ(N)
          CATDEF(N)=CATDEF(N)-EXCESS
          ENDIF

        IF(CATDEF(N) .LT. 0.) THEN
          EXC2(N)=CATDEF(N)
          RUNSRF(N)=RUNSRF(N)-CATDEF(N)
          CATDEF(N)=0.
        ELSE
          EXC2(N)=0.
          ENDIF
c
cc jpw 9/7/99 ----------------------------------------------------------
c  compute remaining elements of the A matrix for covariance forecasting
c ----------------------------------------------------------------------
c
          if (assim) then
            a11(n)=1.-a21(n)
            a22(n)=1.+a12(n)-a32(n)
            a33(n)=1.-a23(n)
            a12(n)=-a12(n)
            a23(n)=-a23(n)
            a32(n)=-a32(n)
          end if
c
cc ---------------------------------------------------------------------
c

        ENDDO

      RETURN
      END
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
C**** [ BEGIN WUPDAT ]
C****
      SUBROUTINE WUPDAT (
     I                     NCH,   ITYP, DTSTEP,
     I                     EVAP, SATCAP, TC, RA, RC,
     I                     RX11,RX21,RX12,RX22,RX13,RX23,RX14,RX24,
     I                     AR1,AR2,AR3,AR4,CDCR1,
     I                     RX1, RX2,EIRFRC,RZEQ,srfmn,WPWET,VGWMAX,
     i                     assim,
     U                     CAPAC, RZEXC, CATDEF, SRFEXC,
     O                     EINT, ESOI, EVEG,
     o                     q11,q22,q33,
     o                     es,ev,et
     &                    )
C****
C**** THIS SUBROUTINE ALLOWS EVAPOTRANSPIRATION TO ADJUST THE WATER
C**** CONTENTS OF THE INTERCEPTION RESERVOIR AND THE SOIL LAYERS.
C****
      IMPLICIT NONE
C****
C****
C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MEMFAC
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON

      PARAMETER (NTYPS = 10, FRSTCH = 1, MEMFAC = 5)
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

C****
      INTEGER NCH
      INTEGER ITYP(NCH), CHNO
      REAL    EVAP(NCH),  SATCAP(NCH),
     &          TC(NCH),      RA(NCH),     RC(NCH),
     &       CAPAC(NCH),  CATDEF(NCH),    RX1(NCH),  SRFEXC(NCH),
     &         RX2(NCH),  EIRFRC(NCH),    RZEQ(NCH), srfmn(nch),
     &        RX11(NCH),RX21(NCH),RX12(NCH),RX22(NCH),RX13(NCH),
     &        RX23(NCH),RX14(NCH),RX24(NCH),AR1(NCH),AR2(NCH),AR3(NCH),
     &         AR4(NCH),CDCR1(NCH),WPWET(NCH),VGWMAX(NCH)
      REAL EINT(NCH), ESOI(NCH), EVEG(NCH), RZEXC(NCH)
      REAL  DTSTEP, EGRO, CNDSAT, CNDUNS, ESATFR, cndv, cnds, WILT
c
cc jpw 9/15/99 ---------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
      real q11(nch),q22(nch),q33(nch)
      logical*1 assim
c
      real es(nch),ev(nch),et(nch)
c
cc ---------------------------------------------------------------------
c
C****
C**** -----------------------------------------------------------------

      DO 100 CHNO = 1, NCH

C**** COMPUTE EFFECTIVE SURFACE CONDUCTANCES IN SATURATED AND UNSATURATED
C**** AREAS:

      CNDSAT=(AR1(CHNO)/RX11(CHNO)) + (AR1(CHNO)/RX21(CHNO))
      CNDUNS=(AR2(CHNO)/RX12(CHNO)) + (AR2(CHNO)/RX22(CHNO)) +
     &       (AR4(CHNO)/RX14(CHNO)) + (AR4(CHNO)/RX24(CHNO))
     
      ESATFR=CNDSAT/(CNDSAT+CNDUNS)

C****
C**** PARTITION EVAP BETWEEN INTERCEPTION AND GROUND RESERVOIRS.
C****

      EINT(CHNO)=EIRFRC(CHNO)*EVAP(CHNO)*DTSTEP
      EGRO = EVAP(CHNO)*DTSTEP - EINT(CHNO)

C**** ENSURE THAT INDIVIDUAL CAPACITIES ARE NOT EXCEEDED.

      IF(EINT(CHNO) .GT. CAPAC(CHNO)) THEN
        EGRO=EGRO+EINT(CHNO)-CAPAC(CHNO)
        EINT(CHNO)=CAPAC(CHNO)
        ENDIF

      WILT=WPWET(CHNO)*VGWMAX(CHNO)
      IF(EGRO .GT. RZEXC(CHNO)+RZEQ(CHNO)-WILT) THEN
c 06.02.98: the minimum is designed to prevent truncation errors 
        EINT(CHNO)=AMIN1(CAPAC(CHNO),
     &        EINT(CHNO)+EGRO-(RZEXC(CHNO)+RZEQ(CHNO)-WILT))
        EGRO=RZEXC(CHNO)+RZEQ(CHNO)-WILT
        ENDIF
c RK: the above test ensures in particular that rzexc+rzeq never < 0.

      ESOI(CHNO)=AMIN1(SRFEXC(CHNO)-SRFMN(CHNO),
     &        EGRO*(AR1(CHNO)*RX11(CHNO)/(RX11(CHNO)+RX21(CHNO)+1.E-20)
     &        +AR2(CHNO)*RX12(CHNO)/(RX12(CHNO)+RX22(CHNO)+1.E-20)
     &        +AR4(CHNO)*RX14(CHNO)/(RX14(CHNO)+RX24(CHNO)+1.E-20)))
      EVEG(CHNO)=EGRO-ESOI(CHNO)

C****
C**** SPECIAL CASE FOR CONDENSATION:
      IF(EVAP(CHNO) .LT. 0.) THEN
          EINT(CHNO)=EVAP(CHNO)*DTSTEP
c 05.20.98: to prevent negative throughfall due to truncation errors 
          EINT(CHNO)=AMIN1(0.,EINT(CHNO))
          ESOI(CHNO)=0.
          EVEG(CHNO)=0.
        ENDIF

C****
C**** REMOVE MOISTURE FROM RESERVOIRS:
C****

      IF (CATDEF(CHNO) .LT. CDCR1(CHNO)) THEN
         CAPAC(CHNO) = AMAX1(0., CAPAC(CHNO) - EINT(CHNO))
         RZEXC(CHNO) = RZEXC(CHNO) - EVEG(CHNO)*(1.-ESATFR)
         SRFEXC(CHNO) = SRFEXC(CHNO) - ESOI(CHNO)*(1.-ESATFR)
         CATDEF(CHNO) = CATDEF(CHNO) + (ESOI(CHNO) + EVEG(CHNO))*ESATFR
         es(chno)=esoi(chno)*(1.-esatfr)
         ev(chno)=eveg(chno)*(1.-esatfr)
         et(chno)=(esoi(chno)+eveg(chno))*esatfr
c 05.12.98: first attempt to include bedrock
      ELSE
         CAPAC(CHNO) = AMAX1(0., CAPAC(CHNO) - EINT(CHNO))
         RZEXC(CHNO) = RZEXC(CHNO) -  EVEG(CHNO)
         SRFEXC(CHNO) = SRFEXC(CHNO) - ESOI(CHNO)
         es(chno)=esoi(chno)
         ev(chno)=eveg(chno)
         et(chno)=0.
      ENDIF
c
cc jpw 9/15/99 ---------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
      if (assim) then
        if (catdef(chno) .lt. cdcr1(chno)) then
          q11(chno)=esoi(chno)*(1.-esatfr)
          q22(chno)=eveg(chno)*(1.-esatfr)
          q33(chno)=(esoi(chno)+eveg(chno))*esatfr
        else
          q11(chno)=esoi(chno)
          q22(chno)=eveg(chno)
          q33(chno)=0.
        end if
      end if
c
cc ---------------------------------------------------------------------
c

C****
 100  CONTINUE
C****
      RETURN
      END
C****
C**** [ END WUPDAT ]
C****
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
C**** [ BEGIN INTERC ]
C**** 
      SUBROUTINE INTERC (
     I                   NCH, ITYP, DTSTEP, TRAINL, TRAINC,SMELT,
     I                   SATCAP, CSOIL, SFRAC,BUG,
     U                   TC, CAPAC,GHFLUX,
     O                   THRU
     &                  )
C****
C**** THIS ROUTINE USES THE PRECIPITATION FORCING TO DETERMINE 
C**** CHANGES IN INTERCEPTION AND SOIL MOISTURE STORAGE.
CCCCCCCCC Changes in snowcover are not treated here anymore
C****
C**** NOTE:  IN THIS FORMULATION, RAIN THAT FALLS ON FROZEN GROUND
C**** RUNS-OFF, RATHER THAN FREEZES.
C****
      IMPLICIT NONE
C****
C**** CHIP HEADER FILE
C****
      INTEGER  NTYPS, FRSTCH, MEMFAC
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON

      PARAMETER (NTYPS = 10, FRSTCH = 1, MEMFAC = 5)
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

C****
      LOGICAL*1 BUG
      INTEGER NCH
      INTEGER ITYP(NCH), CHNO
      REAL  TRAINL(NCH), TRAINC(NCH),  SATCAP(NCH),  
     &          TC(NCH),  CSOIL(NCH),   CAPAC(NCH),     
     &       SMELT(NCH),   THRU(NCH),  GHFLUX(NCH),
     &            SFRAC,      WETINT
      REAL DTSTEP, WATADD, CAVAIL, THRUC, TIMFRL, TIMFRC, 
     &     FWETL, FWETC, THRU1, THRU2, THRUL, XTCORR,SMPERS
C****
C      DATA FWET0  /0.30/
      DATA FWETL /1.00/, FWETC /0.20/, TIMFRL/1.00/
      DATA TIMFRC/0.333/
c value for GSWP
c      TIMFRC/0.125/  

C****
C**** ------------------------------------------------------------------
C**** LOOP OVER CHIPS:
      DO 100 CHNO = 1, NCH

C**** =======================================================
C****
C**** LOAD INTERCEPTION RESERVOIR.  STEP 1: LARGE SCALE CONDENSATION.
C****
C**** DETERMINE XTCORR, THE FRACTION OF A STORM THAT FALLS ON A PREVIOUSLY
C**** WET SURFACE DUE TO THE TIME CORRELATION OF PRECIPITATION POSITION.
C**** (TIME SCALE TIMFRL FOR LARGE SCALE STORMS SET TO ONE FOR FWETL=1
C**** TO REFLECT THE EFFECTIVE LOSS OF "POSITION MEMORY" WHEN STORM 
C**** COVERS ENTIRE GRID SQUARE.)

      XTCORR= (1.-TIMFRL) * 
     &      AMIN1( 1.,(CAPAC(CHNO)/SATCAP(CHNO))/(FWETL*SFRAC) )   

C****
C**** FILL INTERCEPTION RESERVOIR WITH PRECIPITATION.
C**** THRU1 IS FIRST CALCULATED AS THE AMOUNT FALLING THROUGH THE 
C****    CANOPY UNDER THE ASSUMPTION THAT ALL RAIN FALLS RANDOMLY.  
C****    ONLY A FRACTION 1-XTCORR FALLS RANDOMLY, THOUGH, SO THE RESULT 
C****    IS MULTIPLIED BY 1-XTCORR.
C****
      WATADD = TRAINL(CHNO)*DTSTEP + SMELT(CHNO)*DTSTEP
      CAVAIL = ( SATCAP(CHNO) - CAPAC(CHNO) ) * (FWETL*SFRAC)
      WETINT = CAPAC(CHNO)/SATCAP(CHNO)
      IF( WATADD*(1.-WETINT) .LT. CAVAIL ) THEN
          THRU1 = WATADD*WETINT
        ELSE
          THRU1 = (WATADD - CAVAIL)
        ENDIF
      THRU1=THRU1*(1.-XTCORR)

C**** THRU2 IS THE AMOUNT THAT FALLS IMMEDIATELY THROUGH THE CANOPY DUE
C**** TO 'POSITION MEMORY'.

      THRU2=XTCORR*WATADD

      THRUL=THRU1+THRU2

      CAPAC(CHNO)=CAPAC(CHNO)+WATADD-THRU1-THRU2

C****
C**** ---------------------------------------------------
C****
C**** STEP 2: MOIST CONVECTIVE PRECIPITATION.
C****
C**** DETERMINE XTCORR, THE FRACTION OF A STORM THAT FALLS ON A PREVIOUSLY
C**** WET SURFACE DUE TO THE TIME CORRELATION OF PRECIPITATION POSITION.
    
      XTCORR= (1.-TIMFRC) * 
     &     AMIN1( 1.,(CAPAC(CHNO)/SATCAP(CHNO))/(FWETC*SFRAC) )

C****
C**** FILL INTERCEPTION RESERVOIR WITH PRECIPITATION.
C**** THRU1 IS FIRST CALCULATED AS THE AMOUNT FALLING THROUGH THE 
C****    CANOPY UNDER THE ASSUMPTION THAT ALL RAIN FALLS RANDOMLY.  
C****    ONLY A FRACTION 1-XTCORR FALLS RANDOMLY, THOUGH, SO THE RESULT 
C****    IS MULTIPLIED BY 1-XTCORR.
C****
      WATADD = TRAINC(CHNO)*DTSTEP
      CAVAIL = ( SATCAP(CHNO) - CAPAC(CHNO) ) * (FWETC*SFRAC)
      WETINT = CAPAC(CHNO)/SATCAP(CHNO)
      IF( WATADD*(1.-WETINT) .LT. CAVAIL ) THEN
          THRU1 = WATADD*WETINT
        ELSE
          THRU1 = (WATADD - CAVAIL)
        ENDIF
      THRU1=THRU1*(1.-XTCORR)

C**** THRU2 IS THE AMOUNT THAT FALLS IMMEDIATELY THROUGH THE CANOPY DUE
C**** TO 'POSITION MEMORY'.

      THRU2=XTCORR*WATADD

      THRUC=THRU1+THRU2
      CAPAC(CHNO)=CAPAC(CHNO)+WATADD-THRU1-THRU2
C****
      IF (THRUL+THRUC .LT. -1e-5) WRITE(*,*) 'THRU= ',THRUL+THRUC
      THRU(CHNO)=AMAX1(0., THRUL+THRUC)

 100  CONTINUE
C****
      RETURN
      END
C****
C**** [ END INTERC ]
C****
C****
C**** -----------------------------------------------------------------
C**** /////////////////////////////////////////////////////////////////
C**** -----------------------------------------------------------------
C****
      SUBROUTINE BASE (
     I                 NCH,DTSTEP,BF1,BF2,BF3,CDCR1,FRICE,COND,GNU,
     I                 CDCR2,
     i                 assim,
     U                 CATDEF,
     u                 a33,
     O                 BFLOW,bflw
     &                )

      IMPLICIT NONE
      INTEGER NCH,N
      REAL BF1(NCH),BF2(NCH),BF3(NCH),CATDEF(NCH)
      REAL COND(NCH), GNU(NCH)
      REAL BFLOW(NCH),CDCR1(NCH),CDCR2(NCH),FRICE(NCH)
      REAL DTSTEP, ZBAR
c
cc jpw 8/4/99 ----------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
      real a33(nch)
      logical*1 assim
c
      real bflw(nch)
c
cc ---------------------------------------------------------------------
c
      DO N=1,NCH
         ZBAR=SQRT(catdef(n)/bf1(n))-bf2(n)
         BFLOW(N)=(1.-FRICE(N))*1000.*
     &        cond(n)*exp(-(bf3(n)-3.)-gnu(n)*zbar)/gnu(n)
c *1000 is to convert from m/s to mm/s
         IF (CATDEF(N) .GE. CDCR1(N)) BFLOW(N)=0.
c
cc jpw 9/7/99 ----------------------------------------------------------
c include base flow component in the A matrix for covariance forecasting
c ----------------------------------------------------------------------
c
         if (assim .and. bflow(n) .ne. 0.) then
           a33(n)=a33(n)-(500.*cond(n)*dtstep*(1.-frice(n))*
     &            exp(3.-bf3(n)-gnu(n)*(sqrt(catdef(n)/bf1(n))-bf2(n)))/
     &            (bf1(n)*sqrt(catdef(n)/bf1(n))))
         end if
c
cc ---------------------------------------------------------------------
c
         bflw(n)=bflow(n)*dtstep
         CATDEF(N)=CATDEF(N)+BFLW(N)
      ENDDO

      RETURN
      END




