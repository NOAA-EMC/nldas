      SUBROUTINE SFCDIF(ZLM,Z0,THZ0,THLM,SFCSPD,CZIL,AKMS,AKHS)

      IMPLICIT NONE

C     ******************************************************************
C     *                                                                *
C     *                        SURFACE LAYER                           *
C     *                                                                *
C     ******************************************************************
C-----------------------------------------------------------------------

      REAL WWST, WWST2, G, VKRM, EXCM, BETA, BTG, ELFC, WOLD, WNEW
      REAL PIHF, EPSU2, EPSUST, EPSIT, EPSA, ZTMIN, ZTMAX, HPBL, SQVISC
      REAL RIC, RRIC, FHNEU, RFC, RFAC, ZZ, PSLMU, PSLMS, PSLHU, PSLHS
      REAL XX, PSPMU, YY, PSPMS, PSPHU, PSPHS, ZLM, Z0, THZ0, THLM
      REAL SFCSPD, CZIL, AKMS, AKHS, ZILFC, ZU, ZT, RDZ, CXCH
      REAL DTHV, DU2, BTGH, WSTAR2, USTAR, ZSLU, ZSLT, RLOGU, RLOGT
      REAL RLMO, ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4
      REAL XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN, RLMA
C      REAL ZTFC

      INTEGER ITRMX, ILECH, ITR

                             PARAMETER
     &(WWST=1.2,WWST2=WWST*WWST,G=9.8,VKRM=0.40,EXCM=0.001
     &,BETA=1./270.,BTG=BETA*G,ELFC=VKRM*BTG
     &,WOLD=.15,WNEW=1.-WOLD,ITRMX=05,PIHF=3.14159265/2.
C-----------------------------------------------------------------------
     &,EPSU2=1.E-4,EPSUST=0.07,EPSIT=1.E-4,EPSA=1.E-8
     &,ZTMIN=-5.,ZTMAX=1.,HPBL=1000.0
     &,SQVISC=258.2)

      PARAMETER
     &(RIC=0.183,RRIC=1.0/RIC,FHNEU=0.8,RFC=0.191
     &,RFAC=RIC/(FHNEU*RFC*RFC))

C-----------------------------------------------------------------------
C **** LECH'S SURFACE FUNCTIONS ****
      PSLMU(ZZ)=-0.96*log(1.0-4.5*ZZ)
      PSLMS(ZZ)=ZZ*RRIC-2.076*(1.-1./(ZZ+1.))
      PSLHU(ZZ)=-0.96*log(1.0-4.5*ZZ)
      PSLHS(ZZ)=ZZ*RFAC-2.076*(1.-1./(ZZ+1.))

C-----------------------------------------------------------------------
C **** PAULSON'S SURFACE FUNCTIONS *****
      PSPMU(XX)=-2.*log((XX+1.)*0.5)-log((XX*XX+1.)*0.5)+2.*ATAN(XX)
     &          -PIHF
      PSPMS(YY)=5.*YY
      PSPHU(XX)=-2.*log((XX*XX+1.)*0.5)
      PSPHS(YY)=5.*YY
C
      ILECH=0
C***********************************************************************
C-----------------------------------------------------------------------
C ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
CC......ZTFC=0.1
C CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
C Commented out to allow argument-passed value
C      CZIL=0.2
      ZILFC=-CZIL*VKRM*SQVISC

C-----------------------------------------------------------------------
      ZU=Z0
CC.......ZT=Z0*ZTFC
C-----------------------------------------------------------------------
      RDZ=1./ZLM
      CXCH=EXCM*RDZ
C-----------------------------------------------------------------------
      DTHV=THLM-THZ0
C
      DU2=MAX(SFCSPD*SFCSPD,EPSU2)
C--------------BELJARS CORRECTION OF USTAR------------------------------
      BTGH=BTG*HPBL
      WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)
      USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
C--------------ZILITINKEVITCH APPROACH FOR ZT--------------------------------
      ZT=EXP(ZILFC*SQRT(USTAR*Z0))*Z0
C-----------------------------------------------------------------------
        ZSLU=ZLM+ZU
        ZSLT=ZLM+ZT
C
        RLOGU=log(ZSLU/ZU)
        RLOGT=log(ZSLT/ZT)
C
        RLMO=ELFC*AKHS*DTHV/USTAR**3
C-----------------------------------------------------------------------
            DO 100 ITR=1,ITRMX
C--------------1./MONIN-OBUKKHOV LENGTH-SCALE---------------------------
        ZETALT=MAX(ZSLT*RLMO,ZTMIN)
        RLMO=ZETALT/ZSLT
        ZETALU=ZSLU*RLMO
C
        ZETAU=ZU*RLMO
        ZETAT=ZT*RLMO
C-----------------------------------------------------------------------
         IF(ILECH.EQ.0) THEN
          IF(RLMO.LT.0.)THEN
            XLU4=1.-16.*ZETALU
            XLT4=1.-16.*ZETALT
            XU4 =1.-16.*ZETAU
            XT4 =1.-16.*ZETAT
C
            XLU=SQRT(SQRT(XLU4))
            XLT=SQRT(SQRT(XLT4))
            XU =SQRT(SQRT(XU4))
            XT =SQRT(SQRT(XT4))
C
            PSMZ=PSPMU(XU)
            SIMM=PSPMU(XLU)-PSMZ+RLOGU
            PSHZ=PSPHU(XT)
            SIMH=PSPHU(XLT)-PSHZ+RLOGT
          ELSE
            ZETALU=MIN(ZETALU,ZTMAX)
            ZETALT=MIN(ZETALT,ZTMAX)
C
            PSMZ=PSPMS(ZETAU)
            SIMM=PSPMS(ZETALU)-PSMZ+RLOGU
            PSHZ=PSPHS(ZETAT)
            SIMH=PSPHS(ZETALT)-PSHZ+RLOGT
          ENDIF
         ELSE
C-----------------------------------------------------------------------
C ***** LECH'S FUNCTIONS ****
          IF(RLMO.LT.0.)THEN
            PSMZ=PSLMU(ZETAU)
            SIMM=PSLMU(ZETALU)-PSMZ+RLOGU
            PSHZ=PSLHU(ZETAT)
            SIMH=PSLHU(ZETALT)-PSHZ+RLOGT
          ELSE
            ZETALU=MIN(ZETALU,ZTMAX)
            ZETALT=MIN(ZETALT,ZTMAX)
C
            PSMZ=PSLMS(ZETAU)
            SIMM=PSLMS(ZETALU)-PSMZ+RLOGU
            PSHZ=PSLHS(ZETAT)
            SIMH=PSLHS(ZETALT)-PSHZ+RLOGT
          ENDIF
         ENDIF
C--------------BELJAARS CORRECTION FOR USTAR----------------------------
        USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
C--------------ZILITINKEVITCH FIX FOR ZT--------------------------------
        ZT=EXP(ZILFC*SQRT(USTAR*Z0))*Z0

        ZSLT=ZLM+ZT
        RLOGT=log(ZSLT/ZT)
C-----------------------------------------------------------------------
        USTARK=USTAR*VKRM
        AKMS=MAX(USTARK/SIMM,CXCH)
        AKHS=MAX(USTARK/SIMH,CXCH)
C-----------------------------------------------------------------------
        WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)
        RLMN=ELFC*AKHS*DTHV/USTAR**3
C-----------------------------------------------------------------------
        RLMA=RLMO*WOLD+RLMN*WNEW
C-----------------------------------------------------------------------
C      IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110
C-----------------------------------------------------------------------
        RLMO=RLMA
C-----------------------------------------------------------------------
  100   CONTINUE
C-----------------------------------------------------------------------
  110   CONTINUE
C-----------------------------------------------------------------------
        RETURN
        END
