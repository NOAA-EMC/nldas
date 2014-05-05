C MEMBER EX19
C  (from old member FCEX19)
C
      SUBROUTINE EX19(PS,CS,PX,TA,RM,PCTS,RSTS,OWE,SWE,OSC,COVER,
C                             LAST UPDATE: 06/22/95.13:27:10 BY $WC30EA
C                           
CVK  NEW SIMULATED SNOW DEPTH TIME SERIES ADDED
CVK     1   PPX,PPCTS,PRM,TALR)
     1   SNOWD,PPX,PPCTS,PRM,TALR)
          
C.......................................
C     THIS IS THE EXECUTION SUBROUTINE FOR THE SNOW MODEL OPERATION
C        BASED ON NOAA TECH. MOMO. NWS HYDRO-17.
C.......................................
C     SUBROUTINE INITIALLY WRITTEN BY...
C        ERIC ANDERSON - HRL   MAY 1980

c        updated 3/10/97 by mbs for uadj mod

C        UPDATED 4/00 BY V. KOREN TO INCLUDE SNOW DEPTH/TEMPERATURE 
C.......................................
      REAL MFMAX,MFMIN,NMF,LIQW,NEGHS,MBASE,MFC
      DIMENSION PS(1),CS(1),PX(1),TA(1),RM(1),PCTS(1),OWE(1),SWE(1)
      DIMENSION OSC(1),COVER(1),TALR(1),PPX(1),PPCTS(1),PRM(1),RSTS(1)
      DIMENSION SNAME(2),LASTDA(12),IDANG(12),CT(17)
CVK
      DIMENSION SNOWD(1)
C
C     COMMON BLOCKS
      COMMON/OUTCTL/IOUTYP
      COMMON/SNPM19/ALAT,SCF,MFMAX,MFMIN,NMF,UADJ,SI,MBASE,PXTEMP,
     1  PLWHC,TIPM,PA,ADC(11),LMFV,SMFV(12),LAEC,NPTAE,AE(2,14)
      COMMON/SNCO19/WE,NEGHS,LIQW,TINDEX,ACCMAX,SB,SBAESC,SBWS,
C
CVK  ADDED TWO MORE STATES: SNDPT, SNTMP
CVK     1  STORGE,AEADJ,NEXLAG,EXLAG(7)
     1  STORGE,AEADJ,NEXLAG,EXLAG(7),SNDPT,SNTMP
          
      COMMON/SUMS19/SPX,SSFALL,SRM,SMELT,SMELTR,SROBG,DSFALL,DRAIN,
     1 DQNET,DRSL,NDRSP
      COMMON/SNUP19/MFC,SFALLX,WINDC,SCTOL,WETOL,SNOF,UADJC
      COMMON/FSNWUP/IUPWE,IUPSC
      COMMON/FSDATA/NSDV,JHSNW(10),WESNW(10),AESNW(10),WEADD(10)
      COMMON/FPXTYP/NRSV,IJHRS(20),LJHRS(20),PCTSV(20)
      COMMON/FMFC19/NMFC,IJHMF(20),LJHMF(20),VMFC(20)
      COMMON/FADJ19/NUADJ,IJHUA(20),LJHUA(20),VUADJ(20)
C
C    ================================= RCS keyword statements ==========
      CHARACTER*68     RCSKW1,RCSKW2
      DATA             RCSKW1,RCSKW2 /                                 '
C    ===================================================================
C
C
C     DATA STATEMENTS
      DATA SNAME/4HEX19,4H    /
      DATA LASTDA/31,28,31,30,31,31,31,31,30,31,30,31/
      DATA IDANG/285,316,345,10,40,71,101,132,163,193,224,254/
      DATA JSNW/4HSNW1/
C.......................................
C     CONSTANT
C     IF SNOWFALL IS LESS THAN HSNOF/HR -- DO NOT
C        LEAVE CURRENT DEPLETION CURVE.
   90 HSNOF=0.2
C.......................................
C     VALUES OF CONTROL VARIABLES.
C        THOSE NEEDED IN ALL CASES FIRST.
      IDT=PS(14)
      FIDT=IDT
      ITPX=PS(10)
      SNOF=HSNOF*ITPX
      NDT=IDT/ITPX
      NEXLAG=5/ITPX+2
      LRM=PS(17)
      LSWE=PS(20)
      IF (LSWE.EQ.0) GO TO 103
      ITSWE=PS(LSWE+3)
      NSWE=24/ITSWE
  103 LCOVER=PS(22)
CVK
CVK      IF (LCOVER.EQ.0) GO TO 104
      IF (LCOVER.EQ.0) GO TO 216
            
      ITSSC=PS(LCOVER+3)
      NSSC=24/ITSSC

CVK -----------------------------
  216 LSSDPT=PS(31)
      IF (LSSDPT.EQ.0) GO TO 104
      ITSN=PS(LSSDPT+3)
      NSSN=24/ITSN
CVK ----------------------------
  
  104 LSUMS=PS(23)
      LPM=PS(25)
      IPRINT=PS(24)
      IF (MAINUM.NE.1) GO TO 102
      IF (IPRSNW.EQ.0) GO TO 102
      IF (IPRSNW.EQ.1) IPRINT=1
      IF (IPRSNW.EQ.-1) IPRINT=0
  102 IF(NOPROT.EQ.1) IPRINT=0
      IF(NOSNOW.EQ.1) IPRINT=0
C     THE COMMENT SHOULD BE REMOVED FROM THE NEXT LINE WHEN
C     THE SNOW MODEL DISPLAY IS ADDED TO THE ICP PROGRAM.
      IF (IOUTYP.EQ.2) IPRINT=0
      NPS=PS(16)
      IF(NOSNOW.EQ.1) GO TO 101
C
C     THOSE ONLY NEEDED IF NOSNOW=0
      LPCTS=PS(18)
      LOWE=PS(19)
      IF (LOWE.EQ.0) GO TO 108
      ITOWE=PS(LOWE+3)
      NOWE=24/ITOWE
  108 LOSC=PS(21)
      IF (LOSC.EQ.0) GO TO 109
      ITOSC=PS(LOSC+3)
      NOSC=24/ITOSC
  109 LADC=PS(26)
      LTAPM=PS(27)
      LUPPM=PS(28)
      LMFV=PS(29)
      LAEC=PS(30)
C.......................................
C     GET TEMPERATURE INFORMATION AND LAPSE RATES FOR EACH IDT PERIOD
C        IF NEEDED.
  100 IF(NOSNOW.EQ.1) GO TO 110
      IF(LTAPM.EQ.0) GO TO 200
      TAELEV=PS(LTAPM)
      ELEV=PS(LPM+1)
      EDIFF=(TAELEV-ELEV)*0.01
      J=24/IDT
      TALMAX=PS(LTAPM+1)
      TALMIN=PS(LTAPM+2)
      DIFF=TALMAX-TALMIN
      FLOCAL=LOCAL
      DO 105 I=1,J
      FI=I
      TI=(FI-1.0)*FIDT+0.5*FIDT
      TL=TI+FLOCAL
      IF(TL.GT.24.0) TL=TL-24.0
      IF(TL.LT.15.0) GO TO 106
      TALR(I)=TALMAX-((TL-15.0)/15.0)*DIFF
      GO TO 105
  106 IF(TL.GT.6.0) GO TO 107
      TALR(I)=TALMIN+((6.0-TL)/15.0)*DIFF
      GO TO 105
  107 TALR(I)=TALMIN+((TL-6.0)/9.0)*DIFF
  105 CONTINUE
C.......................................
C     RAIN-SNOW ELEVATION PARAMETERS
  200 IF (LAEC.EQ.0) GO TO 110
      NPTAE=PS(LAEC)
      DO 205 I=1,NPTAE
      J=LAEC+5+(I-1)*2
      AE(1,I)=PS(J)
      AE(2,I)=PS(J+1)
  205 CONTINUE
C.......................................
C     PARAMETER VALUES
  110 PXADJ=PS(LPM)
      IF(NOSNOW.EQ.1) GO TO 115
      ELEV=PS(LPM+1)*0.01
      PA=29.9-0.335*ELEV+0.00022*(ELEV**2.4)
      PA=33.86*PA
      SCF=PS(LPM+2)
      MFMAX=PS(LPM+3)
      MFMAX=(MFMAX*FIDT)/6.0
      MFMIN=PS(LPM+4)
      MFMIN=(MFMIN*FIDT)/6.0
      NMF=PS(LPM+7)
      NMF=(NMF*FIDT)/6.0
      UADJ=PS(LPM+5)
      UADJ=(UADJ*ITPX)/6.0
      SI=PS(LPM+6)
      GM=PS(LPM+12)
      PGM=(GM*FIDT)/24.0
      MBASE=PS(LPM+9)
      PXTEMP=PS(LPM+10)
      PLWHC=PS(LPM+11)
      TIPM=PS(LPM+8)
      TIPM=1.0-((1.0-TIPM)**(FIDT/6.0))
      ALAT=PS(LPM+13)
C
C     AREAL DEPLETION CURVE
      ADC(1)=0.05
      DO 111 I=2,10
      J=LADC+I-2
  111 ADC(I)=PS(J)
      ADC(11)=1.0
      IF(LMFV.EQ.0) GO TO 32
      DO 31 I=1,12
       J=LMFV+I-1
   31  SMFV(I)=PS(J)
C
C     PASS PARAMETER VALUES TO GRAPHICS INTERFACE IF REQUESTED.
   32 IF (IOUTYP.EQ.0) GO TO 30
      IF (IDA.NE.IDARUN) GO TO 30
      CALL ICPI19(OPNAME,IDT,ITPX,PXADJ,ELEV,LTAPM,TAELEV,SCF,
     1  MFMAX,MFMIN,NMF,UADJ,SI,GM,MBASE,PXTEMP,PLWHC,TIPM,ALAT,
     2  ADC,LMFV,SMFV,LAEC,NPTAE,AE)
C.......................................
C     UPDATING PARAMETERS
   30 IF(MAINUM.GT.2) GO TO 112
C
C     OPERATIONAL PROGRAMS - ALL PARAMETERS DEFINED.
      WETOL=PS(LUPPM)
      SCTOL=PS(LUPPM+1)
      MFC=1.0
      
c        next statement by mbs 3/13/97 for uadj mod
         uadjc=1.0
      
      SFALLX=PS(LUPPM+5)
      WINDC=PS(LUPPM+6)
      GO TO 115
C
C     CALIBRATION PROGRAMS
  112 IUPWE=0
      IUPSC=0
      IF(LUPPM.EQ.0) GO TO 114
      IF(LOWE.GT.0) IUPWE=1
      IF(LOSC.GT.0) IUPSC=1
      WETOL=PS(LUPPM)
      SCTOL=PS(LUPPM+1)
  114 MFC=1.0
  
c        next statement by mbs 3/13/97 for uadj mod
         uadjc=1.0
       
      SFALLX=1.0
      WINDC=1.0
C.......................................
C     CARRYOVER VALUES.
  115 WE=CS(1)
      NEGHS=CS(2)
      LIQW=CS(3)
      TINDEX=CS(4)
      ACCMAX=CS(5)
      SB=CS(6)
      SBAESC=CS(7)
      SBWS=CS(8)
      STORGE=CS(9)
      AEADJ=CS(10)
      TEX=0.0
      DO 116 I=1,NEXLAG
      EXLAG(I)=CS(10+I)
  116 TEX=TEX+EXLAG(I)

CVK  ADDED TWO MORE NEW STATES
      SNDPT=CS(11+NEXLAG)
      SNTMP=CS(12+NEXLAG)

C
C     IF NOSNOW=1 SET SNOW CARRYOVER TO ZERO.
      IF(NOSNOW.EQ.0) GO TO 117
      CALL ZERO19
      TWE1=0.0
      CHGWE=0.0
      GO TO 118
C.......................................
C     INITIAL VALUES
  117 TWE=WE+LIQW+TEX+STORGE
      TWE1=TWE
      CHGWE=0.0
C
C     CHECK IF INITIAL WATER-EQUIVALENT OR SNOW COVER ARE CHANGED
C        BY MOD INPUT.  OPERATIONAL ONLY.
      IF(MAINUM.NE.1) GO TO 119
      JH=(IDA-1)*24+IHR-IDT
      IF(NSDV.EQ.0) GO TO 119
      DO 10 I=1,NSDV
      N=I
      IF(JH.EQ.JHSNW(I)) GO TO 11
  10  CONTINUE
      GO TO 119
  11  POWE=WESNW(N)
      POSC=AESNW(N)
      IF (WEADD(N).LT.0.0) GO TO 14
      IF (WESNW(N).LT.0.0) GO TO 13
      POWE=POWE+WEADD(N)
      GO TO 14
   13 POWE=TWE+WEADD(N)
   14 IF((POWE.LT.0.0).AND.(POSC.LT.0.0)) GO TO 119
C
C     COMPUTE AREAL EXTENT BASED ON CONDITIONS AT THE START OF THE RUN.
      AESC=0.0
      IF(TWE.EQ.0.0) GO TO 12
      CALL AESC19(WE,LIQW,ACCMAX,SB,SBAESC,SBWS,SI,ADC,AEADJ,AESC)
  12  PCOVER=AESC
C
C     MAKE UPDATE.
      CWE=TWE
      CAESC=PCOVER
      CALL UPDT19(POWE,POSC,TWE,PCOVER,IUPWE,IUPSC,WETOL,
     1SCTOL,SI,ADC)
      IF ((TWE.EQ.CWE).AND.(PCOVER.EQ.CAESC)) GO TO 119
      CHGWE=CHGWE+(TWE-CWE)
C
C     PRINT UPDATE IF PRINTSNOW ON.
      IF(IPRINT.EQ.0) GO TO 119
      WRITE(IPR,904)OPNAME,ISEG
 904  FORMAT(1H0,51HINITIAL SNOW CONDITIONS CHANGED--(SNOW-17 OPERATION,
     11X,2A4,1X,7HSEGMENT,1X,2A4,1H))
      IF(METRIC.EQ.0) GO TO 15
      WRITE(IPR,905)CWE,TWE,CAESC,PCOVER
 905  FORMAT(1H ,10X,12HWE(MM)  OLD=,F5.0,2X,4HNEW=,F5.0,5X,
     1   16HSNOW COVER  OLD=,F4.2,2X,4HNEW=,F4.2)
      GO TO 119
  15  ECWE=CWE/25.4
      ETWE=TWE/25.4
      WRITE(IPR,906)ECWE,ETWE,CAESC,PCOVER
 906  FORMAT(1H ,10X,12HWE(IN)  OLD=,F6.2,2X,4HNEW=,F6.2,5X,
     1  16HSNOW COVER  OLD=,F4.2,2X,4HNEW=,F4.2)
  119 DSFALL=0.0
      DRAIN=0.0
      DQNET=0.0
      DRSL=0.0
      NDRSP=0
      PTWE=TWE
  118 SPX=0.0
      SSFALL=0.0
      SRM=0.0
      SMELT=0.0
      SMELTR=0.0
      SROBG=0.0
      ITITLE=0
      IC=1
      IRS=1
      IMF=1
C     FOLLOWING STATEMENT ADDED BY MBS 3/10/97 FOR UADJ MOD
      IUADJ=1
C.......................................
C     INITIAL TIMING VALUES
      KDA=IDA
      KHR=IHR
      LAST=0
      IF (MAINUM.NE.1) GO TO 310
      IF (IFFG.EQ.0) GO TO 310
C     CALCULATE THE COMPUTATIONAL PERIOD AT OR JUST BEFORE LSTCMPDY
      KDAFFG=LDACPD
      KHRFFG=(LHRCPD/IDT)*IDT
      IF (KHRFFG.GT.0) GO TO 310
      KDAFFG=KDAFFG-1
      KHRFFG=24
  310 IF (NOSNOW.EQ.1) GO TO 120
      CALL MDYH1(KDA,KHR,I,J,L,N,100,0,TZ)
      IDN=IDANG(I)+J
      IMN=I
      IYR=L
C.......................................
C     GET INPUT DATA NEEDED FOR COMPUTATIONAL PERIOD.
  120 KOFF=KDA-IDADAT
C
C     PRECIPITATION DATA - ALSO APPLY ADJUSTMENT
      LPX=KOFF*(24/ITPX)+KHR/ITPX
      DO 121 I=1,NDT
      J=LPX-(NDT-I)
  121 PPX(I)=PX(J)*PXADJ
      IF(NOSNOW.EQ.1) GO TO 150
C
C     PERCENT SNOWFALL DATA
      DO 122 I=1,NDT
      IF(LPCTS.EQ.0) GO TO 123
      J=LPX-(NDT-I)
      PPCTS(I)=PCTS(J)
      GO TO 122
  123 PPCTS(I)=-999.0
  122 CONTINUE
C
C     RAINSNOW, uadj,  AND MFC MODS, OPERATIONAL PROGRAM ONLY
      IF (MAINUM.NE.1) GO TO 124
C     OVERRIDE FORM OF PRECIPITATION, IF MOD INPUT.
      IF(NRSV.EQ.0) GO TO 225
 128  IF(IRS.GT.NRSV) GO TO 225
      JH2=JH+IDT
      IF(JH2.LT.IJHRS(IRS)) GO TO 225
      JH1=JH2-IDT+ITPX
      IF(JH1.LT.IJHRS(IRS)) JH1=IJHRS(IRS)
      IF(JH2.GT.LJHRS(IRS)) JH2=LJHRS(IRS)
      J1=(JH1-JH)/ITPX
      J2=(JH2-JH)/ITPX
      DO 129 I=J1,J2
 129  PPCTS(I)=PCTSV(IRS)
      IF(JH2.LT.LJHRS(IRS)) GO TO 225
      IRS=IRS+1
      GO TO 128
C     DETERMINE MFC BASED ON MOD INPUT
  225 IF (NMFC.EQ.0) GO TO 1124
      IF (IMF.GT.NMFC) GO TO 229
      JH2=JH+IDT
      IF (JH2.LT.IJHMF(IMF)) GO TO 229
      MFC=VMFC(IMF)
      IF (JH2.LT.LJHMF(IMF)) GO TO 1124
      IMF=IMF+1
      GO TO 1124
  229 MFC=1.0
C ....................................................
C     UADJ MOD BY MIKE SMITH 3/10/97
     
C     DETERMINE UADJC BASED ON MOD INPUT
C      UADJC IS MULTIPLIER APPLIED TO UADJ PARAMETER

 
 1124 IF(NUADJ.EQ.0)      GO TO 124
      IF(IUADJ.GT.NUADJ)  GO TO 339
       JH2=JH+IDT
       IF(JH2.LT.IJHUA(IUADJ)) GO TO 339
        UADJC=VUADJ(IUADJ)
        IF(JH2.LT.LJHUA(IUADJ)) GO TO 124
        IUADJ=IUADJ+1
        GO TO 124
  339  UADJC=1.0
C.......................................................
      
      
C
C     AIR TEMPERATURE DATA - CORRECT FOR LAPSE RATE.
  124 L=KOFF*(24/IDT)+KHR/IDT
      PTA=TA(L)
      IF(LTAPM.EQ.0) GO TO 210
      I=KHR/IDT
      PTALR=TALR(I)
      PTA=PTA+EDIFF*PTALR
C
C     RAIN-SNOW ELEVATION DATA
  210 IF (LAEC.GT.0) GO TO 215
      PRSL=-999.0
      GO TO 125
  215 PRSL=RSTS(L)
C
C     OBSERVED WATER-EQUIVALENT AND SNOW COVER DATA.
C        CHECK FOR RUN-TIME DATA FIRST - OPERATIONAL ONLY.
  125 POWE=-99.0
      POSC=-99.0
      IF(MAINUM.NE.1) GO TO 130
      JH=JH+IDT
      IF(NSDV.EQ.0) GO TO 130
      DO 126 I=1,NSDV
      N=I
      IF (JH.EQ.JHSNW(I)) GO TO 127
  126 CONTINUE
      GO TO 130
  127 POWE=WESNW(N)
      POSC=AESNW(N)
      IF (WEADD(N).LT.0.0) GO TO 130
      IF (WESNW(N).LT.0.0) GO TO 131
      POWE=POWE+WEADD(N)
      GO TO 130
  131 POWE=TWE+WEADD(N)
C
C     WATER-EQUIVALENT.                                    00000010
 130  IF(POWE.NE.-99.0) GO TO 137
      IF(LOWE.EQ.0) GO TO 133
      IF((KHR/ITOWE)*ITOWE.EQ.KHR) GO TO 134
 133  POWE=-999.0
      GO TO 137
 134  L=KOFF*NOWE+KHR/ITOWE
      POWE=OWE(L)
C
C     AREAL EXTENT
 137  IF(POSC.NE.-99.0) GO TO 140
      IF(LOSC.EQ.0) GO TO 138
      IF((KHR/ITOSC)*ITOSC.EQ.KHR) GO TO 139
 138  POSC=-999.0
      GO TO 140
  139 L=KOFF*NOSC+KHR/ITOSC
      POSC=OSC(L)
C.......................................
C     DETERMINE IF COMPUTATIONAL PERIOD IS IN THE FUTURE.
  140 IFUT=0
      IF(KDA.LT.LDACPD) GO TO 145
      IF((KDA.EQ.LDACPD).AND.(KHR.LE.LHRCPD)) GO TO 145
      IFUT=1
C.......................................
C     PERFORM SNOW MODEL COMPUTATIONS FOR THE COMPUTATIONAL PERIOD.
C
  145 continue
      CALL PACK19(KDA,KHR,NDT,PTA,PPX,PPCTS,PRSL,POWE,POSC,PGM,PRM,TWE,
     1PCOVER,CWE,CAESC,IFUT,IDT,IBUG,IDN,IMN,IYR,IOUTYP,OPNAME)
      IF (TWE.NE.CWE) CHGWE=CHGWE+(TWE-CWE)
C
C.......................................
C     STORE RESULTS - RAIN+MELT FIRST.
      IF (LRM.EQ.0) GO TO 155
      DO 151 I=1,NDT
      J=LPX-(NDT-I)
  151 RM(J)=PRM(I)
      GO TO 155
C
C     NOSNOW=1 RESULTS.
  150 DO 152 I=1,NDT
      SPX=SPX+PPX(I)
      SRM=SRM+PPX(I)
      SROBG=SROBG+PPX(I)
      TWE=0.0
      PCOVER=0.0
      IF (LRM.EQ.0) GO TO 152
      J=LPX-(NDT-I)
      RM(J)=PPX(I)
  152 CONTINUE
C
C     SIMULATED WATER-EQUIVALENT
  155 IF(LSWE.EQ.0) GO TO 160
      IF((KHR/ITSWE)*ITSWE.NE.KHR) GO TO 160
      L=KOFF*NSWE+KHR/ITSWE
      SWE(L)=TWE
C
C     SIMULATED AREAL SNOW COVER
CVK  CHANGES
CVK  160 IF(LCOVER.EQ.0) GO TO 170
CVK      IF((KHR/ITSSC)*ITSSC.NE.KHR) GO TO 170
  160 IF(LCOVER.EQ.0) GO TO 217
      IF((KHR/ITSSC)*ITSSC.NE.KHR) GO TO 217
            
      L=KOFF*NSSC+KHR/ITSSC
      COVER(L)=PCOVER

CVK  NEW SIMULATED SNOW DEPTH TIME SERIES  ----
  217 IF(LSSDPT .EQ. 0) GOTO 170
      IF((KHR/ITSN)*ITSN .NE. KHR) GOTO 170
      L=KOFF*NSSN+KHR/ITSN
      SNOWD(L)=SNDPT
CVK ------------------------------------------

C.......................................
C     PRINT RESULTS IF REQUESTED.
  170 IF(IPRINT.EQ.0) GO TO 175
      CALL PRSN19(KDA,KHR,IPRINT,ITITLE,PS(2),TWE,CWE,PTWE,POWE,PCOVER,
     1   CAESC,POSC,LAEC,MOPR,IDPR,IYPR,IHPR)
C.......................................
C     STORE FFG CARRYOVER IF REQUESTED.
  175 IF (IFILLC.EQ.0) GO TO 180
      IF (MAINUM.NE.1) GO TO 180
      IF (IFFG.EQ.0) GO TO 180
      IF ((KDA.EQ.KDAFFG).AND.(KHR.EQ.KHRFFG)) GO TO 320
      GO TO 180
  320 CALL CSAV19(CS)
C.......................................
C     STORE CARRYOVER IF REQUESTED.
  180 IF(IFILLC.EQ.0) GO TO 190
      IF(NCSTOR.EQ.0) GO TO 190
      IF(IC.GT.NCSTOR) GO TO 190
      IF((KDA.EQ.ICDAY(IC)).AND.(KHR.EQ.ICHOUR(IC))) GO TO 185
      GO TO 190
C
C     SAVE CARRYOVER
  185 CALL CSAV19(CT)

CVK  two new states added
CVK      NCO=10+NEXLAG
      NCO=12+NEXLAG
      
      CALL FCWTCO(KDA,KHR,CT,NCO)
      IC=IC+1
C.......................................
C     CHECK FOR END OF EXECUTION PERIOD
  190 IF((KDA.EQ.LDA).AND.(KHR.EQ.LHR)) GO TO 195
C
C     INCREMENT TIMES FOR NEXT COMPUTATIONAL PERIOD.
      KHR=KHR+IDT
      IF(KHR.LE.24) GO TO 191
      KHR=IDT
      KDA=KDA+1
      IF (NOSNOW.EQ.0) IDN=IDN+1
      DSFALL=0.0
      DRAIN=0.0
      DQNET=0.0
      DRSL=0.0
      NDRSP=0
  191 IF(ITITLE.EQ.0) GO TO 120
      IF(LAST.GT.0) GO TO 192
      LAST=LASTDA(MOPR)
      IF((MOPR.EQ.2).AND.((IYPR/4)*4.EQ.IYPR)) LAST=LAST+1
  192 IHPR=IHPR+IDT
      IF(IHPR.LE.24) GO TO 120
      IHPR=IHPR-24
      IDPR=IDPR+1
      IF(IDPR.LE.LAST) GO TO 120
      IDPR=1
      MOPR=MOPR+1
      IF(MOPR.LE.12) GO TO 193
      MOPR=1
      IYPR=IYPR+1
  193 LAST=LASTDA(MOPR)
      IF((MOPR.EQ.2).AND.((IYPR/4)*4.EQ.IYPR)) LAST=LAST+1
      GO TO 120
C.......................................
C     END OF EXECUTION PERIOD.
C     COMPUTE WATER BALANCE FOR THE EXECUTION PERIOD.
  195 SBAL=SPX-SRM-(TWE-TWE1)+CHGWE
      IF(ABS(SBAL).LE.1.0) GO TO 196
      WRITE(IPR,911) SBAL
  911 FORMAT(1H0,47H**WARNING** SNOW BALANCE RESIDUAL EXCEEDS 1 MM., 3X,
     19HRESIDUAL=,F7.2)
      CALL WARN
C
C     STORE SUMS IF REQUESTED.
  196 IF(LSUMS.EQ.0) GO TO 197
      PS(LSUMS)=SPX
      PS(LSUMS+1)=SSFALL
      PS(LSUMS+2)=SRM
      PS(LSUMS+3)=SMELT
      PS(LSUMS+4)=SMELTR
      PS(LSUMS+5)=SROBG
      PS(LSUMS+6)=SBAL
C
C     STORE CARRYOVER IF REQUESTED.
  197 IF(IFILLC.EQ.0) GO TO 198
      IF (MAINUM.NE.1) GO TO 330
      IF (IFFG.NE.0) GO TO 198
  330 CALL CSAV19(CS)
C.......................................
C     PRINT DEBUG OUTPUT IF REQUESTED.
  198 RETURN
      END