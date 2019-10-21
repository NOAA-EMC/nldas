      PROGRAM  MAIN     

      IMPLICIT NONE
      
      INCLUDE 'MAIN.h'

      CHARACTER*1 TRP_FLAG
      CHARACTER*11 ENVVAR
      INTEGER IP(NX,NY),MADTT(NX,NY)

      REAL val0(NX,NY)   ! last forcing value (for cases N & C)
      REAL val1(NX,NY)   ! current forcing value (all cases)
      REAL val2(NX,NY)   ! next forcing value (all cases but default)
      REAL val3(NX,NY)   ! ueber-next forcing value (for cases C & L)
      REAL BTIME(NX,NY), ETIME(NX,NY),MTIME(NX,NY), IFLAG(NX,NY)
      REAL WEIGHT1(NX,NY), WEIGHT2(NX,NY)
C SOLAR ZENITH ANGLE

      REAL GMT(NX,NY)
      REAL LHOUR(NX,NY)
      REAL CZMODEL(NX,NY)
      INTEGER ZONE(NX,NY)
      REAL LSTSNW1(NX,NY)	
C --------------------------------------------------------------------
      
C     READ CONTROL FILE
      WRITE(*,*) 'READ CONTROLFILE'
      CNTRFL = 'controlfile'      
      CALL READCNTL(CNTRFL,ICE1,DT1,ZLVL1)

      DO J=1,NY
      DO I=1,NX
      LAT(I,J)=25.0625+(J-1)*.125
      LON(I,J)=235.0625+(I-1)*.125
      IF(LON(I,J).GT.180)LON(I,J)=LON(I,J)-360.
      ICE(I,J)=ICE1
      DT(I,J)=DT1
      ZLVL(I,J)=ZLVL1
      IP(I,J)=0
      MADTT(I,J)=4
      IFLAG(I,J)=1
      ENDDO
      ENDDO
                                            
C     READ MODEL GEOMETRY / LAND SURFACE CHARACTERISTICS
      WRITE(*,*) 'READ MODEL GEOMETRY / LSC'
      CALL READ_LSC(NX,NY,NSOLD,NMONTH,NSOIL,SOILTYP,
     &     LAND_SEA,VEGTYP,SLOPETYP,SOILDEPTH,ALBEDO,SHDFAC,
     &     TBOT,MAXSNOWALB,PRCP_MASK)

      
C     CHECK MODEL GEOMETRY / LAND SURFACE CHARACTERISTICS
      WRITE(*,*) 'CHECK MODEL GEOMETRY / LSC'
      CALL CHECK_LSC(NX,NY,NSOLD,NMONTH,NSOIL,SOILTYP,
     &     LAND_SEA,VEGTYP,SLOPETYP,SOILDEPTH,ALBEDO,SHDFAC,
     &     TBOT,MAXSNOWALB,PRCP_MASK)

C     DEFINE THE OUTPUT MASK FOR GRIB FILES.
      LDASMASK = .TRUE.
      WHERE (LAND_SEA .EQ. 0)
         LDASMASK = .FALSE.
      END WHERE

C     CREATE MAXIMUN/MINIMAL GREENNESS DATA
      SHDMAX = -999.0
      SHDMIN = 999.0
      DO J = 1,NY
         DO I = 1,NX
            IF (LAND_SEA(I,J) .EQ. 1) THEN
               DO K = 1, NMONTH
                  SHDMAX(I,J) = max(SHDMAX(I,J),SHDFAC(K,I,J))
                  SHDMIN(I,J) = min(SHDMIN(I,J),SHDFAC(K,I,J))
               END DO
            END IF
         END DO
      END DO

C     READ INITIAL CONDITIONS
      WRITE(*,*) 'READ INITIAL CONDITIONS'
      CALL READ_INITIAL(NX,NY,NSOLD,T1,STC,SMC,SH2O,CMC,
     &                  SNOWH,SNEQV,CH,CM,LSTSNW1)

C     READ FORCING DATA FOR THE DAY AND CHECK PRECIPITATION
      WRITE(*,*) 'READ FORCING'
C     DO I = 1, (NHOUR-1)
      DO I = 1, NHOUR
         WRITE(HOUR_CH,'(I2.2)') I
CJZ         FILENAME = 'fort.'//HOUR_CH//char(0)
         ENVVAR = 'FORT'//HOUR_CH//char(0)
         CALL GET_ENVIRONMENT_VARIABLE(ENVVAR,FILENAME)
         WRITE(*,*) 'reading atmospheric data ', FILENAME
         WRITE(*,*) 'reading atmospheric data on unit', I
  
         CALL READ_FORCING(NLDAS,F,TAIR(1,1,I),SPFH(1,1,I),PSFC(1,1,I),
     &        UWIND(1,1,I),VWIND(1,1,I),LWDN(1,1,I),EDASPREC,
     &        CAPE,CONVPREC,PRCP(1,1,I),SOLDN(1,1,I),LB,FILENAME,I)
      END DO

      CALL CHECK_FORCING(PRCP,NX,NY,NHOUR,LAND_SEA)

C     READ JULIAN DAY, CALCULATE WEIGTHS AND INTERPOLATE FROM 
C     MONTHLY TO DAILY
      OPEN(30,FILE = 'julday')
      READ(30,*) JULDAY
      CLOSE(30)
      IF (JULDAY .EQ. 0) JULDAY = 365
      CALL CALC_WEIGTHS(JULDAY,W1,W2,M1,M2)
      DO J = 1,NY
         DO I = 1,NX
            ALBEDO_D(I,J) = W1 * ALBEDO(M1,I,J) + W2 * ALBEDO(M2,I,J)
            SHDFAC_D(I,J) = W1 * SHDFAC(M1,I,J) + W2 * SHDFAC(M2,I,J)
         END DO
      END DO

      OPEN(30, file = 'yesterday')
      READ(30,'(A)') YESTERDAY
      CLOSE(30)
      
      OPEN(30, file = 'today')
      READ(30,'(A)') TODAY
      CLOSE(30)


      TRP_FLAG='I'


C     RUN TIME AND SPATIAL LOOP

      DO NT = 1,(NHOUR-1)

         HOUR_NSWRS    = 0.0
         HOUR_LWDN_NET = 0.0
         HOUR_ETA      = 0.0
         HOUR_H        = 0.0
         HOUR_S        = 0.0
         DelColdCont   = 0.0
         HOUR_DSWRF    = 0.0
         HOUR_LWDN     = 0.0

         HOUR_ASNOW    = 0.0
         HOUR_ARAIN    = 0.0
         HOUR_RUNOFF1  = 0.0
         HOUR_RUNOFF2  = 0.0
         HOUR_SNOMLT   = 0.0
         HOUR_CMC      = 0.0

         HOUR_EVP      = 0.0
         HOUR_ETP      = 0.0
         HOUR_EDIR     = 0.0
         HOUR_EC       = 0.0
         HOUR_ETT      = 0.0
         HOUR_ESNOW    = 0.0
         HOUR_ACOND    = 0.0
         HOUR_RC       = 0.0
         HOUR_RCS      = 0.0
         HOUR_RCT      = 0.0
         HOUR_RCQ      = 0.0
         HOUR_RCSOIL   = 0.0
         HOUR_FX       = 0.0
         HOUR_RSMIN    = 0.0
         HOUR_LAI      = 0.0
         HOUR_GVF      = 0.0

         SOILROOT      = 0.0
         SOIL1M        = 0.0
         MSTAV         = 0.0

         DO SUB_DT = 1,4
            WRITE(*,*) NT, SUB_DT
C$OMP PARALLEL DO
            DO J = 1,NY
               DO I = 1,NX
                  IF (LAND_SEA(I,J) .EQ. 1) THEN

c               GMT(I,J)=(SUB_DT-1)*DT(I,J)/3600.+NT-1
                GMT(I,J)=SUB_DT*DT(I,J)/3600.+NT-1
                IF(GMT(I,J).GE.24)GMT(I,J)=GMT(I,J)-24

            CALL LOCALTIME (GMT(I,J),LON(I,J),LHOUR(I,J),ZONE(I,J))
            CALL COSZENITH(LON(I,J),LAT(I,J),LHOUR(I,J),ZONE(I,J),
     &                     JULDAY,CZMODEL(I,J))

C     Temporal Interpolation

         CALL FINTERP(IP(I,J),TRP_FLAG,VAL0(I,J),TAIR(I,J,NT),
     &        TAIR(I,J,NT+1),VAL3(I,J),MADTT(I,J),SUB_DT,DT_TAIR(I,J))
                     
         CALL FINTERP(IP(I,J),TRP_FLAG,VAL0(I,J),SPFH(I,J,NT),
     &        SPFH(I,J,NT+1),VAL3(I,J),MADTT(I,J),SUB_DT,DT_SPFH(I,J))

         CALL FINTERP(IP(I,J),TRP_FLAG,VAL0(I,J),PSFC(I,J,NT),
     &        PSFC(I,J,NT+1),VAL3(I,J),MADTT(I,J),SUB_DT,DT_PSFC(I,J))

         CALL FINTERP(IP(I,J),TRP_FLAG,VAL0(I,J),UWIND(I,J,NT),
     &        UWIND(I,J,NT+1),VAL3(I,J),MADTT(I,J),SUB_DT,DT_UWIND(I,J))

         CALL FINTERP(IP(I,J),TRP_FLAG,VAL0(I,J),VWIND(I,J,NT),
     &        VWIND(I,J,NT+1),VAL3(I,J),MADTT(I,J),SUB_DT,DT_VWIND(I,J))

         CALL FINTERP(IP(I,J),TRP_FLAG,VAL0(I,J),LWDN(I,J,NT),
     &        LWDN(I,J,NT+1),VAL3(I,J),MADTT(I,J),SUB_DT,DT_LWDN(I,J))

         CALL FINTERP(IP(I,J),TRP_FLAG,VAL0(I,J),SOLDN(I,J,NT), 
     &        SOLDN(I,J,NT+1),VAL3(I,J),MADTT(I,J),SUB_DT,DT_SOLDNB(I,J))

          BTIME(I,J)=NT-1
          ETIME(I,J)=NT
          IF(ETIME(I,J).GE.24)ETIME(I,J)=ETIME(I,J)-24
          MTIME(I,J)=GMT(I,J)
          
          CALL ZTERP (IFLAG(I,J),LAT(I,J),LON(I,J),BTIME(I,J),ETIME(I,J),
     +    MTIME(I,J),JULDAY,WEIGHT1(I,J),WEIGHT2(I,J))


       DT_SOLDN(I,J) = SOLDN(I,J,NT)*WEIGHT1(I,J)+ SOLDN(I,J,NT+1)*WEIGHT2(I,J)

C  --------- In cases of small cos(zenith) angles, use linear weighting -----
C  --------- to avoid overly large weights ----------------------------------

          IF(DT_SOLDN(I,J).GT.SOLDN(I,J,NT).AND.DT_SOLDN(I,J).GT.SOLDN(I,J,NT+1))THEN
           IF(CZMODEL(I,J).LT.0.1) THEN  
          DT_SOLDN(I,J) = DT_SOLDNB(I,J)
             ENDIF
             ENDIF
C  ----- Linear interpolation -----------------------------------------------
                    
          DT_PRCP(I,J)  = PRCP(I,J,NT+1)*W_PRECIP(SUB_DT) / DT(I,J)

C     CALCULATE specific humidities
                     CALL QDATAP(DT_TAIR(I,J),DT_PSFC(I,J),DT_SPFH_SAT(I,J))
                     
                     IF (DT_SPFH(I,J) .GE. DT_SPFH_SAT(I,J)) THEN
                        DT_SPFH(I,J) = DT_SPFH_SAT(I,J)
                     END IF
                     
C     CALCULATE SLOPE OF SAT SPECIFIC HUMIDITY CURVE FOR PENMAN: DQSDT2
                     DQSDT2(I,J) = DQSDT(DT_TAIR(I,J), DT_PSFC(I,J))

                     TH2(I,J) = DT_TAIR(I,J) + 0.0098 * ZLVL(I,J)
                     SFCSPD(I,J) = SQRT(DT_UWIND(I,J)*DT_UWIND(I,J)+
     D                             DT_VWIND(I,J)*DT_VWIND(I,J))
                     CALL SFLX (
     C                    ICE(I,J),DT(I,J),ZLVL(I,J),NSOIL(I,J),SOILDEPTH(1,I,J),
     F                    DT_LWDN(I,J),DT_SOLDN(I,J),SOLNET,DT_PSFC(I,J),DT_PRCP(I,J),
     F                    DT_TAIR(I,J),DT_SPFH(I,J),SFCSPD(I,J),
     I                    TH2(I,J),DT_SPFH_SAT(I,J),DQSDT2(I,J),
     S                    VEGTYP(I,J),SOILTYP(I,J),SLOPETYP(I,J),
     S                    SHDFAC_D(I,J),SHDMAX(I,J),SHDMIN(I,J),
     S                    PTU(I,J),ALBEDO_D(I,J),MAXSNOWALB(I,J),TBOT(I,J),
     H                    CMC(I,J),T1(I,J),STC(1,I,J),SMC(1,I,J),
     H                    SH2O(1,I,J),SNOWH(I,J),SNEQV(I,J),
     H                    ALBSNO(I,J),CH(I,J),CM(I,J),
     O                    ETA(I,J),H(I,J),
     O                    EC(I,J),EDIR(I,J),ET(1,I,J),ETT(I,J),
     O                    ESNOW(I,J),DRIP(I,J),DEW(I,J),
     O                    BETA(I,J),ETP(I,J),S(I,J),
     O                    FLX1(I,J),FLX2(I,J),FLX3(I,J),
     O                    SNOMLT(I,J),SNCOVR(I,J),
     O                    RUNOFF1(I,J),RUNOFF2(I,J),RUNOFF3(I,J),
     O                    RC(I,J),PC(I,J),RSMIN(I,J),XLAI(I,J),RCS(I,J),
     O                    RCT(I,J),RCQ(I,J),RCSOIL(I,J),FX(I,J),
     D                    SOILW(I,J),SOILM(I,J),RGL(I,J),HS(I,J),
     P                    SMCWLT(I,J),SMCDRY(I,J),SMCREF(I,J),SMCMAX(I,J),NROOT(I,J),
     P                    CZMODEL(I,J),LSTSNW1(I,J))

C     FORCING VARIABLES
C     Downward shortwave
                     HOUR_DSWRF(I,J) = HOUR_DSWRF(I,J) +
     &                    DT_SOLDN(I,J) / 4.0
C     Downward longwave
                     HOUR_LWDN(I,J) = HOUR_LWDN(I,J)+DT_LWDN(I,J)/4.0

C     Energy balance output variables
C     Net shortwave
                     HOUR_NSWRS(I,J)= HOUR_NSWRS(I,J) - (1.0-
     &                    ALBSNO(I,J)) * DT_SOLDN(I,J) / 4.0
C     Net longwave (when snow emissivity = 0.9)
C                     HOUR_LWDN_NET(I,J) = HOUR_LWDN_NET(I,J) + 
C     &                    ((1-SNCOVR(I,J)*0.1)*5.67E-8*T1(I,J)**4.0/4.0 - DT_LWDN(I,J)/4.0)


C     Net longwave (when snow emissivity = 1.0)
                     HOUR_LWDN_NET(I,J) = HOUR_LWDN_NET(I,J) + 
     &                    ((1-SNCOVR(I,J)*0.0)*5.67E-8*T1(I,J)**4.0/4.0 - DT_LWDN(I,J)/4.0)


C     Actual Evapotranspiration
                     HOUR_ETA(I,J) = HOUR_ETA(I,J) + ETA(I,J) / 4.0
C     Sensible heat flux
                     HOUR_H(I,J)   = HOUR_H(I,J)   + H(I,J) / 4.0
C     Ground heat flux
                     HOUR_S(I,J)   = HOUR_S(I,J) + S(I,J) / 4.0
C     Snow phase change
                     DelColdCont(I,J) = DelColdCont(I,J) -
     &                    (FLX1(I,J)+FLX2(I,J)+FLX3(I,J)) / 4.0

C     Water balance output variables
                     IF (DT_TAIR(I,J) .LE. TFREEZ) THEN
                        HOUR_ASNOW(I,J) = HOUR_ASNOW(I,J)+DT_PRCP(I,J)*DT(I,J)
                     ELSE
                        HOUR_ARAIN(I,J) = HOUR_ARAIN(I,J)+DT_PRCP(I,J)*DT(I,J)
                     END IF
                     HOUR_RUNOFF1(I,J) = HOUR_RUNOFF1(I,J) + 
     &                    RUNOFF1(I,J) * 1000.0 * DT(I,J)
                     HOUR_RUNOFF2(I,J) = HOUR_RUNOFF2(I,J) + 
     &                    (RUNOFF2(I,J) + RUNOFF3(I,J))*1000.0*DT(I,J)
                     HOUR_SNOMLT(I,J) = HOUR_SNOMLT(I,J) + 
     &                    SNOMLT(I,J) * 1000.0

C     Evaporation Components
CYoulong
                     IF(ETP(I,J).LE.0.0) THEN
                     HOUR_EVP(I,J)= HOUR_EVP(I,J)+ETP(I,J)/
     &               4.0/((1.0-SNCOVR(I,J))*2.501E+6+SNCOVR(I,J)*2.83E+6)
                     ELSE
                     HOUR_EVP(I,J)=HOUR_EVP(I,J)+((EDIR(I,J)+EC(I,J)+ETT(I,J))
     &               /(4.0*2.501E+6)+ESNOW(I,J)/(4.0*2.83E+6))
                     ENDIF
CYend    
                     HOUR_ETP(I,J) = HOUR_ETP(I,J) + ETP(I,J)/4.0
                     HOUR_EDIR(I,J) = HOUR_EDIR(I,J) + EDIR(I,J)/4.0
                     HOUR_EC(I,J)   = HOUR_EC(I,J) + EC(I,J)/4.0
                     HOUR_ETT(I,J)  = HOUR_ETT(I,J) + ETT(I,J)/4.0
                     HOUR_ESNOW(I,J) = HOUR_ESNOW(I,J) + ESNOW(I,J) / 4.0
                     HOUR_ACOND(I,J) = HOUR_ACOND(I,J)+CH(I,J)/4.0
                     HOUR_RC(I,J) = HOUR_RC(I,J)+RC(I,J)/4.0
                     HOUR_RCS(I,J) = HOUR_RCS(I,J)+RCS(I,J)/4.0
                     HOUR_RCT(I,J) = HOUR_RCT(I,J)+RCT(I,J)/4.0
                     HOUR_RCQ(I,J) = HOUR_RCQ(I,J)+RCQ(I,J)/4.0
                     HOUR_RCSOIL(I,J) = HOUR_RCSOIL(I,J)+RCSOIL(I,J)/4.0
                     HOUR_FX(I,J) = HOUR_FX(I,J)+FX(I,J)/4.0
                     HOUR_RSMIN(I,J) = HOUR_RSMIN(I,J)+RSMIN(I,J)/4.0
                     HOUR_LAI(I,J) = HOUR_LAI(I,J)+XLAI(I,J)/4.0
                     HOUR_GVF(I,J) = HOUR_GVF(I,J)+SHDFAC_D(I,J)/4.0
                  END IF
                  
               END DO
            END DO
         END DO



C     Subsurface state vaiables
C     Total column soil moisture 
         DO J = 1,NY
            DO I = 1,NX
               IF (LAND_SEA(I,J) .EQ. 1) THEN
                  DO K = 1, NROOT(I,J)
                     SOILROOT(I,J)=SOILROOT(I,J)+SMC(K,I,J)*
     &                    SOILDEPTH(K,I,J)*1000.0
                  END DO
                  DO K = 1, 3
                     SOIL1M(I,J)=SOIL1M(I,J)+SMC(K,I,J)*
     &                    SOILDEPTH(K,I,J)*1000.0
                  END DO
                  DO K = 1, NSOIL(I,J)
                     MSTAV(I,J) = MSTAV(I,J)+SMCMAX(I,J)*
     &                    SOILDEPTH(K,I,J)*1000.0
                     SOILM_OUT(K,I,J) = SMC(K,I,J)*SOILDEPTH(K,I,J)
     &                                  *1000.0
                     SOILL_OUT(K,I,J) = SH2O(K,I,J)*SOILDEPTH(K,I,J)
     &                                  *1000.0
                  END DO
                  SOILM(I,J) = SOILM(I,J) * 1000.0
                  HOUR_CMC(I,J) = CMC(I,J) * 1000.0
                  MSTAV(I,J) = SOILM(I,J) / MSTAV(I,J)
                  HOUR_EVP(I,J) = HOUR_EVP(I,J)*4.0*DT(I,J)
                  END IF

            END DO
         END DO

         
         
         WRITE(*,*) 'writing data'
         CALL GRIBOUT(NX,NY,LDASMASK,NSOLD,NT,YESTERDAY,TODAY,
     &        HOUR_NSWRS,HOUR_LWDN_NET,HOUR_ETA,HOUR_H,HOUR_S,
     &        DelColdCont,HOUR_DSWRF,HOUR_LWDN,HOUR_ASNOW,HOUR_ARAIN,
     &        HOUR_EVP,HOUR_RUNOFF1,HOUR_RUNOFF2,HOUR_SNOMLT,
     &        T1,ALBSNO,SNEQV,HOUR_CMC,STC,SOILM,SOILROOT,SOIL1M,
     &        SOILM_OUT,SOILL_OUT,MSTAV,SOILW,HOUR_EC,HOUR_ETT,
     &        HOUR_EDIR,HOUR_ESNOW,HOUR_ETP,HOUR_ACOND,SNOWH,SNCOVR,
     &        HOUR_RC,HOUR_RCS,HOUR_RCT,HOUR_RCQ,HOUR_RCSOIL,HOUR_FX,
     &        HOUR_RSMIN,HOUR_LAI,HOUR_GVF,DUMMY)
      END DO


      WRITE(*,*) 'writing restart files'
      WRITE(80) ((T1(I,J),I=1,NX),J=1,NY)
      WRITE(81) (((STC(K,I,J),K=1,NSOLD),I=1,NX),J=1,NY)
      WRITE(82) (((SMC(K,I,J),K=1,NSOLD),I=1,NX),J=1,NY)
      WRITE(83) (((SH2O(K,I,J),K=1,NSOLD),I=1,NX),J=1,NY)
      WRITE(84) ((CMC(I,J),I=1,NX),J=1,NY)
      WRITE(85) ((CH(I,J),I=1,NX),J=1,NY)
      WRITE(86) ((CM(I,J),I=1,NX),J=1,NY)
      WRITE(87) ((SNOWH(I,J),I=1,NX),J=1,NY)
      WRITE(88) ((SNEQV(I,J),I=1,NX),J=1,NY)
      WRITE(90) ((LSTSNW1(I,J),I=1,NX),J=1,NY)
      STOP
      END

