      PROGRAM  MAIN_SAC
      
C     LDAS Driver for NWSRFS SAC-SMA model
C     Qingyun Duan, Dag Lohmann Feb. 2000
C     Youlong Xia, modified on June 2013   
C     Youlong Xia, modified on October 2013; Using Noah grib2 output 
 
      IMPLICIT NONE
      
      INCLUDE 'main_sac.h'         
     
C     DRIVER STEP 1 
C     READ CONTROL FILE
      WRITE(*,*) 'READ CONTROLFILE'
      CNTRFL = 'cntl_sac'      
      CALL READCNTL(CNTRFL,ICE,DT,Z,SNOALB,INOAHETP)
      
C     DRIVER STEP 2
C     READ MODEL GEOMETRY / MODEL PARAMETERS
      WRITE(*,*) 'READ MODEL GEOMETRY'
      CALL READ_PAR(NX,NY,NMONTH,SOILTYP,LAND_SEA,VEGTYP,
     &              SLOPETYP,SOILDEPTH,ALBEDO,SHDFAC,ELEV_2d,
C                SAC Parameters
     &              PECLIM,PERATIO,PEADJ,
     &              UZTWM_2d,UZFWM_2d,UZK_2d,PCTIM_2d,ADIMP_2d,
     &              RIVA_2d,ZPERC_2d,REXP_2d,
     &              LZTWM_2d,LZFSM_2d,LZFPM_2d,LZSK_2d,LZPK_2d,
     &              PFREE_2d,SIDE_2d,RSERV_2d,
C                SNOW17 Parameters
     &              SCF_2d,MFMAX_2d,MFMIN_2d,UADJ_2d,SI_2d,
     &              NMF_2d,TIPM_2d,MBASE_2d,PXTEMP_2d,
     &              PLWHC_2d,DAYGM_2d,ADC_2d,NSOLD,PRCP_MASK)

C     DEFINE THE OUTPUT MASK FOR GRIB FILES.
      
      DO J = 1,NY
         DO I = 1,NX
            IF (LAND_SEA(I,J) .EQ. 0) THEN
               LDASMASK(I,J) = .FALSE.
            ELSE
               LDASMASK(I,J) = .TRUE.
            END IF
         END DO
      END DO

C     DRIVER STEP 3
C     READ INITIAL CONDITIONS
      WRITE(*,*) 'READ INITIAL CONDITIONS'
      CALL READ_INITIAL(NX,NY,NSOLD,SMC,SNOWCO)
C      CHECK INITIAL CONDITIONS AND MODEL GEOMETRY
      CALL CHECK_INITIAL(NX,NY,NSOLD,LAND_SEA,SMC,
     &     SOILTYP,VEGTYP,SLOPETYP,ALBEDO,SHDFAC,
     &     NSOIL,NMONTH,SOILDEPTH)

C     DRIVER STEP 4
C      READ FORCING DATA FOR THE DAY
      write(*,*) ' READ FORCING DATA'
      DO I = 1, NHOUR
         WRITE(HOUR_CH,'(I2.2)') I
CYX         FILENAME = 'fort.'//HOUR_CH//char(0)
        ENVVAR = 'FORT'//HOUR_CH//char(0)
        CALL GET_ENVIRONMENT_VARIABLE(ENVVAR,FILENAME)
         WRITE(*,*) 'reading atmospheric data ', FILENAME
         CALL READ_FORCING(NLDAS,F,TAIR(1,1,I),SPFH(1,1,I),PSFC(1,1,I),
     &        UWIND(1,1,I),VWIND(1,1,I),LWDN(1,1,I),EDASPREC,
     &        CAPE,PEVAP,PRCP(1,1,I),SOLDN(1,1,I),LB,FILENAME,
     &        I)
    
      END DO
      CALL CHECK_FORCING(PRCP,NX,NY,NHOUR,LAND_SEA)

C     READ NOAH POTENTIAL EVAPOTRANSPIRATION

      IF (INOAHETP .EQ. 1) THEN
         write(*,*) 'Read NOAH PE'
         DO I = 1, NHOUR
            WRITE(HOUR_CH,'(I2.2)') I+25
CYX            FILENAME = 'fort.'//HOUR_CH//char(0)
            ENVVAR2 = 'FORT'//HOUR_CH//char(0)
            CALL GET_ENVIRONMENT_VARIABLE(ENVVAR2,FILENAME) 
CYX            CALL READ_NOAH_GRIB(nldas,NOAHETP(1,1,I),FILENAME,lb,I+24)
        CALL READ_NOAH_GRIB2(nldas,NOAHETP(1,1,I),FILENAME,I+24)
         END DO
      END IF
C     DRIVER STEP 5
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
            PECLIM_D(I,J) = W1 * PECLIM(M1,I,J) + W2 * PECLIM(M2,I,J)
            PEADJ_D(I,J)  = W1 * PEADJ(M1,I,J)  + W2 * PEADJ(M2,I,J)
         END DO
      END DO
      
      OPEN(30, file = 'yesterday')
      READ(30,'(A)') YESTERDAY
      CLOSE(30)
      CYR = YESTERDAY(1:4)
      CMO = YESTERDAY(5:6)
      CDA = YESTERDAY(7:8)
      read (CYR(1:4),'(i4)') IYEAR
      read (CMO(1:2),'(i2)') IMON
      read (CDA(1:2),'(i2)') IDAY
	
      OPEN(30, file = 'today')
      READ(30,'(A)') TODAY
      CLOSE(30)

C     DRIVER STEP 5
C     RUN TIME AND SPATIAL LOOP

C      open(98,file='test.out',status='unknown')

      DO NT = 2, NHOUR
C        WRITE(*,*) NT,YESTERDAY
         DO J = 1,NY
            DO I = 1,NX
C         DO J = 80,80
C            DO I = 208,208
C               print*, I, J, LAND_SEA(I,J)
               IF (LAND_SEA(I,J) .EQ. 1) THEN
                  DT_TAIR  = TAIR(I,J,NT)
                  DT_SPFH  = SPFH(I,J,NT)
                  DT_PSFC  = PSFC(I,J,NT)
                  DT_UWIND = UWIND(I,J,NT)
                  DT_VWIND = VWIND(I,J,NT)
                  DT_LWDN  = LWDN(I,J,NT)
                  DT_SOLDN = SOLDN(I,J,NT)
                  DT_PRCP  = PRCP(I,J,NT)
                     
C     DRIVER STEP 6 
C     CALCULATE CH (EXCHANGE COEFFICIENT)
                     
                  CALL QDATAP(DT_TAIR,DT_PSFC,DT_SPFH_SAT)
                     
                  IF (DT_SPFH .GE. DT_SPFH_SAT) THEN
                     DT_SPFH = DT_SPFH_SAT
                  END IF
                  RH = DT_SPFH / DT_SPFH_SAT 

C     ASSIGN POTENTIAL EVAPORATION VALUE

                  IF (INOAHETP .EQ. 1) THEN
                     EDMND = NOAHETP(I,J,NT) / (2.501E+6/3600.0)
C  ADJUST PE BY THE SAC_NOAH PE RATIO (ADDED on 10/23/2003):
                     EDMND = EDMND * PERATIO(I,J)
                     write(98,9) YESTERDAY,NT,I,J,EDMND,NOAHETP(I,J,NT)
    9                format(A8,3I4,2f10.4)
                  END IF
                  IF (INOAHETP .EQ. 2) THEN
                     XLAT = 25.0625 + (NY-1) * 0.125
                     IT = NT + 6
                     IF (IT .GT. 24) IT = IT - 24
                     CALL ZENANGL(24,NT,JULDAY,1.,XLAT,
     &                            ZEN,ZENAVG,DAYLIGHT)
                     EDMND = PECLIM_D(I,J) / DAYLIGHT * (ZEN / ZENAVG)

C                     print*,'PECLIM_D(I,J),EDMND,ZEN,ZENAVG,DAYLITE',
C     &               I,J,PECLIM_D(I,J),EDMND,ZEN,ZENAVG,DAYLIGHT

                  END IF
                  IF (INOAHETP .EQ. 0) THEN
                     CALL PENMAN(DT_TAIR,DT_PSFC,DT_SOLDN,DT_LWDN,
     &                           DT_UWIND,DT_VWIND,RH,ALBEDO_D(I,J),
     &                           Z,VEGTYP(I,J),ET_PENMAN)
                     EDMND = ET_PENMAN * DT
C                     print*, I, J, 'EDMND = ', EDMND 
                  END IF

C     ADJUST EDMND BY PE ADJUSTMENT FACT - PEADJ_D
                  EDMND = EDMND * PEADJ_D(I,J)

C                  print*, 'INPUT 1', NT, 
C     &            DT_TAIR, DT_SPFH, DT_PSFC,DT_UWIND,
C     &            DT_VWIND, DT_LWDN, DT_SOLDN, DT_PRCP,
C     &            NOAHETP(I,J,NT),EDMND

C     INITIALIZE PARAMETERS FOR EACH GRID CELL

                  UZTWM = UZTWM_2d(I,J)
                  UZFWM = UZFWM_2d(I,J)
                  UZK   = UZK_2d(I,J) 
                  PCTIM = PCTIM_2d(I,J)
                  ADIMP = ADIMP_2d(I,J)
                  RIVA  = RIVA_2d(I,J)
                  ZPERC = ZPERC_2d(I,J)
                  REXP  = REXP_2d(I,J)
                  LZTWM = LZTWM_2d(I,J)
                  LZFSM = LZFSM_2d(I,J)
                  LZFPM = LZFPM_2d(I,J)
                  LZSK  = LZSK_2d(I,J)
                  LZPK  = LZPK_2d(I,J)
                  PFREE = PFREE_2d(I,J)
                  SIDE  = SIDE_2d(I,J)
                  RSERV = RSERV_2d(I,J)

                  ALAT  = 25.0625 + (NY-1) * 0.125
                  SCF   = SCF_2d(I,J)
                  MFMAX = MFMAX_2d(I,J)
                  MFMIN = MFMIN_2d(I,J)
                  UADJ  = UADJ_2d(I,J)
                  SI    = SI_2d(I,J)
                  NMF   = NMF_2d(I,J)
                  TIPM  = TIPM_2d(I,J)
                  MBASE = MBASE_2d(I,J)
                  PXTEMP= PXTEMP_2d(I,J)
                  PLWHC = PLWHC_2d(I,J)
                  DAYGM = DAYGM_2d(I,J)
                  PA    = PSFC(I,J,NT)/100.
                  ELEV  = ELEV_2d(I,J)
                  DO K = 1, 11
                     ADC(K) = ADC_2d(K,I,J)
                  END DO


      PAREA = 1. - PCTIM - ADIMP 
                  
C     DRIVER STEP 7 

C     INITIALIZE SNOW17
                  DO K = 1, 19
                     CS(K) = SNOWCO(K,I,J)
                  END DO

C     CALL SNOW19
                  IDTS = INT(DT)
                  IDT = INT(DT/3600)
                  TMP = DT_TAIR - 273.15
                  ISNOW = 0
                  IF (NT .EQ. 1) TPREV(I,J) = TMP
                  IF (DT_PRCP .GT. 0.0 .AND. TMP .LT. PXTEMP) ISNOW = 1
                  IF (CS(18) .GT. 0.0) ISNOW = 1
                  IF (ISNOW .EQ. 1) THEN
                     CALL EXSNOW19(IDTS,IDT,IDAY,IMON,IYEAR,
C     SNOW17 INPUT/OUTPUT AND 
     &                    DT_PRCP,TMP,RAIM,SNEQV(I,J),SNOW,SNOWH(I,J),
C     SNOW17 PARAMETERS
     &                    ALAT,SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,
     &                    PXTEMP,PLWHC,DAYGM,ELEV,PA,ADC(1),
C     SNOW17 CARRYOVER VARIABLES
     &                    CS(1),TPREV(I,J))

C                     print*, I,J,NT,IDAY,IMON,IYEAR,DT_PRCP,
C     &                    TMP,RAIM,SNEQV(I,J),SNOW,SNOWH(I,J),CS
C     SNOW17 PARAMETERS
C                     print*,'SNOW19 Parameters - ',
C     &                    ALAT,SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,
C     &                    PXTEMP,PLWHC,DAYGM,ELEV,PA,ADC
                  END IF
C     UPDATE TPREV(I,J)
                  TPREV(I,J) = TMP

C     UPDATE SNOW17 CARRYOVER VARIABLES
                  DO K = 1, 19
                     SNOWCO(K,I,J) = CS(K)
                  END DO

C     CALL SAC-SMA
                  IF (TMP .LE. 0.0) THEN
                     LIQUID_PREC(I,J) = 0.0
                     FROZEN_PREC(I,J) = DT_PRCP
                  ELSE
                     LIQUID_PREC(I,J) = DT_PRCP
                     FROZEN_PREC(I,J) = 0.0
                  END IF
                  IF (ISNOW .EQ. 1) THEN
                     SNOWMELT(I,J) = RAIM - LIQUID_PREC(I,J)
                     RAIM = RAIM
                  ELSE
                     SNOWMELT(I,J) = 0.0
                     RAIM = DT_PRCP
                  END IF
                  
C     INITIALIZE SAC Model
                  
                  UZTWC = SMC(1,I,J)
                  UZFWC = SMC(2,I,J)
                  LZTWC = SMC(3,I,J)
                  LZFSC = SMC(4,I,J)
                  LZFPC = SMC(5,I,J)
                  ADIMC = SMC(6,I,J)
                  
C                  print*, 'WRITE BEFORE MODEL CALL'
C                  print*, YESTERDAY,NT,DT_PRCP,RAIM,TMP,EDMND,
C     &                    'SAC STATES AND OUTPUTS  ',
C     &                    QS,QG,EVP,
C     &                    'SAC Parameters  ',
C     &                    UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,
C     &                    REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,
C     &                    SIDE,RSERV,
C     &                    'SAC State variables  ',
C     &                    UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC

                  IF (EDMND .LT. 0.0) EDMND = 0.0
                  CALL EXSAC(NSOLD,DT,RAIM,TMP,EDMND,
C     SAC PARAMETERS
     &                       UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,
     &                       REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,
     &                       SIDE,RSERV,
C     SAC State variables  ',
     &                       UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC,
C     SAC OUTPUTS
     &                       QS,QG,Q,EVP)

C                  print*, 'WRITE AFTER MODEL CALL'
C                  print*, YESTERDAY,NT,DT_PRCP,TMP,EDMND,
C     &                    'SAC STATES AND OUTPUTS  ',
C     &                    QS,QG,EVP,
C     &                    'SAC Parameters  ',
C     &                    UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,
C     &                    REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,
C     &                    SIDE,RSERV,
C     &                    'SAC State variables  ',
C     &                    UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC

                  ETA(I,J) = EVP
                  ETP(I,J) = EDMND * (2.501E+6/3600.0) 
                  RUNOFF1(I,J) = QS
                  RUNOFF2(I,J) = QG
                  
                  SMC(1,I,J) = UZTWC
                  SMC(2,I,J) = UZFWC
                  SMC(3,I,J) = LZTWC
                  SMC(4,I,J) = LZFSC
                  SMC(5,I,J) = LZFPC
                  SMC(6,I,J) = ADIMC
                  

                  SOILM(I,J) = (UZTWC + UZFWC + LZTWC + LZFSC + LZFPC) * PAREA + ADIMC * ADIMP 
C this is a bug
C                  SOILM(I,J) = 0.0
C                  DO K = 1,6
C                     SOILM(I,J) = SOILM(I,J) +  SMC(K,I,J)
C                  END DO
                  
C                  write(98,100) yesterday,NT,DT_PRCP,TMP,EDMND,EVP,QS,
C     &                 QG,(SMC(K,I,J),K=1,6),ISNOW,SNOW,SNOWMELT(I,J),
C     &                 (CS(K),K=1,19)
C100               format(A8,i2.2,6f8.4,6f8.2,i2,21f8.2)

               END IF
            END DO
         END DO

         SH2O = SMC
         
         LUGB = NHOUR + 1

         CALL GRIBOUT(NX,NY,LDASMASK,NSOLD,NT-1,YESTERDAY,TODAY,
     &        LIQUID_PREC, FROZEN_PREC, ETA, ETP,
     &        RUNOFF1,RUNOFF2, SNOWMELT, SNEQV, SNOWH, SMC, 
     &        SH2O,SOILM,DUMMY,LUGB)
         
      END DO
      
      WRITE(*,*) 'writing restart files'
      WRITE(77) (((SMC(K,I,J),K=1,NSOLD),I=1,NX),J=1,NY)
      WRITE(78) (((SNOWCO(K,I,J),K=1,19),I=1,NX),J=1,NY)
      
      END

