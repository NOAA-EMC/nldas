CCCCCCCCC  SUBROUTINE SFLX AND ASSOCIATED PHYSICS SUBROUTINE  CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PURPOSE:  SUB-DRIVER FOR "NOAH/OSU LSM" FAMILY OF PHYSICS SUBROUTINES
C           FOR A SOIL/VEG/SNOWPACK LAND-SURFACE MODEL TO UPDATE SOIL
C           MOISTURE, SOIL ICE, SOIL TEMPERATURE, SKIN TEMPERATURE,
C           SNOWPACK WATER CONTENT, SNOWDEPTH, AND ALL TERMS
C           OF THE SURFACE ENERGY BALANCE AND SURFACE WATER
C           BALANCE (EXCLUDING INPUT ATMOSPHERIC FORCINGS OF
C           DOWNWARD RADIATION AND PRECIP)
C
C  28 Apr 2002: Kristi Arsenault; Added NOAH LSM 2.5 code to LDAS
C  15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new 
C                GRIB libraries
C  28 May 2002: Kristi Arsenault; Updated NOAH code to 2.5.1, and
C               corrected problem with FX in DEVAP function.
C  12 Jun 2002: Kristi Arsenault; Updated NOAH code to 2.5.2 
C  04 Nov 2002: Kristi Arsenault; Incorporated new TBOT fields
C  24 Jun 2003: Kristi Arsenault; Updated Noah LSM to 2.6 version
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC THE FOLLOWING SUBROUTINES ARE IN ALPHABETICAL ORDER CCCCCCCCCC

CCCCC 1. PHYSICS SUBROUTINE ==>  SUBROUTINE ALCALC  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ALCALC (ALB,SNOALB,SHDFAC,SHDMIN,SNCOVR,TSNOW,ALBEDO)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C CALCULATE ALBEDO INCLUDING SNOW EFFECT (0 -> 1)
C   ALB     SNOWFREE ALBEDO
C   SNOALB  MAXIMUM (DEEP) SNOW ALBEDO
C   SHDFAC    AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C   SHDMIN    MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C   SNCOVR  FRACTIONAL SNOW COVER
C   ALBEDO  SURFACE ALBEDO INCLUDING SNOW EFFECT
C   TSNOW   SNOW SURFACE TEMPERATURE (K)
C ----------------------------------------------------------------------
      REAL ALB, SNOALB, SHDFAC, SHDMIN, SNCOVR, ALBEDO, TSNOW
      
C ----------------------------------------------------------------------
C SNOALB IS ARGUMENT REPRESENTING MAXIMUM ALBEDO OVER DEEP SNOW,
C AS PASSED INTO SFLX, AND ADAPTED FROM THE SATELLITE-BASED MAXIMUM 
C SNOW ALBEDO FIELDS PROVIDED BY D. ROBINSON AND G. KUKLA 
C (1985, JCAM, VOL 24, 402-411)
C ----------------------------------------------------------------------
C         changed in version 2.6 on June 2nd 2003
C          ALBEDO = ALB + (1.0-(SHDFAC-SHDMIN))*SNCOVR*(SNOALB-ALB) 
          ALBEDO = ALB + SNCOVR*(SNOALB-ALB)
          IF (ALBEDO .GT. SNOALB) ALBEDO=SNOALB

C     BASE FORMULATION (DICKINSON ET AL., 1986, COGLEY ET AL., 1990)
C          IF (TSNOW.LE.263.16) THEN
C            ALBEDO=SNOALB
C          ELSE
C            IF (TSNOW.LT.273.16) THEN
C              TM=0.1*(TSNOW-263.16)
C              ALBEDO=0.5*((0.9-0.2*(TM**3))+(0.8-0.16*(TM**3)))
C            ELSE
C              ALBEDO=0.67
C            ENDIF
C          ENDIF

C     ISBA FORMULATION (VERSEGHY, 1991; BAKER ET AL., 1990)
C          IF (TSNOW.LT.273.16) THEN
C            ALBEDO=SNOALB-0.008*DT/86400
C          ELSE
C            ALBEDO=(SNOALB-0.5)*EXP(-0.24*DT/86400)+0.5
C          ENDIF

C ----------------------------------------------------------------------
C END SUBROUTINE ALCALC
C ----------------------------------------------------------------------
      RETURN
      END

CCCC 2. PHYSICS SUBROUTINE ==>  SUBROUTINE CANRES  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CANRES (SOLAR,CH,SFCTMP,Q2,SFCPRS,SMC,ZSOIL,NSOIL,
     &                   SMCWLT,SMCREF,RSMIN,RC,PC,NROOT,Q2SAT,DQSDT2, 
     &                   TOPT,RSMAX,RGL,HS,XLAI,
     &                   RCS,RCT,RCQ,RCSOIL)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE CANRES                    
C ----------------------------------------------------------------------
C CALCULATE CANOPY RESISTANCE WHICH DEPENDS ON INCOMING SOLAR RADIATION,
C AIR TEMPERATURE, ATMOSPHERIC WATER VAPOR PRESSURE DEFICIT AT THE
C LOWEST MODEL LEVEL, AND SOIL MOISTURE (PREFERABLY UNFROZEN SOIL
C MOISTURE RATHER THAN TOTAL)
C ----------------------------------------------------------------------
C SOURCE:  JARVIS (1976), NOILHAN AND PLANTON (1989, MWR), JACQUEMIN AND
C NOILHAN (1990, BLM)
C SEE ALSO:  CHEN ET AL (1996, JGR, VOL 101(D3), 7251-7268), EQNS 12-14
C AND TABLE 2 OF SEC. 3.1.2         
C ----------------------------------------------------------------------
C INPUT:
C   SOLAR   INCOMING SOLAR RADIATION
C   CH      SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE
C   SFCTMP  AIR TEMPERATURE AT 1ST LEVEL ABOVE GROUND
C   Q2      AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND
C   Q2SAT   SATURATION AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND
C   DQSDT2  SLOPE OF SATURATION HUMIDITY FUNCTION WRT TEMP
C   SFCPRS  SURFACE PRESSURE
C   SMC     VOLUMETRIC SOIL MOISTURE 
C   ZSOIL   SOIL DEPTH (NEGATIVE SIGN, AS IT IS BELOW GROUND)
C   NSOIL   NO. OF SOIL LAYERS
C   NROOT   NO. OF SOIL LAYERS IN ROOT ZONE (1.LE.NROOT.LE.NSOIL)
C   XLAI    LEAF AREA INDEX
C   SMCWLT  WILTING POINT
C   SMCREF  REFERENCE SOIL MOISTURE (WHERE SOIL WATER DEFICIT STRESS
C             SETS IN)
C RSMIN, RSMAX, TOPT, RGL, HS ARE CANOPY STRESS PARAMETERS SET IN
C   SURBOUTINE REDPRM
C OUTPUT:
C   PC  PLANT COEFFICIENT
C   RC  CANOPY RESISTANCE
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER K
      INTEGER NROOT
      INTEGER NSOIL

      REAL CH
      REAL CP
      REAL DELTA
      REAL DQSDT2
      REAL FF
      REAL GX
      REAL HS
      REAL P
      REAL PART(NSOLD) 
      REAL PC
      REAL Q2
      REAL Q2SAT
      REAL RC
      REAL RSMIN
      REAL RCQ
      REAL RCS
      REAL RCSOIL
      REAL RCT
      REAL RD
      REAL RGL
      REAL RR
      REAL RSMAX
      REAL SFCPRS
      REAL SFCTMP
      REAL SIGMA
      REAL SLV
      REAL SMC(NSOIL)
      REAL SMCREF
      REAL SMCWLT
      REAL SOLAR
      REAL TOPT
      REAL SLVCP
      REAL ST1
      REAL TAIR4
      REAL XLAI
      REAL ZSOIL(NSOIL)

      PARAMETER(CP = 1004.5)
      PARAMETER(RD = 287.04)
      PARAMETER(SIGMA = 5.67E-8)
      PARAMETER(SLV = 2.501000E6)

C ----------------------------------------------------------------------
C INITIALIZE CANOPY RESISTANCE MULTIPLIER TERMS.
C ----------------------------------------------------------------------
      RCS = 0.0
      RCT = 0.0
      RCQ = 0.0
      RCSOIL = 0.0
      RC = 0.0

C ----------------------------------------------------------------------
C CONTRIBUTION DUE TO INCOMING SOLAR RADIATION
C ----------------------------------------------------------------------
      FF = 0.55*2.0*SOLAR/(RGL*XLAI)
      RCS = (FF + RSMIN/RSMAX) / (1.0 + FF)
      RCS = MAX(RCS,0.0001)

C ----------------------------------------------------------------------
C CONTRIBUTION DUE TO AIR TEMPERATURE AT FIRST MODEL LEVEL ABOVE GROUND
C RCT EXPRESSION FROM NOILHAN AND PLANTON (1989, MWR).
C ----------------------------------------------------------------------
      RCT = 1.0 - 0.0016*((TOPT-SFCTMP)**2.0)
      RCT = MAX(RCT,0.0001)

C ----------------------------------------------------------------------
C CONTRIBUTION DUE TO VAPOR PRESSURE DEFICIT AT FIRST MODEL LEVEL.
C RCQ EXPRESSION FROM SSIB 
C ----------------------------------------------------------------------
      RCQ = 1.0/(1.0+HS*(Q2SAT-Q2))
      RCQ = MAX(RCQ,0.01)

C ----------------------------------------------------------------------
C CONTRIBUTION DUE TO SOIL MOISTURE AVAILABILITY.
C DETERMINE CONTRIBUTION FROM EACH SOIL LAYER, THEN ADD THEM UP.
C ----------------------------------------------------------------------
      GX = (SMC(1) - SMCWLT) / (SMCREF - SMCWLT)
      IF (GX .GT. 1.) GX = 1.
      IF (GX .LT. 0.) GX = 0.

C ----------------------------------------------------------------------
C USE SOIL DEPTH AS WEIGHTING FACTOR
C ----------------------------------------------------------------------
      PART(1) = (ZSOIL(1)/ZSOIL(NROOT)) * GX
C ----------------------------------------------------------------------
C USE ROOT DISTRIBUTION AS WEIGHTING FACTOR
C      PART(1) = RTDIS(1) * GX
C ----------------------------------------------------------------------
      DO K = 2,NROOT
        GX = (SMC(K) - SMCWLT) / (SMCREF - SMCWLT)
        IF (GX .GT. 1.) GX = 1.
        IF (GX .LT. 0.) GX = 0.
C ----------------------------------------------------------------------
C USE SOIL DEPTH AS WEIGHTING FACTOR        
C ----------------------------------------------------------------------
        PART(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT)) * GX
C ----------------------------------------------------------------------
C USE ROOT DISTRIBUTION AS WEIGHTING FACTOR
C        PART(K) = RTDIS(K) * GX 
C ----------------------------------------------------------------------
      END DO

      DO K = 1,NROOT
        RCSOIL = RCSOIL+PART(K)
      END DO
      RCSOIL = MAX(RCSOIL,0.0001)

C ----------------------------------------------------------------------
C DETERMINE CANOPY RESISTANCE DUE TO ALL FACTORS.  CONVERT CANOPY
C RESISTANCE (RC) TO PLANT COEFFICIENT (PC) TO BE USED WITH POTENTIAL
C EVAP IN DETERMINING ACTUAL EVAP.  PC IS DETERMINED BY:
C   PC * LINERIZED PENMAN POTENTIAL EVAP =
C   PENMAN-MONTEITH ACTUAL EVAPORATION (CONTAINING RC TERM).
C ----------------------------------------------------------------------
      RC = RSMIN/(XLAI*RCS*RCT*RCQ*RCSOIL)

c      TAIR4 = SFCTMP**4.
c      ST1 = (4.*SIGMA*RD)/CP
c      SLVCP = SLV/CP
c      RR = ST1*TAIR4/(SFCPRS*CH) + 1.0
      RR = (4.*SIGMA*RD/CP)*(SFCTMP**4.)/(SFCPRS*CH) + 1.0
      DELTA = (SLV/CP)*DQSDT2

      PC = (RR+DELTA)/(RR*(1.+RC*CH)+DELTA)

C ----------------------------------------------------------------------
C END SUBROUTINE CANRES
C ----------------------------------------------------------------------
      RETURN
      END


      FUNCTION CSNOW (DSNOW)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C FUNCTION CSNOW
C ----------------------------------------------------------------------
C CALCULATE SNOW TERMAL CONDUCTIVITY
C ----------------------------------------------------------------------
      REAL C
      REAL DSNOW
      REAL CSNOW
      REAL UNIT

      PARAMETER(UNIT = 0.11631) 
                                         
C ----------------------------------------------------------------------
C CSNOW IN UNITS OF CAL/(CM*HR*C), RETURNED IN W/(M*C)
C BASIC VERSION IS DYACHKOVA EQUATION (1960), FOR RANGE 0.1-0.4
C ----------------------------------------------------------------------
      C=0.328*10**(2.25*DSNOW)
      CSNOW=UNIT*C

C ----------------------------------------------------------------------
C DE VAUX EQUATION (1933), IN RANGE 0.1-0.6
C ----------------------------------------------------------------------
C      CSNOW=0.0293*(1.+100.*DSNOW**2)
      
C ----------------------------------------------------------------------
C E. ANDERSEN FROM FLERCHINGER
C ----------------------------------------------------------------------
C      CSNOW=0.021+2.51*DSNOW**2        
      
C ----------------------------------------------------------------------
C END FUNCTION CSNOW
C ----------------------------------------------------------------------
      RETURN                                                      
      END


      FUNCTION DEVAP (ETP1,SMC,ZSOIL,SHDFAC,SMCMAX,BEXP,
     &                DKSAT,DWSAT,SMCDRY,SMCREF,SMCWLT,FXEXP)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C FUNCTION DEVAP
C ----------------------------------------------------------------------
C CALCULATE DIRECT SOIL EVAPORATION
C ----------------------------------------------------------------------
      REAL BEXP
      REAL DEVAP
      REAL DKSAT
      REAL DWSAT
      REAL ETP1
      REAL FX
      REAL FXEXP
      REAL SHDFAC
      REAL SMC
      REAL SMCDRY
      REAL SMCMAX
      REAL ZSOIL
      REAL SMCREF
      REAL SMCWLT

C ----------------------------------------------------------------------
C DIRECT EVAP A FUNCTION OF RELATIVE SOIL MOISTURE AVAILABILITY, LINEAR
C WHEN FXEXP=1.
C ----------------------------------------------------------------------
      FX = ( (SMC - SMCDRY) / (SMCMAX - SMCDRY) )**FXEXP

C ----------------------------------------------------------------------
C FX > 1 REPRESENTS DEMAND CONTROL
C FX < 1 REPRESENTS FLUX CONTROL
C ----------------------------------------------------------------------
      FX = MAX ( MIN ( FX, 1. ) ,0. )

C ----------------------------------------------------------------------
C ALLOW FOR THE DIRECT-EVAP-REDUCING EFFECT OF SHADE
C ----------------------------------------------------------------------
      DEVAP = FX * ( 1.0 - SHDFAC ) * ETP1

C ----------------------------------------------------------------------
C END FUNCTION DEVAP
C ----------------------------------------------------------------------
      RETURN
      END


CCCCC 3. PHYSICS SUBROUTINE ==>  SUBROUTINE HRT  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE HRT (RHSTS,STC,SMC,SMCMAX,NSOIL,ZSOIL,YY,ZZ1,
     &                TBOT,ZBOT,PSISAT,SH2O,DT,BEXP,
     &                F1,DF1,QUARTZ,CSOIL,AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE HRT
C ----------------------------------------------------------------------
C CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
C THERMAL DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX
C COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      LOGICAL*1 ITAVG

      INTEGER I
      INTEGER K
      INTEGER NSOIL

C ----------------------------------------------------------------------
C DECLARE WORK ARRAYS NEEDED IN TRI-DIAGONAL IMPLICIT SOLVER
C ----------------------------------------------------------------------
      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

C ----------------------------------------------------------------------
C DECLARATIONS
C ----------------------------------------------------------------------
      REAL BEXP
      REAL CAIR
      REAL CH2O
      REAL CICE
      REAL CSOIL
      REAL DDZ
      REAL DDZ2
      REAL DENOM
      REAL DF1
      REAL DF1N
      REAL DF1K
      REAL DT
      REAL DTSDZ
      REAL DTSDZ2
      REAL F1
      REAL HCPCT
      REAL PSISAT
      REAL QUARTZ
      REAL QTOT
      REAL RHSTS(NSOIL)
      REAL SSOIL
      REAL SICE
      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SMCMAX
      REAL SNKSRC
      REAL STC(NSOIL)
      REAL T0
      REAL TAVG
      REAL TBK
      REAL TBK1
      REAL TBOT
      REAL ZBOT
      REAL TSNSR
      REAL TSURF
      REAL YY
      REAL ZSOIL(NSOIL)
      REAL ZZ1

      PARAMETER(T0 = 273.15)

C ----------------------------------------------------------------------
C SET SPECIFIC HEAT CAPACITIES OF AIR, WATER, ICE, SOIL MINERAL       
C ----------------------------------------------------------------------
      PARAMETER(CAIR = 1004.0)
      PARAMETER(CH2O = 4.2E6)
      PARAMETER(CICE = 2.106E6)
C NOTE: CSOIL NOW SET IN ROUTINE REDPRM AND PASSED IN
C      PARAMETER(CSOIL = 1.26E6)

C ----------------------------------------------------------------------
C INITIALIZE LOGICAL FOR SOIL LAYER TEMPERATURE AVERAGING.
C ----------------------------------------------------------------------
      ITAVG = .TRUE.
C      ITAVG = .FALSE.

C ----------------------------------------------------------------------
C BEGIN SECTION FOR TOP SOIL LAYER
C ----------------------------------------------------------------------
C CALC THE HEAT CAPACITY OF THE TOP SOIL LAYER
C ----------------------------------------------------------------------
      HCPCT = SH2O(1)*CH2O + (1.0-SMCMAX)*CSOIL + (SMCMAX-SMC(1))*CAIR
     &        + ( SMC(1) - SH2O(1) )*CICE

C ----------------------------------------------------------------------
C CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
C ----------------------------------------------------------------------
      DDZ = 1.0 / ( -0.5 * ZSOIL(2) )
      AI(1) = 0.0
      CI(1) = (DF1 * DDZ) / (ZSOIL(1) * HCPCT)
      BI(1) = -CI(1) + DF1 / (0.5 * ZSOIL(1) * ZSOIL(1)*HCPCT*ZZ1)

C ----------------------------------------------------------------------
C CALCULATE THE VERTICAL SOIL TEMP GRADIENT BTWN THE 1ST AND 2ND SOIL
C LAYERS.  THEN CALCULATE THE SUBSURFACE HEAT FLUX. USE THE TEMP
C GRADIENT AND SUBSFC HEAT FLUX TO CALC "RIGHT-HAND SIDE TENDENCY
C TERMS", OR "RHSTS", FOR TOP SOIL LAYER.
C ----------------------------------------------------------------------
      DTSDZ = (STC(1) - STC(2)) / (-0.5 * ZSOIL(2))
      SSOIL = DF1 * (STC(1) - YY) / (0.5 * ZSOIL(1) * ZZ1)
      RHSTS(1) = (DF1 * DTSDZ - SSOIL) / (ZSOIL(1) * HCPCT)

C ----------------------------------------------------------------------
C NEXT CAPTURE THE VERTICAL DIFFERENCE OF THE HEAT FLUX AT TOP AND
C BOTTOM OF FIRST SOIL LAYER FOR USE IN HEAT FLUX CONSTRAINT APPLIED TO
C POTENTIAL SOIL FREEZING/THAWING IN ROUTINE SNKSRC.
C ----------------------------------------------------------------------
      QTOT = SSOIL - DF1*DTSDZ

C ----------------------------------------------------------------------
C IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):
C SET TEMP "TSURF" AT TOP OF SOIL COLUMN (FOR USE IN FREEZING SOIL
C PHYSICS LATER IN FUNCTION SUBROUTINE SNKSRC).  IF SNOWPACK CONTENT IS
C ZERO, THEN TSURF EXPRESSION BELOW GIVES TSURF = SKIN TEMP.  IF
C SNOWPACK IS NONZERO (HENCE ARGUMENT ZZ1=1), THEN TSURF EXPRESSION
C BELOW YIELDS SOIL COLUMN TOP TEMPERATURE UNDER SNOWPACK.  THEN
C CALCULATE TEMPERATURE AT BOTTOM INTERFACE OF 1ST SOIL LAYER FOR USE
C LATER IN FUNCTION SUBROUTINE SNKSRC
C ----------------------------------------------------------------------
      IF (ITAVG) THEN 
        TSURF = (YY + (ZZ1-1) * STC(1)) / ZZ1
        CALL TBND (STC(1),STC(2),ZSOIL,ZBOT,1,NSOIL,TBK)
      ENDIF

C ----------------------------------------------------------------------
C CALCULATE FROZEN WATER CONTENT IN 1ST SOIL LAYER. 
C ----------------------------------------------------------------------
      SICE = SMC(1) - SH2O(1)

C ----------------------------------------------------------------------
C IF FROZEN WATER PRESENT OR ANY OF LAYER-1 MID-POINT OR BOUNDING
C INTERFACE TEMPERATURES BELOW FREEZING, THEN CALL SNKSRC TO
C COMPUTE HEAT SOURCE/SINK (AND CHANGE IN FROZEN WATER CONTENT)
C DUE TO POSSIBLE SOIL WATER PHASE CHANGE
C ----------------------------------------------------------------------
      IF ( (SICE   .GT. 0.) .OR. (TSURF .LT. T0) .OR.
     &     (STC(1) .LT. T0) .OR. (TBK   .LT. T0) ) THEN

        IF (ITAVG) THEN 
          CALL TMPAVG(TAVG,TSURF,STC(1),TBK,ZSOIL,NSOIL,1)
        ELSE
          TAVG = STC(1)
        ENDIF
        TSNSR = SNKSRC (TAVG,SMC(1),SH2O(1), 
     &    ZSOIL,NSOIL,SMCMAX,PSISAT,BEXP,DT,1,QTOT)

        RHSTS(1) = RHSTS(1) - TSNSR / ( ZSOIL(1) * HCPCT )
      ENDIF
 
C ----------------------------------------------------------------------
C THIS ENDS SECTION FOR TOP SOIL LAYER.
C ----------------------------------------------------------------------
C INITIALIZE DDZ2
C ----------------------------------------------------------------------
      DDZ2 = 0.0

C ----------------------------------------------------------------------
C LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS
C (EXCEPT SUBSFC OR "GROUND" HEAT FLUX NOT REPEATED IN LOWER LAYERS)
C ----------------------------------------------------------------------
      DF1K = DF1
      DO K = 2,NSOIL

C ----------------------------------------------------------------------
C CALCULATE HEAT CAPACITY FOR THIS SOIL LAYER.
C ----------------------------------------------------------------------
        HCPCT = SH2O(K)*CH2O +(1.0-SMCMAX)*CSOIL +(SMCMAX-SMC(K))*CAIR
     &        + ( SMC(K) - SH2O(K) )*CICE

        IF (K .NE. NSOIL) THEN
C ----------------------------------------------------------------------
C THIS SECTION FOR LAYER 2 OR GREATER, BUT NOT LAST LAYER.
C ----------------------------------------------------------------------
C CALCULATE THERMAL DIFFUSIVITY FOR THIS LAYER.
C ----------------------------------------------------------------------
          CALL TDFCND (DF1N,SMC(K),QUARTZ,SMCMAX,SH2O(K))

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER
C ----------------------------------------------------------------------
          DENOM = 0.5 * ( ZSOIL(K-1) - ZSOIL(K+1) )
          DTSDZ2 = ( STC(K) - STC(K+1) ) / DENOM

C ----------------------------------------------------------------------
C CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
C ----------------------------------------------------------------------
          DDZ2 = 2. / (ZSOIL(K-1) - ZSOIL(K+1))
          CI(K) = -DF1N * DDZ2 / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)

C ----------------------------------------------------------------------
C IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):  CALCULATE
C TEMP AT BOTTOM OF LAYER.
C ----------------------------------------------------------------------
          IF (ITAVG) THEN 
            CALL TBND (STC(K),STC(K+1),ZSOIL,ZBOT,K,NSOIL,TBK1)
          ENDIF
        ELSE

C ----------------------------------------------------------------------
C SPECIAL CASE OF BOTTOM SOIL LAYER:  CALCULATE THERMAL DIFFUSIVITY FOR
C BOTTOM LAYER.
C ----------------------------------------------------------------------
          CALL TDFCND (DF1N,SMC(K),QUARTZ,SMCMAX,SH2O(K))

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT THRU BOTTOM LAYER.
C ----------------------------------------------------------------------
          DENOM = .5 * (ZSOIL(K-1) + ZSOIL(K)) - ZBOT
          DTSDZ2 = (STC(K)-TBOT) / DENOM

C ----------------------------------------------------------------------
C SET MATRIX COEF, CI TO ZERO IF BOTTOM LAYER.
C ----------------------------------------------------------------------
          CI(K) = 0.

C ----------------------------------------------------------------------
C IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):  CALCULATE
C TEMP AT BOTTOM OF LAST LAYER.
C ----------------------------------------------------------------------
          IF (ITAVG) THEN 
            CALL TBND (STC(K),TBOT,ZSOIL,ZBOT,K,NSOIL,TBK1)
          ENDIF 

        ENDIF
C ----------------------------------------------------------------------
C THIS ENDS SPECIAL LOOP FOR BOTTOM LAYER.
C ----------------------------------------------------------------------
C CALCULATE RHSTS FOR THIS LAYER AFTER CALC'NG A PARTIAL PRODUCT.
C ----------------------------------------------------------------------
        DENOM = ( ZSOIL(K) - ZSOIL(K-1) ) * HCPCT
        RHSTS(K) = ( DF1N * DTSDZ2 - DF1K * DTSDZ ) / DENOM
        QTOT = -1.0*DENOM*RHSTS(K)
        SICE = SMC(K) - SH2O(K)

        IF ( (SICE .GT. 0.) .OR. (TBK .LT. T0) .OR.
     &     (STC(K) .LT. T0) .OR. (TBK1 .LT. T0) ) THEN

          IF (ITAVG) THEN 
            CALL TMPAVG(TAVG,TBK,STC(K),TBK1,ZSOIL,NSOIL,K)
          ELSE
            TAVG = STC(K)
          ENDIF
          TSNSR = SNKSRC(TAVG,SMC(K),SH2O(K),ZSOIL,NSOIL,
     &                   SMCMAX,PSISAT,BEXP,DT,K,QTOT)
          RHSTS(K) = RHSTS(K) - TSNSR / DENOM
        ENDIF 

C ----------------------------------------------------------------------
C CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER.
C ----------------------------------------------------------------------
        AI(K) = - DF1 * DDZ / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        BI(K) = -(AI(K) + CI(K))

C ----------------------------------------------------------------------
C RESET VALUES OF DF1, DTSDZ, DDZ, AND TBK FOR LOOP TO NEXT SOIL LAYER.
C ----------------------------------------------------------------------
        TBK   = TBK1
        DF1K  = DF1N
        DTSDZ = DTSDZ2
        DDZ   = DDZ2
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE HRT
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 4. PHYSICS SUBROUTINE ==>  SUBROUTINE HRTICE  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE HRTICE (RHSTS,STC,NSOIL,ZSOIL,YY,ZZ1,DF1,AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE HRTICE
C ----------------------------------------------------------------------
C CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
C THERMAL DIFFUSION EQUATION IN THE CASE OF SEA-ICE PACK.  ALSO TO
C COMPUTE (PREPARE) THE MATRIX COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX
C OF THE IMPLICIT TIME SCHEME.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER K
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

      REAL DDZ
      REAL DDZ2
      REAL DENOM
      REAL DF1
      REAL DTSDZ
      REAL DTSDZ2
      REAL HCPCT
      REAL RHSTS(NSOIL)
      REAL SSOIL
      REAL STC(NSOIL)
      REAL TBOT
      REAL YY
      REAL ZBOT
      REAL ZSOIL(NSOIL)
      REAL ZZ1

      DATA TBOT /271.16/

C ----------------------------------------------------------------------
C SET A NOMINAL UNIVERSAL VALUE OF THE SEA-ICE SPECIFIC HEAT CAPACITY,
C HCPCT = 1880.0*917.0.
C ----------------------------------------------------------------------
      PARAMETER(HCPCT = 1.72396E+6)

C ----------------------------------------------------------------------
C THE INPUT ARGUMENT DF1 IS A UNIVERSALLY CONSTANT VALUE OF SEA-ICE
C THERMAL DIFFUSIVITY, SET IN ROUTINE SNOPAC AS DF1 = 2.2.
C ----------------------------------------------------------------------
C SET ICE PACK DEPTH.  USE TBOT AS ICE PACK LOWER BOUNDARY TEMPERATURE
C (THAT OF UNFROZEN SEA WATER AT BOTTOM OF SEA ICE PACK).  ASSUME ICE
C PACK IS OF N=NSOIL LAYERS SPANNING A UNIFORM CONSTANT ICE PACK
C THICKNESS AS DEFINED BY ZSOIL(NSOIL) IN ROUTINE SFLX.
C ----------------------------------------------------------------------
      ZBOT = ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
C ----------------------------------------------------------------------
      DDZ = 1.0 / ( -0.5 * ZSOIL(2) )
      AI(1) = 0.0
      CI(1) = (DF1 * DDZ) / (ZSOIL(1) * HCPCT)
      BI(1) = -CI(1) + DF1/(0.5 * ZSOIL(1) * ZSOIL(1) * HCPCT * ZZ1)

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT BTWN THE TOP AND 2ND SOIL LAYERS.
C RECALC/ADJUST THE SOIL HEAT FLUX.  USE THE GRADIENT AND FLUX TO CALC
C RHSTS FOR THE TOP SOIL LAYER.
C ----------------------------------------------------------------------
      DTSDZ = ( STC(1) - STC(2) ) / ( -0.5 * ZSOIL(2) )
      SSOIL = DF1 * ( STC(1) - YY ) / ( 0.5 * ZSOIL(1) * ZZ1 )
      RHSTS(1) = ( DF1 * DTSDZ - SSOIL ) / ( ZSOIL(1) * HCPCT )

C ----------------------------------------------------------------------
C INITIALIZE DDZ2
C ----------------------------------------------------------------------
      DDZ2 = 0.0

C ----------------------------------------------------------------------
C LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
        IF (K .NE. NSOIL) THEN

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER.
C ----------------------------------------------------------------------
          DENOM = 0.5 * ( ZSOIL(K-1) - ZSOIL(K+1) )
          DTSDZ2 = ( STC(K) - STC(K+1) ) / DENOM

C ----------------------------------------------------------------------
C CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT.
C ----------------------------------------------------------------------
          DDZ2 = 2. / (ZSOIL(K-1) - ZSOIL(K+1))
          CI(K) = -DF1 * DDZ2 / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        ELSE

C ----------------------------------------------------------------------
C CALC THE VERTICAL SOIL TEMP GRADIENT THRU THE LOWEST LAYER.
C ----------------------------------------------------------------------
          DTSDZ2 = (STC(K)-TBOT)/(.5 * (ZSOIL(K-1) + ZSOIL(K))-ZBOT)

C ----------------------------------------------------------------------
C SET MATRIX COEF, CI TO ZERO.
C ----------------------------------------------------------------------
          CI(K) = 0.
        ENDIF

C ----------------------------------------------------------------------
C CALC RHSTS FOR THIS LAYER AFTER CALC'NG A PARTIAL PRODUCT.
C ----------------------------------------------------------------------
        DENOM = ( ZSOIL(K) - ZSOIL(K-1) ) * HCPCT
        RHSTS(K) = ( DF1 * DTSDZ2 - DF1 * DTSDZ ) / DENOM

C ----------------------------------------------------------------------
C CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER.
C ----------------------------------------------------------------------
        AI(K) = - DF1 * DDZ / ((ZSOIL(K-1) - ZSOIL(K)) * HCPCT)
        BI(K) = -(AI(K) + CI(K))

C ----------------------------------------------------------------------
C RESET VALUES OF DTSDZ AND DDZ FOR LOOP TO NEXT SOIL LYR.
C ----------------------------------------------------------------------
        DTSDZ = DTSDZ2
        DDZ   = DDZ2

      END DO
C ----------------------------------------------------------------------
C END SUBROUTINE HRTICE
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 5. PHYSICS SUBROUTINE ==>  SUBROUTINE HSTEP  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE HSTEP (STCOUT,STCIN,RHSTS,DT,NSOIL,AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE HSTEP
C ----------------------------------------------------------------------
C CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER K
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)
      REAL CIin(NSOLD)
      REAL DT
      REAL RHSTS(NSOIL)
      REAL RHSTSin(NSOIL)
      REAL STCIN(NSOIL)
      REAL STCOUT(NSOIL)

C ----------------------------------------------------------------------
C CREATE FINITE DIFFERENCE VALUES FOR USE IN ROSR12 ROUTINE
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        RHSTS(K) = RHSTS(K) * DT
        AI(K) = AI(K) * DT
        BI(K) = 1. + BI(K) * DT
        CI(K) = CI(K) * DT
      END DO

C ----------------------------------------------------------------------
C COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
         RHSTSin(K) = RHSTS(K)
      END DO
      DO K = 1,NSOLD
        CIin(K) = CI(K)
      END DO

C ----------------------------------------------------------------------
C SOLVE THE TRI-DIAGONAL MATRIX EQUATION
C ----------------------------------------------------------------------
      CALL ROSR12(CI,AI,BI,CIin,RHSTSin,RHSTS,NSOIL)

C ----------------------------------------------------------------------
C CALC/UPDATE THE SOIL TEMPS USING MATRIX SOLUTION
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        STCOUT(K) = STCIN(K) + CI(K)
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE HSTEP
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 6. PHYSICS SUBROUTINE ==>  SUBROUTINE NOPAC  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NOPAC(ETP,ETA,PRCP,SMC,SMCMAX,SMCWLT,
     &                 SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,SHDFAC,
     &                 SBETA,Q2,T1,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,
     &                 STC,EPSCA,BEXP,PC,RCH,RR,CFACTR, 
     &                 SH2O,SLOPE,KDT,FRZFACT,PSISAT,ZSOIL,
     &                 DKSAT,DWSAT,TBOT,ZBOT,RUNOFF1,RUNOFF2,
     &                 RUNOFF3,EDIR,EC,ET,ETT,NROOT,ICE,RTDIS,
     &                 QUARTZ,FXEXP,CSOIL,
     &                 BETA,DRIP,DEW,FLX1,FLX2,FLX3)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE NOPAC
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES AND UPDATE SOIL MOISTURE
C CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN NO SNOW PACK IS
C PRESENT.
C ----------------------------------------------------------------------
      INTEGER ICE
      INTEGER NROOT
      INTEGER NSOIL

      REAL BEXP
      REAL BETA
      REAL CFACTR
      REAL CMC
      REAL CMCMAX
      REAL CP
      REAL CSOIL
      REAL DEW
      REAL DF1
      REAL DKSAT
      REAL DRIP
      REAL DT
      REAL DWSAT
      REAL EC
      REAL EDIR
      REAL EPSCA
      REAL ETA
      REAL ETA1
      REAL ETP
      REAL ETP1
      REAL ET(NSOIL)
      REAL ETT
      REAL FDOWN
      REAL F1
      REAL FXEXP
      REAL FLX1
      REAL FLX2
      REAL FLX3
      REAL FRZFACT
      REAL KDT
      REAL PC
      REAL PRCP
      REAL PRCP1
      REAL PSISAT
      REAL Q2
      REAL QUARTZ
      REAL RCH
      REAL RR
      REAL RTDIS(NSOIL)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RUNOFF3
      REAL SSOIL
      REAL SBETA
      REAL SFCTMP
      REAL SHDFAC
      REAL SH2O(NSOIL)
      REAL SIGMA
      REAL SLOPE
      REAL SMC(NSOIL)
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL STC(NSOIL)
      REAL T1
      REAL T24
      REAL TBOT
      REAL TH2
      REAL YY
      REAL YYNUM
      REAL ZBOT
      REAL ZSOIL(NSOIL)
      REAL ZZ1

      PARAMETER(CP = 1004.5)
      PARAMETER(SIGMA = 5.67E-8)

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE:
C CONVERT ETP FROM KG M-2 S-1 TO MS-1 AND INITIALIZE DEW.
C ----------------------------------------------------------------------
      PRCP1 = PRCP * 0.001
      ETP1 = ETP * 0.001
      DEW = 0.0

      IF (ETP .GT. 0.0) THEN

C ----------------------------------------------------------------------
C CONVERT PRCP FROM 'KG M-2 S-1' TO 'M S-1'.
C ----------------------------------------------------------------------
           CALL SMFLX (ETA1,SMC,NSOIL,CMC,ETP1,DT,PRCP1,ZSOIL,
     &                 SH2O,SLOPE,KDT,FRZFACT,
     &                 SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
     &                 SMCREF,SHDFAC,CMCMAX,
     &                 SMCDRY,CFACTR,RUNOFF1,RUNOFF2,RUNOFF3, 
     &                 EDIR,EC,ET,ETT,SFCTMP,Q2,NROOT,RTDIS,FXEXP,
     &                 DRIP)

C ----------------------------------------------------------------------
C       CONVERT MODELED EVAPOTRANSPIRATION FM  M S-1  TO  KG M-2 S-1
C ----------------------------------------------------------------------
        ETA = ETA1 * 1000.0

      ELSE

C ----------------------------------------------------------------------
C IF ETP < 0, ASSUME DEW FORMS (TRANSFORM ETP1 INTO DEW AND REINITIALIZE
C ETP1 TO ZERO).
C ----------------------------------------------------------------------
        DEW = -ETP1
        ETP1 = 0.0

C ----------------------------------------------------------------------
C CONVERT PRCP FROM 'KG M-2 S-1' TO 'M S-1' AND ADD DEW AMOUNT.
C ----------------------------------------------------------------------
        PRCP1 = PRCP1 + DEW
C
      CALL SMFLX (ETA1,SMC,NSOIL,CMC,ETP1,DT,PRCP1,ZSOIL,
     &            SH2O,SLOPE,KDT,FRZFACT,
     &            SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
     &            SMCREF,SHDFAC,CMCMAX,
     &            SMCDRY,CFACTR,RUNOFF1,RUNOFF2,RUNOFF3, 
     &            EDIR,EC,ET,ETT,SFCTMP,Q2,NROOT,RTDIS,FXEXP,
     &            DRIP)

C ----------------------------------------------------------------------
C CONVERT MODELED EVAPOTRANSPIRATION FROM 'M S-1' TO 'KG M-2 S-1'.
C ----------------------------------------------------------------------
        ETA = ETA1 * 1000.0

      ENDIF

C ----------------------------------------------------------------------
C BASED ON ETP AND E VALUES, DETERMINE BETA
C ----------------------------------------------------------------------
      IF ( ETP .LE. 0.0 ) THEN
        BETA = 0.0
        IF ( ETP .LT. 0.0 ) THEN
          BETA = 1.0
          ETA = ETP
        ENDIF
      ELSE
        BETA = ETA / ETP
      ENDIF

C ----------------------------------------------------------------------
C GET SOIL THERMAL DIFFUXIVITY/CONDUCTIVITY FOR TOP SOIL LYR,
C CALC. ADJUSTED TOP LYR SOIL TEMP AND ADJUSTED SOIL FLUX, THEN
C CALL SHFLX TO COMPUTE/UPDATE SOIL HEAT FLUX AND SOIL TEMPS.
C ----------------------------------------------------------------------
      CALL TDFCND (DF1,SMC(1),QUARTZ,SMCMAX,SH2O(1))

C ----------------------------------------------------------------------
C VEGETATION GREENNESS FRACTION REDUCTION IN SUBSURFACE HEAT FLUX 
C VIA REDUCTION FACTOR, WHICH IS CONVENIENT TO APPLY HERE TO THERMAL 
C DIFFUSIVITY THAT IS LATER USED IN HRT TO COMPUTE SUB SFC HEAT FLUX
C (SEE ADDITIONAL COMMENTS ON VEG EFFECT SUB-SFC HEAT FLX IN 
C ROUTINE SFLX)
C ----------------------------------------------------------------------
      DF1 = DF1 * EXP(SBETA*SHDFAC)

C ----------------------------------------------------------------------
C COMPUTE INTERMEDIATE TERMS PASSED TO ROUTINE HRT (VIA ROUTINE 
C SHFLX BELOW) FOR USE IN COMPUTING SUBSURFACE HEAT FLUX IN HRT
C ----------------------------------------------------------------------
      YYNUM = FDOWN - SIGMA * T24
      YY = SFCTMP + (YYNUM/RCH+TH2-SFCTMP-BETA*EPSCA) / RR
      ZZ1 = DF1 / ( -0.5 * ZSOIL(1) * RCH * RR ) + 1.0

      CALL SHFLX (SSOIL,STC,SMC,SMCMAX,NSOIL,T1,DT,YY,ZZ1,ZSOIL,
     &            TBOT,ZBOT,SMCWLT,PSISAT,SH2O,BEXP,F1,DF1,ICE,
     &            QUARTZ,CSOIL)

C ----------------------------------------------------------------------
C SET FLX1 AND FLX3 (SNOPACK PHASE CHANGE HEAT FLUXES) TO ZERO SINCE
C THEY ARE NOT USED HERE IN SNOPAC.  FLX2 (FREEZING RAIN HEAT FLUX) WAS
C SIMILARLY INITIALIZED IN THE PENMAN ROUTINE.
C ----------------------------------------------------------------------
      FLX1 = 0.0
      FLX3 = 0.0

C ----------------------------------------------------------------------
C END SUBROUTINE NOPAC
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 7. PHYSICS SUBROUTINE ==>  SUBROUTINE PENMAN  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PENMAN (SFCTMP,SFCPRS,CH,T2V,TH2,PRCP,FDOWN,T24,SSOIL,
     &                   Q2,Q2SAT,ETP,RCH,EPSCA,RR,SNOWNG,FRZGRA,
     &                   DQSDT2,FLX2)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE PENMAN
C ----------------------------------------------------------------------
C CALCULATE POTENTIAL EVAPORATION FOR THE CURRENT POINT.  VARIOUS
C PARTIAL SUMS/PRODUCTS ARE ALSO CALCULATED AND PASSED BACK TO THE
C CALLING ROUTINE FOR LATER USE.
C ----------------------------------------------------------------------
      LOGICAL*1 SNOWNG
      LOGICAL*1 FRZGRA

      REAL A
      REAL BETA
      REAL CH
      REAL CP
      REAL CPH2O
      REAL CPICE
      REAL DELTA
      REAL DQSDT2
      REAL ELCP
      REAL EPSCA
      REAL ETP
      REAL FDOWN
      REAL FLX2
      REAL FNET
      REAL LSUBC
      REAL LSUBF
      REAL PRCP
      REAL Q2
      REAL Q2SAT
      REAL R
      REAL RAD
      REAL RCH
      REAL RHO
      REAL RR
      REAL SSOIL
      REAL SFCPRS
      REAL SFCTMP
      REAL SIGMA
      REAL T24
      REAL T2V
      REAL TH2

      PARAMETER(CP = 1004.6)
      PARAMETER(CPH2O = 4.218E+3)
      PARAMETER(CPICE = 2.106E+3)
      PARAMETER(R = 287.04)
      PARAMETER(ELCP = 2.4888E+3)
      PARAMETER(LSUBF = 3.335E+5)
      PARAMETER(LSUBC = 2.501000E+6)
      PARAMETER(SIGMA = 5.67E-8)

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE:
C ----------------------------------------------------------------------
      FLX2 = 0.0

C ----------------------------------------------------------------------
C PREPARE PARTIAL QUANTITIES FOR PENMAN EQUATION.
C ----------------------------------------------------------------------
      DELTA = ELCP * DQSDT2
      T24 = SFCTMP * SFCTMP * SFCTMP * SFCTMP
      RR = T24 * 6.48E-8 /(SFCPRS * CH) + 1.0
      RHO = SFCPRS / (R * T2V)
      RCH = RHO * CP * CH

C ----------------------------------------------------------------------
C ADJUST THE PARTIAL SUMS / PRODUCTS WITH THE LATENT HEAT
C EFFECTS CAUSED BY FALLING PRECIPITATION.
C ----------------------------------------------------------------------
      IF (.NOT. SNOWNG) THEN
        IF (PRCP .GT. 0.0) RR = RR + CPH2O*PRCP/RCH
      ELSE
        RR = RR + CPICE*PRCP/RCH
      ENDIF

      FNET = FDOWN - SIGMA*T24 - SSOIL

C ----------------------------------------------------------------------
C INCLUDE THE LATENT HEAT EFFECTS OF FRZNG RAIN CONVERTING TO ICE ON
C IMPACT IN THE CALCULATION OF FLX2 AND FNET.
C ----------------------------------------------------------------------
      IF (FRZGRA) THEN
        FLX2 = -LSUBF * PRCP
        FNET = FNET - FLX2
      ENDIF

C ----------------------------------------------------------------------
C FINISH PENMAN EQUATION CALCULATIONS.
C ----------------------------------------------------------------------
      RAD = FNET/RCH + TH2 - SFCTMP
      A = ELCP * (Q2SAT - Q2)
      EPSCA = (A*RR + RAD*DELTA) / (DELTA + RR)
      ETP = EPSCA * RCH / LSUBC

C ----------------------------------------------------------------------
C END SUBROUTINE PENMAN
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 8. PHYSICS SUBROUTINE ==>  SUBROUTINE ROSR12  CCCCCCCCCCCCCCCCCCC

      SUBROUTINE ROSR12 (P,A,B,C,D,DELTA,NSOIL)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE ROSR12
C ----------------------------------------------------------------------
C INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
C ###                                            ### ###  ###   ###  ###
C #B(1), C(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
C #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
C # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
C # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
C # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
C # .                                          .   # #  .   # = #   .  #
C # .                                          .   # #  .   #   #   .  #
C # .                                          .   # #  .   #   #   .  #
C # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
C # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
C # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
C ###                                            ### ###  ###   ###  ###
C ----------------------------------------------------------------------
      INTEGER K
      INTEGER KK
      INTEGER NSOIL
      
      REAL A(NSOIL)
      REAL B(NSOIL)
      REAL C(NSOIL)
      REAL D(NSOIL)
      REAL DELTA(NSOIL)
      REAL P(NSOIL)
      
C ----------------------------------------------------------------------
C INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER
C ----------------------------------------------------------------------
      C(NSOIL) = 0.0

C ----------------------------------------------------------------------
C SOLVE THE COEFS FOR THE 1ST SOIL LAYER
C ----------------------------------------------------------------------
      P(1) = -C(1) / B(1)
      DELTA(1) = D(1) / B(1)

C ----------------------------------------------------------------------
C SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
        P(K) = -C(K) * ( 1.0 / (B(K) + A (K) * P(K-1)) )
        DELTA(K) = (D(K)-A(K)*DELTA(K-1))*(1.0/(B(K)+A(K)*P(K-1)))
      END DO

C ----------------------------------------------------------------------
C SET P TO DELTA FOR LOWEST SOIL LAYER
C ----------------------------------------------------------------------
      P(NSOIL) = DELTA(NSOIL)

C ----------------------------------------------------------------------
C ADJUST P FOR SOIL LAYERS 2 THRU NSOIL
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
         KK = NSOIL - K + 1
         P(KK) = P(KK) * P(KK+1) + DELTA(KK)
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE ROSR12
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 9. PHYSICS SUBROUTINE ==>  SUBROUTINE SFCDIF  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SFCDIF (ZLM,Z0,THZ0,THLM,SFCSPD,CZIL,AKMS,AKHS)
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SFCDIF
C ----------------------------------------------------------------------
C CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS.
C SEE CHEN ET AL (1997, BLM)
C ----------------------------------------------------------------------
      
      REAL WWST, WWST2, G, VKRM, EXCM, BETA, BTG, ELFC, WOLD, WNEW
      REAL PIHF, EPSU2, EPSUST, EPSIT, EPSA, ZTMIN, ZTMAX, HPBL, SQVISC
      REAL RIC, RRIC, FHNEU, RFC, RFAC, ZZ, PSLMU, PSLMS, PSLHU, PSLHS
      REAL XX, PSPMU, YY, PSPMS, PSPHU, PSPHS, ZLM, Z0, THZ0, THLM
      REAL SFCSPD, CZIL, AKMS, AKHS, ZILFC, ZU, ZT, RDZ, CXCH
      REAL DTHV, DU2, BTGH, WSTAR2, USTAR, ZSLU, ZSLT, RLOGU, RLOGT
      REAL RLMO, ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4
      REAL XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN, RLMA
CCC   ......REAL ZTFC
      
      INTEGER ITRMX, ILECH, ITR
      
      PARAMETER
     &     (WWST=1.2,WWST2=WWST*WWST,G=9.8,VKRM=0.40,EXCM=0.001
     &     ,BETA=1./270.,BTG=BETA*G,ELFC=VKRM*BTG
     &     ,WOLD=.15,WNEW=1.-WOLD,ITRMX=05,PIHF=3.14159265/2.)
C ----------------------------------------------------------------------
      PARAMETER
     &     (EPSU2=1.E-4,EPSUST=0.07,EPSIT=1.E-4,EPSA=1.E-8
     &     ,ZTMIN=-5.,ZTMAX=1.,HPBL=1000.0
     &     ,SQVISC=258.2)
C ----------------------------------------------------------------------
      PARAMETER
     &     (RIC=0.183,RRIC=1.0/RIC,FHNEU=0.8,RFC=0.191
     &     ,RFAC=RIC/(FHNEU*RFC*RFC))
      
C ----------------------------------------------------------------------
C NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
C ----------------------------------------------------------------------
C LECH'S SURFACE FUNCTIONS
C ----------------------------------------------------------------------
      PSLMU(ZZ)=-0.96*log(1.0-4.5*ZZ)
      PSLMS(ZZ)=ZZ*RRIC-2.076*(1.-1./(ZZ+1.))
      PSLHU(ZZ)=-0.96*log(1.0-4.5*ZZ)
      PSLHS(ZZ)=ZZ*RFAC-2.076*(1.-1./(ZZ+1.))
      
C ----------------------------------------------------------------------
C PAULSON'S SURFACE FUNCTIONS
C ----------------------------------------------------------------------
      PSPMU(XX)=-2.*log((XX+1.)*0.5)-log((XX*XX+1.)*0.5)+2.*ATAN(XX)
     &     -PIHF
      PSPMS(YY)=5.*YY
      PSPHU(XX)=-2.*log((XX*XX+1.)*0.5)
      PSPHS(YY)=5.*YY

C ----------------------------------------------------------------------
C THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND
C OVER SOLID SURFACE (LAND, SEA-ICE).  
C ----------------------------------------------------------------------
      ILECH=0
      
C ----------------------------------------------------------------------
C     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
C     C......ZTFC=0.1
C     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
C ----------------------------------------------------------------------
      ZILFC=-CZIL*VKRM*SQVISC
      
C ----------------------------------------------------------------------
      ZU=Z0
C     C.......ZT=Z0*ZTFC
      RDZ=1./ZLM
      CXCH=EXCM*RDZ
      DTHV=THLM-THZ0     
      DU2=MAX(SFCSPD*SFCSPD,EPSU2)

C ----------------------------------------------------------------------
C BELJARS CORRECTION OF USTAR
C ----------------------------------------------------------------------
      BTGH=BTG*HPBL
ccc   If statements to avoid TANGENT LINEAR problems near zero
      IF (BTGH*AKHS*DTHV .NE. 0.0) THEN
         WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)
      ELSE
         WSTAR2=0.0
      ENDIF
      USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)

C ----------------------------------------------------------------------
C ZILITINKEVITCH APPROACH FOR ZT
C ----------------------------------------------------------------------
      ZT=EXP(ZILFC*SQRT(USTAR*Z0))*Z0

C ----------------------------------------------------------------------
      ZSLU=ZLM+ZU
      ZSLT=ZLM+ZT
C     PRINT*,'ZSLT=',ZSLT
C     PRINT*,'ZLM=',ZLM
C     PRINT*,'ZT=',ZT
C     
      RLOGU=log(ZSLU/ZU)
      RLOGT=log(ZSLT/ZT)
C     
      RLMO=ELFC*AKHS*DTHV/USTAR**3
C     PRINT*,'RLMO=',RLMO
C     PRINT*,'ELFC=',ELFC
C     PRINT*,'AKHS=',AKHS
C     PRINT*,'DTHV=',DTHV
C     PRINT*,'USTAR=',USTAR
      
      DO ITR=1,ITRMX
C ----------------------------------------------------------------------
C 1./MONIN-OBUKKHOV LENGTH-SCALE
C ----------------------------------------------------------------------
         ZETALT=MAX(ZSLT*RLMO,ZTMIN)
         RLMO=ZETALT/ZSLT
         ZETALU=ZSLU*RLMO
         ZETAU=ZU*RLMO
         ZETAT=ZT*RLMO
         
         IF(ILECH.EQ.0) THEN
            IF(RLMO.LT.0.)THEN
               XLU4=1.-16.*ZETALU
               XLT4=1.-16.*ZETALT
               XU4 =1.-16.*ZETAU
               XT4 =1.-16.*ZETAT
               
               XLU=SQRT(SQRT(XLU4))
               XLT=SQRT(SQRT(XLT4))
               XU =SQRT(SQRT(XU4))
               XT =SQRT(SQRT(XT4))
               
               PSMZ=PSPMU(XU)
C     PRINT*,'-----------1------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSPMU(ZETAU)=',PSPMU(ZETAU)
C     PRINT*,'XU=',XU
C     PRINT*,'------------------------'
               SIMM=PSPMU(XLU)-PSMZ+RLOGU
               PSHZ=PSPHU(XT)
               SIMH=PSPHU(XLT)-PSHZ+RLOGT
            ELSE
               ZETALU=MIN(ZETALU,ZTMAX)
               ZETALT=MIN(ZETALT,ZTMAX)
               PSMZ=PSPMS(ZETAU)
C     PRINT*,'-----------2------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSPMS(ZETAU)=',PSPMS(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
               SIMM=PSPMS(ZETALU)-PSMZ+RLOGU
               PSHZ=PSPHS(ZETAT)
               SIMH=PSPHS(ZETALT)-PSHZ+RLOGT
            ENDIF
         ELSE
C ----------------------------------------------------------------------
C LECH'S FUNCTIONS
C ----------------------------------------------------------------------
            IF(RLMO.LT.0.)THEN
               PSMZ=PSLMU(ZETAU)
C     PRINT*,'-----------3------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSLMU(ZETAU)=',PSLMU(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
               SIMM=PSLMU(ZETALU)-PSMZ+RLOGU
               PSHZ=PSLHU(ZETAT)
               SIMH=PSLHU(ZETALT)-PSHZ+RLOGT
            ELSE
               ZETALU=MIN(ZETALU,ZTMAX)
               ZETALT=MIN(ZETALT,ZTMAX)
C     
               PSMZ=PSLMS(ZETAU)
C     PRINT*,'-----------4------------'
C     PRINT*,'PSMZ=',PSMZ
C     PRINT*,'PSLMS(ZETAU)=',PSLMS(ZETAU)
C     PRINT*,'ZETAU=',ZETAU
C     PRINT*,'------------------------'
               SIMM=PSLMS(ZETALU)-PSMZ+RLOGU
               PSHZ=PSLHS(ZETAT)
               SIMH=PSLHS(ZETALT)-PSHZ+RLOGT
            ENDIF
         ENDIF
C ----------------------------------------------------------------------
C BELJAARS CORRECTION FOR USTAR
C ----------------------------------------------------------------------
         USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)

C ----------------------------------------------------------------------
C ZILITINKEVITCH FIX FOR ZT
C ----------------------------------------------------------------------
         ZT=EXP(ZILFC*SQRT(USTAR*Z0))*Z0
         
         ZSLT=ZLM+ZT
         RLOGT=log(ZSLT/ZT)
C-----------------------------------------------------------------------
         USTARK=USTAR*VKRM
         AKMS=MAX(USTARK/SIMM,CXCH)
         AKHS=MAX(USTARK/SIMH,CXCH)
C-----------------------------------------------------------------------
C IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
C-----------------------------------------------------------------------
         IF (BTGH*AKHS*DTHV .NE. 0.0) THEN
            WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)
         ELSE
            WSTAR2=0.0
         ENDIF
         RLMN=ELFC*AKHS*DTHV/USTAR**3
C-----------------------------------------------------------------------
         RLMA=RLMO*WOLD+RLMN*WNEW
C-----------------------------------------------------------------------
C     IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110
C-----------------------------------------------------------------------
         RLMO=RLMA
C-----------------------------------------------------------------------
      END DO

C     PRINT*,'----------------------------'
C     PRINT*,'SFCDIF OUTPUT !  ! ! ! ! ! ! ! !  !   !    !'
C     PRINT*,'ZLM=',ZLM
C     PRINT*,'Z0=',Z0
C     PRINT*,'THZ0=',THZ0
C     PRINT*,'THLM=',THLM
C     PRINT*,'SFCSPD=',SFCSPD
C     PRINT*,'CZIL=',CZIL
C     PRINT*,'AKMS=',AKMS
C     PRINT*,'AKHS=',AKHS
C     PRINT*,'----------------------------'
C     
C ----------------------------------------------------------------------
C END SUBROUTINE SFCDIF
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 10. PHYSICS SUBROUTINE ==>  SUBROUTINE SFLX  CCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SFLX (ntl,
     C  ICE,DT,ZLVL,NSOIL,SLDPTH,
     F  LWDN,SOLDN,SFCPRS,PRCP,SFCTMP,Q2,SFCSPD,
     I  TH2,Q2SAT,DQSDT2,
     S  SLOPE,SHDFAC,SHDMIN,PTU,ALB,SNOALB,
     S  RSMIN,RGL,HS,SNUP,Z0,XLAI,NROOT,
     S  PSISAT,BEXP,DKSAT,SMCMAX,QUARTZ,DWSAT,
     S  SMCWLT,SMCREF,SMCDRY,F1,KDT,FRZX,FRZFACT,TBOT,
     H  CMC,T1,STC,SMC,SH2O,SNOWH,SNEQV,ALBEDO,CH,CM,
     O  EVP,ETA,SHEAT,
C ----------------------------------------------------------------------
C OUTPUTS, DIAGNOSTICS, PARAMETERS BELOW GENERALLY NOT NECESSARY WHEN
C COUPLED WITH E.G. A NWP MODEL (SUCH AS THE NOAA/NWS/NCEP MESOSCALE ETA
C MODEL).  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES. 
C ----------------------------------------------------------------------
     O  EC,EDIR,ET,ETT,ESNOW,DRIP,DEW,
     O  BETA,ETP,SSOIL,
     O  FLX1,FLX2,FLX3,
     O  SNOMLT,SNCOVR,
     O  RUNOFF1,RUNOFF2,RUNOFF3,
     O  RC,PC,RCS,RCT,RCQ,RCSOIL,
     D  SOILW,SOILT,SOILM,
     M  LDAS,TILE)
C     P  SMCWLT,SMCDRY,SMCREF,SMCMAX,NROOT)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables

      IMPLICIT NONE

      type (ldasdec) LDAS
      type (tiledec) TILE

C ----------------------------------------------------------------------
C SUBROUTINE SFLX - VERSION 2.6 - June 2nd 2003
C ----------------------------------------------------------------------
C SUB-DRIVER FOR "NOAH/OSU LSM" FAMILY OF PHYSICS SUBROUTINES FOR A
C SOIL/VEG/SNOWPACK LAND-SURFACE MODEL TO UPDATE SOIL MOISTURE, SOIL
C ICE, SOIL TEMPERATURE, SKIN TEMPERATURE, SNOWPACK WATER CONTENT,
C SNOWDEPTH, AND ALL TERMS OF THE SURFACE ENERGY BALANCE AND SURFACE
C WATER BALANCE (EXCLUDING INPUT ATMOSPHERIC FORCINGS OF DOWNWARD
C RADIATION AND PRECIP)
C ----------------------------------------------------------------------
C SFLX ARGUMENT LIST KEY:
C ----------------------------------------------------------------------
C  C  CONFIGURATION INFORMATION
C  F  FORCING DATA
C  I  OTHER (INPUT) FORCING DATA
C  S  SURFACE CHARACTERISTICS
C  H  HISTORY (STATE) VARIABLES
C  O  OUTPUT VARIABLES
C  D  DIAGNOSTIC OUTPUT
C ----------------------------------------------------------------------
C 1. CONFIGURATION INFORMATION (C):
C ----------------------------------------------------------------------
C   ICE	       SEA-ICE FLAG  (=1: SEA-ICE, =0: LAND)
C   DT	       TIMESTEP (SEC) (DT SHOULD NOT EXCEED 3600 SECS, RECOMMEND
C                1800 SECS OR LESS)
C   ZLVL       HEIGHT (M) ABOVE GROUND OF ATMOSPHERIC FORCING VARIABLES
C   NSOIL      NUMBER OF SOIL LAYERS (AT LEAST 2, AND NOT GREATER THAN
C                PARAMETER NSOLD SET BELOW)
C   SLDPTH     THE THICKNESS OF EACH SOIL LAYER (M)
C ----------------------------------------------------------------------
C 2. FORCING DATA (F):
C ----------------------------------------------------------------------
C   LWDN       LW DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET LONGWAVE)
C   SOLDN      SOLAR DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET SOLAR)
C   SFCPRS     PRESSURE AT HEIGHT ZLVL ABOVE GROUND (PASCALS)
C   PRCP       PRECIP RATE (KG M-2 S-1) (NOTE, THIS IS A RATE)
C   SFCTMP     AIR TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND
C   TH2        AIR POTENTIAL TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND
C   Q2         MIXING RATIO AT HEIGHT ZLVL ABOVE GROUND (KG KG-1)
C ----------------------------------------------------------------------
C 3. OTHER FORCING (INPUT) DATA (I):
C ----------------------------------------------------------------------
C   SFCSPD     WIND SPEED (M S-1) AT HEIGHT ZLVL ABOVE GROUND
C   Q2SAT      SAT MIXING RATIO AT HEIGHT ZLVL ABOVE GROUND (KG KG-1)
C   DQSDT2     SLOPE OF SAT SPECIFIC HUMIDITY CURVE AT T=SFCTMP
C                (KG KG-1 K-1)
C ----------------------------------------------------------------------
C 4. CANOPY/SOIL CHARACTERISTICS (S):
C ----------------------------------------------------------------------
C   VEGTYP     VEGETATION TYPE (INTEGER INDEX)
C   SOILTYP    SOIL TYPE (INTEGER INDEX)
C   SLOPETYP   CLASS OF SFC SLOPE (INTEGER INDEX)
C   SHDFAC     AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C                (FRACTION= 0.0-1.0)
C   SHDMIN     MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
C                (FRACTION= 0.0-1.0) <= SHDFAC
C   PTU        PHOTO THERMAL UNIT (PLANT PHENOLOGY FOR ANNUALS/CROPS)
C                (NOT YET USED, BUT PASSED TO REDPRM FOR FUTURE USE IN
C                VEG PARMS)
C   ALB        BACKROUND SNOW-FREE SURFACE ALBEDO (FRACTION), FOR JULIAN
C                DAY OF YEAR (USUALLY FROM TEMPORAL INTERPOLATION OF
C                MONTHLY MEAN VALUES' CALLING PROG MAY OR MAY NOT
C                INCLUDE DIURNAL SUN ANGLE EFFECT)
C   SNOALB     UPPER BOUND ON MAXIMUM ALBEDO OVER DEEP SNOW (E.G. FROM
C                ROBINSON AND KUKLA, 1985, J. CLIM. & APPL. METEOR.)
C   TBOT       BOTTOM SOIL TEMPERATURE (LOCAL YEARLY-MEAN SFC AIR
C                TEMPERATURE)
C ----------------------------------------------------------------------
C 5. HISTORY (STATE) VARIABLES (H):
C ----------------------------------------------------------------------
C  CMC         CANOPY MOISTURE CONTENT (M)
C  T1          GROUND/CANOPY/SNOWPACK) EFFECTIVE SKIN TEMPERATURE (K)
C  STC(NSOIL)  SOIL TEMP (K)
C  SMC(NSOIL)  TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
C  SH2O(NSOIL) UNFROZEN SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
C                NOTE: FROZEN SOIL MOISTURE = SMC - SH2O
C  SNOWH       ACTUAL SNOW DEPTH (M)
C  SNEQV       LIQUID WATER-EQUIVALENT SNOW DEPTH (M)
C                NOTE: SNOW DENSITY = SNEQV/SNOWH
C  ALBEDO      SURFACE ALBEDO INCLUDING SNOW EFFECT (UNITLESS FRACTION)
C                =SNOW-FREE ALBEDO (ALB) WHEN SNEQV=0, OR
C                =FCT(MSNOALB,ALB,VEGTYP,SHDFAC,SHDMIN) WHEN SNEQV>0
C  CH          SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE
C                (M S-1); NOTE: CH IS TECHNICALLY A CONDUCTANCE SINCE
C                IT HAS BEEN MULTIPLIED BY WIND SPEED.
C  CM          SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM (M S-1); NOTE:
C                CM IS TECHNICALLY A CONDUCTANCE SINCE IT HAS BEEN
C                MULTIPLIED BY WIND SPEED.  CM IS NOT NEEDED IN SFLX
C ----------------------------------------------------------------------
C 6. OUTPUT (O):
C ----------------------------------------------------------------------
C OUTPUT VARIABLES NECESSARY FOR A COUPLED NUMERICAL WEATHER PREDICTION
C MODEL, E.G. NOAA/NWS/NCEP MESOSCALE ETA MODEL.  FOR THIS APPLICATION,
C THE REMAINING OUTPUT/DIAGNOSTIC/PARAMETER BLOCKS BELOW ARE NOT
C NECESSARY.  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES.
C   ETA        ACTUAL LATENT HEAT FLUX (KG M-2 S-1: NEGATIVE, IF UP FROM
C	         SURFACE)
C   SHEAT      SENSIBLE HEAT FLUX (W M-2: NEGATIVE, IF UPWARD FROM
C	         SURFACE)
C ----------------------------------------------------------------------
C   EC         CANOPY WATER EVAPORATION (M S-1)
C   EDIR       DIRECT SOIL EVAPORATION (M S-1)
C   ET(NSOIL)  PLANT TRANSPIRATION FROM A PARTICULAR ROOT (SOIL) LAYER
C                 (M S-1)
C   ETT        TOTAL PLANT TRANSPIRATION (M S-1)
C   ESNOW      SUBLIMATION FROM (OR DEPOSITION TO IF <0) SNOWPACK
C                (M S-1)
C   DRIP       THROUGH-FALL OF PRECIP AND/OR DEW IN EXCESS OF CANOPY
C                WATER-HOLDING CAPACITY (M)
C   DEW        DEWFALL (OR FROSTFALL FOR T<273.15) (M)
C ----------------------------------------------------------------------
C   BETA       RATIO OF ACTUAL/POTENTIAL EVAP (DIMENSIONLESS)
C   ETP        POTENTIAL EVAPORATION (KG M-2 S-1)
C   SSOIL      SOIL HEAT FLUX (W M-2: NEGATIVE IF DOWNWARD FROM SURFACE)
C ----------------------------------------------------------------------
C   FLX1       PRECIP-SNOW SFC (W M-2)
C   FLX2       FREEZING RAIN LATENT HEAT FLUX (W M-2)
C   FLX3       PHASE-CHANGE HEAT FLUX FROM SNOWMELT (W M-2)
C ----------------------------------------------------------------------
C   SNOMLT     SNOW MELT (M) (WATER EQUIVALENT)
C   SNCOVR     FRACTIONAL SNOW COVER (UNITLESS FRACTION, 0-1)
C ----------------------------------------------------------------------
C   RUNOFF1    SURFACE RUNOFF (M S-1), NOT INFILTRATING THE SURFACE
C   RUNOFF2    SUBSURFACE RUNOFF (M S-1), DRAINAGE OUT BOTTOM OF LAST
C                SOIL LAYER
C   RUNOFF3    NUMERICAL TRUNCTATION IN EXCESS OF POROSITY (SMCMAX)
C                FOR A GIVEN SOIL LAYER AT THE END OF A TIME STEP
C ----------------------------------------------------------------------
C   RC         CANOPY RESISTANCE (S M-1)
C   PC         PLANT COEFFICIENT (UNITLESS FRACTION, 0-1) WHERE PC*ETP
C                = ACTUAL TRANSP
C   XLAI       LEAF AREA INDEX (DIMENSIONLESS)
C   RSMIN      MINIMUM CANOPY RESISTANCE (S M-1)
C   RCS        INCOMING SOLAR RC FACTOR (DIMENSIONLESS)
C   RCT        AIR TEMPERATURE RC FACTOR (DIMENSIONLESS)
C   RCQ        ATMOS VAPOR PRESSURE DEFICIT RC FACTOR (DIMENSIONLESS)
C   RCSOIL     SOIL MOISTURE RC FACTOR (DIMENSIONLESS)
C ----------------------------------------------------------------------
C 7. DIAGNOSTIC OUTPUT (D):
C ----------------------------------------------------------------------
C   SOILW      AVAILABLE SOIL MOISTURE IN ROOT ZONE (UNITLESS FRACTION
C	         BETWEEN SMCWLT AND SMCMAX)
C   SOILM      TOTAL SOIL COLUMN MOISTURE CONTENT (FROZEN+UNFROZEN) (M) 
C ----------------------------------------------------------------------
C 8. PARAMETERS (P):
C ----------------------------------------------------------------------
C   SMCWLT     WILTING POINT (VOLUMETRIC)
C   SMCDRY     DRY SOIL MOISTURE THRESHOLD WHERE DIRECT EVAP FRM TOP
C                LAYER ENDS (VOLUMETRIC)
C   SMCREF     SOIL MOISTURE THRESHOLD WHERE TRANSPIRATION BEGINS TO
C                STRESS (VOLUMETRIC)
C   SMCMAX     POROSITY, I.E. SATURATED VALUE OF SOIL MOISTURE
C                (VOLUMETRIC)
C   NROOT      NUMBER OF ROOT LAYERS, A FUNCTION OF VEG TYPE, DETERMINED
C              IN SUBROUTINE REDPRM.
C ----------------------------------------------------------------------
      INTEGER, PARAMETER:: NSOLD=20
      INTEGER NTL

C ----------------------------------------------------------------------
C DECLARATIONS - LOGICAL
C ----------------------------------------------------------------------
      LOGICAL*1 FRZGRA
      LOGICAL*1 SATURATED
      LOGICAL*1 SNOWNG

C ----------------------------------------------------------------------
C DECLARATIONS - INTEGER
C ----------------------------------------------------------------------
      INTEGER ICE
      INTEGER I,K
      INTEGER KZ
      INTEGER NSOIL
      INTEGER NROOT

C ----------------------------------------------------------------------
C DECLARATIONS - REAL
C ----------------------------------------------------------------------
      REAL ALBEDO
      REAL ALB
      REAL BEXP
      REAL BETA
      REAL CFACTR
      REAL CH
      REAL CM
      REAL CMC
      REAL CMCMAX
      REAL CP
      REAL CSNOW
      REAL CSOIL
      REAL CZIL
      REAL DEW
      REAL DF1
      REAL DF1H
      REAL DF1A
      REAL DKSAT
      REAL DT
      REAL DWSAT
      REAL DQSDT2
      REAL DSOIL
      REAL DTOT
      REAL DRIP
      REAL EC
      REAL EDIR
      REAL ESNOW
      REAL ET(NSOIL)
      REAL ETT
      REAL FRCSNO
      REAL FRCSOI
      REAL EPSCA
      REAL ETA
      REAL ETP
      REAL EVP
      REAL FDOWN
      REAL F1
      REAL FLX1
      REAL FLX2
      REAL FLX3
      REAL FXEXP
      REAL FRZX
      REAL FRZFACT
      REAL SHEAT
      REAL HS
      REAL KDT
      REAL LWDN
      REAL LVH2O
      REAL PC
      REAL PRCP
      REAL PTU
      REAL PRCP1
      REAL PSISAT
      REAL Q2
      REAL Q2SAT
      REAL QUARTZ
      REAL R
      REAL RCH
      REAL REFKDT
      REAL RR
      REAL RTDIS(NSOLD)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RGL
      REAL RUNOFF3
      REAL RSMAX
      REAL RC
      REAL RSMIN
      REAL RCQ
      REAL RCS
      REAL RCSOIL
      REAL RCT
      REAL RSNOW
      REAL SNDENS
      REAL SNCOND 
      REAL SSOIL
      REAL SBETA
      REAL SFCPRS
      REAL SFCSPD
      REAL SFCTMP
      REAL SHDFAC
      REAL SHDMIN
      REAL SH2O(NSOIL)
      REAL SLDPTH(NSOIL)
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL SMC(NSOIL)
      REAL SNEQV
      REAL SNCOVR
      REAL SNOWH
      REAL SN_NEW
      REAL SLOPE
      REAL SNUP
      REAL SALP
      REAL SNOALB
      REAL STC(NSOIL)
      REAL SNOMLT
      REAL SOLDN
      REAL SOILM
      REAL SOILW
      REAL SOILT
      REAL SOILWM
      REAL SOILWW
      REAL SOILTM
      REAL SOILTW     
      REAL T1
      REAL T1V
      REAL T24
      REAL T2V
      REAL TBOT
      REAL TH2
      REAL TH2V
      REAL TOPT
      REAL TFREEZ
      REAL TSNOW
      REAL XLAI
      REAL ZLVL
      REAL ZBOT
      REAL Z0
      REAL ZSOIL(NSOLD)

C ----------------------------------------------------------------------
C DECLARATIONS - PARAMETERS
C ----------------------------------------------------------------------
      PARAMETER(TFREEZ = 273.15)
      PARAMETER(LVH2O = 2.501E+6)
      PARAMETER(R = 287.04)
      PARAMETER(CP = 1004.5)

      PARAMETER(CFACTR = 0.5)    ! CANOPY WATER PARAMETER
      PARAMETER(CMCMAX = 0.5E-3) ! CANOPY WATER PARAMETER
      PARAMETER(RSMAX = 5000.0)  ! MAX. STOMATAL RESISTANCE
      PARAMETER(TOPT = 298.0)    ! OPTIMUM TRANSPIR AIR TEMP
      PARAMETER(SBETA = -2.0)    ! TO CALC VEG EFFECT ON SHFLX

      REAL, PARAMETER:: FRZK=0.15 ! Ice content threshold in soil

      REAL, PARAMETER:: REFDK = 2.0E-6  ! Reference
      PARAMETER(REFKDT = 3.0)

      PARAMETER(FXEXP = 2.0)     ! BARE SOIL EVAP EXP USED IN DEVAP
      PARAMETER(CSOIL = 2.00E+6) ! SOIL HEAT CAPACITY [J M-3 K-1]

C-- SPECIFY DEPTH[M] OF LOWER BOUNDARY SOIL TEMPERATURE.
C     PARAMETER(ZBOT = -3.0)
      PARAMETER(ZBOT = -8.0)

C-- SPECIFY SNOW DISTRIBUTION SHAPE PARAMETER SALP - SHAPE PARAMETER
C    OF DISTRIBUTION FUNCTION OF SNOW COVER. (from ANDERSON, HYDRO-17)
C    BEST FIT IS WHEN SALP = 2.6
C      PARAMETER(SALP = 2.6)
C - Changed for version 2.6 June 2nd 2003 *
      PARAMETER(SALP = 4.0)

C--  PARAMETER USED TO CALCULATE ROUGHNESS LENGTH OF HEAT.
C      PARAMETER(CZIL = 0.2)
C - Changed for version 2.6 June 2nd 2003 *
      PARAMETER(CZIL = 0.075)

C ----------------------------------------------------------------------
C   INITIALIZATION
C ----------------------------------------------------------------------
      RUNOFF1 = 0.0
      RUNOFF2 = 0.0
      RUNOFF3 = 0.0
      SNOMLT = 0.0

C ----------------------------------------------------------------------
C  THE VARIABLE "ICE" IS A FLAG DENOTING SEA-ICE CASE 
C ----------------------------------------------------------------------
      IF (ICE .EQ. 1) THEN
C ----------------------------------------------------------------------
C SEA-ICE LAYERS ARE EQUAL THICKNESS AND SUM TO 3 METERS
C ----------------------------------------------------------------------
        DO KZ = 1,NSOIL
          ZSOIL(KZ) = -3.*FLOAT(KZ)/FLOAT(NSOIL)
        END DO

      ELSE
C ----------------------------------------------------------------------
C CALCULATE DEPTH (NEGATIVE) BELOW GROUND FROM TOP SKIN SFC TO BOTTOM OF
C   EACH SOIL LAYER.  NOTE:  SIGN OF ZSOIL IS NEGATIVE (DENOTING BELOW
C   GROUND)
C ----------------------------------------------------------------------
        ZSOIL(1) = -SLDPTH(1)
        DO KZ = 2,NSOIL
          ZSOIL(KZ) = -SLDPTH(KZ)+ZSOIL(KZ-1)
        END DO

      ENDIF
         
C ----------------------------------------------------------------------
C CALCULATE ROOT DISTRIBUTION.  PRESENT VERSION ASSUMES UNIFORM
C DISTRIBUTION BASED ON SOIL LAYER DEPTHS.
C ----------------------------------------------------------------------
      DO I = 1,NROOT
        RTDIS(I) = -SLDPTH(I)/ZSOIL(NROOT)
      END DO

C ----------------------------------------------------------------------
C  INITIALIZE PRECIPITATION LOGICALS.
C ----------------------------------------------------------------------
      SNOWNG = .FALSE.
      FRZGRA = .FALSE.

C ----------------------------------------------------------------------
C IF SEA-ICE CASE, ASSIGN DEFAULT WATER-EQUIV SNOW ON TOP
C ----------------------------------------------------------------------
      IF (ICE .EQ. 1) THEN
        SNEQV = 0.01
        SNOWH = 0.05
      ENDIF

C ----------------------------------------------------------------------
C IF INPUT SNOWPACK IS NONZERO, THEN COMPUTE SNOW DENSITY "SNDENS" AND
C   SNOW THERMAL CONDUCTIVITY "SNCOND" (NOTE THAT CSNOW IS A FUNCTION
C   SUBROUTINE)
C ----------------------------------------------------------------------
      IF (SNEQV .EQ. 0.0) THEN
        SNDENS = 0.0
        SNOWH = 0.0
        SNCOND = 1.0
      ELSE
        SNDENS = SNEQV/SNOWH
        SNCOND = CSNOW(SNDENS) 
      ENDIF

C ----------------------------------------------------------------------
C DETERMINE IF IT'S PRECIPITATING AND WHAT KIND OF PRECIP IT IS.
C IF IT'S PRCPING AND THE AIR TEMP IS COLDER THAN 0 C, IT'S SNOWING!
C IF IT'S PRCPING AND THE AIR TEMP IS WARMER THAN 0 C, BUT THE GRND
C TEMP IS COLDER THAN 0 C, FREEZING RAIN IS PRESUMED TO BE FALLING.
C ----------------------------------------------------------------------
      IF (PRCP .GT. 0.0) THEN
        IF (SFCTMP .LE. TFREEZ) THEN
          SNOWNG = .TRUE.
        ELSE
          IF (T1 .LE. TFREEZ) FRZGRA = .TRUE.
        ENDIF
      ENDIF

C ----------------------------------------------------------------------
C IF EITHER PRCP FLAG IS SET, DETERMINE NEW SNOWFALL (CONVERTING PRCP
C RATE FROM KG M-2 S-1 TO A LIQUID EQUIV SNOW DEPTH IN METERS) AND ADD
C IT TO THE EXISTING SNOWPACK.
C NOTE THAT SINCE ALL PRECIP IS ADDED TO SNOWPACK, NO PRECIP INFILTRATES
C INTO THE SOIL SO THAT PRCP1 IS SET TO ZERO.
C ----------------------------------------------------------------------
      IF ( (SNOWNG) .OR. (FRZGRA) ) THEN
        SN_NEW = PRCP * DT * 0.001
        SNEQV = SNEQV + SN_NEW
        PRCP1 = 0.0
C ----------------------------------------------------------------------
C UPDATE SNOW DENSITY BASED ON NEW SNOWFALL, USING OLD AND NEW SNOW.
C UPDATE SNOW THERMAL CONDUCTIVITY
C ----------------------------------------------------------------------
        CALL SNOW_NEW (SFCTMP,SN_NEW,SNOWH,SNDENS)
        SNCOND = CSNOW (SNDENS) 

      ELSE
C ----------------------------------------------------------------------
C PRECIP IS LIQUID (RAIN), HENCE SAVE IN THE PRECIP VARIABLE THAT
C LATER CAN WHOLELY OR PARTIALLY INFILTRATE THE SOIL (ALONG WITH 
C ANY CANOPY "DRIP" ADDED TO THIS LATER)
C ----------------------------------------------------------------------
        PRCP1 = PRCP

      ENDIF

C ----------------------------------------------------------------------
C DETERMINE SNOWCOVER AND ALBEDO OVER LAND.
C ----------------------------------------------------------------------
      IF (ICE .EQ. 0) THEN

C ----------------------------------------------------------------------
C IF SNOW DEPTH=0, SET SNOW FRACTION=0, ALBEDO=SNOW FREE ALBEDO.
C ----------------------------------------------------------------------
        IF (SNEQV .EQ. 0.0) THEN
          SNCOVR = 0.0
          ALBEDO = ALB

        ELSE
C ----------------------------------------------------------------------
C DETERMINE SNOW FRACTIONAL COVERAGE.
C DETERMINE SURFACE ALBEDO MODIFICATION DUE TO SNOWDEPTH STATE.
C ----------------------------------------------------------------------
          CALL SNFRAC (SNEQV,SNUP,SALP,SNOWH,SNCOVR)
          CALL ALCALC (ALB,SNOALB,SHDFAC,SHDMIN,SNCOVR,TSNOW,ALBEDO)
        ENDIF

      ELSE
C ----------------------------------------------------------------------
C SNOW COVER, ALBEDO OVER SEA-ICE
C ----------------------------------------------------------------------
        SNCOVR = 1.0
C   changed in version 2.6 on June 2nd 2003
        ALBEDO = 0.60
        ALBEDO = 0.65
      ENDIF

C ----------------------------------------------------------------------
C THERMAL CONDUCTIVITY FOR SEA-ICE CASE
C ----------------------------------------------------------------------
      IF (ICE .EQ. 1) THEN
        DF1 = 2.2

      ELSE
C ----------------------------------------------------------------------
C NEXT CALCULATE THE SUBSURFACE HEAT FLUX, WHICH FIRST REQUIRES
C CALCULATION OF THE THERMAL DIFFUSIVITY.  TREATMENT OF THE
C LATTER FOLLOWS THAT ON PAGES 148-149 FROM "HEAT TRANSFER IN 
C COLD CLIMATES", BY V. J. LUNARDINI (PUBLISHED IN 1981 
C BY VAN NOSTRAND REINHOLD CO.) I.E. TREATMENT OF TWO CONTIGUOUS 
C "PLANE PARALLEL" MEDIUMS (NAMELY HERE THE FIRST SOIL LAYER 
C AND THE SNOWPACK LAYER, IF ANY). THIS DIFFUSIVITY TREATMENT 
C BEHAVES WELL FOR BOTH ZERO AND NONZERO SNOWPACK, INCLUDING THE 
C LIMIT OF VERY THIN SNOWPACK.  THIS TREATMENT ALSO ELIMINATES
C THE NEED TO IMPOSE AN ARBITRARY UPPER BOUND ON SUBSURFACE 
C HEAT FLUX WHEN THE SNOWPACK BECOMES EXTREMELY THIN.
C ----------------------------------------------------------------------
C FIRST CALCULATE THERMAL DIFFUSIVITY OF TOP SOIL LAYER, USING
C BOTH THE FROZEN AND LIQUID SOIL MOISTURE, FOLLOWING THE 
C SOIL THERMAL DIFFUSIVITY FUNCTION OF PETERS-LIDARD ET AL.
C (1998,JAS, VOL 55, 1209-1224), WHICH REQUIRES THE SPECIFYING
C THE QUARTZ CONTENT OF THE GIVEN SOIL CLASS (SEE ROUTINE REDPRM)
C ----------------------------------------------------------------------
        CALL TDFCND (DF1,SMC(1),QUARTZ,SMCMAX,SH2O(1))

C ----------------------------------------------------------------------
C NEXT ADD SUBSURFACE HEAT FLUX REDUCTION EFFECT FROM THE 
C OVERLYING GREEN CANOPY, ADAPTED FROM SECTION 2.1.2 OF 
C PETERS-LIDARD ET AL. (1997, JGR, VOL 102(D4))
C ----------------------------------------------------------------------
        DF1 = DF1 * EXP(SBETA*SHDFAC)
      ENDIF

C ----------------------------------------------------------------------
C FINALLY "PLANE PARALLEL" SNOWPACK EFFECT FOLLOWING 
C V.J. LINARDINI REFERENCE CITED ABOVE. NOTE THAT DTOT IS
C COMBINED DEPTH OF SNOWDEPTH AND THICKNESS OF FIRST SOIL LAYER
C ----------------------------------------------------------------------
      DSOIL = -(0.5 * ZSOIL(1))

      IF (SNEQV .EQ. 0.) THEN
        SSOIL = DF1 * (T1 - STC(1) ) / DSOIL
      ELSE
        DTOT = SNOWH + DSOIL
        FRCSNO = SNOWH/DTOT
        FRCSOI = DSOIL/DTOT

C 1. HARMONIC MEAN (SERIES FLOW)
C        DF1 = (SNCOND*DF1)/(FRCSOI*SNCOND+FRCSNO*DF1)
        DF1H = (SNCOND*DF1)/(FRCSOI*SNCOND+FRCSNO*DF1)

C 2. ARITHMETIC MEAN (PARALLEL FLOW)
C        DF1 = FRCSNO*SNCOND + FRCSOI*DF1
        DF1A = FRCSNO*SNCOND + FRCSOI*DF1
C
C 3. GEOMETRIC MEAN (INTERMEDIATE BETWEEN HARMONIC AND ARITHMETIC MEAN)
C        DF1 = (SNCOND**FRCSNO)*(DF1**FRCSOI)
C TEST - MBEK, 10 Jan 2002
C weigh DF by snow fraction
c        DF1 = DF1H*SNCOVR + DF1A*(1.0-SNCOVR)
c        DF1 = DF1H*SNCOVR + DF1*(1.0-SNCOVR)
        DF1 = DF1A*SNCOVR + DF1*(1.0-SNCOVR)

C ----------------------------------------------------------------------
C CALCULATE SUBSURFACE HEAT FLUX, SSOIL, FROM FINAL THERMAL DIFFUSIVITY
C OF SURFACE MEDIUMS, DF1 ABOVE, AND SKIN TEMPERATURE AND TOP 
C MID-LAYER SOIL TEMPERATURE
C ----------------------------------------------------------------------
        SSOIL = DF1 * (T1 - STC(1) ) / DTOT
      ENDIF

C ----------------------------------------------------------------------
C DETERMINE SURFACE ROUGHNESS OVER SNOWPACK USING SNOW CONDITION FROM
C THE PREVIOUS TIMESTEP.
C ----------------------------------------------------------------------
      IF (SNCOVR .GT. 0.) THEN
        CALL SNOWZ0 (SNCOVR,Z0)
      ENDIF

C ----------------------------------------------------------------------
C NEXT CALL ROUTINE SFCDIF TO CALCULATE THE SFC EXCHANGE COEF (CH) FOR
C HEAT AND MOISTURE.
C
C NOTE !!!
C COMMENT OUT CALL SFCDIF, IF SFCDIF ALREADY CALLED IN CALLING PROGRAM
C (SUCH AS IN COUPLED ATMOSPHERIC MODEL).
C
C NOTE !!!
C DO NOT CALL SFCDIF UNTIL AFTER ABOVE CALL TO REDPRM, IN CASE
C ALTERNATIVE VALUES OF ROUGHNESS LENGTH (Z0) AND ZILINTINKEVICH COEF
C (CZIL) ARE SET THERE VIA NAMELIST I/O.
C
C NOTE !!!
C ROUTINE SFCDIF RETURNS A CH THAT REPRESENTS THE WIND SPD TIMES THE
C "ORIGINAL" NONDIMENSIONAL "Ch" TYPICAL IN LITERATURE.  HENCE THE CH
C RETURNED FROM SFCDIF HAS UNITS OF M/S.  THE IMPORTANT COMPANION
C COEFFICIENT OF CH, CARRIED HERE AS "RCH", IS THE CH FROM SFCDIF TIMES
C AIR DENSITY AND PARAMETER "CP".  "RCH" IS COMPUTED IN "CALL PENMAN".
C RCH RATHER THAN CH IS THE COEFF USUALLY INVOKED LATER IN EQNS.
C
C NOTE !!!
C SFCDIF ALSO RETURNS THE SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM, CM,
C ALSO KNOWN AS THE SURFACE DRAGE COEFFICIENT, BUT CM IS NOT USED HERE.
C ----------------------------------------------------------------------
C CALC VIRTUAL TEMPS AND VIRTUAL POTENTIAL TEMPS NEEDED BY SUBROUTINES
C SFCDIF AND PENMAN.
C ----------------------------------------------------------------------
      T2V = SFCTMP * (1.0 + 0.61 * Q2 )
C ----------------------------------------------------------------------
C COMMENT OUT BELOW 2 LINES IF CALL SFCDIF IS COMMENTED OUT, I.E. IN THE
C COUPLED MODEL.
C ----------------------------------------------------------------------
      T1V = T1 * (1.0 + 0.61 * Q2)
      TH2V = TH2 * (1.0 + 0.61 * Q2)

      CALL SFCDIF (ZLVL,Z0,T1V,TH2V,SFCSPD,CZIL,CM,CH)

C ----------------------------------------------------------------------
C CALCULATE TOTAL DOWNWARD RADIATION (SOLAR PLUS LONGWAVE) NEEDED IN
C PENMAN EP SUBROUTINE THAT FOLLOWS
C ----------------------------------------------------------------------
      FDOWN = SOLDN*(1.0-ALBEDO) + LWDN

C ----------------------------------------------------------------------
C CALL PENMAN SUBROUTINE TO CALCULATE POTENTIAL EVAPORATION (ETP), AND
C OTHER PARTIAL PRODUCTS AND SUMS SAVE IN COMMON/RITE FOR LATER
C CALCULATIONS.
C ----------------------------------------------------------------------
       CALL PENMAN (SFCTMP,SFCPRS,CH,T2V,TH2,PRCP,FDOWN,T24,SSOIL,
     &              Q2,Q2SAT,ETP,RCH,EPSCA,RR,SNOWNG,FRZGRA,
     &              DQSDT2,FLX2)

C ----------------------------------------------------------------------
C CALL CANRES TO CALCULATE THE CANOPY RESISTANCE AND CONVERT IT INTO PC
C IF NONZERO GREENNESS FRACTION
C ----------------------------------------------------------------------
      IF (SHDFAC .GT. 0.) THEN
      
C ----------------------------------------------------------------------
C  FROZEN GROUND EXTENSION: TOTAL SOIL WATER "SMC" WAS REPLACED 
C  BY UNFROZEN SOIL WATER "SH2O" IN CALL TO CANRES BELOW
C ----------------------------------------------------------------------
        CALL CANRES (SOLDN,CH,SFCTMP,Q2,SFCPRS,SH2O,ZSOIL,NSOIL,
     &               SMCWLT,SMCREF,RSMIN,RC,PC,NROOT,Q2SAT,DQSDT2,
     &               TOPT,RSMAX,RGL,HS,XLAI,
     &               RCS,RCT,RCQ,RCSOIL)

      ENDIF

C ----------------------------------------------------------------------
C NOW DECIDE MAJOR PATHWAY BRANCH TO TAKE DEPENDING ON WHETHER SNOWPACK
C EXISTS OR NOT:
C ----------------------------------------------------------------------
      ESNOW = 0.0
      IF (SNEQV .EQ. 0.0) THEN
        CALL NOPAC (ETP,ETA,PRCP,SMC,SMCMAX,SMCWLT,
     &     	    SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,SHDFAC,
     &     	    SBETA,Q2,T1,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,
     &     	    STC,EPSCA,BEXP,PC,RCH,RR,CFACTR,
     &     	    SH2O,SLOPE,KDT,FRZX,PSISAT,ZSOIL,
     &     	    DKSAT,DWSAT,TBOT,ZBOT,RUNOFF1,RUNOFF2,
     &     	    RUNOFF3,EDIR,EC,ET,ETT,NROOT,ICE,RTDIS,
     &     	    QUARTZ,FXEXP,CSOIL,
     &     	    BETA,DRIP,DEW,FLX1,FLX2,FLX3)
      ELSE
        CALL SNOPAC (ETP,ETA,PRCP,PRCP1,SNOWNG,SMC,SMCMAX,SMCWLT,
     &               SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,
     &               SBETA,DF1,
     &               Q2,T1,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,STC,EPSCA,
     &               SFCPRS,BEXP,PC,RCH,RR,CFACTR,SNCOVR,SNEQV,SNDENS,
     &               SNOWH,SH2O,SLOPE,KDT,FRZX,PSISAT,SNUP,
     &               ZSOIL,DWSAT,DKSAT,TBOT,ZBOT,SHDFAC,RUNOFF1,
     &               RUNOFF2,RUNOFF3,EDIR,EC,ET,ETT,NROOT,SNOMLT,
     &               ICE,RTDIS,QUARTZ,FXEXP,CSOIL,
     &               BETA,DRIP,DEW,FLX1,FLX2,FLX3)
        ESNOW = ETA
      ENDIF

C ----------------------------------------------------------------------
C   PREPARE SENSIBLE HEAT (H) FOR RETURN TO PARENT MODEL
C ----------------------------------------------------------------------
      SHEAT = -(CH * CP * SFCPRS)/(R * T2V) * ( TH2 - T1 )
          
C ----------------------------------------------------------------------
C  CONVERT UNITS AND/OR SIGN OF TOTAL EVAP (ETA), POTENTIAL EVAP (ETP),
C  SUBSURFACE HEAT FLUX (S), AND RUNOFFS FOR WHAT PARENT MODEL EXPECTS
C  CONVERT ETA FROM KG M-2 S-1 TO W M-2
C ----------------------------------------------------------------------
      ETA = ETA*LVH2O
      ETA = ETA*LVH2O
      ETP = ETP*LVH2O

C ----------------------------------------------------------------------
C CONVERT THE SIGN OF SOIL HEAT FLUX SO THAT:
C   SSOIL>0: WARM THE SURFACE  (NIGHT TIME)
C   SSOIL<0: COOL THE SURFACE  (DAY TIME)
C ----------------------------------------------------------------------
      SSOIL = -1.0*SSOIL      

C ----------------------------------------------------------------------
C  CONVERT RUNOFF3 (INTERNAL LAYER RUNOFF FROM SUPERSAT) FROM M TO M S-1
C  AND ADD TO SUBSURFACE RUNOFF/DRAINAGE/BASEFLOW
C ----------------------------------------------------------------------
      RUNOFF3 = RUNOFF3/DT
      RUNOFF2 = RUNOFF2+RUNOFF3

C ----------------------------------------------------------------------
C TOTAL COLUMN SOIL MOISTURE IN METERS (SOILM)  ROOT-ZONE 
C ----------------------------------------------------------------------
      SOILM = -1.0*SMC(1)*ZSOIL(1)
      DO K = 2,NSOIL
        SOILM = SOILM+SMC(K)*(ZSOIL(K-1)-ZSOIL(K))
      END DO

C ----------------------------------------------------------------------
C ROOT-ZONE SOIL MOISTURE AVAILABILITY (FRACTION) RELATIVE
C TO POROSITY/SATURATION (SOILW; aka, MSTAVRZ)
C ----------------------------------------------------------------------
      SOILWM = -1.0*(SMCMAX-SMCWLT)*ZSOIL(1)
      SOILWW = -1.0*(SMC(1)-SMCWLT)*ZSOIL(1)
      DO K = 2,NROOT
        SOILWM = SOILWM+(SMCMAX-SMCWLT)*(ZSOIL(K-1)-ZSOIL(K))
        SOILWW = SOILWW+(SMC(K)-SMCWLT)*(ZSOIL(K-1)-ZSOIL(K))
      END DO
      SOILW = SOILWW/SOILWM

C ----------------------------------------------------------------------
C TOTAL COL SOIL MOISTURE AVAIL RELATIVE TO POROSITY/SATURATION (SOILT) 
C  (aka, MSTAVTOT)
C ----------------------------------------------------------------------
      SOILTM = -1.0*(SMCMAX-SMCWLT)*ZSOIL(1)
      SOILTW = -1.0*(SMC(1)-SMCWLT)*ZSOIL(1)
      DO K = 2,NSOIL
        SOILTM = SOILTM+(SMCMAX-SMCWLT)*(ZSOIL(K-1)-ZSOIL(K))
        SOILTW = SOILTW+(SMC(K)-SMCWLT)*(ZSOIL(K-1)-ZSOIL(K))
      END DO
      SOILT = SOILTW/SOILTM

C ----------------------------------------------------------------------
C END SUBROUTINE SFLX
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 11. PHYSICS SUBROUTINE ==>  SUBROUTINE SHFLX  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SHFLX (SSOIL,STC,SMC,SMCMAX,NSOIL,T1,DT,YY,ZZ1,ZSOIL,
     &                  TBOT,ZBOT,SMCWLT,PSISAT,SH2O,BEXP,F1,DF1,ICE,
     &                  QUARTZ,CSOIL)
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SHFLX
C ----------------------------------------------------------------------
C UPDATE THE TEMPERATURE STATE OF THE SOIL COLUMN BASED ON THE THERMAL
C DIFFUSION EQUATION AND UPDATE THE FROZEN SOIL MOISTURE CONTENT BASED
C ON THE TEMPERATURE.
C ----------------------------------------------------------------------
      INTEGER, PARAMETER:: NSOLD=20

      INTEGER I
      INTEGER ICE
      INTEGER IFRZ
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

      REAL BEXP
      REAL CSOIL
      REAL DF1
      REAL DT
      REAL F1
      REAL PSISAT
      REAL QUARTZ
      REAL RHSTS(NSOLD)
      REAL SSOIL
      REAL SH2O(NSOIL)
      REAL SMC(NSOIL)
      REAL SMCMAX
      REAL SMCWLT
      REAL STC(NSOIL)
      REAL STCF(NSOLD)
      REAL T0
      REAL T1
      REAL TBOT
      REAL YY
      REAL ZBOT
      REAL ZSOIL(NSOIL)
      REAL ZZ1

      PARAMETER(T0 = 273.15)

C ----------------------------------------------------------------------
C HRT ROUTINE CALCS THE RIGHT HAND SIDE OF THE SOIL TEMP DIF EQN
C ----------------------------------------------------------------------
      IF (ICE.EQ.1) THEN

C ----------------------------------------------------------------------
C SEA-ICE CASE
C ----------------------------------------------------------------------
         CALL HRTICE (RHSTS,STC,NSOIL,ZSOIL,YY,ZZ1,DF1,AI,BI,CI)

         CALL HSTEP (STCF,STC,RHSTS,DT,NSOIL,AI,BI,CI)
         
      ELSE

C ----------------------------------------------------------------------
C LAND-MASS CASE
C ----------------------------------------------------------------------
         CALL HRT (RHSTS,STC,SMC,SMCMAX,NSOIL,ZSOIL,YY,ZZ1,TBOT,
     &             ZBOT,PSISAT,SH2O,DT,
     &             BEXP,F1,DF1,QUARTZ,CSOIL,AI,BI,CI)
         
         CALL HSTEP (STCF,STC,RHSTS,DT,NSOIL,AI,BI,CI)

      ENDIF

      DO I = 1,NSOIL
         STC(I) = STCF(I)
      END DO
      
C ----------------------------------------------------------------------
C IN THE NO SNOWPACK CASE (VIA ROUTINE NOPAC BRANCH,) UPDATE THE GRND
C (SKIN) TEMPERATURE HERE IN RESPONSE TO THE UPDATED SOIL TEMPERATURE 
C PROFILE ABOVE.  (NOTE: INSPECTION OF ROUTINE SNOPAC SHOWS THAT T1
C BELOW IS A DUMMY VARIABLE ONLY, AS SKIN TEMPERATURE IS UPDATED
C DIFFERENTLY IN ROUTINE SNOPAC) 
C ----------------------------------------------------------------------
      T1 = (YY + (ZZ1 - 1.0) * STC(1)) / ZZ1

C ----------------------------------------------------------------------
C CALCULATE SURFACE SOIL HEAT FLUX
C ----------------------------------------------------------------------
      SSOIL = DF1 * (STC(1) - T1) / (0.5 * ZSOIL(1))

C ----------------------------------------------------------------------
C END SUBROUTINE SHFLX
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 12. PHYSICS SUBROUTINE ==>  SUBROUTINE SMFLX  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SMFLX (ETA1,SMC,NSOIL,CMC,ETP1,DT,PRCP1,ZSOIL,
     &                  SH2O,SLOPE,KDT,FRZFACT,
     &                  SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
     &                  SMCREF,SHDFAC,CMCMAX,
     &                  SMCDRY,CFACTR,RUNOFF1,RUNOFF2,RUNOFF3,
     &                  EDIR,EC,ET,ETT,SFCTMP,Q2,NROOT,RTDIS,FXEXP,
     &                  DRIP)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SMFLX
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE FLUX.  THE SOIL MOISTURE CONTENT (SMC - A PER
C UNIT VOLUME MEASUREMENT) IS A DEPENDENT VARIABLE THAT IS UPDATED WITH
C PROGNOSTIC EQNS. THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED.
C FROZEN GROUND VERSION:  NEW STATES ADDED: SH2O, AND FROZEN GROUND
C CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER I
      INTEGER K
      INTEGER NSOIL
      INTEGER NROOT

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)

      REAL BEXP
      REAL BETA
      REAL CFACTR
      REAL CMC
      REAL CMC2MS
      REAL CMCMAX
      REAL DEVAP
      REAL DKSAT
      REAL DRIP
      REAL DT
      REAL DUMMY
      REAL DWSAT
      REAL EC
      REAL EDIR
      REAL ET(NSOIL)
      REAL ETA1
      REAL ETP1
      REAL ETT
      REAL EXCESS
      REAL FRZFACT
      REAL FXEXP
      REAL KDT
      REAL PC
      REAL PCPDRP
      REAL PRCP1
      REAL Q2
      REAL RHSCT
      REAL RHSTT(NSOLD)
      REAL RTDIS(NSOIL)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RUNOFF3
      REAL SFCTMP
      REAL SHDFAC
      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SICE(NSOLD)
      REAL SH2OA(NSOLD)
      REAL SH2OFG(NSOLD)
      REAL SLOPE
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL TRHSCT
      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C TEMPERATURE CRITERIA FOR SNOWFALL TFREEZ SHOULD HAVE SAME VALUE AS IN
C SFLX.F (MAIN SUBROUTINE)
C ----------------------------------------------------------------------
      REAL TFREEZ
      PARAMETER(TFREEZ = 273.15)

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE IF THE POTENTIAL EVAPOTRANSPIRATION IS
C GREATER THAN ZERO.
C ----------------------------------------------------------------------
      DUMMY = 0.
      EDIR = 0.
      EC = 0.
      ETT = 0.
      DO K = 1,NSOIL
        ET(K) = 0.
      END DO

      IF (ETP1 .GT. 0.0) THEN

C ----------------------------------------------------------------------
C RETRIEVE DIRECT EVAPORATION FROM SOIL SURFACE.  CALL THIS FUNCTION
C ONLY IF VEG COVER NOT COMPLETE.
C FROZEN GROUND VERSION:  SH2O STATES REPLACE SMC STATES.
C ----------------------------------------------------------------------
        IF (SHDFAC .LT. 1.) THEN
          EDIR = DEVAP(ETP1,SH2O(1),ZSOIL(1),SHDFAC,SMCMAX,
     &                 BEXP,DKSAT,DWSAT,SMCDRY,SMCREF,SMCWLT,FXEXP)
        ENDIF

C ----------------------------------------------------------------------
C INITIALIZE PLANT TOTAL TRANSPIRATION, RETRIEVE PLANT TRANSPIRATION,
C AND ACCUMULATE IT FOR ALL SOIL LAYERS.
C ----------------------------------------------------------------------
        IF (SHDFAC.GT.0.0) THEN

          CALL TRANSP (ET,NSOIL,ETP1,SH2O,CMC,ZSOIL,SHDFAC,SMCWLT,
     &                 CMCMAX,PC,CFACTR,SMCREF,SFCTMP,Q2,NROOT,RTDIS)

          DO K = 1,NSOIL
            ETT = ETT + ET ( K )
          END DO

C ----------------------------------------------------------------------
C CALCULATE CANOPY EVAPORATION.
C IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR CMC=0.0.
C ----------------------------------------------------------------------
          IF (CMC .GT. 0.0) THEN
            EC = SHDFAC * ( ( CMC / CMCMAX ) ** CFACTR ) * ETP1
          ELSE
            EC = 0.0
          ENDIF

C ----------------------------------------------------------------------
C EC SHOULD BE LIMITED BY THE TOTAL AMOUNT OF AVAILABLE WATER ON THE
C CANOPY.  -F.CHEN, 18-OCT-1994
C ----------------------------------------------------------------------
          CMC2MS = CMC / DT
          EC = MIN ( CMC2MS, EC )
        ENDIF
      ENDIF

C ----------------------------------------------------------------------
C TOTAL UP EVAP AND TRANSP TYPES TO OBTAIN ACTUAL EVAPOTRANSP
C ----------------------------------------------------------------------
      ETA1 = EDIR + ETT + EC
      
C ----------------------------------------------------------------------
C COMPUTE THE RIGHT HAND SIDE OF THE CANOPY EQN TERM ( RHSCT )
C ----------------------------------------------------------------------
      RHSCT = SHDFAC * PRCP1 - EC

C ----------------------------------------------------------------------
C CONVERT RHSCT (A RATE) TO TRHSCT (AN AMOUNT) AND ADD IT TO EXISTING
C CMC.  IF RESULTING AMT EXCEEDS MAX CAPACITY, IT BECOMES DRIP AND WILL
C FALL TO THE GRND.
C ----------------------------------------------------------------------
      DRIP = 0.
      TRHSCT = DT * RHSCT
      EXCESS = CMC + TRHSCT
      IF (EXCESS .GT. CMCMAX) DRIP = EXCESS - CMCMAX

C ----------------------------------------------------------------------
C PCPDRP IS THE COMBINED PRCP1 AND DRIP (FROM CMC) THAT GOES INTO THE
C SOIL
C ----------------------------------------------------------------------
      PCPDRP = (1. - SHDFAC) * PRCP1 + DRIP / DT

C ----------------------------------------------------------------------
C STORE ICE CONTENT AT EACH SOIL LAYER BEFORE CALLING SRT & SSTEP
C ----------------------------------------------------------------------
      DO I = 1,NSOIL
        SICE(I) = SMC(I) - SH2O(I)
      END DO
            
C ----------------------------------------------------------------------
C CALL SUBROUTINES SRT AND SSTEP TO SOLVE THE SOIL MOISTURE
C TENDENCY EQUATIONS. 
C
C IF THE INFILTRATING PRECIP RATE IS NONTRIVIAL,
C   (WE CONSIDER NONTRIVIAL TO BE A PRECIP TOTAL OVER THE TIME STEP 
C    EXCEEDING ONE ONE-THOUSANDTH OF THE WATER HOLDING CAPACITY OF 
C    THE FIRST SOIL LAYER)
C THEN CALL THE SRT/SSTEP SUBROUTINE PAIR TWICE IN THE MANNER OF 
C   TIME SCHEME "F" (IMPLICIT STATE, AVERAGED COEFFICIENT)
C   OF SECTION 2 OF KALNAY AND KANAMITSU (1988, MWR, VOL 116, 
C   PAGES 1945-1958)TO MINIMIZE 2-DELTA-T OSCILLATIONS IN THE 
C   SOIL MOISTURE VALUE OF THE TOP SOIL LAYER THAT CAN ARISE BECAUSE
C   OF THE EXTREME NONLINEAR DEPENDENCE OF THE SOIL HYDRAULIC 
C   DIFFUSIVITY COEFFICIENT AND THE HYDRAULIC CONDUCTIVITY ON THE
C   SOIL MOISTURE STATE
C OTHERWISE CALL THE SRT/SSTEP SUBROUTINE PAIR ONCE IN THE MANNER OF
C   TIME SCHEME "D" (IMPLICIT STATE, EXPLICIT COEFFICIENT) 
C   OF SECTION 2 OF KALNAY AND KANAMITSU
C PCPDRP IS UNITS OF KG/M**2/S OR MM/S, ZSOIL IS NEGATIVE DEPTH IN M 
C ----------------------------------------------------------------------
C      IF ( PCPDRP .GT. 0.0 ) THEN
      IF ( (PCPDRP*DT) .GT. (0.001*1000.0*(-ZSOIL(1))*SMCMAX) ) THEN

C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C SMC STATES REPLACED BY SH2O STATES IN SRT SUBR.  SH2O & SICE STATES
C INCLUDED IN SSTEP SUBR.  FROZEN GROUND CORRECTION FACTOR, FRZFACT
C ADDED.  ALL WATER BALANCE CALCULATIONS USING UNFROZEN WATER
C ----------------------------------------------------------------------
        CALL SRT (RHSTT,EDIR,ET,SH2O,SH2O,NSOIL,PCPDRP,ZSOIL,
     &            DWSAT,DKSAT,SMCMAX,BEXP,RUNOFF1, 
     &            RUNOFF2,DT,SMCWLT,SLOPE,KDT,FRZFACT,SICE,AI,BI,CI)
         
        CALL SSTEP (SH2OFG,SH2O,DUMMY,RHSTT,RHSCT,DT,NSOIL,SMCMAX,
     &              CMCMAX,RUNOFF3,ZSOIL,SMC,SICE,AI,BI,CI)
         
        DO K = 1,NSOIL
          SH2OA(K) = (SH2O(K) + SH2OFG(K)) * 0.5
        END DO
        
        CALL SRT (RHSTT,EDIR,ET,SH2O,SH2OA,NSOIL,PCPDRP,ZSOIL,
     &            DWSAT,DKSAT,SMCMAX,BEXP,RUNOFF1,
     &            RUNOFF2,DT,SMCWLT,SLOPE,KDT,FRZFACT,SICE,AI,BI,CI)
         
        CALL SSTEP (SH2O,SH2O,CMC,RHSTT,RHSCT,DT,NSOIL,SMCMAX,
     &              CMCMAX,RUNOFF3,ZSOIL,SMC,SICE,AI,BI,CI)
         
      ELSE
         
        CALL SRT (RHSTT,EDIR,ET,SH2O,SH2O,NSOIL,PCPDRP,ZSOIL,
     &            DWSAT,DKSAT,SMCMAX,BEXP,RUNOFF1,
     &            RUNOFF2,DT,SMCWLT,SLOPE,KDT,FRZFACT,SICE,AI,BI,CI)

        CALL SSTEP (SH2O,SH2O,CMC,RHSTT,RHSCT,DT,NSOIL,SMCMAX,
     &              CMCMAX,RUNOFF3,ZSOIL,SMC,SICE,AI,BI,CI)
         
      ENDIF
      
c      RUNOF = RUNOFF

C ----------------------------------------------------------------------
C END SUBROUTINE SMFLX
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 13. PHYSICS SUBROUTINE ==>  SUBROUTINE SNFRAC CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SNFRAC (SNEQV,SNUP,SALP,SNOWH,SNCOVR)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SNFRAC
C ----------------------------------------------------------------------
C CALCULATE SNOW FRACTION (0 -> 1)
C SNEQV   SNOW WATER EQUIVALENT (M)
C SNUP    THRESHOLD SNEQV DEPTH ABOVE WHICH SNCOVR=1
C SALP    TUNING PARAMETER
C SNCOVR  FRACTIONAL SNOW COVER
C ----------------------------------------------------------------------
      REAL SNEQV, SNUP, SALP, SNCOVR, RSNOW, Z0N, SNOWH
      
C ----------------------------------------------------------------------
C SNUP IS VEG-CLASS DEPENDENT SNOWDEPTH THRESHHOLD (SET IN ROUTINE
C REDPRM) ABOVE WHICH SNOCVR=1.
C ----------------------------------------------------------------------
          IF (SNEQV .LT. SNUP) THEN
            RSNOW = SNEQV/SNUP
            SNCOVR = 1. - ( EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))
          ELSE
            SNCOVR = 1.0
          ENDIF

          Z0N=0.035 
C     FORMULATION OF DICKINSON ET AL. 1986

C        SNCOVR=SNOWH/(SNOWH + 5*Z0N)

C     FORMULATION OF MARSHALL ET AL. 1994
C        SNCOVR=SNEQV/(SNEQV + 2*Z0N)

C ----------------------------------------------------------------------
C END SUBROUTINE SNFRAC
C ----------------------------------------------------------------------
      RETURN
      END


      FUNCTION SNKSRC (TAVG,SMC,SH2O,ZSOIL,NSOIL,
     &                 SMCMAX,PSISAT,BEXP,DT,K,QTOT) 
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C FUNCTION SNKSRC
C ----------------------------------------------------------------------
C CALCULATE SINK/SOURCE TERM OF THE TERMAL DIFFUSION EQUATION. (SH2O) IS
C AVAILABLE LIQUED WATER.
C ----------------------------------------------------------------------
      INTEGER K
      INTEGER NSOIL
      
      REAL BEXP
      REAL DF
      REAL DH2O
      REAL DT
      REAL DZ
      REAL DZH
      REAL FREE
      REAL FRH2O
      REAL HLICE
      REAL PSISAT
      REAL QTOT
      REAL SH2O
      REAL SMC
      REAL SMCMAX
      REAL SNKSRC
      REAL T0
      REAL TAVG
      REAL TDN
      REAL TM
      REAL TUP
      REAL TZ
      REAL X0
      REAL XDN
      REAL XH2O
      REAL XUP
      REAL ZSOIL (NSOIL)

      PARAMETER(DH2O = 1.0000E3)
      PARAMETER(HLICE = 3.3350E5)
      PARAMETER(T0 = 2.7315E2)
      
      IF (K .EQ. 1) THEN
        DZ = -ZSOIL(1)
      ELSE
        DZ = ZSOIL(K-1)-ZSOIL(K)
      ENDIF

C ----------------------------------------------------------------------
C VIA FUNCTION FRH2O, COMPUTE POTENTIAL OR 'EQUILIBRIUM' UNFROZEN
C SUPERCOOLED FREE WATER FOR GIVEN SOIL TYPE AND SOIL LAYER TEMPERATURE.
C FUNCTION FRH20 INVOKES EQN (17) FROM V. KOREN ET AL (1999, JGR, VOL.
C 104, PG 19573).  (ASIDE:  LATTER EQN IN JOURNAL IN CENTIGRADE UNITS.
C ROUTINE FRH2O USE FORM OF EQN IN KELVIN UNITS.)
C ----------------------------------------------------------------------
      FREE = FRH2O(TAVG,SMC,SH2O,SMCMAX,BEXP,PSISAT)

C ----------------------------------------------------------------------
C IN NEXT BLOCK OF CODE, INVOKE EQN 18 OF V. KOREN ET AL (1999, JGR,
C VOL. 104, PG 19573.)  THAT IS, FIRST ESTIMATE THE NEW AMOUNTOF LIQUID
C WATER, 'XH2O', IMPLIED BY THE SUM OF (1) THE LIQUID WATER AT THE BEGIN
C OF CURRENT TIME STEP, AND (2) THE FREEZE OF THAW CHANGE IN LIQUID
C WATER IMPLIED BY THE HEAT FLUX 'QTOT' PASSED IN FROM ROUTINE HRT.
C SECOND, DETERMINE IF XH2O NEEDS TO BE BOUNDED BY 'FREE' (EQUIL AMT) OR
C IF 'FREE' NEEDS TO BE BOUNDED BY XH2O.
C ----------------------------------------------------------------------
      XH2O = SH2O + QTOT*DT/(DH2O*HLICE*DZ)

C ----------------------------------------------------------------------
C FIRST, IF FREEZING AND REMAINING LIQUID LESS THAN LOWER BOUND, THEN
C REDUCE EXTENT OF FREEZING, THEREBY LETTING SOME OR ALL OF HEAT FLUX
C QTOT COOL THE SOIL TEMP LATER IN ROUTINE HRT.
C ----------------------------------------------------------------------
      IF ( XH2O .LT. SH2O .AND. XH2O .LT. FREE) THEN 
        IF ( FREE .GT. SH2O ) THEN
          XH2O = SH2O
        ELSE
          XH2O = FREE
        ENDIF
      ENDIF
              
C ----------------------------------------------------------------------
C SECOND, IF THAWING AND THE INCREASE IN LIQUID WATER GREATER THAN UPPER
C BOUND, THEN REDUCE EXTENT OF THAW, THEREBY LETTING SOME OR ALL OF HEAT
C FLUX QTOT WARM THE SOIL TEMP LATER IN ROUTINE HRT.
C ----------------------------------------------------------------------
      IF ( XH2O .GT. SH2O .AND. XH2O .GT. FREE )  THEN
        IF ( FREE .LT. SH2O ) THEN
          XH2O = SH2O
        ELSE
          XH2O = FREE
        ENDIF
      ENDIF 

      IF (XH2O .LT. 0.) XH2O = 0.
      IF (XH2O .GT. SMC) XH2O = SMC

C ----------------------------------------------------------------------
C CALCULATE PHASE-CHANGE HEAT SOURCE/SINK TERM FOR USE IN ROUTINE HRT
C AND UPDATE LIQUID WATER TO REFLCET FINAL FREEZE/THAW INCREMENT.
C ----------------------------------------------------------------------
      SNKSRC = -DH2O*HLICE*DZ*(XH2O-SH2O)/DT
      SH2O = XH2O
      
C ----------------------------------------------------------------------
C END FUNCTION SNKSRC
C ----------------------------------------------------------------------
77    RETURN
      END

CCCCC 14. PHYSICS SUBROUTINE ==>  SUBROUTINE SNOPAC CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SNOPAC (ETP,ETA,PRCP,PRCP1,SNOWNG,SMC,SMCMAX,SMCWLT,
     &                   SMCREF,SMCDRY,CMC,CMCMAX,NSOIL,DT,
     &                   SBETA,DF1,
     &                   Q2,T1,SFCTMP,T24,TH2,FDOWN,F1,SSOIL,STC,EPSCA,
     &                   SFCPRS,BEXP,PC,RCH,RR,CFACTR,SNCOVR,ESD,SNDENS,
     &                   SNOWH,SH2O,SLOPE,KDT,FRZFACT,PSISAT,SNUP,
     &                   ZSOIL,DWSAT,DKSAT,TBOT,ZBOT,SHDFAC,RUNOFF1,
     &                   RUNOFF2,RUNOFF3,EDIR,EC,ET,ETT,NROOT,SNOMLT,
     &                   ICE,RTDIS,QUARTZ,FXEXP,CSOIL,
     &                   BETA,DRIP,DEW,FLX1,FLX2,FLX3)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SNOPAC
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES & UPDATE SOIL MOISTURE
C CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN A SNOW PACK IS
C PRESENT.
C ----------------------------------------------------------------------
      INTEGER ICE
      INTEGER NROOT
      INTEGER NSOIL

      LOGICAL*1 SNOWNG

      REAL BEXP
      REAL BETA
      REAL CFACTR
      REAL CMC
      REAL CMCMAX
      REAL CP
      REAL CPH2O
      REAL CPICE
      REAL CSOIL
      REAL DENOM
      REAL DEW
      REAL DF1
      REAL DKSAT
      REAL DRIP
      REAL DSOIL
      REAL DTOT
      REAL DT
      REAL DWSAT
      REAL EC
      REAL EDIR
      REAL EPSCA
      REAL ESD
      REAL ESDMIN
      REAL EXPSNO
      REAL EXPSOI
      REAL ETA
      REAL ETA1
      REAL ETP
      REAL ETP1
      REAL ETP2
      REAL ET(NSOIL)
      REAL ETT
      REAL EX
      REAL EXPFAC
      REAL FDOWN
      REAL FXEXP
      REAL FLX1
      REAL FLX2
      REAL FLX3
      REAL F1
      REAL KDT
      REAL LSUBF
      REAL LSUBC
      REAL LSUBS
      REAL PC
      REAL PRCP
      REAL PRCP1
      REAL Q2
      REAL RCH
      REAL RR
      REAL RTDIS(NSOIL)
      REAL SSOIL
      REAL SBETA
      REAL SSOIL1
      REAL SFCTMP
      REAL SHDFAC
      REAL SIGMA
      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL SNOMLT
      REAL SNOWH
      REAL STC(NSOIL)
      REAL T1
      REAL T11
      REAL T12
      REAL T12A
      REAL T12B
      REAL T24
      REAL TBOT
      REAL ZBOT
      REAL TH2
      REAL YY
      REAL ZSOIL(NSOIL)
      REAL ZZ1
      REAL TFREEZ
      REAL SALP
      REAL SFCPRS
      REAL SLOPE
      REAL FRZFACT
      REAL PSISAT
      REAL SNUP
      REAL RUNOFF1
      REAL RUNOFF2
      REAL RUNOFF3
      REAL QUARTZ
      REAL SNDENS
      REAL SNCOND
      REAL RSNOW
      REAL SNCOVR
      REAL QSAT
      REAL ETP3
      REAL SEH
      REAL T14
      REAL CSNOW

      PARAMETER(CP = 1004.5)
      PARAMETER(CPH2O = 4.218E+3)
      PARAMETER(CPICE = 2.106E+3)
      PARAMETER(ESDMIN = 1.E-6)
      PARAMETER(LSUBF = 3.335E+5)
      PARAMETER(LSUBC = 2.501000E+6)
      PARAMETER(LSUBS = 2.83E+6)
      PARAMETER(SIGMA = 5.67E-8)
      PARAMETER(TFREEZ = 273.15)

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE:
C CONVERT POTENTIAL EVAP (ETP) FROM KG M-2 S-1 TO M S-1 AND THEN TO AN
C AMOUNT (M) GIVEN TIMESTEP (DT) AND CALL IT AN EFFECTIVE SNOWPACK
C REDUCTION AMOUNT, ETP2 (M).  THIS IS THE AMOUNT THE SNOWPACK WOULD BE
C REDUCED DUE TO EVAPORATION FROM THE SNOW SFC DURING THE TIMESTEP.
C EVAPORATION WILL PROCEED AT THE POTENTIAL RATE UNLESS THE SNOW DEPTH
C IS LESS THAN THE EXPECTED SNOWPACK REDUCTION.
C IF SEAICE (ICE=1), BETA REMAINS=1.
C ----------------------------------------------------------------------
      PRCP1 = PRCP1*0.001

      ETP2 = ETP * 0.001 * DT
      BETA = 1.0
      IF (ICE .NE. 1) THEN
        IF (ESD .LT. ETP2) THEN
          BETA = ESD / ETP2
        ENDIF
      ENDIF

C ----------------------------------------------------------------------
C IF ETP<0 (DOWNWARD) THEN DEWFALL (=FROSTFALL IN THIS CASE).
C ----------------------------------------------------------------------
      DEW = 0.0
      IF (ETP .LT. 0.0) THEN
        DEW = -ETP * 0.001
      ENDIF

C ----------------------------------------------------------------------
C IF PRECIP IS FALLING, CALCULATE HEAT FLUX FROM SNOW SFC TO NEWLY
C ACCUMULATING PRECIP.  NOTE THAT THIS REFLECTS THE FLUX APPROPRIATE FOR
C THE NOT-YET-UPDATED SKIN TEMPERATURE (T1).  ASSUMES TEMPERATURE OF THE
C SNOWFALL STRIKING THE GOUND IS =SFCTMP (LOWEST MODEL LEVEL AIR TEMP).
C ----------------------------------------------------------------------
      FLX1 = 0.0
      IF (SNOWNG) THEN
        FLX1 = CPICE * PRCP * (T1 - SFCTMP)
      ELSE
        IF (PRCP .GT. 0.0) FLX1 = CPH2O * PRCP * (T1 - SFCTMP)
      ENDIF

C ----------------------------------------------------------------------
C CALCULATE AN 'EFFECTIVE SNOW-GRND SFC TEMP' (T12) BASED ON HEAT FLUXES
C BETWEEN THE SNOW PACK AND THE SOIL AND ON NET RADIATION.
C INCLUDE FLX1 (PRECIP-SNOW SFC) AND FLX2 (FREEZING RAIN LATENT HEAT)
C FLUXES.
C FLX2 REFLECTS FREEZING RAIN LATENT HEAT FLUX USING T1 CALCULATED IN
C PENMAN.
C ----------------------------------------------------------------------
      DSOIL = -(0.5 * ZSOIL(1))
      DTOT = SNOWH + DSOIL
      DENOM = 1.0 + DF1 / (DTOT * RR * RCH)
      T12A = ( (FDOWN-FLX1-FLX2-SIGMA*T24)/RCH
     &       + TH2 - SFCTMP - BETA*EPSCA ) / RR
      T12B = DF1 * STC(1) / (DTOT * RR * RCH)
      T12 = (SFCTMP + T12A + T12B) / DENOM      

C ----------------------------------------------------------------------
C IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS AT OR BELOW FREEZING, NO SNOW
C MELT WILL OCCUR.  SET THE SKIN TEMP TO THIS EFFECTIVE TEMP.  REDUCE
C (BY SUBLIMINATION ) OR INCREASE (BY FROST) THE DEPTH OF THE SNOWPACK,
C DEPENDING ON SIGN OF ETP.
C UPDATE SOIL HEAT FLUX (SSOIL) USING NEW SKIN TEMPERATURE (T1)
C SINCE NO SNOWMELT, SET ACCUMULATED SNOWMELT TO ZERO, SET 'EFFECTIVE'
C PRECIP FROM SNOWMELT TO ZERO, SET PHASE-CHANGE HEAT FLUX FROM SNOWMELT
C TO ZERO.
C ----------------------------------------------------------------------
      IF (T12 .LE. TFREEZ) THEN
        T1 = T12
        SSOIL = DF1 * (T1 - STC(1)) / DTOT
        ESD = MAX(0.0, ESD-ETP2)
        FLX3 = 0.0
        EX = 0.0
        SNOMLT = 0.0

      ELSE
C ----------------------------------------------------------------------
C IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS ABOVE FREEZING, SNOW MELT
C WILL OCCUR.  CALL THE SNOW MELT RATE,EX AND AMT, SNOMLT.  REVISE THE
C EFFECTIVE SNOW DEPTH.  REVISE THE SKIN TEMP BECAUSE IT WOULD HAVE CHGD
C DUE TO THE LATENT HEAT RELEASED BY THE MELTING. CALC THE LATENT HEAT
C RELEASED, FLX3. SET THE EFFECTIVE PRECIP, PRCP1 TO THE SNOW MELT RATE,
C EX FOR USE IN SMFLX.  ADJUSTMENT TO T1 TO ACCOUNT FOR SNOW PATCHES.
C CALCULATE QSAT VALID AT FREEZING POINT.  NOTE THAT ESAT (SATURATION
C VAPOR PRESSURE) VALUE OF 6.11E+2 USED HERE IS THAT VALID AT FRZZING
C POINT.  NOTE THAT ETP FROM CALL PENMAN IN SFLX IS IGNORED HERE IN
C FAVOR OF BULK ETP OVER 'OPEN WATER' AT FREEZING TEMP.
C UPDATE SOIL HEAT FLUX (S) USING NEW SKIN TEMPERATURE (T1)
C ----------------------------------------------------------------------
        T1 = TFREEZ * SNCOVR + T12 * (1.0 - SNCOVR)
        QSAT = (0.622*6.11E2)/(SFCPRS-0.378*6.11E2)
        ETP = RCH*(QSAT-Q2)/CP
        ETP2 = ETP*0.001*DT
        BETA = 1.0
        SSOIL = DF1 * (T1 - STC(1)) / DTOT
	
C ----------------------------------------------------------------------
C IF POTENTIAL EVAP (SUBLIMATION) GREATER THAN DEPTH OF SNOWPACK.
C BETA<1
C SNOWPACK HAS SUBLIMATED AWAY, SET DEPTH TO ZERO.
C ----------------------------------------------------------------------
        IF (ESD .LE. ETP2) THEN
          BETA = ESD / ETP2
          ESD = 0.0
          EX = 0.0
          SNOMLT = 0.0
	  
        ELSE
C ----------------------------------------------------------------------
C POTENTIAL EVAP (SUBLIMATION) LESS THAN DEPTH OF SNOWPACK, RETAIN
C   BETA=1.
C SNOWPACK (ESD) REDUCED BY POTENTIAL EVAP RATE
C ETP3 (CONVERT TO FLUX)
C ----------------------------------------------------------------------
          ESD = ESD-ETP2
          ETP3 = ETP*LSUBC
          SEH = RCH*(T1-TH2)
          T14 = T1*T1
          T14 = T14*T14
          FLX3 = FDOWN - FLX1 - FLX2 - SIGMA*T14 - SSOIL - SEH - ETP3
          IF (FLX3 .LE .0.0) FLX3 = 0.0
          EX = FLX3*0.001/LSUBF
C ----------------------------------------------------------------------
C SNOWMELT REDUCTION DEPENDING ON SNOW COVER
C IF SNOW COVER LESS THAN 5% NO SNOWMELT REDUCTION
C ***NOTE:  DOES 'IF' BELOW FAIL TO MATCH THE MELT WATER WITH THE MELT
C           ENERGY?
C ----------------------------------------------------------------------
          IF (SNCOVR .GT. 0.05) EX = EX * SNCOVR
          SNOMLT = EX * DT

C ----------------------------------------------------------------------
C ESDMIN REPRESENTS A SNOWPACK DEPTH THRESHOLD VALUE BELOW WHICH WE
C CHOOSE NOT TO RETAIN ANY SNOWPACK, AND INSTEAD INCLUDE IT IN SNOWMELT.
C ----------------------------------------------------------------------
          IF (ESD-SNOMLT .GE. ESDMIN) THEN
            ESD = ESD - SNOMLT

          ELSE
C ----------------------------------------------------------------------
C SNOWMELT EXCEEDS SNOW DEPTH
C ----------------------------------------------------------------------
            EX = ESD/DT
            FLX3 = EX*1000.0*LSUBF
            SNOMLT = ESD
            ESD = 0.0

          ENDIF

C ----------------------------------------------------------------------
C END OF 'ESD .LE. ETP2' IF-BLOCK
C ----------------------------------------------------------------------
        ENDIF

        PRCP1 = PRCP1 + EX

C ----------------------------------------------------------------------
C END OF 'T12 .LE. TFREEZ' IF-BLOCK
C ----------------------------------------------------------------------
      ENDIF
         
C ----------------------------------------------------------------------
C FINAL BETA NOW IN HAND, SO COMPUTE EVAPORATION.  EVAP EQUALS ETP
C UNLESS BETA<1.
C ----------------------------------------------------------------------
      ETA = BETA*ETP

C ----------------------------------------------------------------------
C SET THE EFFECTIVE POTNL EVAPOTRANSP (ETP1) TO ZERO SINCE THIS IS SNOW
C CASE, SO SURFACE EVAP NOT CALCULATED FROM EDIR, EC, OR ETT IN SMFLX
C (BELOW).
C IF SEAICE (ICE=1) SKIP CALL TO SMFLX.
C SMFLX RETURNS UPDATED SOIL MOISTURE VALUES.  IN THIS, THE SNOW PACK
C CASE, ETA1 IS NOT USED IN CALCULATION OF EVAP.
C ----------------------------------------------------------------------
      ETP1 = 0.0
      IF (ICE .NE. 1) THEN
        CALL SMFLX (ETA1,SMC,NSOIL,CMC,ETP1,DT,PRCP1,ZSOIL,
     &              SH2O,SLOPE,KDT,FRZFACT,
     &              SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
     &              SMCREF,SHDFAC,CMCMAX,
     &              SMCDRY,CFACTR,RUNOFF1,RUNOFF2,RUNOFF3,
     &              EDIR,EC,ET,ETT,SFCTMP,Q2,NROOT,RTDIS,FXEXP,
     &              DRIP)

      ENDIF

C ----------------------------------------------------------------------
C BEFORE CALL SHFLX IN THIS SNOWPACK CASE, SET ZZ1 AND YY ARGUMENTS TO
C SPECIAL VALUES THAT ENSURE THAT GROUND HEAT FLUX CALCULATED IN SHFLX
C MATCHES THAT ALREADY COMPUTER FOR BELOW THE SNOWPACK, THUS THE SFC
C HEAT FLUX TO BE COMPUTED IN SHFLX WILL EFFECTIVELY BE THE FLUX AT THE
C SNOW TOP SURFACE.  T11 IS A DUMMY ARGUEMENT SO WE WILL NOT USE THE
C SKIN TEMP VALUE AS REVISED BY SHFLX.
C ----------------------------------------------------------------------
      ZZ1 = 1.0
      YY = STC(1)-0.5*SSOIL*ZSOIL(1)*ZZ1/DF1
      T11 = T1

C ----------------------------------------------------------------------
C SHFLX WILL CALC/UPDATE THE SOIL TEMPS.  NOTE:  THE SUB-SFC HEAT FLUX 
C (SSOIL1) AND THE SKIN TEMP (T11) OUTPUT FROM THIS SHFLX CALL ARE NOT
C USED  IN ANY SUBSEQUENT CALCULATIONS. RATHER, THEY ARE DUMMY VARIABLES
C HERE IN THE SNOPAC CASE, SINCE THE SKIN TEMP AND SUB-SFC HEAT FLUX ARE
C UPDATED INSTEAD NEAR THE BEGINNING OF THE CALL TO SNOPAC.
C ----------------------------------------------------------------------
      CALL SHFLX (SSOIL1,STC,SMC,SMCMAX,NSOIL,T11,DT,YY,ZZ1,ZSOIL,
     &            TBOT,ZBOT,SMCWLT,PSISAT,SH2O,BEXP,F1,DF1,ICE,
     &            QUARTZ,CSOIL)
      
C ----------------------------------------------------------------------
C SNOW DEPTH AND DENSITY ADJUSTMENT BASED ON SNOW COMPACTION.  YY IS
C ASSUMED TO BE THE SOIL TEMPERTURE AT THE TOP OF THE SOIL COLUMN.
C ----------------------------------------------------------------------
      IF (ESD .GT. 0.) THEN
        CALL SNOWPACK (ESD,DT,SNOWH,SNDENS,T1,YY)
      ELSE
        ESD = 0.
        SNOWH = 0.
        SNDENS = 0.
        SNCOND = 1.
      ENDIF

C ----------------------------------------------------------------------
C END SUBROUTINE SNOPAC
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 15. PHYSICS SUBROUTINE ==>  SUBROUTINE SNOWPACK CCCCCCCCCCCCCCCCCC
      SUBROUTINE SNOWPACK (ESD,DTSEC,SNOWH,SNDENS,TSNOW,TSOIL)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SNOWPACK
C ----------------------------------------------------------------------
C CALCULATE COMPACTION OF SNOWPACK UNDER CONDITIONS OF INCREASING SNOW
C DENSITY, AS OBTAINED FROM AN APPROXIMATE SOLUTION OF E. ANDERSON'S
C DIFFERENTIAL EQUATION (3.29), NOAA TECHNICAL REPORT NWS 19, BY VICTOR
C KOREN, 03/25/95.
C ----------------------------------------------------------------------
C ESD     WATER EQUIVALENT OF SNOW (M)
C DTSEC   TIME STEP (SEC)
C SNOWH   SNOW DEPTH (M)
C SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
C TSNOW   SNOW SURFACE TEMPERATURE (K)
C TSOIL   SOIL SURFACE TEMPERATURE (K)
C
C SUBROUTINE WILL RETURN NEW VALUES OF SNOWH AND SNDENS
C ----------------------------------------------------------------------
      INTEGER IPOL, J

      REAL BFAC,C1,C2,SNDENS,DSX,DTHR,DTSEC,DW,SNOWHC,SNOWH,PEXP,TAVGC,
     &     TSNOW,TSNOWC,TSOIL,TSOILC,ESD,ESDC,ESDCX,G,KN

      PARAMETER(C1 = 0.01, C2=21.0, G=9.81, KN=4000.0)

C ----------------------------------------------------------------------
C CONVERSION INTO SIMULATION UNITS
C ----------------------------------------------------------------------
      SNOWHC = SNOWH*100.
      ESDC = ESD*100.
      DTHR = DTSEC/3600.
      TSNOWC = TSNOW-273.15
      TSOILC = TSOIL-273.15

C ----------------------------------------------------------------------
C CALCULATING OF AVERAGE TEMPERATURE OF SNOW PACK
C ----------------------------------------------------------------------
      TAVGC = 0.5*(TSNOWC+TSOILC)                                    

C ----------------------------------------------------------------------
C CALCULATING OF SNOW DEPTH AND DENSITY AS A RESULT OF COMPACTION
C  SNDENS=DS0*(EXP(BFAC*ESD)-1.)/(BFAC*ESD)
C  BFAC=DTHR*C1*EXP(0.08*TAVGC-C2*DS0)
C NOTE: BFAC*ESD IN SNDENS EQN ABOVE HAS TO BE CAREFULLY TREATED
C NUMERICALLY BELOW:
C   C1 IS THE FRACTIONAL INCREASE IN DENSITY (1/(CM*HR)) 
C   C2 IS A CONSTANT (CM3/G) KOJIMA ESTIMATED AS 21 CMS/G
C ----------------------------------------------------------------------
      IF (ESDC .GT. 1.E-2) THEN
        ESDCX = ESDC
      ELSE
        ESDCX = 1.E-2
      ENDIF
      BFAC = DTHR*C1*EXP(0.08*TAVGC-C2*SNDENS)

C      DSX = SNDENS*((DEXP(BFAC*ESDC)-1.)/(BFAC*ESDC))
C ----------------------------------------------------------------------
C THE FUNCTION OF THE FORM (e**x-1)/x IMBEDDED IN ABOVE EXPRESSION
C FOR DSX WAS CAUSING NUMERICAL DIFFICULTIES WHEN THE DENOMINATOR "x"
C (I.E. BFAC*ESDC) BECAME ZERO OR APPROACHED ZERO (DESPITE THE FACT THAT
C THE ANALYTICAL FUNCTION (e**x-1)/x HAS A WELL DEFINED LIMIT AS 
C "x" APPROACHES ZERO), HENCE BELOW WE REPLACE THE (e**x-1)/x 
C EXPRESSION WITH AN EQUIVALENT, NUMERICALLY WELL-BEHAVED 
C POLYNOMIAL EXPANSION.
C
C NUMBER OF TERMS OF POLYNOMIAL EXPANSION, AND HENCE ITS ACCURACY, 
C IS GOVERNED BY ITERATION LIMIT "IPOL".
C      IPOL GREATER THAN 9 ONLY MAKES A DIFFERENCE ON DOUBLE
C            PRECISION (RELATIVE ERRORS GIVEN IN PERCENT %).
C       IPOL=9, FOR REL.ERROR <~ 1.6 E-6 % (8 SIGNIFICANT DIGITS)
C       IPOL=8, FOR REL.ERROR <~ 1.8 E-5 % (7 SIGNIFICANT DIGITS)
C       IPOL=7, FOR REL.ERROR <~ 1.8 E-4 % ...
C ----------------------------------------------------------------------
      IPOL = 4
      PEXP = 0.
      DO J = IPOL,1,-1
C        PEXP = (1. + PEXP)*BFAC*ESDC/REAL(J+1) 
        PEXP = (1. + PEXP)*BFAC*ESDCX/REAL(J+1) 
      END DO
      PEXP = PEXP + 1.

      DSX = SNDENS*(PEXP)
C ----------------------------------------------------------------------
C ABOVE LINE ENDS POLYNOMIAL SUBSTITUTION
C ----------------------------------------------------------------------
C     END OF KOREAN FORMULATION

C     BASE FORMULATION (COGLEY ET AL., 1990)
C     CONVERT DENSITY FROM G/CM3 TO KG/M3
C       DSM=SNDENS*1000.0
 
C       DSX=DSM+DTSEC*0.5*DSM*G*ESD/
C    &      (1E7*EXP(-0.02*DSM+KN/(TAVGC+273.16)-14.643))
 
C     CONVERT DENSITY FROM KG/M3 TO G/CM3
C       DSX=DSX/1000.0

C     END OF COGLEY ET AL. FORMULATION 

C ----------------------------------------------------------------------
C SET UPPER/LOWER LIMIT ON SNOW DENSITY
C ----------------------------------------------------------------------
      IF (DSX .GT. 0.40) DSX = 0.40
      IF (DSX .LT. 0.05) DSX = 0.05
      SNDENS = DSX
C ----------------------------------------------------------------------
C UPDATE OF SNOW DEPTH AND DENSITY DEPENDING ON LIQUID WATER DURING
C SNOWMELT.  ASSUMED THAT 13% OF LIQUID WATER CAN BE STORED IN SNOW PER
C DAY DURING SNOWMELT TILL SNOW DENSITY 0.40.
C ----------------------------------------------------------------------
      IF (TSNOWC .GE. 0.) THEN
        DW = 0.13*DTHR/24.
        SNDENS = SNDENS*(1.-DW)+DW
        IF (SNDENS .GT. 0.40) SNDENS = 0.40
      ENDIF

C ----------------------------------------------------------------------
C CALCULATE SNOW DEPTH (CM) FROM SNOW WATER EQUIVALENT AND SNOW DENSITY.
C CHANGE SNOW DEPTH UNITS TO METERS
C ----------------------------------------------------------------------
      SNOWHC = ESDC/SNDENS
      SNOWH = SNOWHC*0.01

C ----------------------------------------------------------------------
C END SUBROUTINE SNOWPACK
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 16. PHYSICS SUBROUTINE ==>  SUBROUTINE SNOWZ0 CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SNOWZ0 (SNCOVR,Z0)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SNOWZ0
C ----------------------------------------------------------------------
C CALCULATE TOTAL ROUGHNESS LENGTH OVER SNOW
C SNCOVR  FRACTIONAL SNOW COVER
C Z0      ROUGHNESS LENGTH (m)
C Z0S     SNOW ROUGHNESS LENGTH:=0.001 (m)
C ----------------------------------------------------------------------
      REAL SNCOVR, Z0, Z0S
c      PARAMETER (Z0S=0.001)
      
C CURRENT NOAH LSM CONDITION - MBEK, 09-OCT-2001
      Z0S = Z0
C
      Z0 = (1-SNCOVR)*Z0 + SNCOVR*Z0S
C ----------------------------------------------------------------------
C END SUBROUTINE SNOWZ0
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 17. PHYSICS SUBROUTINE ==>  SUBROUTINE SNOW_NEW CCCCCCCCCCCCCCCCCC

      SUBROUTINE SNOW_NEW (TEMP,NEWSN,SNOWH,SNDENS)

      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE SNOW_NEW
C ----------------------------------------------------------------------
C CALCULATE SNOW DEPTH AND DENSITITY TO ACCOUNT FOR THE NEW SNOWFALL.
C NEW VALUES OF SNOW DEPTH & DENSITY RETURNED.
C
C TEMP    AIR TEMPERATURE (K)
C NEWSN   NEW SNOWFALL (M)
C SNOWH   SNOW DEPTH (M)
C SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
C ----------------------------------------------------------------------
      REAL SNDENS
      REAL DSNEW
      REAL SNOWHC
      REAL HNEWC
      REAL SNOWH
      REAL NEWSN
      REAL NEWSNC
      REAL TEMP 
      REAL TEMPC
      
C ----------------------------------------------------------------------
C CONVERSION INTO SIMULATION UNITS      
C ----------------------------------------------------------------------
      SNOWHC = SNOWH*100.
      NEWSNC = NEWSN*100.
      TEMPC = TEMP-273.15
      
C ----------------------------------------------------------------------
C CALCULATING NEW SNOWFALL DENSITY DEPENDING ON TEMPERATURE
C EQUATION FROM GOTTLIB L. 'A GENERAL RUNOFF MODEL FOR SNOWCOVERED
C AND GLACIERIZED BASIN', 6TH NORDIC HYDROLOGICAL CONFERENCE,
C VEMADOLEN, SWEDEN, 1980, 172-177PP.
C-----------------------------------------------------------------------
      IF (TEMPC .LE. -15.) THEN
        DSNEW = 0.05
      ELSE                                                      
        DSNEW = 0.05+0.0017*(TEMPC+15.)**1.5
      ENDIF
      
C ----------------------------------------------------------------------
C ADJUSTMENT OF SNOW DENSITY DEPENDING ON NEW SNOWFALL      
C ----------------------------------------------------------------------
      HNEWC = NEWSNC/DSNEW
      SNDENS = (SNOWHC*SNDENS+HNEWC*DSNEW)/(SNOWHC+HNEWC)
      SNOWHC = SNOWHC+HNEWC
      SNOWH = SNOWHC*0.01
      
C ----------------------------------------------------------------------
C END SUBROUTINE SNOW_NEW
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 18. PHYSICS SUBROUTINE ==>  SUBROUTINE SRT  CCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SRT (RHSTT,EDIR,ET,SH2O,SH2OA,NSOIL,PCPDRP,
     &                ZSOIL,DWSAT,DKSAT,SMCMAX,BEXP,RUNOFF1, 
     &                RUNOFF2,DT,SMCWLT,SLOPE,KDT,FRZX,SICE,AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SRT
C ----------------------------------------------------------------------
C CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
C WATER DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX
C COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER CVFRZ      
      INTEGER IALP1
      INTEGER IOHINF
      INTEGER J
      INTEGER JJ      
      INTEGER K
      INTEGER KS
      INTEGER NSOIL

      REAL ACRT
      REAL AI(NSOLD)
      REAL BEXP
      REAL BI(NSOLD)
      REAL CI(NSOLD)
      REAL DD
      REAL DDT
      REAL DDZ
      REAL DDZ2
      REAL DENOM
      REAL DENOM2
      REAL DICE
      REAL DKSAT
      REAL DMAX(NSOLD)
      REAL DSMDZ
      REAL DSMDZ2
      REAL DT
      REAL DT1
      REAL DWSAT
      REAL EDIR
      REAL ET(NSOIL)
      REAL FCR
      REAL FRZX
      REAL INFMAX
      REAL KDT
      REAL MXSMC
      REAL MXSMC2
      REAL NUMER
      REAL PCPDRP
      REAL PDDUM
      REAL PX
      REAL RHSTT(NSOIL)
      REAL RUNOFF1
      REAL RUNOFF2
      REAL SH2O(NSOIL)
      REAL SH2OA(NSOIL)
      REAL SICE(NSOIL)
      REAL SICEMAX
      REAL SLOPE
      REAL SLOPX
      REAL SMCAV
      REAL SMCMAX
      REAL SMCWLT
      REAL SSTT
      REAL SUM
      REAL VAL
      REAL WCND
      REAL WCND2
      REAL WDF
      REAL WDF2
      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C REFERENCE FROZEN GROUND PARAMETER, CVFRZ, IS A SHAPE PARAMETER OF
C AREAL DISTRIBUTION FUNCTION OF SOIL ICE CONTENT WHICH EQUALS 1/CV.
C CV IS A COEFFICIENT OF SPATIAL VARIATION OF SOIL ICE CONTENT.  BASED
C ON FIELD DATA CV DEPENDS ON AREAL MEAN OF FROZEN DEPTH, AND IT CLOSE
C TO CONSTANT = 0.6 IF AREAL MEAN FROZEN DEPTH IS ABOVE 20 CM.  THAT IS
C WHY PARAMETER CVFRZ = 3 (INT{1/0.6*0.6}).
C CURRENT LOGIC DOESN'T ALLOW CVFRZ BE BIGGER THAN 3
C ----------------------------------------------------------------------
        PARAMETER(CVFRZ = 3)
        
C ----------------------------------------------------------------------
C DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF.  INCLUDE THE
C INFILTRATION FORMULE FROM SCHAAKE AND KOREN MODEL.
C MODIFIED BY Q DUAN
C ----------------------------------------------------------------------
      IOHINF=1

C ----------------------------------------------------------------------
C LET SICEMAX BE THE GREATEST, IF ANY, FROZEN WATER CONTENT WITHIN SOIL
C LAYERS.
C ----------------------------------------------------------------------
      SICEMAX = 0.0
      DO KS=1,NSOIL
       IF (SICE(KS) .GT. SICEMAX) SICEMAX = SICE(KS)
      END DO

C ----------------------------------------------------------------------
C DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF
C ----------------------------------------------------------------------
      PDDUM = PCPDRP
      RUNOFF1 = 0.0
      IF (PCPDRP .NE. 0.0) THEN

C ----------------------------------------------------------------------
C MODIFIED BY Q. DUAN, 5/16/94
C ----------------------------------------------------------------------
C        IF (IOHINF .EQ. 1) THEN

        DT1 = DT/86400.
        SMCAV = SMCMAX - SMCWLT
        DMAX(1)=-ZSOIL(1)*SMCAV

C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C ----------------------------------------------------------------------
        DICE = -ZSOIL(1) * SICE(1)
          
        DMAX(1)=DMAX(1)*(1.0 - (SH2OA(1)+SICE(1)-SMCWLT)/SMCAV)
        DD=DMAX(1)

        DO KS=2,NSOIL
          
C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C ----------------------------------------------------------------------
          DICE = DICE + ( ZSOIL(KS-1) - ZSOIL(KS) ) * SICE(KS)
         
          DMAX(KS) = (ZSOIL(KS-1)-ZSOIL(KS))*SMCAV
          DMAX(KS) = DMAX(KS)*(1.0 - (SH2OA(KS)+SICE(KS)-SMCWLT)/SMCAV)
          DD = DD+DMAX(KS)
        END DO

C ----------------------------------------------------------------------
C VAL = (1.-EXP(-KDT*SQRT(DT1)))
C IN BELOW, REMOVE THE SQRT IN ABOVE
C ----------------------------------------------------------------------
        VAL = (1.-EXP(-KDT*DT1))
        DDT = DD*VAL
        PX = PCPDRP*DT  
        IF (PX .LT. 0.0) PX = 0.0
        INFMAX = (PX*(DDT/(PX+DDT)))/DT
          
C ----------------------------------------------------------------------
C FROZEN GROUND VERSION:
C REDUCTION OF INFILTRATION BASED ON FROZEN GROUND PARAMETERS
C ----------------------------------------------------------------------
        FCR = 1. 
        IF (DICE .GT. 1.E-2) THEN 
          ACRT = CVFRZ * FRZX / DICE 
          SUM = 1.
          IALP1 = CVFRZ - 1 
          DO J = 1,IALP1
            K = 1
            DO JJ = J+1,IALP1
              K = K * JJ
            END DO
            SUM = SUM + (ACRT ** ( CVFRZ-J)) / FLOAT (K) 
          END DO
          FCR = 1. - EXP(-ACRT) * SUM 
        ENDIF 
        INFMAX = INFMAX * FCR

C ----------------------------------------------------------------------
C CORRECTION OF INFILTRATION LIMITATION:
C IF INFMAX .LE. HYDROLIC CONDUCTIVITY ASSIGN INFMAX THE VALUE OF
C HYDROLIC CONDUCTIVITY
C ----------------------------------------------------------------------
C         MXSMC = MAX ( SH2OA(1), SH2OA(2) ) 
        MXSMC = SH2OA(1)

        CALL WDFCND (WDF,WCND,MXSMC,SMCMAX,BEXP,DKSAT,DWSAT,
     &               SICEMAX)

        INFMAX = MAX(INFMAX,WCND)
        INFMAX = MIN(INFMAX,PX)

        IF (PCPDRP .GT. INFMAX) THEN
          RUNOFF1 = PCPDRP - INFMAX
          PDDUM = INFMAX
        ENDIF

      ENDIF

C ----------------------------------------------------------------------
C TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN LINE
C BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
C 'MXSMC = MAX(SH2OA(1), SH2OA(2))'
C ----------------------------------------------------------------------
      MXSMC = SH2OA(1)

      CALL WDFCND (WDF,WCND,MXSMC,SMCMAX,BEXP,DKSAT,DWSAT,
     &             SICEMAX)
 
C ----------------------------------------------------------------------
C CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
C ----------------------------------------------------------------------
      DDZ = 1. / ( -.5 * ZSOIL(2) )
      AI(1) = 0.0
      BI(1) = WDF * DDZ / ( -ZSOIL(1) )
      CI(1) = -BI(1)

C ----------------------------------------------------------------------
C CALC RHSTT FOR THE TOP LAYER AFTER CALC'NG THE VERTICAL SOIL MOISTURE
C GRADIENT BTWN THE TOP AND NEXT TO TOP LAYERS.
C ----------------------------------------------------------------------
      DSMDZ = ( SH2O(1) - SH2O(2) ) / ( -.5 * ZSOIL(2) )
      RHSTT(1) = (WDF * DSMDZ + WCND - PDDUM + EDIR + ET(1))/ZSOIL(1)
      SSTT = WDF * DSMDZ + WCND + EDIR + ET(1)

C ----------------------------------------------------------------------
C INITIALIZE DDZ2
C ----------------------------------------------------------------------
      DDZ2 = 0.0

C ----------------------------------------------------------------------
C LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABV PROCESS
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
        DENOM2 = (ZSOIL(K-1) - ZSOIL(K))
        IF (K .NE. NSOIL) THEN
          SLOPX = 1.

C ----------------------------------------------------------------------
C AGAIN, TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN
C LINE BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
C 'MXSMC2 = MAX (SH2OA(K), SH2OA(K+1))'
C ----------------------------------------------------------------------
          MXSMC2 = SH2OA(K)

          CALL WDFCND (WDF2,WCND2,MXSMC2,SMCMAX,BEXP,DKSAT,DWSAT,
     &                 SICEMAX)

C ----------------------------------------------------------------------
C CALC SOME PARTIAL PRODUCTS FOR LATER USE IN CALC'NG RHSTT
C ----------------------------------------------------------------------
          DENOM = (ZSOIL(K-1) - ZSOIL(K+1))
          DSMDZ2 = (SH2O(K) - SH2O(K+1)) / (DENOM * 0.5)

C ----------------------------------------------------------------------
C CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
C ----------------------------------------------------------------------
          DDZ2 = 2.0 / DENOM
          CI(K) = -WDF2 * DDZ2 / DENOM2
        ELSE

C ----------------------------------------------------------------------
C SLOPE OF BOTTOM LAYER IS INTRODUCED
C ----------------------------------------------------------------------
          SLOPX = SLOPE

C ----------------------------------------------------------------------
C RETRIEVE THE SOIL WATER DIFFUSIVITY AND HYDRAULIC CONDUCTIVITY FOR
C THIS LAYER
C ----------------------------------------------------------------------
          CALL WDFCND (WDF2,WCND2,SH2OA(NSOIL),SMCMAX,BEXP,DKSAT,DWSAT,
     &                 SICEMAX)

C ----------------------------------------------------------------------
C CALC A PARTIAL PRODUCT FOR LATER USE IN CALC'NG RHSTT
C ----------------------------------------------------------------------
          DSMDZ2 = 0.0

C ----------------------------------------------------------------------
C SET MATRIX COEF CI TO ZERO
C ----------------------------------------------------------------------
          CI(K) = 0.0
        ENDIF

C ----------------------------------------------------------------------
C CALC RHSTT FOR THIS LAYER AFTER CALC'NG ITS NUMERATOR
C ----------------------------------------------------------------------
        NUMER = (WDF2 * DSMDZ2) + SLOPX * WCND2 - (WDF * DSMDZ) 
     &    - WCND + ET(K)
        RHSTT(K) = NUMER / (-DENOM2)

C ----------------------------------------------------------------------
C CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER
C ----------------------------------------------------------------------
        AI(K) = -WDF * DDZ / DENOM2
        BI(K) = -( AI(K) + CI(K) )

C ----------------------------------------------------------------------
C RESET VALUES OF WDF, WCND, DSMDZ, AND DDZ FOR LOOP TO NEXT LYR
C RUNOFF2:  SUB-SURFACE OR BASEFLOW RUNOFF
C ----------------------------------------------------------------------
        IF (K .EQ. NSOIL) THEN
          RUNOFF2 = SLOPX * WCND2
        ENDIF

        IF (K .NE. NSOIL) THEN
          WDF = WDF2
          WCND = WCND2
          DSMDZ = DSMDZ2
          DDZ = DDZ2
        ENDIF
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE SRT
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 19. PHYSICS SUBROUTINE ==>  SUBROUTINE SSTEP  CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SSTEP (SH2OOUT,SH2OIN,CMC,RHSTT,RHSCT,DT,
     &                  NSOIL,SMCMAX,CMCMAX,RUNOFF3,ZSOIL,SMC,SICE,
     &                  AI,BI,CI)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SSTEP
C ----------------------------------------------------------------------
C CALCULATE/UPDATE SOIL MOISTURE CONTENT VALUES AND CANOPY MOISTURE
C CONTENT VALUES.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER I
      INTEGER K 
      INTEGER KK11
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)
      REAL CIin(NSOLD)
      REAL CMC
      REAL CMCMAX
      REAL DDZ
      REAL DT
      REAL RHSCT
      REAL RHSTT(NSOIL)
      REAL RHSTTin(NSOIL)
      REAL RUNOFF3
      REAL SH2OIN(NSOIL)
      REAL SH2OOUT(NSOIL)
      REAL SICE(NSOIL)
      REAL SMC(NSOIL)
      REAL SMCMAX
      REAL STOT
      REAL WPLUS
      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C CREATE 'AMOUNT' VALUES OF VARIABLES TO BE INPUT TO THE
C TRI-DIAGONAL MATRIX ROUTINE.
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        RHSTT(K) = RHSTT(K) * DT
        AI(K) = AI(K) * DT
        BI(K) = 1. + BI(K) * DT
        CI(K) = CI(K) * DT
      END DO

C ----------------------------------------------------------------------
C COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        RHSTTin(K) = RHSTT(K)
      END DO
      DO K = 1,NSOLD
        CIin(K) = CI(K)
      END DO

C ----------------------------------------------------------------------
C CALL ROSR12 TO SOLVE THE TRI-DIAGONAL MATRIX
C ----------------------------------------------------------------------
      CALL ROSR12 (CI,AI,BI,CIin,RHSTTin,RHSTT,NSOIL)

C ----------------------------------------------------------------------
C SUM THE PREVIOUS SMC VALUE AND THE MATRIX SOLUTION TO GET A
C NEW VALUE.  MIN ALLOWABLE VALUE OF SMC WILL BE 0.02.
C RUNOFF3: RUNOFF WITHIN SOIL LAYERS
C ----------------------------------------------------------------------
      WPLUS = 0.0
      RUNOFF3 = 0.
      DDZ = -ZSOIL(1)
      
      DO K = 1,NSOIL
        IF (K .NE. 1) DDZ = ZSOIL(K - 1) - ZSOIL(K)
        SH2OOUT(K) = SH2OIN(K) + CI(K) + WPLUS / DDZ

        STOT = SH2OOUT(K) + SICE(K)
        IF (STOT .GT. SMCMAX) THEN
          IF (K .EQ. 1) THEN
            DDZ = -ZSOIL(1)
          ELSE
            KK11 = K - 1
            DDZ = -ZSOIL(K) + ZSOIL(KK11)
          ENDIF
          WPLUS = (STOT-SMCMAX) * DDZ
        ELSE
          WPLUS = 0.
        ENDIF
        SMC(K) = MAX ( MIN(STOT,SMCMAX),0.02 )
        SH2OOUT(K) = MAX((SMC(K)-SICE(K)),0.0)
      END DO

      RUNOFF3 = WPLUS

C ----------------------------------------------------------------------
C UPDATE CANOPY WATER CONTENT/INTERCEPTION (CMC).  CONVERT RHSCT TO 
C AN 'AMOUNT' VALUE AND ADD TO PREVIOUS CMC VALUE TO GET NEW CMC.
C ----------------------------------------------------------------------
      CMC = CMC + DT * RHSCT
      IF (CMC .LT. 1.E-20) CMC=0.0
      CMC = MIN(CMC,CMCMAX)

C ----------------------------------------------------------------------
C END SUBROUTINE SSTEP
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 20. PHYSICS SUBROUTINE ==>  SUBROUTINE TBND CCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TBND (TU,TB,ZSOIL,ZBOT,K,NSOIL,TBND1)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE TBND
C ----------------------------------------------------------------------
C CALCULATE TEMPERATURE ON THE BOUNDARY OF THE LAYER BY INTERPOLATION OF
C THE MIDDLE LAYER TEMPERATURES
C ----------------------------------------------------------------------
      INTEGER NSOIL
      INTEGER K

      REAL TBND1
      REAL T0
      REAL TU
      REAL TB
      REAL ZB
      REAL ZBOT
      REAL ZUP
      REAL ZSOIL (NSOIL)

      PARAMETER(T0 = 273.15)

C ----------------------------------------------------------------------
C USE SURFACE TEMPERATURE ON THE TOP OF THE FIRST LAYER
C ----------------------------------------------------------------------
      IF (K .EQ. 1) THEN
        ZUP = 0.
      ELSE
        ZUP = ZSOIL(K-1)
      ENDIF

C ----------------------------------------------------------------------
C USE DEPTH OF THE CONSTANT BOTTOM TEMPERATURE WHEN INTERPOLATE
C TEMPERATURE INTO THE LAST LAYER BOUNDARY
C ----------------------------------------------------------------------
      IF (K .EQ. NSOIL) THEN
        ZB = 2.*ZBOT-ZSOIL(K)
      ELSE
        ZB = ZSOIL(K+1)
      ENDIF

C ----------------------------------------------------------------------
C LINEAR INTERPOLATION BETWEEN THE AVERAGE LAYER TEMPERATURES
C ----------------------------------------------------------------------
      TBND1 = TU+(TB-TU)*(ZUP-ZSOIL(K))/(ZUP-ZB)
      
C ----------------------------------------------------------------------
C END SUBROUTINE TBND
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 21. PHYSICS SUBROUTINE ==>  SUBROUTINE TDFCND CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TDFCND ( DF, SMC, QZ,  SMCMAX, SH2O)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE TDFCND
C ----------------------------------------------------------------------
C CALCULATE THERMAL DIFFUSIVITY AND CONDUCTIVITY OF THE SOIL FOR A GIVEN
C POINT AND TIME.
C ----------------------------------------------------------------------
C PETERS-LIDARD APPROACH (PETERS-LIDARD et al., 1998)
C June 2001 CHANGES: FROZEN SOIL CONDITION.
C ----------------------------------------------------------------------
       REAL DF
       REAL GAMMD
       REAL THKDRY
       REAL AKE
       REAL THKICE
       REAL THKO
       REAL THKQTZ
       REAL THKSAT
       REAL THKS
       REAL THKW
       REAL QZ
       REAL SATRATIO
       REAL SH2O
       REAL SMC
       REAL SMCMAX
       REAL XU
       REAL XUNFROZ

C ----------------------------------------------------------------------
C WE NOW GET QUARTZ AS AN INPUT ARGUMENT (SET IN ROUTINE REDPRM):
C      DATA QUARTZ /0.82, 0.10, 0.25, 0.60, 0.52, 
C     &             0.35, 0.60, 0.40, 0.82/
C ----------------------------------------------------------------------
C IF THE SOIL HAS ANY MOISTURE CONTENT COMPUTE A PARTIAL SUM/PRODUCT
C OTHERWISE USE A CONSTANT VALUE WHICH WORKS WELL WITH MOST SOILS
C ----------------------------------------------------------------------
C  THKW ......WATER THERMAL CONDUCTIVITY
C  THKQTZ ....THERMAL CONDUCTIVITY FOR QUARTZ
C  THKO ......THERMAL CONDUCTIVITY FOR OTHER SOIL COMPONENTS
C  THKS ......THERMAL CONDUCTIVITY FOR THE SOLIDS COMBINED(QUARTZ+OTHER)
C  THKICE ....ICE THERMAL CONDUCTIVITY
C  SMCMAX ....POROSITY (= SMCMAX)
C  QZ .........QUARTZ CONTENT (SOIL TYPE DEPENDENT)
C ----------------------------------------------------------------------
C USE AS IN PETERS-LIDARD, 1998 (MODIF. FROM JOHANSEN, 1975).
C
C                                  PABLO GRUNMANN, 08/17/98
C REFS.:
C      FAROUKI, O.T.,1986: THERMAL PROPERTIES OF SOILS. SERIES ON ROCK 
C              AND SOIL MECHANICS, VOL. 11, TRANS TECH, 136 PP.
C      JOHANSEN, O., 1975: THERMAL CONDUCTIVITY OF SOILS. PH.D. THESIS,
C              UNIVERSITY OF TRONDHEIM,
C      PETERS-LIDARD, C. D., ET AL., 1998: THE EFFECT OF SOIL THERMAL 
C              CONDUCTIVITY PARAMETERIZATION ON SURFACE ENERGY FLUXES
C              AND TEMPERATURES. JOURNAL OF THE ATMOSPHERIC SCIENCES,
C              VOL. 55, PP. 1209-1224.
C ----------------------------------------------------------------------
C NEEDS PARAMETERS
C POROSITY(SOIL TYPE):
C      POROS = SMCMAX
C SATURATION RATIO:
      SATRATIO = SMC/SMCMAX

C PARAMETERS  W/(M.K)
      THKICE = 2.2
      THKW = 0.57
      THKO = 2.0
C      IF (QZ .LE. 0.2) THKO = 3.0
      THKQTZ = 7.7
C SOLIDS' CONDUCTIVITY      
      THKS = (THKQTZ**QZ)*(THKO**(1.- QZ))

C UNFROZEN FRACTION (FROM 1., i.e., 100%LIQUID, TO 0. (100% FROZEN))
      XUNFROZ = (SH2O + 1.E-9) / (SMC + 1.E-9)

C UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
      XU=XUNFROZ*SMCMAX 
C SATURATED THERMAL CONDUCTIVITY
      THKSAT = THKS**(1.-SMCMAX)*THKICE**(SMCMAX-XU)*THKW**(XU)

C DRY DENSITY IN KG/M3
      GAMMD = (1. - SMCMAX)*2700.

C DRY THERMAL CONDUCTIVITY IN W.M-1.K-1
      THKDRY = (0.135*GAMMD + 64.7)/(2700. - 0.947*GAMMD)

      IF ( (SH2O + 0.0005) .LT. SMC ) THEN
C FROZEN
              AKE = SATRATIO
      ELSE
C UNFROZEN
C RANGE OF VALIDITY FOR THE KERSTEN NUMBER (AKE)
          IF ( SATRATIO .GT. 0.1 ) THEN

C KERSTEN NUMBER (USING "FINE" FORMULA, VALID FOR SOILS CONTAINING AT 
C LEAST 5% OF PARTICLES WITH DIAMETER LESS THAN 2.E-6 METERS.)
C (FOR "COARSE" FORMULA, SEE PETERS-LIDARD ET AL., 1998).

              AKE = LOG10(SATRATIO) + 1.0

          ELSE

C USE K = KDRY
              AKE = 0.0

          ENDIF
      ENDIF

C  THERMAL CONDUCTIVITY

       DF = AKE*(THKSAT - THKDRY) + THKDRY

C ----------------------------------------------------------------------
C END SUBROUTINE TDFCND
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 22. PHYSICS SUBROUTINE ==>  SUBROUTINE TMPAVG CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TMPAVG (TAVG,TUP,TM,TDN,ZSOIL,NSOIL,K) 
      
      IMPLICIT NONE
      
C ----------------------------------------------------------------------
C SUBROUTINE TMPAVG
C ----------------------------------------------------------------------
C CALCULATE SOIL LAYER AVERAGE TEMPERATURE (TAVG) IN FREEZING/THAWING
C LAYER USING UP, DOWN, AND MIDDLE LAYER TEMPERATURES (TUP, TDN, TM),
C WHERE TUP IS AT TOP BOUNDARY OF LAYER, TDN IS AT BOTTOM BOUNDARY OF
C LAYER.  TM IS LAYER PROGNOSTIC STATE TEMPERATURE.
C ----------------------------------------------------------------------
      INTEGER K
      INTEGER NSOIL

      REAL DZ
      REAL DZH
      REAL T0
      REAL TAVG
      REAL TDN
      REAL TM
      REAL TUP
      REAL X0
      REAL XDN
      REAL XUP
      REAL ZSOIL (NSOIL)

      PARAMETER(T0 = 2.7315E2)

C ----------------------------------------------------------------------
      IF (K .EQ. 1) THEN
        DZ = -ZSOIL(1)
      ELSE
        DZ = ZSOIL(K-1)-ZSOIL(K)
      ENDIF

      DZH=DZ*0.5

      IF (TUP .LT. T0) THEN
        IF (TM .LT. T0) THEN
          IF (TDN .LT. T0) THEN
C ----------------------------------------------------------------------
C TUP, TM, TDN < T0
C ----------------------------------------------------------------------
            TAVG = (TUP + 2.0*TM + TDN)/ 4.0            
          ELSE
C ----------------------------------------------------------------------
C TUP & TM < T0,  TDN >= T0
C ----------------------------------------------------------------------
            X0 = (T0 - TM) * DZH / (TDN - TM)
            TAVG = 0.5 * (TUP*DZH+TM*(DZH+X0)+T0*(2.*DZH-X0)) / DZ
          ENDIF      
        ELSE
          IF (TDN .LT. T0) THEN
C ----------------------------------------------------------------------
C TUP < T0, TM >= T0, TDN < T0
C ----------------------------------------------------------------------
            XUP  = (T0-TUP) * DZH / (TM-TUP)
            XDN  = DZH - (T0-TM) * DZH / (TDN-TM)
            TAVG = 0.5 * (TUP*XUP+T0*(2.*DZ-XUP-XDN)+TDN*XDN) / DZ
          ELSE
C ----------------------------------------------------------------------
C TUP < T0, TM >= T0, TDN >= T0
C ----------------------------------------------------------------------
            XUP  = (T0-TUP) * DZH / (TM-TUP)
            TAVG = 0.5 * (TUP*XUP+T0*(2.*DZ-XUP)) / DZ
          ENDIF   
        ENDIF
      ELSE
        IF (TM .LT. T0) THEN
          IF (TDN .LT. T0) THEN
C ----------------------------------------------------------------------
C TUP >= T0, TM < T0, TDN < T0
C ----------------------------------------------------------------------
            XUP  = DZH - (T0-TUP) * DZH / (TM-TUP)
            TAVG = 0.5 * (T0*(DZ-XUP)+TM*(DZH+XUP)+TDN*DZH) / DZ
          ELSE
C ----------------------------------------------------------------------
C TUP >= T0, TM < T0, TDN >= T0
C ----------------------------------------------------------------------
            XUP  = DZH - (T0-TUP) * DZH / (TM-TUP)
            XDN  = (T0-TM) * DZH / (TDN-TM)
            TAVG = 0.5 * (T0*(2.*DZ-XUP-XDN)+TM*(XUP+XDN)) / DZ
          ENDIF   
        ELSE
          IF (TDN .LT. T0) THEN
C ----------------------------------------------------------------------
C TUP >= T0, TM >= T0, TDN < T0
C ----------------------------------------------------------------------
            XDN  = DZH - (T0-TM) * DZH / (TDN-TM)
            TAVG = (T0*(DZ-XDN)+0.5*(T0+TDN)*XDN) / DZ                 
          ELSE
C ----------------------------------------------------------------------
C TUP >= T0, TM >= T0, TDN >= T0
C ----------------------------------------------------------------------
            TAVG = (TUP + 2.0*TM + TDN) / 4.0
          ENDIF
        ENDIF
      ENDIF
C ----------------------------------------------------------------------
C END SUBROUTINE TMPAVG
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 23. PHYSICS SUBROUTINE ==>  SUBROUTINE TRANSP CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TRANSP (ET,NSOIL,ETP1,SMC,CMC,ZSOIL,SHDFAC,SMCWLT,
     &                   CMCMAX,PC,CFACTR,SMCREF,SFCTMP,Q2,NROOT,RTDIS)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE TRANSP
C ----------------------------------------------------------------------
C CALCULATE TRANSPIRATION FOR THE VEG CLASS.
C ----------------------------------------------------------------------
      INTEGER I
      INTEGER K
      INTEGER NSOIL
      INTEGER NROOT

      REAL CFACTR
      REAL CMC
      REAL CMCMAX
      REAL DENOM
      REAL ET(NSOIL)
      REAL ETP1
      REAL ETP1A
      REAL GX (7)
C.....REAL PART(NSOIL)
      REAL PC
      REAL Q2
      REAL RTDIS(NSOIL)
      REAL RTX
      REAL SFCTMP
      REAL SGX
      REAL SHDFAC
      REAL SMC(NSOIL)
      REAL SMCREF
      REAL SMCWLT
      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C INITIALIZE PLANT TRANSP TO ZERO FOR ALL SOIL LAYERS.
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        ET(K) = 0.
      END DO

C ----------------------------------------------------------------------
C CALCULATE AN 'ADJUSTED' POTENTIAL TRANSPIRATION
C IF STATEMENT BELOW TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
C NOTE: GX AND OTHER TERMS BELOW REDISTRIBUTE TRANSPIRATION BY LAYER,
C ET(K), AS A FUNCTION OF SOIL MOISTURE AVAILABILITY, WHILE PRESERVING
C TOTAL ETP1A.
C ----------------------------------------------------------------------
      IF (CMC .NE. 0.0) THEN
        ETP1A = SHDFAC * PC * ETP1 * (1.0 - (CMC /CMCMAX) ** CFACTR)
      ELSE
        ETP1A = SHDFAC * PC * ETP1
      ENDIF
      
      SGX = 0.0
      DO I = 1,NROOT
        GX(I) = ( SMC(I) - SMCWLT ) / ( SMCREF - SMCWLT )
        GX(I) = MAX ( MIN ( GX(I), 1. ), 0. )
        SGX = SGX + GX (I)
      END DO
      SGX = SGX / NROOT
      
      DENOM = 0.
      DO I = 1,NROOT
        RTX = RTDIS(I) + GX(I) - SGX
        GX(I) = GX(I) * MAX ( RTX, 0. )
        DENOM = DENOM + GX(I)
      END DO
      IF (DENOM .LE. 0.0) DENOM = 1.

      DO I = 1,NROOT
        ET(I) = ETP1A * GX(I) / DENOM
      END DO

C ----------------------------------------------------------------------
C ABOVE CODE ASSUMES A VERTICALLY UNIFORM ROOT DISTRIBUTION
C CODE BELOW TESTS A VARIABLE ROOT DISTRIBUTION
C ----------------------------------------------------------------------
C      ET(1) = ( ZSOIL(1) / ZSOIL(NROOT) ) * GX * ETP1A
C      ET(1) = ( ZSOIL(1) / ZSOIL(NROOT) ) * ETP1A
C ----------------------------------------------------------------------
C USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
C ----------------------------------------------------------------------
C      ET(1) = RTDIS(1) * ETP1A
C      ET(1) = ETP1A * PART(1)
C ----------------------------------------------------------------------
C LOOP DOWN THRU THE SOIL LAYERS REPEATING THE OPERATION ABOVE,
C BUT USING THE THICKNESS OF THE SOIL LAYER (RATHER THAN THE
C ABSOLUTE DEPTH OF EACH LAYER) IN THE FINAL CALCULATION.
C ----------------------------------------------------------------------
C      DO K = 2,NROOT
C        GX = ( SMC(K) - SMCWLT ) / ( SMCREF - SMCWLT )
C        GX = MAX ( MIN ( GX, 1. ), 0. )
C TEST CANOPY RESISTANCE
C        GX = 1.0
C        ET(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT))*GX*ETP1A
C        ET(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT))*ETP1A
C ----------------------------------------------------------------------
C USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
C ----------------------------------------------------------------------
C        ET(K) = RTDIS(K) * ETP1A
C        ET(K) = ETP1A*PART(K)
C      END DO      
C ----------------------------------------------------------------------
C END SUBROUTINE TRANSP
C ----------------------------------------------------------------------
      RETURN
      END

CCCCC 24. PHYSICS SUBROUTINE ==>  SUBROUTINE WDFCND CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE WDFCND (WDF,WCND,SMC,SMCMAX,BEXP,DKSAT,DWSAT,
     &                   SICEMAX)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE WDFCND
C ----------------------------------------------------------------------
C CALCULATE SOIL WATER DIFFUSIVITY AND SOIL HYDRAULIC CONDUCTIVITY.
C ----------------------------------------------------------------------
      REAL BEXP
      REAL DKSAT
      REAL DWSAT
      REAL EXPON
      REAL FACTR1
      REAL FACTR2
      REAL SICEMAX
      REAL SMC
      REAL SMCMAX
      REAL VKwgt
      REAL WCND
      REAL WDF

C ----------------------------------------------------------------------
C     CALC THE RATIO OF THE ACTUAL TO THE MAX PSBL SOIL H2O CONTENT
C ----------------------------------------------------------------------
      SMC = SMC
      SMCMAX = SMCMAX
      FACTR1 = 0.2 / SMCMAX
      FACTR2 = SMC / SMCMAX

C ----------------------------------------------------------------------
C PREP AN EXPNTL COEF AND CALC THE SOIL WATER DIFFUSIVITY
C ----------------------------------------------------------------------
      EXPON = BEXP + 2.0
      WDF = DWSAT * FACTR2 ** EXPON

C ----------------------------------------------------------------------
C FROZEN SOIL HYDRAULIC DIFFUSIVITY.  VERY SENSITIVE TO THE VERTICAL
C GRADIENT OF UNFROZEN WATER. THE LATTER GRADIENT CAN BECOME VERY
C EXTREME IN FREEZING/THAWING SITUATIONS, AND GIVEN THE RELATIVELY 
C FEW AND THICK SOIL LAYERS, THIS GRADIENT SUFFERES SERIOUS 
C TRUNCTION ERRORS YIELDING ERRONEOUSLY HIGH VERTICAL TRANSPORTS OF
C UNFROZEN WATER IN BOTH DIRECTIONS FROM HUGE HYDRAULIC DIFFUSIVITY.  
C THEREFORE, WE FOUND WE HAD TO ARBITRARILY CONSTRAIN WDF 
C --
C VERSION D_10CM: ........  FACTR1 = 0.2/SMCMAX
C WEIGHTED APPROACH...................... PABLO GRUNMANN, 28_SEP_1999.
C ----------------------------------------------------------------------
      IF (SICEMAX .GT. 0.0)  THEN
        VKWGT = 1./(1.+(500.*SICEMAX)**3.)
        WDF = VKWGT*WDF + (1.- VKWGT)*DWSAT*FACTR1**EXPON
      ENDIF

C ----------------------------------------------------------------------
C RESET THE EXPNTL COEF AND CALC THE HYDRAULIC CONDUCTIVITY
C ----------------------------------------------------------------------
      EXPON = (2.0 * BEXP) + 3.0
      WCND = DKSAT * FACTR2 ** EXPON

C ----------------------------------------------------------------------
C END SUBROUTINE WDFCND
C ----------------------------------------------------------------------
      RETURN
      END
