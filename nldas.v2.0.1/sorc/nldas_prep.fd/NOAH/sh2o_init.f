!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! sh2oinit.f: 
!
! DESCRIPTION:
!  To do a 'cold start' initialization of soil liquid water SH2O for 
!   NOAH LSM, though also adaptable for other land-surface models,
!   using either GDAS or Eta forcing data. 
!
! *NOTE:
!   LDAS%STARTCODE = NOAH_IC: 4=model init 
!
! REVISION HISTORY:
!  NCEP;  Original code developed by NCEP for Eta model
!           subroutines to initialize soil liquid water, SH2O   
!  04 Nov 2002: Kristi Arsenault; Modified code to be used with noahrst.f
!                                 to initialize NOAH with NCEP forcing data  
!
!=========================================================================

      SUBROUTINE SH2OINIT(SMC,STC,SMCMAX,PSIS,BETA,SH2O)

      IMPLICIT NONE      

!=== Local Variables =====================================================
      INTEGER :: C,R,T,I,J,L,N ! Loop counters
      INTEGER :: VCLASS,NC,NR,NCH
      INTEGER VEGT  ! Tile veg type

      REAL T1       ! NOAH Skin Temperature (K)
      REAL STC      ! NOAH Soil Layer Temperature (K)
      REAL SMC      ! NOAH Soil Layer Total Moisture (liq+frzn) 
      REAL SH2O     ! NOAH Soil Layer Liquid Moisture

      REAL  SOILMOIST,SOILTEMP
      REAL  POR1,POR2,POR3,POR4
      REAL  SMCWLT                  ! NOAH Wilting point for soil moisture (vol.)
      REAL  SLDPTH1                 ! Soil layer thicknesses (m)
      REAL  SLDPTH2, SLDPTH3, SLDPTH4
      REAL  PSIS                    ! Saturated soil potential
      REAL  BETA                    ! B-parameter 
      REAL  SMCMAX                  ! Max soil moisture content (porosity)
      REAL  BX
      REAL  FK
      REAL  FRH2O 

      REAL :: RHOICE=917.0          ! Density of ice
      REAL :: WT1,WT2               ! Weights for soil wetness initialization
      REAL :: HLICE=3.335E5         ! Ice parameter
      REAL :: GRAV=9.81             ! Gravity (m s-1) 
      REAL :: T0=273.15             ! Freezing point of water
      REAL :: BLIM=5.5              ! B-parameter upper limit

!=== End Variable Definition =============================================

C ----------------------------------------------------------------------
C COLD START:  determine liquid soil water content (SH2O)
C NSOIL number of soil layers
C ----------------------------------------------------------------------
C  SH2O <= SMC for T < 273.149K (-0.001C)
        IF (STC .LT. 273.149) THEN
C ----------------------------------------------------------------------
C first guess following explicit solution for Flerchinger Eqn from Koren 
C et al, JGR, 1999, Eqn 17 (KCOUNT=0 in FUNCTION FRH2O). 
           BX = BETA
           IF ( BETA .GT. BLIM )  BX = BLIM

!       print *, 'SMC,STC,SMCMAX,PSIS: ',SMC,STC,SMCMAX,PSIS
           FK=(((HLICE/(GRAV*(-PSIS)))*
     .        ((STC-T0)/STC))**(-1/BX))*SMCMAX
           IF (FK .LT. 0.02) FK = 0.02
           SH2O = MIN ( FK, SMC )
C ----------------------------------------------------------------------
C now use iterative solution for liquid soil water content using 
C FUNCTION FRH2O with the initial guess for SH2O from above explicit
C first guess.
           SH2O = FRH2O(STC,SMC,SH2O,SMCMAX,BETA,PSIS)

        ELSE
C ----------------------------------------------------------------------
C  SH2O = SMC for T => 273.149K (-0.001C)
           SH2O=SMC
C ----------------------------------------------------------------------
        ENDIF

      RETURN 
      END 

!===========================================================================

       FUNCTION FRH2O(TKELV,SMC,SH2O,SMCMAX,BEXP,PSIS)

       IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  PURPOSE:  CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL WATER CONTENT 
C  IF TEMPERATURE IS BELOW 273.15K (T0).  REQUIRES NEWTON-TYPE ITERATION 
C  TO SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF 
C  KOREN ET AL. (1999, JGR, VOL 104(D16), 19569-19585).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c New version (JUNE 2001): much faster and more accurate newton iteration 
c achieved by first taking log of eqn cited above -- less than 4
c (typically 1 or 2) iterations achieves convergence.  Also, explicit  
c 1-step solution option for special case of parameter Ck=0, which reduces 
c the original implicit equation to a simpler explicit form, known as the 
c ""Flerchinger Eqn". Improved handling of solution in the limit of  
c freezing point temperature T0. 
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT:
C   TKELV.........Temperature (Kelvin)
C   SMC...........Total soil moisture content (volumetric)
C   SH2O..........Liquid soil moisture content (volumetric)
C   SMCMAX........Saturation soil moisture content (from REDPRM)
C   B.............Soil type "B" parameter (from REDPRM)
C   PSIS..........Saturated soil matric potential (from REDPRM)
C
C OUTPUT:
C   FRH2O.........supercooled liquid water content.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       REAL BEXP
       REAL BLIM
       REAL BX
       REAL CK
       REAL DENOM
       REAL DF
       REAL DH2O
       REAL DICE
       REAL DSWL
       REAL ERROR
       REAL FK
       REAL FRH2O
       REAL GS
       REAL HLICE
       REAL PSIS
       REAL SH2O
       REAL SMC
       REAL SMCMAX
       REAL SWL
       REAL SWLK
       REAL TKELV
       REAL T0

       INTEGER NLOG
       INTEGER KCOUNT

       PARAMETER (CK=8.0)
C      PARAMETER (CK=0.0)
       PARAMETER (BLIM=5.5)
C      PARAMETER (BLIM=7.0)
       PARAMETER (ERROR=0.005)

       PARAMETER (HLICE=3.335E5)
       PARAMETER (GS = 9.81)
       PARAMETER (DICE=920.0)
       PARAMETER (DH2O=1000.0)
       PARAMETER (T0=273.15)

C  ###   LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)  ####
C  ###   SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT  ####
C  ###   IS NON-REALISTICALLY HIGH AT VERY LOW TEMPERATURES    ####
C##################################################################

      BX = BEXP
      IF ( BEXP .GT. BLIM ) BX = BLIM
C------------------------------------------------------------------

C INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
      NLOG=0
      KCOUNT=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (T0), SH2O = SMC 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF (TKELV .GT. (T0 - 1.E-3)) THEN
         FRH2O=SMC

      ELSE

        IF (CK .NE. 0.0) THEN

C -------------------------------------------------------------
C OPTION 1: ITERATED SOLUTION FOR NONZERO CK
C IN KOREN ET AL, JGR, 1999, EQN 17
C -------------------------------------------------------------
C INITIAL GUESS FOR SWL (frozen content)
         SWL = SMC-SH2O

C KEEP WITHIN BOUNDS.
          IF (SWL .GT. (SMC-0.02)) SWL=SMC-0.02
          IF(SWL .LT. 0.) SWL=0. 
C--------------------------------------------------------------
C  START OF ITERATIONS 
C--------------------------------------------------------------
         DO WHILE (NLOG .LT. 10 .AND. KCOUNT .EQ. 0)
          NLOG = NLOG+1
          DF = ALOG(( PSIS*GS/HLICE ) * ( ( 1.+CK*SWL )**2. ) *
     &      ( SMCMAX/(SMC-SWL) )**BX) - ALOG(-(TKELV-T0)/TKELV)
          DENOM = 2. * CK / ( 1.+CK*SWL ) + BX / ( SMC - SWL )
          SWLK = SWL - DF/DENOM
C BOUNDS USEFUL FOR MATHEMATICAL SOLUTION.
          IF (SWLK .GT. (SMC-0.02)) SWLK = SMC - 0.02
          IF(SWLK .LT. 0.) SWLK = 0.
C MATHEMATICAL SOLUTION BOUNDS APPLIED.
          DSWL=ABS(SWLK-SWL)
          SWL=SWLK 

C---------------------------------------------------------------
C IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)  
C WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED. 
C---------------------------------------------------------------
          IF ( DSWL .LE. ERROR )  THEN
            KCOUNT=KCOUNT+1
          END IF
         END DO 
C---------------------------------------------------------------
C  END OF ITERATIONS 
C---------------------------------------------------------------
C BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION.
         FRH2O = SMC - SWL

CCCCCCCCCCCCCCCCCCCCCCCC END OPTION 1 CCCCCCCCCCCCCCCCCCCCCCCCCCC

        ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0
C IN KOREN ET AL., JGR, 1999, EQN 17
C----------------------------------------------------------------
        IF (KCOUNT .EQ. 0) THEN
c      Print*,'Flerchinger used in NEW version. Iterations=',NLOG
          FK=(((HLICE/(GS*(-PSIS)))*((TKELV-T0)/TKELV))**
     .      (-1/BX))*SMCMAX
C  APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION
          IF (FK .LT. 0.02) FK = 0.02
          FRH2O = MIN ( FK, SMC )

CCCCCCCCCCCCCCCCCCCCCCCCC END OPTION 2 CCCCCCCCCCCCCCCCCCCCCCCCCC

        ENDIF

      ENDIF

      RETURN
      END

