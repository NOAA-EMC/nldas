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
! mos_tile.f:
!
! DESCRIPTION:
!
! REVISION HISTORY:
!  18  May 2000: Brian Cosgrove; Fixed a divide by zero error in the RSURFP
!                subroutine.  If denominator is 0, set to 1E-37 instead.
!  22 Aug. 2000: Brian Cosgrove; Modified code for output of
!                standard LDAS output variables.  LDAS and MOS modules 
!                are now passed into MOSTILE
!  28 Sep. 2000: Brian Cosgrove; Modified condition statement for output
!                of canopy conductance...only compute if canopy resistance
!                is greater than 100 (previous cutoff was 0)
!  17 Nov. 2000: Brian Cosgrove; Fixed MNDY variable--data was being 
!                written in a shape that didn't correspond to the array
!                shape using old style FORTRAN methods
!  06 Feb. 2001: Brian Cosgrove; removed ALHX from function m_esat declaration. 
!                This function declaration now matches the call in mos_main call.
!  24 May 2001: Brian Cosgrove; Fixed divide by zero bug in SMFAC subroutine
!  02 Jul 2001: Urszula Jambor: Applied Brian's SWET1 fix in M_INTERC call.
!  13 Jul 2001: Matt Rodell; Modified RSOIL calculation for SOILCO==>zero
!  05 Sep 2001: Brian Cosgrove; Added DUMMY2 variable to calls,
!               CSOIL0 was changed by someone, sometime from 70000 to 175000
!  21 Feb 2001: Brian Cosgrove; Altered ASTRO subroutine so that it
!               no longer uses the irun variable to determine how many
!               lats/lons to process.  Now it only process 1 lat/lon. 
!               Before this, it was taking in a scaler variable (tile%lat)
!               but was expecting an array of dimension 'irun' which even
!               if irun set to 1 could be bad programming
!  15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new 
!                GRIB libraries
!  17 Nov 2003: Brian Cosgrove; Changed CSOIL back to value of 70000 which
!               is Kosters original value.  The PSAS value of 175000
!               leads to bad ground heat flux values when the PSAS
!               scheme is not used...which is most of the time.
!=========================================================================
 

         SUBROUTINE MOSTILE (RESET,
     I            DTSTEP, TRAINL,TRAINC, TSNOW,  UM,
     I            ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,
     I            TM,     QM,     CD,     SUNANG, PARDIR,  PARDIF,
     I            SWNET,  HLWDWN, PSUR,    ZLAI,   GREEN,  Z2,
     I            SQSCAT, RSOIL1, RSOIL2,   RDC,    U2FAC,
     I            QSATTC, DQSDTC, ALWRAD, BLWRAD,
     U            DTCANAL, TC, TD, QA, SWET1, SWET2, SWET3, CAPAC, SNOW,
     O            EVAP, SHFLUX, RUNOFF, BOMB,
     O            EINT,   ESOI,   EVEG,   ESNO, SMELT, HLATN,
     O            HLWUP,GDRAIN,RUNSRF,FWSOIL,
     O            STRDG1, STRDG2, STRDG3, STRDG4,
     O            WSMAX1, WSMAX2, WSMAX3, CONDRY, GHFLUX,
     &            POR1,POR2,POR3,water1,water2,water3,
     &         	  VGBEEX,VGPSAX,VGCSAX,VGZDEX,
     &    	  VGSLOX,VGPH1X,VGPH2X,VGRPLX,
     &     	  VGWMAX,SNWMID,
     &     	  VGCHIL,VGZMEW,VGRST1,VGRST2,VGRST3, 
     &     	  VGRDRS,VGZ2,VGROTD,
     &     	  VGDFAC,VGTLL,VGTU,
     &            VGTCF1,VGTCF2,VGTCF3,III,LDAS,TILE)

! Declare modules and data structures

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      IMPLICIT NONE
      type (ldasdec) LDAS
      type (tiledec) TILE

C   SCCS VERSION @(#)lsm.f	1.2 11/3/92 

CCCC *
CCCC *       This subroutine computes the outgoing fluxes of sensible and
CCCC * latent heat from the land surface and updates the surface prognostic
CCCC * variables.
CCCC *


C
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER   FRSTCH, MemFac,ISINDX
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE, outrc, outra
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1, MemFac = 5)
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

C       SOIL PARAMETERS PASSED INTO TILE
        REAL VGBEEX             !Soil paramter B related to pore size
				!  distribution index
        REAL VGPSAX             !Soil moisture potential of a saturated
                                !  soil (m)
        REAL VGCSAX             !Hydraulic conductivity of the soil
                                !  at saturation (m/s)
        REAL VGZDEX(NLAY)       !Soil layer thickness (m) of layer i
	REAL DELZ12		!Distance between the centers of the
				!  first and second soil layers
	REAL DELZ23		!Distance between the centers of the
				!  second and third soil layers
	REAL VGSLOX		!Cosine of theta
	REAL VGPH1X		!Soil moisture potential above which
				!  vegetation is not moisture-stressed (m)
	REAL VGPH2X		!Soil moisture potential below which
                                !  transpiration ceases due to wilting (m)
	REAL VGRPLX		!Average resistance to moisture transport
                                !  within the vegetation itself (s)
	REAL VGWMAX(NLAY)	!Moisture holding capacity of soil
				!  layer i (kg/m2)
	REAL SNWMID		!
	REAL VGCHIL		!Parameter describing departure of leaf
				!  angles from a spherical distribution
	REAL VGZMEW		!
	REAL VGRST1		!Stomatal resistance parameter a
        REAL VGRST2             !Stomatal resistance parameter b
        REAL VGRST3             !Stomatal resistance parameter c
        REAL VGRDRS             !Resistance to moisture transport per unit
                                !  root length used to derive RSOIL1
        REAL VGZ2               !Height of the canopy, Z2 is set to this,
                                !  and it's used to derive U2FAC
        REAL VGROTD             !Rooting depth parameter used to derive
                                !  RSOIL1 and RSOIL2
        REAL VGDFAC             !Parameter controlling vapor pressure
                                !  deficit stress (1/mb)
	REAL VGTLL		!Temperature below which temperature
				!  stress prevents transpiration (K)
	REAL VGTU		!Temperature above which temperature
                                !  stress prevents transpiration (K)
	REAL VGTCF1		!Coefficient in temperature stress
				!  equation (K-4)
        REAL VGTCF2             !Coefficient in temperature stress
                                !  equation (K-3)
        REAL VGTCF3             !Coefficient in temperature stress
                                !  equation (K-2)
	REAL GREEN		!Current greenness fraction of vegetation
	REAL ZLAI		!Current leaf area index of vegetation
        REAL VGZ0           !Month Parameter roughness length
                                !  Table 7 in Mosaic manual
        REAL VGROTL        !Month parameter variation of root length density
                                !  Table 9 in Mosaic manual
        REAL VGRDC         !Month parameter variation of vegetation specific
                                !  constant used to determin subcanopy
                                !  aerodynamic resistance.
                                !  Table 8 in Mosaic manual
        REAL VGDD          !Month parameter variation of zero plane
                                !  displacement height
                                !  Table 10 in Mosaic manual
        REAL VGROCA             !Average cross sectional area of root,


c  This is for testing
        integer ival,jval


C       VEGETATION PARAMETERS PASSED INTO TILE

CCCC *
C      INTEGER ITYP
      REAL  DTSTEP,TRAINL,TRAINC,TSNOW,UM,
     &       ETURB,  DEDQA, HSTURB, DHSDTC,
     &         TM ,     CD, SUNANG, DHSDQA,
     &         QM , PARDIR, PARDIF,  SWNET,
     &      HLWDWN,   PSUR, 
     &          Z2, SQSCAT,  DEDTC
      REAL  RSOIL1, RSOIL2,    RDC,  U2FAC,
     &     QSATTC, DQSDTC,  ALWRAD, BLWRAD,
     &     DTCANAL,    TC,      TD,     QA,   BOMB,
     &       SWET1,  SWET2,  SWET3,  CAPAC, 
     &        SNOW,   EVAP, SHFLUX, RUNOFF
      REAL   EINT,    ESOI,   EVEG,   ESNO,
     &     STRDG1,  STRDG2, STRDG3, STRDG4,
     &      SMELT,   HLATN,  HLWUP, GDRAIN,
     &     RUNSRF,  FWSOIL, WSMAX1, WSMAX2,
     &     WSMAX3,  CONDRY, GHFLUX
      REAL water1,water2,water3
      INTEGER RESET
      REAL SWET(NLAY),DELTC,DELEA,CSOIL0,WSOI12     

      REAL PHLAY(NLAY), AKLAY(NLAY), SWET12,
     &      CSOIL,
     &       RCUN,     VPDSTR,   ESATTX,
     &       VPDSTX
      REAL EMAXRT,  FTEMP,
     &        PHR,      SOILCO,    RC,
     &        EAX,          TX,   RCX,
     &     DRCDTC,      SATCAP,   PAR,
     &       PDIR,       DUMMY, DUMMY2
      REAL FTEMPX,  DRCDEA
C The following are now passed in and only of unit or layer # length
      REAL  DEDEA,  DHSDEA,     EM, ESATTC,
     &     DESDTC,      EA,     RA,   ALHX,
     &     WETEQ1, WETEQ2,
     &        RX1,     RX2, SNWFRC, POTFRC,
     &     ESNFRC,  EIRFRC,   FCAN,  EPFRC,
     &     DEFCIT, EADJST, RTBS
      real cmpbug,POR1,POR2,POR3
	integer III    !the actual tile number this can be removed
CC The following are just passed into here to be passed
CC on to M_RCUNST


CCCC * - - - - - - - -
CCCC * THE FOLLOWING 3 VARIABLES MUST MAINTAIN CONSISTENT VALUES.  ZDEP12(N) =
CCCC *  (VGZDEX(1,N)+VGZDEX(2,N))/2. ; SIMILARLY FOR ZDEP23(N).
CCCC *
C
C  Due to the variability of the layer depths
c  the value for DELZ12 and DELZ23 will
c  be computed in the program below

C    DATA DELZ12 / 0.75, 0.75, 0.75, 0.245, 0.245, 0.095, 0.0092,
C     8        0.0092, 0.095, 0.0092 /
C      DATA DELZ23 / 1.74, 1.74, 1.74, 0.735, 0.735, 0.585, 0.155,
C     8        0.155, 0.235, 0.155 /
CCCC * - - - - - - - -
CCCC *
      DATA DELTC /0.01/, DELEA /0.001/
C      DATA CSOIL0 /175000./
      DATA CSOIL0 /70000./

CCCC *
CCCC * --------------------------------------------------------------------- 
CCCC *
CFPP$ EXPAND (M_TMPFAC, M_RCANOP, SOIL, M_WUPDAT, SMFAC)
CCCC *
CCCC * Expand data as specified by ITYP into arrays of size NCH.
CCCC *
        DELZ12=( (VGZDEX(1)+VGZDEX(2)) /2.0)
	  DELZ23=( (VGZDEX(2)+VGZDEX(3)) /2.0)


       BOMB = 0.

        POR1=VGWMAX(SFCLY)/
     1             (VGZDEX(SFCLY)*1000.)
        POR2=VGWMAX(ROOTLY)/
     1             (VGZDEX(ROOTLY)*1000.)
        POR3=VGWMAX(RECHLY)/
     1             (VGZDEX(RECHLY)*1000.)
       IF(RESET.eq.1)THEN  !convert ETA volumetric to 
                           !%saturation soil moisture
        SWET1=SWET1/POR1       
        IF(SWET1.gt.1.0)SWET1=0.95
        IF(SWET1.lt.0.0)SWET1=0.05
        SWET2=SWET2/POR2      
        IF(SWET2.gt.1.0)SWET2=0.95
        IF(SWET2.lt.0.0)SWET2=0.05
        SWET3=SWET3/POR3       
        IF(SWET3.gt.1.0)SWET3=0.95
        IF(SWET3.lt.0.0)SWET3=0.05
       ENDIF


CCCC *
CCCC * Pre-process input arrays as necessary:

      DEDQA  = AMAX1(  DEDQA, 500./ALHE )
      DEDTC  = AMAX1(  DEDTC,   0. )
      DHSDQA = AMAX1( DHSDQA,   0. )
      DHSDTC = AMAX1( DHSDTC, -10. )
       
      EM     = QM * PSUR / EPSILON
      EA     = QA * PSUR / EPSILON
      ESATTC = QSATTC * PSUR / EPSILON

      DESDTC = DQSDTC * PSUR / EPSILON
C      DEDEA  = DEDQA * EPSILON / ( PSUR * ALHX )
  
      DEDEA  = DEDQA * EPSILON / PSUR
      DHSDEA = DHSDQA * EPSILON / PSUR
C      DEDTC  = DEDTC / ALHX
C      ETURB  = ETURB / ALHX
      PAR    = (PARDIR + PARDIF + 1.E-20)
      PDIR   = PARDIR / PAR
      RA     = ONE / ( CD * UM )
      SATCAP = 0.2 * ZLAI
      CSOIL  = CSOIL0
      SWET(SFCLY ) = amax1(amin1(SWET1,1.),0.)
      SWET(ROOTLY) = amax1(amin1(SWET2,1.),0.)
      SWET(RECHLY) = amax1(amin1(SWET3,1.),0.)
      CAPAC = amax1(amin1(CAPAC,SATCAP),0.)
CCCC *
     
      SNWFRC = SNOW / ( SNOW + SNWMID )
      FCAN = AMIN1( 1., AMAX1(0.,CAPAC/SATCAP) )
      POTFRC=1.-(1.-SNWFRC)*(1.-FCAN)


      WSOI12=SWET(SFCLY ) * VGWMAX(SFCLY) +
     &       SWET(ROOTLY) * VGWMAX(ROOTLY)
      SWET12 = WSOI12 /
     &          (VGWMAX(SFCLY) + VGWMAX(ROOTLY))

      EMAXRT = ( SNOW + CAPAC + WSOI12 ) / DTSTEP
CCCC *
      RUNOFF = 0.
      RUNSRF = 0.
      GHFLUX = 0.
CCCC * (SMELT is zeroed in FLUXES)
C      SMELT  = 0.
      FWSOIL = 0.
 100  CONTINUE

CCCC *
CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CCCC *       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
CCCC *          (M_RCUNST  computes the unstressed resistance,
CCCC *           M_VPDFAC  computes the vapor pressure deficit stress term,
CCCC *           M_TMPFAC  computes the temperature stress term,
CCCC *           SOIL    computes soil moisture potentials and conds.,
CCCC *           SMFAC   computes rc given a leaf water potential stress,
CCCC *           RSURFP  computes rc given a parallel resist. from the
CCCC *                   surface,
CCCC *           M_RCANOP  computes rc corrected for snow and interception.)
CCCC *

      CALL M_RCUNST (SUNANG, SQSCAT, PDIR,
     I               PAR, ZLAI, GREEN,
     O               RCUN, III,
     &               VGCHIL,VGZMEW,VGRST1,VGRST2,VGRST3)

      CALL M_VPDFAC ( ESATTC, EA,   
     O               VPDSTR, VGDFAC
     &              )
 
      CALL M_TMPFAC (TC, FTEMP,
     &      VGTLL, VGTU,
     &      VGTCF1, VGTCF2, VGTCF3, III) 

      CALL SOIL ( SWET12, VGBEEX, VGPSAX, VGCSAX, DELZ12,
     O           PHR, SOILCO, DUMMY
     &          )
      cmpbug=0.

      CALL SMFAC ( ESATTC, EA,   PHR, SOILCO, RCUN,
     I            VPDSTR, FTEMP,  TC,     PSUR, Z2,  RSOIL1, RSOIL2,
     I            VGPH1X, VGPH2X, VGRPLX,
     O            RC, III
     &           )

      CALL RSURFP ( UM, U2FAC, Z2, RDC, SWET(SFCLY),   
     I             ESATTC, EA,
     U             RC,
     O             RX1, RX2, III
     &            )

      CALL M_RCANOP ( CAPAC, SNOW, SATCAP, RA, ETURB, SNWFRC,
     I               POTFRC,
     U               RC
     &              )
CCCC *
CCCC * -    -    -    -    -    -    -    -    -    -    -    -    -    -
CCCC *
CCCC * Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
CCCC *

      TX = TC + DELTC
      ESATTX = ESATTC + DESDTC * DELTC
      EAX = EA + DELEA
  120 CONTINUE

CCCC *
CCCC * temperature:
      CALL M_VPDFAC (ESATTX, EA, VPDSTX, VGDFAC)
      CALL M_TMPFAC (TX, FTEMPX,
     &            VGTLL, VGTU,
     &            VGTCF1, VGTCF2, VGTCF3, III)
      CALL SMFAC (  ESATTX, EA,   PHR, SOILCO, RCUN,
     I            VPDSTX, FTEMPX, TX,     PSUR, Z2,  RSOIL1, RSOIL2,
     I            VGPH1X, VGPH2X, VGRPLX,
     O            RCX, III
     &           )
      CALL RSURFP (UM, U2FAC, Z2, RDC, SWET(SFCLY),
     I             ESATTX, EA,
     &             RCX,
     O             DUMMY,DUMMY2, III
     &            )
      CALL M_RCANOP (CAPAC,SNOW,SATCAP,RA,ETURB,SNWFRC,
     &               POTFRC, RCX)
CCCC *

      DRCDTC = (RCX - RC) / DELTC
  125 CONTINUE
CCCC *
       
CCCC * vapor pressure:

      CALL M_VPDFAC (ESATTC, EAX, VPDSTX, VGDFAC)
      CALL SMFAC ( ESATTC, EAX,  PHR, SOILCO, RCUN,
     I            VPDSTX, FTEMP, TC,     PSUR, Z2,  RSOIL1, RSOIL2,
     I            VGPH1X, VGPH2X, VGRPLX,
     O            RCX, III
     &           )
      CALL RSURFP (UM, U2FAC, Z2, RDC, SWET(SFCLY),
     I             ESATTC, EAX,
     &             RCX,
     O             DUMMY,DUMMY2, III
     &            )
      CALL M_RCANOP (CAPAC,SNOW,SATCAP,RA,ETURB,SNWFRC,
     &               POTFRC, RCX)
CCCC *

      DRCDEA = (RCX - RC) / DELEA
  150 CONTINUE
CCCC *

CCCC *
CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CCCC *       STEP 2: Solve the energy balance at the surface.
CCCC *
CCCC *
CCCC * Determine effective latent heat of vaporization based on
CCCC * assumed fractional coverage of snow.  This requires the
CCCC * determination of the fractions of moisture taken from the various
CCCC * reservoirs:
CCCC *         ESNFRC=fraction of total evap. coming from snow
CCCC *         EIRFRC=fraction of total evap. coming from interception res.
      
      RTBS=RX1*RX2/(RX1+RX2)
      EPFRC=POTFRC * ( RA + RTBS ) /
     &                    ( RA + POTFRC*RTBS )
      ESNFRC=EPFRC*SNWFRC/
     &                      (SNWFRC+FCAN+1.E-20)
      EIRFRC=EPFRC*FCAN/(SNWFRC+FCAN+1.E-20)
      ALHX = (1.-ESNFRC)*ALHE + ESNFRC*ALHS

      CONDRY = 1./(RTBS+1.E-20)

  200 CONTINUE
      CALL M_FLUXES (DTSTEP, ESATTC, DESDTC,   ALHX, 
     I               ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC, 
     I               RC, DRCDEA, DRCDTC, 
     I               SWNET, HLWDWN, ALWRAD, BLWRAD, DTCANAL,ESNFRC,
     I               TM,     EM,  CSOIL,   PSUR, EMAXRT, VGWMAX,              
     U               TC,     TD,     EA,   SWET,   SNOW,
     O               RUNOFF,   EVAP, SHFLUX,  SMELT,  HLWUP, BOMB,
     O               GHFLUX,III
     &              )
       
CCCC *
CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CCCC *       STEP 3: Update the water quantities in the soil and in the
CCCC *          interception reservoir.
CCCC *          (M_WUPDAT  reduces moisture contents using calculated
CCCC *           evap.,
CCCC *           GWATER  computes darcian flux of water between soil
CCCC *           layers.)
CCCC *

      CALL M_WUPDAT (DTSTEP,
     I               EVAP, SATCAP, VGWMAX, TC, RA, RC,
     I               RX1, RX2,ESNFRC,EIRFRC,
     U               CAPAC, SNOW, SWET, RUNOFF,
     O               EINT, ESOI, EVEG, ESNO, RUNSRF,FWSOIL,
     &            III  )

CCCC * 
CCCC * Correct any energy balance inconsistency.
CCCC *

      IF(EVAP .GT. 0.) THEN
        EADJST=(ESNO/DTSTEP) - ESNFRC*EVAP
        SHFLUX=SHFLUX-EADJST*(ALHS-ALHE)
        ENDIF
  300 CONTINUE

      CALL SOIL (SWET (SFCLY), VGBEEX,
     I             VGPSAX,       VGCSAX,      DELZ12,
     O             PHLAY(SFCLY), AKLAY(SFCLY),  DUMMY
     &            )

      CALL SOIL (  SWET(ROOTLY), VGBEEX,    
     I             VGPSAX,       VGCSAX,      DELZ12,
     O             PHLAY(ROOTLY), AKLAY(ROOTLY),   WETEQ1
     &            )

      CALL SOIL ( SWET(RECHLY), VGBEEX,
     I             VGPSAX,       VGCSAX,      DELZ23,
     O             PHLAY(RECHLY), AKLAY(RECHLY),   WETEQ2
     &            )

      CALL GWATER (VGWMAX, PHLAY, AKLAY, TC, DTSTEP,
     I             VGZDEX, VGSLOX, WETEQ1, WETEQ2,
     I             VGPSAX, VGCSAX,
     U             SWET, RUNOFF, GDRAIN,
     &            III)

CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CCCC *
CCCC *       STEP 4: Allow precipitation to fill interception reservoir
CCCC *          and top soil layer
CCCC *

      CALL SOIL ( SWET(ROOTLY), VGBEEX,
     I             VGPSAX,       VGCSAX,      DELZ12,
     O             PHLAY(ROOTLY), AKLAY(ROOTLY),   WETEQ1
     &            )
       
      SWET1=SWET(SFCLY)
 
      CALL M_INTERC ( DTSTEP, TRAINL, TRAINC, TSNOW, SATCAP,
     I               VGWMAX, CSOIL, WETEQ1,
     U               TC, CAPAC, SNOW, SWET1, RUNOFF, RUNSRF,
c     U               TC, CAPAC, SNOW, SWET, RUNOFF, RUNSRF,
     &               SMELT, FWSOIL, GHFLUX
     &              )
 
      
CCCC *
CCCC *
CCCC * Process data for return to GCM:


      QA = EA * EPSILON / PSUR
C      DEDTC = DEDTC * ALHX
C      ETURB = ETURB * ALHX
      HLATN  = EVAP  * ALHX
c      SWET1 = SWET(SFCLY)
      SWET2 = SWET(ROOTLY)
      SWET3 = SWET(RECHLY)


      IF(SWET1.LT.1.E-10) SWET1 = 0.0
      IF(SWET2.LT.1.E-10) SWET2 = 0.0
      IF(SWET3.LT.1.E-10) SWET3 = 0.0
      IF(CAPAC.LT.1.E-10) CAPAC = 0.0
      IF(SNOW .LT.1.E-10) SNOW  = 0.0
      IF(RUNOFF.LT.1.E-10) RUNOFF = 0.0
      EINT = EINT * ALHE / DTSTEP
      ESOI = ESOI * ALHE / DTSTEP
      EVEG = EVEG * ALHE / DTSTEP
      ESNO = ESNO * ALHS / DTSTEP
      DEFCIT = EVAP * ( RC + RA )
      STRDG1 = DEFCIT / RA
      STRDG2 = DEFCIT / ( RA + RCUN )
      STRDG3 = DEFCIT / ( RA + RCUN/VPDSTR )
      STRDG4 = DEFCIT / ( RA +
     *                   RCUN/(VPDSTR*FTEMP) )
      WSMAX1 = VGWMAX(1)
      WSMAX2 = VGWMAX(2)
      WSMAX3 = VGWMAX(3)

c       Create variables that hold the depth of water in each layer
c       do this by multiplying max depth of water in each layer
c       (WSMAX) by the degree of saturation in each layer
c       (SWET)
        water1=SWET1*WSMAX1
        water2=SWET2*WSMAX2
        water3=SWET3*WSMAX3

 2000 CONTINUE
CCCC *
	if (rc.gt.100.0) then
	TILE%CC=1.0/rc 
	else
c	if conapy resistance is 0, then make 1/rc undefined..occurs
c	in negative evaporation processes etc...
	TILE%CC=ldas%udef
	endif
	TILE%AC=CD*UM
      RETURN
      END
CCCC *
CCCC * [ END CHIP ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN M_INTERC ]
CCCC * 
      SUBROUTINE M_INTERC ( DTSTEP, TRAINL, TRAINC,
     I                     TSNOW, SATCAP, WMAX, CSOIL, WETEQ1,
     U                     TC, CAPAC, SNOW,  SWET1,  RUNOFF, RUNSRF,
     U                     SMELT, FWSOIL, GHFLUX
     &                    )
CCCC *
CCCC * This routine uses the precipitation forcing to determine 
CCCC * changes in interception, snowcover, and soil moisture storage.
CCCC *
CCCC * Note:  In this formulation, rain that falls on frozen ground
CCCC * runs-off, rather than freezes.
CCCC *
      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER   FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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
CCCC *
      REAL  TRAINL, TRAINC, TSNOW,  SATCAP,  
     &      WMAX(NLAY), TC,   CSOIL,   CAPAC,     
     &      SNOW, SWET1,  RUNOFF,  SMELT,
     &     RUNSRF, FWSOIL, GHFLUX, WETEQ1, WETINT
      REAL DTSTEP, SNOWM, WATADD, CAVAIL, THRUC, WRUNC, WRUNL,
     &     TIMFRL, TIMFRC, FWETL, FWETC, THRU1, THRU2, THRUL, XTCORR,
     &     WETFRC, SMPERS
CCCC *
C      DATA FWET0  /0.30/
      DATA FWETL /1.00/, FWETC /0.20/, TIMFRL/1.00/, TIMFRC/0.125/
CCCC *
CCCC * ------------------------------------------------------------------

CCCC *
CCCC * Add to snow cover.  Melt snow if necessary.
      SNOW = SNOW + TSNOW*DTSTEP
      SNOWM = 0.
      IF( SNOW.GT.0. .AND. TC.GT.TF ) THEN
        SNOWM = AMIN1( SNOW,
     &       AMAX1( 0., (TC-TF)*CSOIL/ALHM )  )
        IF( SNOWM .EQ. SNOW ) THEN
            TC = TC - SNOWM * ALHM / CSOIL
            SNOW=0.
          ELSE
            TC=TF
            SNOW = SNOW - SNOWM
          ENDIF
        SMPERS=SNOWM/DTSTEP
        SMELT=SMELT+SMPERS
        GHFLUX=GHFLUX-SMPERS*ALHM
        ENDIF
CCCC *
CCCC * =======================================================
CCCC *
CCCC * Load interception reservoir.  STEP 1: Large scale condensation.
CCCC *
CCCC * Determine XTCORR, the fraction of a storm that falls on a previously
CCCC * wet surface due to the time correlation of precipitation position.
CCCC * (Time scale TIMFRL for large scale storms set to one for FWETL=1
CCCC * to reflect the effective loss of "position memory" when storm 
CCCC * covers entire grid square.)

      XTCORR= (1.-TIMFRL) * AMIN1( 1.,(CAPAC/SATCAP)/FWETL )    
CCCC *
CCCC * Fill interception reservoir with precipitation.
CCCC * THRU1 is first calculated as the amount falling through the 
CCCC *    canopy under the assumption that all rain falls randomly.  
CCCC *    only a fraction 1-XTCORR falls randomly, though, so the result 
CCCC *    is multiplied by 1-XTCORR.
CCCC *
      WATADD = TRAINL*DTSTEP + SMELT*DTSTEP
      CAVAIL = ( SATCAP - CAPAC ) * FWETL
      WETINT = CAPAC/SATCAP
      IF( WATADD*(1.-WETINT) .LT. CAVAIL ) THEN
          THRU1 = WATADD*WETINT
        ELSE
          THRU1 = (WATADD - CAVAIL)
        ENDIF
      THRU1=THRU1*(1.-XTCORR)

CCCC * THRU2 is the amount that falls immediately through the canopy due
CCCC * to 'position memory'.

      THRU2=XTCORR*WATADD

      THRUL=THRU1+THRU2
      CAPAC=CAPAC+WATADD-THRU1-THRU2
CCCC *
CCCC * ---------------------------------------------------
CCCC *
CCCC * STEP 2: Moist convective precipitation.
CCCC *
CCCC * Determine XTCORR, the fraction of a storm that falls on a previously
CCCC * wet surface due to the time correlation of precipitation position.

      XTCORR= (1.-TIMFRC) * AMIN1( 1.,(CAPAC/SATCAP)/FWETC )    

CCCC *
CCCC * Fill interception reservoir with precipitation.
CCCC * THRU1 is first calculated as the amount falling through the 
CCCC *    canopy under the assumption that all rain falls randomly.  
CCCC *    only a fraction 1-XTCORR falls randomly, though, so the result 
CCCC *    is multiplied by 1-XTCORR.
CCCC *
      WATADD = TRAINC*DTSTEP
      CAVAIL = ( SATCAP - CAPAC ) * FWETC
      WETINT = CAPAC/SATCAP
      IF( WATADD*(1.-WETINT) .LT. CAVAIL ) THEN
          THRU1 = WATADD*WETINT
        ELSE
          THRU1 = (WATADD - CAVAIL)
        ENDIF
      THRU1=THRU1*(1.-XTCORR)

CCCC * THRU2 is the amount that falls immediately through the canopy due
CCCC * to 'position memory'.

      THRU2=XTCORR*WATADD

      THRUC=THRU1+THRU2
      CAPAC=CAPAC+WATADD-THRU1-THRU2
CCCC *
CCCC ** =================================================================
CCCC *
CCCC * Add precipitation moisture to soil, generate runoff.  if 
CCCC * temperature is below freezing, all precipitation runs off.
CCCC *
CCCC * STEP 1: Large scale condensation:
CCCC *
      IF(SWET1.GT.WETEQ1 .AND. WETEQ1.NE.1.) THEN
          WETFRC=(SWET1-WETEQ1)/(1.-WETEQ1)
        ELSE
          WETFRC=0.
        ENDIF

      CAVAIL = ( 1.-WETFRC)*FWETL*WMAX(1)*(1-WETEQ1)
      IF ( THRUL * (1-WETFRC) .LT. CAVAIL ) THEN
          WRUNL = THRUL * WETFRC
          SWET1 = SWET1  
     *    + THRUL * ( 1. - WETFRC) / WMAX(1)
        ELSE
          WRUNL = THRUL - CAVAIL
          SWET1 = SWET1 + CAVAIL / WMAX(1)
        ENDIF
CCCC *

CCCC * STEP 2: Moist convective precipitation:
CCCC *
      IF(SWET1.GT.WETEQ1 .AND. WETEQ1.NE.1.) THEN
          WETFRC=(SWET1-WETEQ1)/(1.-WETEQ1)
        ELSE
          WETFRC=0.
        ENDIF

      CAVAIL = (1.-WETFRC)*FWETC*WMAX(1)*(1-WETEQ1)
      IF ( THRUC * (1-WETFRC) .LT. CAVAIL ) THEN
          WRUNC = THRUC * WETFRC
          SWET1 = SWET1  
     *    + THRUC * ( 1. - WETFRC) / WMAX(1)
        ELSE
          WRUNC = THRUC - CAVAIL
          SWET1 = SWET1 + CAVAIL / WMAX(1)
        ENDIF
CCCC *
      RUNOFF = RUNOFF + (WRUNC+WRUNL)/DTSTEP
      RUNSRF = RUNSRF + (WRUNC+WRUNL)/DTSTEP
      FWSOIL = FWSOIL + (THRUC+THRUL-WRUNC-WRUNL)/DTSTEP
CCCC *
 100  CONTINUE
CCCC *
      RETURN
      END
CCCC *
CCCC * [ END M_INTERC ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN M_RCUNST ]
CCCC *
      SUBROUTINE M_RCUNST ( SUNANG, SQSCAT, PDIR,
     I                     PAR, ZLAI, GREEN,
     O                     RCUN, III,
     &                     VGCHIL,VGZMEW,VGRST1,VGRST2,VGRST3)
CCCC *
CCCC *     This subroutine calculates the unstressed canopy resistance.
CCCC * (p. 1353, Sellers 1985.)  Extinction coefficients are computed first.
CCCC *
      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER   FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY
      Integer III  !actual tile #
      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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
CCCC *

      REAL  SUNANG,   PDIR,   PAR,   ZLAI,
     &      SQSCAT,  GREEN,  RCUN

      REAL   VGCHIL, VGZMEW,         
     &       VGRST1, VGRST2, VGRST3

      REAL            RHO4,         EXTK1,         EXTK2,
     &               RCINV,         GAMMA,          EKAT,    DUM1,
     &                DUM2,          DUM3,            AA,      BB,
     &                  ZK,            CC,          PAR0





 
CCCC * First compute optical parameters.
CCCC * (Note: CHIL is constrained to be >= 0.01, as in SiB calcs.)

      AA = 0.5 - (0.633 + 0.330*VGCHIL)*VGCHIL
      BB = 0.877 * ( ONE - 2.*AA )
      CC =  ( AA + BB*SUNANG ) / SUNANG

      EXTK1 =  CC * SQSCAT
      EXTK2 = (ONE / VGZMEW) * SQSCAT

      DUM1 =      PDIR  *   CC
      DUM2 = (ONE-PDIR) * ( BB*(ONE/3.+PIE/4.) + AA*1.5 )

CCCC * Bound extinction coefficient by 50./ZLAI:

      ZK =     PDIR *AMIN1( EXTK1, 50./ZLAI ) +
     &    (ONE-PDIR)*AMIN1( EXTK2, 50./ZLAI )

CCCC * Now compute unstressed canopy resistance:

      GAMMA = VGRST1 / VGRST3 + VGRST2

      EKAT = EXP( ZK*ZLAI )
      PAR0=AMAX1( PAR, .000001 )
      RHO4 = GAMMA / (PAR0 * (DUM1 + DUM2))

      DUM1 = (VGRST2 - GAMMA) / (GAMMA + 1.E-20)
      DUM2 = (RHO4 * EKAT + ONE) / (RHO4 + ONE)
	
      DUM3 = ZK * VGRST3
      RCINV = ( DUM1*ALOG(DUM2) + ZK*ZLAI ) / DUM3

      RCUN = ONE / (RCINV * GREEN + 1.E-10)

 100  CONTINUE


      RETURN
      END
CCCC *
CCCC * [ END M_RCUNST ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN SOIL ]
CCCC *
      SUBROUTINE SOIL ( WET, VGBEEX, VGPSAX, VGCSAX, DELZ,
     O         PHR, SOILCO, WETEQ
     &                )
CCCC *
CCCC * This subroutine returns soil moisture potential and conductivity.

      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER  FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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
CCCC *

      REAL    WET,    PHR, SOILCO, VGPSAX,
     &    VGCSAX, VGBEEX,  DELZ,  WETEQ,
     &     WEXPB,  WET0,  PHEQ

      WET0  = AMAX1(WET,0.01)
      WEXPB = WET0**VGBEEX

      PHR   = VGPSAX / WEXPB
      SOILCO = VGCSAX * WEXPB * WEXPB * WET0 * WET0 * WET0

      PHEQ = PHR - DELZ
      WETEQ = ( PHEQ/VGPSAX ) ** ( -1/VGBEEX )

C 100  CONTINUE


      RETURN
      END
CCCC *
CCCC * [ END SOIL ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN M_FLUXES ]
CCCC *
      SUBROUTINE M_FLUXES (DTSTEP, ESATTC, DESDTC,   ALHX,
     I                    ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,
     I                       RC, DRCDEA, DRCDTC,
     I                    SWNET, HLWDWN, ALWRAD, BLWRAD, DTCANAL,ESNFRC,
     I                       TM,     EM,  CSOIL,   PSUR, EMAXRT,   WMAX,
     U                       TC,     TD,     EA,  SWET1,   SNOW,
     O                   RUNOFF,   EVAP, SHFLUX,  SMELT,  HLWUP, BOMB,
     O                   GHFLUX,III
     &                  )
CCCC *
CCCC * This subroutine computes the fluxes of latent and sensible heat
CCCC * from the surface through an energy balance calculation.
CCCC *
      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER  FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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

CCCC *

      REAL       DTSTEP, ESATTC, DESDTC,   ALHX,
     &       ETURB,  DEDEA,  DEDTC,
     &      HSTURB, DHSDEA, DHSDTC,
     &          RC, DRCDEA, DRCDTC,
     &       SWNET, HLWDWN, ALWRAD, BLWRAD,DTCANAL,
     &          TM,     EM,  CSOIL,   PSUR,
     &      EMAXRT,         WMAX(NLAY),
     &          TC,     TD,     EA, ESNFRC,
     &               SWET1(NLAY),   SNOW,
     &      RUNOFF,   EVAP, SHFLUX,  SMELT,
     &       HLWUP,   BOMB, GHFLUX
      INTEGER III       !the actual tile #, can be removed
      REAL     HLWTC,  CDEEPS,      Q0,  RHOAIR,   CONST,  DHLWTC,
     &        EPLANT,     A11,     A12,     A21,     A22,      F0,
     &           DEA,     DTC,  SNLEFT,     Q0X,  Q0SNOW,
     &         EANEW,  ESATNW,  EHARMN, DETERM, DENOM

      LOGICAL*1 DEBUG, CHOKE
      DATA DEBUG /.FALSE./
CCCC *
CCCC * -------------------------------------------------------------------
      
CCCC *
      HLWTC = ALWRAD + BLWRAD * TC
      CDEEPS = PIE * CSOIL * 2. / 86400.
      RHOAIR = PSUR * 100. / (RGAS * TC)
      CONST = RHOAIR * EPSILON / PSUR
      DHLWTC = BLWRAD
CCCC *
CCCC * Compute matrix elements A11, A22, AND Q0 (energy balance equation).
CCCC *
	     
      A11 = CSOIL/DTSTEP +
     &        DHLWTC +
     &        DHSDTC +
     &        ALHX*DEDTC +
     &        CDEEPS
	    
      A12 = DHSDEA + ALHX * DEDEA
	 
	     Q0 =  SWNET +
     &        HLWDWN -
     &        HLWTC -
     &        HSTURB -
     &        ALHX * ETURB -
     &        CDEEPS * (TC - TD)+
     &        CSOIL/DTSTEP*DTCANAL
     
CCCC *
CCCC * Compute matrix elements A21, A22, and F0 (canopy water budget  
CCCC * equation) and solve for fluxes.  Three cases are considered:
CCCC *
CCCC * 1. Standard case: RC>0.
CCCC * 2. RC = 0.  Can only occur if CIR is full or ETURB is negative.
CCCC *
      CHOKE = .TRUE.

      IF( RC .GT. 0.) THEN 
          EPLANT = CONST * (ESATTC - EA) / RC
          IF(EPLANT*ETURB.GT.0.) THEN
              EHARMN = 2.*EPLANT*ETURB / (EPLANT + ETURB)
            ELSE
              EHARMN=0.
            ENDIF
CCCC *
CCCC *            Some limitations to A21 and A22 are applied:
CCCC *            we assume that the increase in plant evaporation
CCCC *            due to an increase in either TC or EA balances 
CCCC *            or outweighs any decrease due to RC changes.
CCCC *

          A21 =  -DEDTC*RC +
     &      amax1(0., CONST*DESDTC - EHARMN*DRCDTC )
          A22 = -( RC*DEDEA +
     &               amax1( 0., CONST + EHARMN*DRCDEA )   )

          F0 = RC * (ETURB - EPLANT)
          DETERM = AMIN1( A12*A21/(A11*A22) - 1., -0.1 )
          DEA = ( Q0*A21 - A11*F0 ) / ( DETERM * A11*A22 )
          DTC = ( Q0 - A12*DEA ) / A11
          EVAP = ETURB + DEDEA*DEA + DEDTC*DTC
          SHFLUX = HSTURB + DHSDEA*DEA + DHSDTC*DTC
          DENOM = DETERM * A11*A22
        ELSE
          CHOKE = .FALSE.
          A21 = -DESDTC
          A22 = 1.
          F0 = ESATTC - EA
          DEA = ( Q0*A21 - A11*F0 ) / ( A12*A21 - A11*A22 )
          DTC = ( Q0 - A12*DEA ) / A11
          EVAP = ETURB + DEDEA*DEA + DEDTC*DTC
          SHFLUX = HSTURB + DHSDEA*DEA + DHSDTC*DTC
          DENOM = A12 * A21 - A11*A22
        ENDIF


CCCC * - - - - - - - - - - - - - - - - - - - -
CCCC * Account for snowmelt, if necessary.
CCCC * - - - - - - - - - - - - - - - - - - - -
      SMELT=0.
      SNLEFT=SNOW-EVAP*DTSTEP*ESNFRC
      IF(TC+DTC .GT. TF .AND. SNLEFT.GT.0.) THEN

CCCC *      First re-calculate energy balance under assumption that all
CCCC *      snow is melted.

        Q0X=Q0-ALHM*SNLEFT/DTSTEP
        SMELT=SNLEFT/DTSTEP
      
            DEA = ( Q0X*A21 - A11*F0 ) / DENOM
            DTC = ( Q0X - A12*DEA ) / A11

CCCC *      If TC+DTC is now less than TF, too much snow has melted.  Now
CCCC *      solve for balance assuming only some of the snow has melted.
CCCC *      Set TC to TF.

        IF(TC+DTC .LT. TF) THEN
          DTC=TF-TC
	
          DEA=(F0-A21*DTC)/A22
          Q0SNOW=A11*DTC+A12*DEA
          SMELT=(Q0-Q0SNOW)/ALHM
          ENDIF

        EVAP = ETURB + DEDEA*DEA +
     &                   DEDTC*DTC
        SHFLUX = HSTURB + DHSDEA*DEA +
     &                   DHSDTC*DTC

        ENDIF

CCCC * - - - - - - - - - - - - - - - - - - - - - - -
CCCC * Adjustments

CCCC * 1. Adjust deltas and fluxes if all available water evaporates
CCCC *    during time step:
CCCC *
      IF( EVAP .GT. EMAXRT ) THEN
        CHOKE = .FALSE.
        DEA = EM - EA
        DTC = 
     &   (Q0 + ALHX*(ETURB-EMAXRT) - DHSDEA*DEA) 
     &          /  ( A11 - ALHX*DEDTC )
        EVAP = EMAXRT
        SHFLUX = HSTURB + DHSDEA*DEA +
     &                    DHSDTC*DTC
        SMELT=0.
        ENDIF
CCCC *
CCCC * Adjust DEA and DTC if solutions were pathological:
CCCC *
      ESATNW = ESATTC+DESDTC*DTC
      EANEW = EA + DEA



CCCC * 2. Pathological cases. 

CCCC * Case 1: EVAP is positive in presence of negative gradient.
CCCC * Case 2: EVAP and ETURB have opposite sign, implying that 
CCCC * "virtual effects" derivatives are meaningless and thus that we
CCCC * don't know the proper tendency terms.
CCCC * In both cases, assume zero evaporation for the time step.

C      IF( ( EVAP(CHNO) .LT. 0. .AND. EM(CHNO).LT.ESATNW )
C     &    .OR. ( EVAP(CHNO)*ETURB(CHNO) .LT. 0. ) )  THEN 
      IF( EVAP .LT. 0. .AND. EM.LT.ESATNW )  THEN 
        CHOKE = .FALSE.
        DEA = EM - EA
        DTC = ( Q0 + ALHX*ETURB - DHSDEA*DEA ) / 
     &            ( A11 - ALHX*DEDTC )
        EVAP = 0.
        SHFLUX = HSTURB + DHSDEA*DEA +
     &                                            DHSDTC*DTC

        IF(TC+DTC .GT. TF .AND. SNOW.GT.0.) THEN
          Q0X=Q0-ALHM*SNOW/DTSTEP
          SMELT=SNOW/DTSTEP
          DTC = ( Q0X + ALHX*ETURB - DHSDEA*DEA ) / 
     &              ( A11 - ALHX*DEDTC )
          IF(TC+DTC .LT. TF) THEN
            DTC=TF-TC
	
            Q0SNOW=A11*DTC+A12*DEA
            SMELT=(Q0-Q0SNOW)/ALHM
            ENDIF
          SHFLUX = HSTURB + DHSDEA*DEA +
     &                                            DHSDTC*DTC
          ENDIF
        ENDIF



CCCC * 3. Exceesive dea change: apply "choke".

      IF( CHOKE .AND. ABS(DEA) .GT. 0.5*EA ) THEN
        DEA = SIGN(.5*EA,DEA)
        DTC = ( Q0 - A12*DEA ) / A11
        EVAP = ETURB + DEDEA*DEA + DEDTC*DTC
        SHFLUX = HSTURB + DHSDEA*DEA +
     &                                            DHSDTC*DTC
        IF(TC+DTC .GT. TF .AND. SNOW.GT.0.) THEN
          SNLEFT=SNOW-EVAP*DTSTEP*ESNFRC
          Q0X=Q0-ALHM*SNLEFT/DTSTEP
          SMELT=SNLEFT/DTSTEP
          DTC = ( Q0X - A12*DEA ) / A11
          IF(TC+DTC .LT. TF) THEN
            DTC=TF-TC
	  
            Q0SNOW=A11*DTC+A12*DEA
            SMELT=(Q0-Q0SNOW)/ALHM
            ENDIF
          EVAP = ETURB + DEDEA*DEA + DEDTC*DTC
          SHFLUX = HSTURB + DHSDEA*DEA +
     &                                            DHSDTC*DTC
          ENDIF
        ENDIF

CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
      TC = TC + DTC
      EA = EA + DEA
      TD = TD + DTSTEP * CDEEPS * (TC - TD) / (CSOIL*67.73)
      HLWUP = HLWTC + DHLWTC*DTC
      SNOW=SNOW-SMELT*DTSTEP
      HLWTC = ALWRAD + BLWRAD * TC
      GHFLUX=SWNET+HLWDWN-HLWTC-ALHX*EVAP-SHFLUX-SMELT*ALHM

CCCC * Make sure EA remains positive
      EA = AMAX1(EA, 0.0)

  200 CONTINUE



      RETURN
      END
CCCC *
CCCC * [ END M_FLUXES ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN M_VPDFAC ]
CCCC *
      SUBROUTINE M_VPDFAC ( ESATTC, EA,
     O                      VPDSTR,
     &                      VGDFAC        )
CCCC *
CCCC * This subroutine computes the vapor pressure deficit stress.
CCCC *
      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER  FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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
CCCC *
      
      REAL    ESATTC, EA,  VPDSTR
      REAL  VGDFAC
CCCC *
CCCC * -----------------------------------------------------------------

CCCC *
c      VPDSTR = 1. - (ESATTC-EA) * VGDFAC
c      VPDSTR  = AMIN1( 1., AMAX1( VPDSTR, 1.E-10 ) )
      VPDSTR = 1.
CCCC *
 100  CONTINUE
CCCC *
      RETURN
      END
CCCC *
CCCC * [ END M_VPDFAC ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN M_TMPFAC ]
CCCC *
      SUBROUTINE M_TMPFAC ( TC,
     O                   FTEMP,
     &      VGTLL, VGTU,
     &      VGTCF1, VGTCF2, VGTCF3, III)
CCCC *
CCCC * Compute temperature stress factor.
CCCC *
CCCC * Removing all "MemFac" references
      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER  FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY
	INTEGER III   !the tile #

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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
CCCC *
C
      REAL TC, FTEMP
      REAL  VGTLL, VGTU
      REAL  VGTCF1, VGTCF2, VGTCF3

CCCC *
      FTEMP = (TC - VGTLL) *
     &              (TC - VGTU) *
     &                    ( VGTCF1*TC*TC +
     &                      VGTCF2*TC +
     &                      VGTCF3 )
      IF ( TC .LE. VGTLL .OR. TC .GE. VGTU )
     &      FTEMP  = 1.E-10
      FTEMP = AMIN1( 1., AMAX1( FTEMP, 1.E-10 ) )
CCCC *
 100  CONTINUE
CCCC *
      RETURN
      END
CCCC *
CCCC * [ END M_TMPFAC ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN SMFAC ]
CCCC *
      SUBROUTINE SMFAC ( ESATTC, EA, PHR,  SOILCO,
     I                  RCUN,   VPDSTR, FTEMP,  TC, PSUR, Z2,
     I                  RSOIL1, RSOIL2, VGPH1X, VGPH2X, VGRPLX,
     O                  RC, III
     &                 )
CCCC *
CCCC * This subroutine estimates RC after computing the 
CCCC * leaf water potential stress.
CCCC *
      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER  FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY
	Integer III   !the actual tile # we are on

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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
CCCC *
      INTEGER ITYP
      REAL  ESATTC, EA,     PHR,  SOILCO,
     &        RCUN,VPDSTR, FTEMP, TC,
     &        PSUR,  Z2,  RSOIL1, RSOIL2,
     &      VGPH1X, VGPH2X,  VGRPLX, RC
      REAL   RCUNTD,  RHOAIR,   CONST,     DEF,     D12,     DR2,
     &        RSOIL,      R0,    EEST,  RSTFAC
CCCC *
CCCC * -----------------------------------------------------------------

CCCC *
      RCUNTD = RCUN / ( VPDSTR*FTEMP )
      RHOAIR = PSUR * 100. / ( RGAS*TC )

CCCC *
      CONST = RHOAIR * EPSILON / PSUR
      DEF = ( ESATTC - EA ) * CONST
      D12 = VGPH1X - VGPH2X
      DR2 = PHR - Z2 - VGPH2X

!     Modified to prevent error when SOILCO approaches zero.
!      RSOIL = RSOIL1 + RSOIL2/SOILCO
      RSOIL = RSOIL1 + RSOIL2/(SOILCO + 1.0E-20)

      R0 = ( VGRPLX + RSOIL ) / RHOW
      IF (ABS(RCUNTD*D12 + DEF*R0).LT.0.00000001) THEN
        PRINT *,'WARNING WARNING WARNING, problem at tile ',iii
        PRINT *,'RCUNTD*D12 + DEF*R0 too small ',
     &  RCUNTD*D12 + DEF*R0,' set to 0.00000001'
        EEST = DEF*DR2 / 0.00000001
      ELSE
        EEST = DEF*DR2 / ( RCUNTD*D12 + DEF*R0 )
      ENDIF
      RSTFAC = ( DR2 - R0*EEST ) / D12
      RSTFAC = AMIN1( 1., AMAX1( 0.001, RSTFAC ) )
      RC = RCUNTD / RSTFAC
CCCC *
 100  CONTINUE
CCCC *
      RETURN
      END
CCCC *
CCCC * [ END SMFAC ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN RSURFP ]
CCCC *
      SUBROUTINE RSURFP ( UM, U2FAC, Z2, RDC, WET, ESATTC, EA,
     U                   RC,
     O                   RX1, RX2, III
     &                  )
CCCC *
      IMPLICIT NONE

      REAL  UM,   U2FAC,     Z2,    RDC, 
     &      WET,  ESATTC,     EA,
     &      RC,     RX1,    RX2
      REAL  U2, RSURF, HESAT
	Integer III   !the actual tile number we are on
CCCC *
CCCC * -----------------------------------------------------------------

CCCC *
      U2 = UM * U2FAC
c      RSURF = RDC / U2 + 30. / (1.E-20 + WET)
      RSURF = RDC / U2 + 26. + 6. / (1.E-10 + WET)**2
CCCC * Account for subsaturated humidity at soil surface:
CCCC *
      HESAT = ESATTC * MIN( 1., WET*2. )
c	print *,'calculated hesat',hesat
      IF( EA .LT. HESAT ) THEN
c	print *,'in if loop to calc rsurf'
c	print *,'rsurf=',RSURF
c	print *,'ESATTC=',ESATTC,'HESAT=',HESAT
c	print *,'(ESATTC-HESAT)=',(ESATTC-HESAT)
c	print *,'(HESAT-EA)=',(HESAT-EA)
	IF ((HESAT-EA).NE.0) THEN 
c	print *,'hesat-ea not equal 0 so go on'
          RSURF=RSURF*( 1. + (ESATTC-HESAT)/(HESAT-EA) )
C	If denominator is 0, then set to 1E-37 to fix divide by
C	zero error
        ELSEIF ((HESAT-EA).EQ.0) THEN
	print *,'would have been zero, so set to 1E-37'
          RSURF=RSURF*( 1. + (ESATTC-HESAT)/(1E-37) )
	ENDIF
c	print *,'calculated rsurf',rsurf
        ELSE
c	print *,'setting rsurf=1.E10'
          RSURF=1.E10
        ENDIF

      RX1=RC
      RX2=RSURF

      RC = RC * RSURF / ( RC + RSURF )
CCCC *
 100  CONTINUE
CCCC *
      RETURN
      END
CCCC *
CCCC * [ END RSURFP ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN M_RCANOP ]
CCCC *
      SUBROUTINE M_RCANOP ( CAPAC, SNOW, SATCAP, RA, ETURB, SNWFRC,
     I                     POTFRC,
     U                     RC
     &                    )
CCCC *
CCCC * The effective latent heat resistance RC depends on the quantity 
CCCC * of interception reservoir water and the snow cover.  POTFRC
CCCC * is the fraction of the tile from which potential evaporation
CCCC * occurs.
CCCC *
      IMPLICIT NONE


      REAL CAPAC, SNOW, SATCAP, RA,
     &     ETURB, RC, SNWFRC, POTFRC
      REAL ETCRIT,RAMPFC

CCCC * (Note: ETCRIT arbitrarily set to ~-5 W/m2, or -2.e-6 mm/sec.)
      DATA ETCRIT/ -2.E-6 /
CCCC *
CCCC * -----------------------------------------------------------------

 
CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CCCC * Case 1: Vegetation present (SATCAP GT .001).  Potential evap. from
cCCC * both interception reservoir and snow.

      IF(SATCAP.GT..001) THEN
        RC=RC*(1.-POTFRC)/
     &                ( 1.+POTFRC*RC/RA )
        ENDIF

CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CCCC * Case 2: Vegetation absent (SATCAP LE .001).  Potential evap. from
cCCC * snow only.

      IF(SATCAP .LE. .001) THEN
        RC=RC*(1.-SNWFRC)/
     &                  ( 1.+SNWFRC*RC/RA )
        ENDIF

CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CCCC * Assume RC=0 for condensation (dew).
CCCC * RAMPFC is used to ensure continuity in RC.

      RAMPFC=ETURB/ETCRIT
      IF ( RAMPFC .GE. 0. ) RC = RC*(1.-RAMPFC)
      IF ( RAMPFC .GT. 1. ) RC = 0.
CCCC *
 100  CONTINUE
CCCC *
      RETURN
      END
CCCC *
CCCC * [ END M_RCANOP ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC ** [ BEGIN M_WUPDAT ]
CCCC *
      SUBROUTINE M_WUPDAT ( DTSTEP,
     I                     EVAP, SATCAP, VGWMAX, TC, RA, RC,
     I                     RX1, RX2,ESNFRC,EIRFRC,
     U                     CAPAC, SNOW, SWET, RUNOFF,
     O                     EINT, ESOI, EVEG, ESNO, RUNSRF, FWSOIL,
     &                   III )
CCCC *
CCCC * This subroutine allows evapotranspiration to adjust the water
CCCC * contents of the interception reservoir and the soil layers.
CCCC *
      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER III
      INTEGER  FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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
CCCC *

      REAL    EVAP,  SATCAP,  VGWMAX(NLAY),
     &          TC,      RA,             RC,
     &       CAPAC,    SNOW,    SWET(NLAY),
     &      RUNOFF,     RX1,
     &         RX2,  RUNSRF,  FWSOIL,
     &      ESNFRC,  EIRFRC
      REAL EINT, ESOI, EVEG, ESNO
      REAL  DTSTEP, EGRO, FWS, THRU, DEWRUN,
     &     WTOTAL,WLAY1,WLAY2,ELAY1,ELAY2,EGROI
CCCC *
CCCC * -----------------------------------------------------------------

CCCC *
CCCC * Partition evap between interception, snow, and ground reservoirs.
CCCC *
      WLAY1 = SWET(SFCLY) * VGWMAX(SFCLY)
      WLAY2 = SWET(ROOTLY) * VGWMAX(ROOTLY)
      WTOTAL = WLAY1 + WLAY2
CCCC *
      ESNO=ESNFRC*EVAP*DTSTEP
      EINT=EIRFRC*EVAP*DTSTEP
      EGRO = EVAP*DTSTEP - ESNO - EINT

CCCC * Ensure that individual capacities are not exceeded.

      IF(ESNO .GT. SNOW) THEN
        EINT=EINT+(ESNO-SNOW)
        ESNO=SNOW
        ENDIF
      EGROI=EGRO+EINT
      IF(EGROI .GT. CAPAC+WTOTAL) THEN
        ESNO=ESNO+EGROI-(CAPAC+WTOTAL)
        EGROI=CAPAC+WTOTAL
        ENDIF

      EINT=EGROI-EGRO
      IF(EINT .GT. CAPAC) THEN
        EGRO=EGRO+EINT-CAPAC
        EINT=CAPAC
        ENDIF
      IF(EGRO .GT. WTOTAL) THEN
        EINT=EINT+EGRO-WTOTAL
        EGRO=WTOTAL
        ENDIF
CCCC *
CCCC * Separate egro into surface-evaporation/transpiration components:
CCCC *
      IF( RX1+RX2 .NE. 0. ) THEN
          ESOI=EGRO*RX1/(RX1+RX2)
          EVEG=EGRO - ESOI
        ELSE
          ESOI=EGRO/2.
          EVEG=EGRO/2.
        ENDIF
CCCC *
CCCC * Translate ESOI and EVEG into evaporation fluxes from layers 1 and 2:
CCCC *
      IF( WTOTAL .GT. 0. ) THEN
          FWS = EVEG / WTOTAL 
        ELSE
          FWS = 1.
        ENDIF
      ELAY1 = WLAY1*FWS + ESOI
      ELAY2 = WLAY2*FWS
CCCC *
CCCC * Ensure that enough soil water is available in each layer:
CCCC *
      IF( ELAY1 .GT. WLAY1 ) THEN
        ELAY2 = ELAY2 + (ELAY1 - WLAY1)
        ELAY1 = WLAY1
        ENDIF
      IF( ELAY2 .GT. WLAY2 ) THEN
        ELAY1 = ELAY1 + (ELAY2 - WLAY2)
        ELAY2 = WLAY2
        ENDIF

CCCC * Special case for condensation:
      IF(EVAP .LT. 0.) THEN
          EINT=(1.-ESNFRC)*EVAP*DTSTEP
          ELAY1=0.
          ELAY2=0.
          ESOI=0.
          EVEG=0.
          ESNO=ESNFRC*EVAP*DTSTEP
        ENDIF

CCCC *
CCCC * Remove moisture from reservoirs:
CCCC *
      SNOW = SNOW - ESNO
      CAPAC = CAPAC - EINT
      SWET(SFCLY) = (WLAY1 - ELAY1) / VGWMAX(SFCLY)
      SWET(ROOTLY) = (WLAY2 - ELAY2) / VGWMAX(ROOTLY)
CCCC *
CCCC * Ensure against numerical precision problems:
      SWET(SFCLY) = AMIN1( 1., AMAX1( 0., SWET(SFCLY) ) )
      SWET(ROOTLY) = AMIN1( 1., AMAX1( 0., SWET(ROOTLY) ) )
CCCC *
CCCC *
CCCC * -------------------------------------------------
CCCC * DEWFALL:
CCCC * 
CCCC * If dewfall adds to cir, insure that it doesn't fill 
CCCC * beyond capacity.  If resulting throughfall adds to top soil layer, 
CCCC * insure that it also doesn't fill beyond capacity.
CCCC * 
      IF( CAPAC .GT. SATCAP ) THEN
        THRU = CAPAC - SATCAP
        CAPAC = SATCAP
        SWET(SFCLY) = SWET(SFCLY) + 
     &                         THRU / VGWMAX(SFCLY)
        DEWRUN=0.
        IF ( SWET(SFCLY) .GT. 1. ) THEN
          DEWRUN = ( SWET(SFCLY) - 1. ) * VGWMAX(SFCLY)
          SWET(SFCLY) = 1.
          RUNOFF = RUNOFF + DEWRUN/DTSTEP
          RUNSRF = RUNSRF + DEWRUN/DTSTEP
          ENDIF
        FWSOIL = FWSOIL + (THRU-DEWRUN)/DTSTEP     
        ENDIF
CCCC *
 100  CONTINUE
CCCC *
      RETURN
      END
CCCC *
CCCC * [ END M_WUPDAT ]
CCCC *
CCCC * -----------------------------------------------------------------
CCCC * /////////////////////////////////////////////////////////////////
CCCC * -----------------------------------------------------------------
CCCC *
CCCC * [ BEGIN GWATER ]
CCCC *
      SUBROUTINE GWATER (WSMAX,  PHLAY, AKLAY, TC,
     I                     DTSTEP, VGZDEX, VGSLOX, WETEQ1, WETEQ2,
     I                     PHSAT,   AKSAT,
     U                     SWET,   RUNOFF,  GDRAIN,
     &                 III   )
CCCC *
CCCC * This subroutine computes diffusion between soil layers.
CCCC *
      IMPLICIT NONE
CCCC *
CCCC * CHIP HEADER FILE
CCCC *
      INTEGER   FRSTCH
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY
      INTEGER III
      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW

      PARAMETER (FRSTCH = 1)
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
CCCC *

      REAL  VGSLOX, RUNOFF, GDRAIN
      REAL  ZDEP12,   AKAVE,  GWFLUX,  ZDEP23,  HALFMX,  DHDZ,
     &    FAREA

      REAL   WSMAX(NLAY),   PHLAY(NLAY),
     &       AKLAY(NLAY),             TC,
     &                  DTSTEP,    SWET(NLAY),
     &      VGZDEX(NLAY),  WETEQ1,  WETEQ2,
     &      PHSAT, AKSAT

CCCC *
CCCC * ----------------------------------------------------------------
  
CCCC *
CCCC * Diffusion between layer 1 and 2:
      ZDEP12 = VGZDEX(SFCLY) + VGZDEX(ROOTLY)
      IF( SWET(SFCLY) .GT. WETEQ1 ) THEN
          FAREA=(SWET(SFCLY) - WETEQ1) / (1. - WETEQ1)
          DHDZ = 1. + 2.*(PHSAT-PHLAY(ROOTLY))/ZDEP12
          AKAVE = AKSAT
        ELSE
          FAREA = 1.
          DHDZ = 1. + 2.*(PHLAY(SFCLY)-PHLAY(ROOTLY))/ZDEP12
          AKAVE = AKLAY(ROOTLY)
        ENDIF
CCCC *
      GWFLUX = 1000. * DTSTEP * AKAVE * DHDZ * FAREA
CCCC *
CCCC * Test for limits on water holding capacity.
CCCC *
      HALFMX=0.5*ABS( SWET(SFCLY)-WETEQ1 )
     &                               * WSMAX(SFCLY)
      IF (GWFLUX .GE. 0.) THEN
          GWFLUX = AMIN1( GWFLUX, 
     &               SWET(SFCLY) * WSMAX(SFCLY) )
          GWFLUX = AMIN1( GWFLUX,    WSMAX(ROOTLY) -
     &              SWET(ROOTLY) * WSMAX(ROOTLY) )
          GWFLUX = AMIN1( GWFLUX, HALFMX )
        ELSE
          GWFLUX = -AMIN1( ABS(GWFLUX),
     &                SWET(ROOTLY) * WSMAX(ROOTLY) )
          GWFLUX = -AMIN1( ABS(GWFLUX),   WSMAX(SFCLY) -
     &                SWET(SFCLY) * WSMAX(SFCLY) )
          GWFLUX = -AMIN1( ABS(GWFLUX),  HALFMX )
        ENDIF
CCCC *
CCCC * Prevent diffusion when ground is frozen:
      IF (TC .LT. TF) GWFLUX = 0.
CCCC *
       
CCCC * Update water contents
      SWET(SFCLY) = SWET(SFCLY) -
     &                        GWFLUX / WSMAX(SFCLY)
      SWET(ROOTLY) = SWET(ROOTLY) +
     &                        GWFLUX / WSMAX(ROOTLY)
CCCC *
CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CCCC * Diffusion between root and recharge layers:
      ZDEP23 = VGZDEX(ROOTLY) + VGZDEX(RECHLY)
      DHDZ=1. + 2. * (PHLAY(ROOTLY) - PHLAY(RECHLY)) / ZDEP23
      IF(DHDZ.GE.0.) THEN
          AKAVE = AKLAY(ROOTLY)
        ELSE
          AKAVE = AKLAY(RECHLY)
        ENDIF
      GWFLUX = 1000. * DTSTEP * AKAVE * DHDZ
CCCC *
CCCC * Test for limits on water holding capacity:
      HALFMX=0.5*ABS( SWET(ROOTLY)-WETEQ2 )
     &                               * WSMAX(ROOTLY)
      IF (GWFLUX .GE. 0.) THEN
          GWFLUX = AMIN1( GWFLUX,
     &                    SWET(ROOTLY) * WSMAX(ROOTLY) )
          GWFLUX = AMIN1( GWFLUX,    WSMAX(RECHLY) -
     &                    SWET(RECHLY) * WSMAX(RECHLY) )
          GWFLUX = AMIN1( GWFLUX, HALFMX )
        ELSE
          GWFLUX = -AMIN1( ABS(GWFLUX),
     &                    SWET(RECHLY) * WSMAX(RECHLY) )
          GWFLUX = -AMIN1( ABS(GWFLUX), WSMAX(ROOTLY) -
     &                    SWET(ROOTLY) * WSMAX(ROOTLY) )
          GWFLUX = -AMIN1( ABS(GWFLUX), HALFMX )
        ENDIF
CCCC *
CCCC * Prevent diffusion when ground is frozen:
      IF(TC .LE. TF) GWFLUX = 0.
CCCC *
CCCC * Update water contents
      SWET(ROOTLY) = SWET(ROOTLY) -
     &                            GWFLUX / WSMAX(ROOTLY)
      SWET(RECHLY) = SWET(RECHLY) +
     &                            GWFLUX / WSMAX(RECHLY)
CCCC *
CCCC * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CCCC * Flux of moisture out of lower soil layer to water table
CCCC * (approximation to SiB)
CCCC *
      GWFLUX = VGSLOX * AKLAY(RECHLY) * 1000. * DTSTEP
      IF (TC .LE. TF) GWFLUX = 0.
      GWFLUX = AMIN1( GWFLUX,
     &                 SWET(RECHLY) * WSMAX(RECHLY) )
      SWET(RECHLY) = SWET(RECHLY) -
     &                                GWFLUX / WSMAX(RECHLY)
CCCC *
      GDRAIN=GWFLUX/DTSTEP
      RUNOFF  = RUNOFF  + GWFLUX/DTSTEP
CCCC *
 100  CONTINUE
CCCC *
      RETURN
      END
CCCC *
CCCC * [ END GWATER ]







CCCC *
      SUBROUTINE M_PMONTH (JDAY,ALAT,
     I      GREEN, LAI,
     O      Z0, SQSCAT, Z2, DZM,
     O      RSOIL1, RSOIL2, SATCAP, RDC,U2FAC,
     &      VGRDRS,VGZ2,VGROTD,VGROCA,
     &      VGRDC,VGROTL,VGZ0,VGDD,
     &      VGTR12,VGTR11,VGRF11,VGRF12)

      IMPLICIT NONE

cCCC *     This subroutine sets seasonally varying
cCCC *  variables.  NOTE: LAI is stored as LAI (SIB) / VCOVER (SIB).
cCCC *
cCCC *
      INTEGER  JDAY
cCCC *
cCCC *
      REAL GREEN,  LAI,    SQSCAT, Z0,
     &     Z2,     D,      DZM,    RSOIL1,
     &     RSOIL2, SATCAP, RDC,    U2FAC,
     &     ALAT
cCCC *
      INTEGER ChNo, K, KDAY, KLOW
cCCC *
CC These are now passed in so removed is the 12,10 
      REAL VGRDC,  VGROTL,
     &     VGZ0,   VGDD
CC These are now passed in so removed is (10)   
      REAL VGRF11,VGRF12,VGTR11,VGTR12
C  These are now passed in so they are length 1
      REAL  VGROCA,VGRDRS,VGROTD,VGZ2
cCCC *
      REAL  ALPHAF, CUNI, DUM1, DUM2, ROOTL, SCAT,
     &   VKC, VROOT, WHI, WLO, RLENGT, DPARAM
cCCC *
cCCC * --------------------------------------------------------------------

      DATA VKC /0.35/

CCCC *
CCCC *
CCCC * COMPUTE SOME PARAMETERS DIRECTLY (this is already done
CCCC * In the interpolation stuff up at the Sub.driver.f level
        RDC  =  VGRDC 
        ROOTL = VGROTL
        RLENGT=VGZ0	  
        Z0=RLENGT
        DPARAM=VGDD
CCC
CCC
        VROOT = ROOTL * VGROCA 
        DUM1 = ALOG (VROOT / (1 - VROOT))
        DUM2 = 1. / (8. * 3.14159 * ROOTL)
        ALPHAF = DUM2 * (VROOT - 3. -2. * DUM1)
        RSOIL1  =
     &    VGRDRS  / (ROOTL * VGROTD )
        RSOIL2  = ALPHAF / VGROTD
CCCC *
        SCAT = GREEN  *
     &         (VGTR11  + VGRF11 ) +
     &         (1. - GREEN ) *
     &         (VGTR12  + VGRF12 )
        SQSCAT  = SQRT (1. - SCAT)
CCCC *
        DZM  = 600.

CCCC *  
        CUNI = (ALOG (0.025 * DZM  / RLENGT) / VKC) + 8.4
        U2FAC  =
     &    ALOG ((VGZ2  - DPARAM) / RLENGT) /
     &    (CUNI * VKC)
        Z2  = VGZ2 
CCCC *
        SATCAP  = 0.2 * LAI
 20   CONTINUE
CCCC *
      RETURN
      END
CCCC *

* [ END M_PMONTH ]
cCCC *

 


      SUBROUTINE ASTRO (YEAR,MONTH,DAY,SEC,ALAT,ALON,IRUN,COSZ,RA)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  INPUT:
C  ======
C    YEAR,MONTH,DAY,SEC      : CURRENT YYMMDD
C    ALAT:LATITUDES  IN DEGREES.
C    ALON:LONGITUDES IN DEGREES. (0 = GREENWICH, + = EAST).
C    IRUN      : # OF POINTS TO CALCULATE
C
C  OUTPUT:
C  =======
C    COSZ  : COSINE OF ZENITH ANGLE.
C    RA          : EARTH-SUN DISTANCE IN UNITS OF
C                  THE ORBITS SEMI-MAJOR AXIS.
C
C  NOTE:
C  =====
C  THE INSOLATION AT THE TOP OF THE ATMOSPHERE IS:
C
C  S = (SOLAR CONSTANT)*(1/RA**2)*COSZ,
C
C  WHERE:
C  RA AND COSZ ARE THE TWO OUTPUTS OF THIS SUBROUTINE.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C*                  GODDARD LABORATORY FOR ATMOSPHERES               C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

      implicit none

c Input Variables
c ---------------
      integer irun
      real    cosz, alat, alon, ra

c Local Variables
c ---------------
      integer year, day, sec, month, iday, idayp1
      integer dayscy
      integer i,n,nsecf,k,km,kp

      real hc
      real pi, zero, one, two, six, dg2rd, yrlen, eqnx, ob, ecc, per
      real daylen, fac, thm, thp, thnow, zs, zc, sj, cj

      parameter ( pi    = 3.1415926535898)
      parameter ( zero  = 0.0 )
      parameter ( one   = 1.0 )
      parameter ( two   = 2.0 )
      parameter ( six   = 6.0 )
      parameter ( dg2rd = pi/180. )

      parameter ( yrlen  = 365.25  )
      parameter ( dayscy = 365*4+1 )
      parameter ( eqnx   =  80.9028)
      parameter ( ob     =  23.45*dg2rd )
      parameter ( ecc    =   0.0167 )
      parameter ( per    = 102.0*dg2rd)
      parameter ( daylen = 86400.)

C      REAL   TH(DAYSCY),T0,T1,T2,T3,T4,FUN,Y
C      REAL MNDY(12,4)
c     REAL MNDY(48,1)	
	REAL TH(DAYSCY),T0,T1,T2,T3,T4,FUN,Y,MNDY(48,1)

      LOGICAL*1 FIRST
      DATA    FIRST/.TRUE./
      SAVE

      DATA MNDY /0,31,60,91,121,152,182,213,244,274,305,335,366,
     .           397,34*0 /

c      DATA MNDY /0,31,60,91,121,152,182,213,244,274,305,335,366,
c     .           397,22*0 /


      FUN(Y) = (TWO*PI/((ONE-ECC**2)**1.5))*(ONE/YRLEN)
     .       * (ONE - ECC*COS(Y-PER)) ** 2
C      NSECF(N) = N/10000*3600 + MOD(N,10000)/100* 60 + MOD(N,100)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                 COMPUTE DAY-ANGLES FOR 4-YEAR CYCLE              C  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

      IF(FIRST) THEN
           DO 100 I=15,48
           MNDY(I,1) = MNDY(I-12,1) + 365
100        CONTINUE

           KM  = INT(EQNX) + 1
           FAC = KM-EQNX
           T0 = ZERO
           T1 = FUN(T0         )*FAC
           T2 = FUN(ZERO+T1/TWO)*FAC
           T3 = FUN(ZERO+T2/TWO)*FAC
           T4 = FUN(ZERO+T3    )*FAC
           TH(KM) = (T1 + TWO*(T2 + T3) + T4) / SIX

           DO 200 K=2,DAYSCY
           T1 = FUN(TH(KM)       )
           T2 = FUN(TH(KM)+T1/TWO)
           T3 = FUN(TH(KM)+T2/TWO)
           T4 = FUN(TH(KM)+T3    )
           KP = MOD(KM,DAYSCY) + 1
           TH(KP) = TH(KM) + (T1 + TWO*(T2 + T3) + T4) / SIX
           KM = KP
 200       CONTINUE

           FIRST=.FALSE.
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
CCC            COMPUTE EARTH-SUN DISTANCE TO CURRENT SECOND          C  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c      IDAY   = DAY + MNDY(MONTH,MOD(YEAR,4)+1)

      IDAY   = DAY + MNDY(MONTH*(MOD(YEAR,4)+1),1)
      IDAYP1 = MOD( IDAY,DAYSCY) + 1
      THM    = MOD( TH(IDAY)  ,TWO*PI)
      THP    = MOD( TH(IDAYP1),TWO*PI)

      IF(THP.LT.THM) THP = THP + TWO*PI
      FAC   = FLOAT(SEC)/DAYLEN
      THNOW = THM*(ONE-FAC) + THP*FAC

      ZS = SIN(THNOW) * SIN(OB)
      ZC = SQRT(ONE-ZS*ZS)
      RA = (1.-ECC*ECC) / ( ONE-ECC*COS(THNOW-PER) )

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                 COMPUTE COSINE OF THE ZENITH ANGLE               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

      FAC  = FAC*TWO*PI + PI

      HC = COS( FAC+ALON*DG2RD )
      SJ = SIN(ALAT*DG2RD)
      CJ = SQRT(ONE-SJ*SJ)

          COSZ = SJ*ZS + CJ*ZC*HC
      IF( COSZ.LT.ZERO ) COSZ = ZERO

      RETURN
      END






CCCC * ------------------------------------------------------------------
CCCC * //////////////////////////////////////////////////////////////////
CCCC * ------------------------------------------------------------------
CCCC *
      REAL FUNCTION M_QSAT(T,PR,ALHX)
CCCC * 
      IMPLICIT NONE

      REAL  T, PR, ALHX, C1, C2, C3, ALHE, ALHS, ALHM, C1LOG
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      DATA C1/3.797915/, C2/7.93252E-6/, C3/2.166847E-3/,
     &           C1LOG/1.33445/
CCCC *
c     M_QSAT = C1*EXP(ALHX*(C2-C3/T))/PR
      M_QSAT = C1*EXP((ALHX/ALHE)*(21.18123-C1LOG-5418./T))/PR
CCCC *
      RETURN
      END


CCCC * ------------------------------------------------------------------
CCCC * //////////////////////////////////////////////////////////////////
CCCC * ------------------------------------------------------------------
CCCC *

c****
      real function m_esat(t)
c****
      implicit none
      integer  NT, frstch, memfac
      integer  nlay, sfcly, rootly, rechly

      real  zero, one, pie
      real alhe, alhs, alhm, tf, stefan, rgas, shw, shi, rhow, grav
      real epsilon, nosnow

      parameter (NT = 10, frstch = 1, memfac = 5)
      parameter (nlay = 3)
      parameter (sfcly = 1, rootly = sfcly + 1, rechly = rootly + 1)

      parameter (zero = 0., one = 1., pie = 3.14159265)
      parameter (alhe = 2.4548e6, alhs = 2.8368e6, alhm = alhs-alhe)
      parameter (tf = 273.16)
      parameter (stefan = 5.669e-8)
      parameter (rgas = .286*1003.5)
      parameter (shw = 4200., shi = 2060.)
      parameter (rhow = 1000.)
      parameter (grav = 9.81)
      parameter (epsilon = 18.01/28.97)
      parameter (nosnow = 0.)

      real  t
c****
      m_esat = exp(21.18123 - 5418./t) / epsilon
c****
      return
      end

