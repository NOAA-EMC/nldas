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
! mos_main.f: 
!
! DESCRIPTION:
!  MOSAIC MODEL: this is the main 1-D code for the MOSAIC model.  
!
!  In the old LDAS driver there were middle steps (such as "Process" and "chip"
!   that were called for various data manipulation
!   it is my hope that all that can be done here and in an explained fashion
!
!  This Subroutine will be called for EACH tile
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  5  Nov 1999: Paul Houser; Significant F90 Revision
!  4  Apr 2000: Jeff Walker; Altered names of subroutines
! 11  Apr 2000: Brian Cosgrove; Removed duplicate delcarations of variables
! 12  May 2000: Brian Cosgrove/Jared Entin; Added code to prevent 
!               numerical instablility due to small values of wind
!               and humidity
! 22 Aug. 2000: Brian Cosgrove; Modified code for output of
!               standard LDAS output variables.  Added LDAS and MOS
!               modules into call for MOSTILE and also added many
!               new variables in the MOS%RETURN section 
! 25 Aug. 2000: Brian Cosgrove; Fixed snowcover fraction output variable
! 21 Sep. 2000: Brian Cosgrove; Fixed code so that rain output field
!               is set to zero when snow output field is greater than zero
! 27 Sep. 2000: Brian Cosgrove; Fixed code so that 1 meter soil moisture
!               is correctly calculated (FACTOR variable had been
!               incorrectly computed before this fix)
! 23 Mar. 2001: Jon Radakovich; Updated for PSAS temperature assimilation
! 19 Apr. 2001  Updated scheme for calculating albedo, using sibalb module
! 05 Sep. 2001: Brian Cosgrove; Added in volumetric variables soilwm,soilww
!               for use in LDAS output.  Zthick changed from 50 to 10 sometime
!               by someone else.  Added volumetric output at 4 more levels
!               corresponding to OK mesonet levels.  Changed calculation
!               of soil wetness output variables
! 05 Feb. 2002: Brian Cosgrove; Changed Zthick back to 50 from 10 after talk
!               with Jon Radokovich.  10 is correct for observation height, but 
!               was causing fluxes that were too high, so changed back to
!               'incorrect' value of 50
! 22 Apr. 2002: Urszula Jambor; Added conditional to suppress calculation
!               of dewpt temp. if using Aaron Berg's reanalysis ECMWF data.
! 13 Sep. 2002: Urszula Jambor; Reversed suppression of dew point temp.
!               calculation if using Aaron Berg's reanalysis ECMWF data.
! 27 Nov. 2002: Urszula Jambor; Restricted mos%snow to max. of 100000.0
!               to prevent problems in GRIB output snow fields.
! 12 Dec. 2002: Brian Cosgrove; Fixed usage of Wiltpoint variable.  Before,
!               Wiltpoint1 and Wiltpoint2 were used in calculation of 
!               root zone soil moisture availability...now, only Wiltpoint2
!               is used since wiltpoint1 is not the correct wilting point
!               needed for the calculation.
! 13 May 2003:  Brian Cosgrove; Fixed Albedo output variable so that
!               it is now in units of percent, and not a fraction
!=========================================================================

      SUBROUTINE MOS_MAIN(TN,LDAS,TILE,MOS,SIB)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      USE mos_module       ! Mosaic tile variables
      use sibalb_module    ! SiB-based coefficients for albedo calculation.
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (tiledec) TILE   
      type (mosdec) MOS
      type (sibalbdec) sib

!=== Local Variables =====================================================
        integer kkk
        INTEGER FLAG            !Flag for top 1 meter soil mst derivation
                                !1=have value, 0=don't have value
        INTEGER K,J,I		!Loop counters
        INTEGER SEC	        !Current number of seconds into day
        INTEGER RESET		!Flag for TILE SUBROUTINE
				!1=convert ETA volumetric to %saturation
                                !  soil moisture because MOSAIC soil
                                !  moisture is being reset to ETA values
                                !  and unit conversion needs to occur
	INTEGER TN              !which tile we are actually on, can be removed
	REAL DEPTH		!Depth value used in top 1meter soil mst calcs
	REAL FACTOR	        !Factor by which to multiply soil moisture
				!in a layer to arrive at its contribution
				!to soil mst in top 1 meter
	REAL TEMP1,TEMP2        !Factors used in top 1 meter soil mst calcs
	REAL ALBEDO1,ALBEDO2
	REAL ALBEDO3,ALBEDO4
	REAL M_ESAT		!Function which uses DPT2M to
				!  compute EM.
        REAL CTT                ! 
	REAL M_QSAT		!Function 
	REAL SNWMID		!
	REAL ALHE
	REAL ALHS
	REAL ALHM
	REAL ESATTC
	REAL GRAV,RGAS,CP	!Constants
	REAL VKRMN,EPSI         !Constants
	REAL B,C,D,ZTHICK	!Constants
	REAL RHOAIR	        !Density of air
	REAL QM			!Parameter from EM,EPSI,PSUR
	REAL QC                 !Parameter from EA,EPSI,PSUR
        REAL DTV                !
        REAL RI                 !
	  REAL NegRI              !
        REAL DRIDTC             !
        REAL DRIDQC             !
        REAL DUMMY              !
        REAL CN                 !
        REAL T                  !
        REAL FU                 !
        REAL FT                 !
        REAL DFTDRI             !
        REAL CS                 !
        REAL R                  !
        REAL S                  !
        REAL DFTDQC             !
        REAL DFTDTC             !
        REAL HSTURB             !Sensible heat flux computed by GCM
	REAL ETURB		!Evaporation rate computed by GCM
	REAL DEDQA		!Derivative of evaporation rate
				!  w.r.t. specific humidity
	REAL DEDTC		!Derivative of evaporation rate
                                !  w.r.t. surface-canopy system
                                !  temperature
	REAL DHSDQA		!
	REAL DHSDTC		!Derivative of sensible heat flux
                                !  w.r.t. surface-canopy system
                                !  temperature
	REAL DEDEA		!
	REAL DHSDEA		!
	REAL AVISDR		!Visible direct snow albedo
	REAL ANIRDR		!Near infrared direct snow albedo
	REAL AVISDF		!Visible diffuse snow albedo
	REAL ANIRDF		!Near infrared diffuse snow albedo
	REAL ALBAVE		!Average snow albedo
	REAL SWNET		!Net shortwave radiation absorbed
				!  by the surface
	REAL PAR		!Set to 0.5*SWDN
	REAL PDIR		!Set to a value of 0.5
	REAL TOTALB		!Set equal to ALBAVES
	REAL ZLAI		!Set equal to LAI
	REAL QSATTC		!The saturated specific humidity
                                !  based on the surface-canopy 
                                !  system temperature
	REAL DQSDTC		!The derivative of the saturated
	 			!  specific humidity with respect
                                !  to the surface-canopy system
                                !  temperature
	REAL PARDIR		!The direct component of the 
				!  photosynthetically active
                                !  radiation (W/m2)
 	REAL PARDIF		!The diffuse component of the 
                                !  photosynthetically active
                                !  radiation (W/m2)
	REAL CD			!Drag coefficient
	REAL LWATTC		!
	REAL DLWDTC		!
	REAL ALWRAD		!First term in longwave radiation
                                !  linearization (W/m2)
	REAL BLWRAD		!Second term in longwave radiation
                                !  linearization (W/m2K)
	REAL EVAP		!Evaporation rate (W/m2)
	REAL SHFLUX		!Sensible heat flux (W/m2)
	REAL RUNOFF		!Total runoff generated (Sum of
				!  surface runoff and moisture
				!  diffusion flux at the bottom of
                                !  the lowest soil layer (kg/m2sec)
	REAL BOMB		!Flag set to 0 in SUBROUTINE TILE
	REAL ESOI 		!The evaporation rate from bare soil (W/m2)
	REAL EINT		!The interception loss (W/m2)
	REAL ESNO		!The evaporation rate from snowpack (W/m2)
	REAL EVEG		!The transpiration rate (W/m2)
	REAL SMELT		!The rate of snowmelt (kg/m2sec)
	REAL HLATN		!The latent heat flux (W/m2)
	REAL HLWUP		!The outgoing longwave radiation flux (W/m2)
	REAL GDRAIN		!The diffusion of moisture across the
				!  bottom of the root zone (middle soil layer)
	REAL RUNSRF		!The overland flow (kg/m2sec)
	REAL FWSOIL		!The infiltration of rainwater into
				!  the top soil layer (kg/m2sec)
	REAL STRDG1		!Diagnostic not currently used
        REAL STRDG2             !Diagnostic not currently used
        REAL STRDG3             !Diagnostic not currently used
        REAL STRDG4             !Diagnostic not currently used
	REAL WSMAX1		!Parameter derived from VGWMAX
        REAL WSMAX2             !Parameter derived from VGWMAX
        REAL WSMAX3             !Parameter derived from VGWMAX
c	  
	 REAL CONDRY		!
c	 
	 REAL GHFLUX		!Ground heat flux (W/m2) ?
c	 
	 REAL POR1		!Porosity of first soil layer?
        REAL POR2               !Porosity of second soil layer?
        REAL POR3               !Porosity of third soil layer?
	REAL WATER1		!SOWET(1)*WSMAX1
        REAL WATER2             !SOWET(2)*WSMAX2
        REAL WATER3             !SOWET(3)*WSMAX3
	REAL COND		!Set to 1/RA


      REAL  DUMBSNOW     !for the Call to SIBALB
	INTEGER DUMBIRUN    !for the Call to SIBALB

	REAL  XV2WIND     !for calculating the wind speed from vector amounts
	REAL  YV2WIND     !for calculating the wind speed from vector amounts
	
	   Real MLTOP,MLBOT

         Real  Wiltpoint1,Wiltpoint2,Wiltpoint3
	   Real TopLayer,MidLayer,BotLayer
        REAL EPSILON            !
        PARAMETER (ALHE = 2.4548e6)
	PARAMETER (ALHS = 2.8368e6)
	PARAMETER (ALHM = ALHS-ALHE)
      	PARAMETER (EPSILON = 18.01/28.97)


C   MOSAIC MONTHLY VEGETATION PARAMETERS
C   HOWEVER, these have already been interpolated
C   So the value passed in is singular for each variable
C				
        REAL GREEN       ! greenness fraction of vegetation
			        !  Holds 12 months of data for interpolation
        REAL LAI         ! leaf area index of vegetation
                                !  Table 5 in Mosaic manual
	REAL VGZ0           !Monthly roughness length
			        !  Table 7 in Mosaic manual
			        !  Used to derive Z0 and U2FAC
	REAL VGROTL         !Monthly variation of root length density
                                !  Table 9 in Mosaic manual
			        !  Used to derive RSOIL1 and RSOIL2
	REAL VGRDC	        !Monthly variation of vegetation specific
			        !  constant used to determin subcanopy
                                !  aerodynamic resistance.
                                !  Table 8 in Mosaic manual
			        !  Used to derive RDC
	REAL VGDD	        !Monthly variation of zero plane
                                !  displacement height
                                !  Table 10 in Mosaic manual
                                !  Used to derive U2FAC



				
C   MOSAIC STATIC PARAMETERS (Passed into mos_middle in MOS%VEGP array)
C  THESE ARE THE STATIC VEGETATION ONES
        REAL VGRF11             !Vegetation parameter used to derive SQSCAT
        REAL VGRF12             !Vegetation parameter used to derive SQSCAT
        REAL VGTR11             !Vegetation parameter used to derive SQSCAT
        REAL VGTR12             !Vegetation parameter used to derive SQSCAT
        REAL VGZ2               !Height of the canopy, Z2 is set to this,
			              !  and it's used to derive U2FAC
        REAL VGROTD             !Rooting depth parameter used to derive 
                                !  RSOIL1 and RSOIL2
	  REAL VGRDRS	        !Resistance to moisture transport per unit
                                !  root length used to derive RSOIL1
	  REAL VGROCA	        !Average cross sectional area of root,
                                !  Used to derive RSOIL2
        REAL VKC	              !Parameter used to derive U2FAC (same for
                                !  all vegetation types (.35))
        REAL VGPH1X             !Soil moisture potential above which
                                !  vegetation is not moisture-stressed (m)
        REAL VGPH2X             !Soil moisture potential below which
                                !  transpiration ceases due to wilting (m)
        REAL VGRPLX             !Average resistance to moisture transport
                                !  within the vegetation itself (s)
        REAL VGCHIL             !Parameter describing departure of leaf
                                !  angles from a spherical distribution
        REAL VGZMEW             !
        REAL VGRST1             !Stomatal resistance parameter a
        REAL VGRST2             !Stomatal resistance parameter b
        REAL VGRST3             !Stomatal resistance parameter c
        REAL VGDFAC             !Parameter controlling vapor pressure
                                !  deficit stress (1/mb)
        REAL VGTLL              !Temperature below which temperature
                                !  stress prevents transpiration (K)
        REAL VGTU               !Temperature above which temperature
                                !  stress prevents transpiration (K)
        REAL VGTCF1             !Coefficient in temperature stress
                                !  equation (K-4)
        REAL VGTCF2             !Coefficient in temperature stress
                                !  equation (K-3)
        REAL VGTCF3             !Coefficient in temperature stress
                                !  equation (K-2)

C MOSAIC STATIC SOIL PARAMETERS

        REAL VGBEEX             !Soil paramter B related to pore size
                                !  distribution index
        REAL VGPSAX             !Soil moisture potential of a saturated
                                !  soil (m)
        REAL VGCSAX             !Hydraulic conductivity of the soil
                                !  at saturation (m/s)
        REAL VGZDEX(3)          !Soil layer thickness (m) of layer i
        REAL VGSLOX             !Cosine of theta

        REAL VGWMAX(3)          !Moisture holding capacity of soil
                                !  layer i (kg/m2)



C   MOSAIC TEMPORALLY INTERPOLATED OR OTHERWISE DERIVED PARAMETERS
C   THESE VALUES ARE GENERALLY BASED ON THE VALUES CONTAINED IN THE
C   VEGPARAM ARRAY OR ON FORCING VARIABLES
					
C        REAL GREEN              !Current greenness fraction of vegetation
C                                !  Used among other things to derive
C                                !  the SQSCAT parameter
C        REAL LAI                !Current leaf area index of vegetation
        REAL Z0		        !Current roughness parameter (based on VGZ0)
        REAL SQSCAT             !Vegetation parameter derived from
                                !  VGRF11,VGRF12,VGTR11,VGTR12,GREEN
	REAL Z2		        !Height of the canopy, set to VGZ2
	REAL DZM	        !Parameter that is set to 600 inside 'PMONTH'
                                !  also used to derive U2FAC
	REAL RSOIL1	        !Parameter that is derived from
			        !  VGRDRS, VGROTD, VGROTL
	REAL RSOIL2	        !Parameter that is derived from
                                !  VGROCA, VGROTD, VGROTL
        REAL SATCAP             !Parameter that is derived from
                                !  LAI
	REAL RDC	        !Parameter derived from VGRDC
	REAL U2FAC	        !Parameter derived from VGZZ, VGZ0,
                                !  DZM
        REAL EM                 !Value computed from M_ESAT
        REAL EA                 !Value computed from QA, PSUR
                                !  and EPSILON
        REAL PSUR               !Surface Pressure (Mb), from SFCPRS/100
        REAL DPT2M              !Dew Point (K), derived from HUMID
        REAL EAI                !Variable used to derive DPT2M
        REAL ESA                !Variable used to derive DPT2M
        REAL SUNANG             !Cosine of the Solar Zenith angle computed
                                !  in the ASTRO subroutine
        REAL RA                 !Earth-Sun distance in units of the
                                !  orbits semi-major axis
        REAL TRAINC             !Convective rainfall rate (kg/m2sec)
                                !  derived from CPCP and TPCP
        REAL TRAINL             !Large-scale rainfall rate (kg/m2sec)
                                !  derived from TPCP and TRAINC
        REAL TSNOW              !The snowfall rate (kg/m2sec)
	REAL SNWFRC		!
	REAL ALBED		!Albedo on 0 to 1 scale derived from
                                !  the AT variable FROM EDAS
	


C MOSAIC FORCING VARIABLES (IN TILE SPACE)
	REAL WND         	!Wind Speed (m/s)
	REAL TMP2M	     	!2 Meter Temperature (K)
	REAL HUMID		!2 Meter Humidity (kg/kg) (is corrected for
                                !numerical instability by setting low 
                                !values to .00001
	REAL HUMIDORIG          !2 Meter Humidity (uncorrected version 
                                !of HUMID)
	REAL CPCP		!Convective Precipitation (kg/m2sec)
	REAL TPCP		!Total Precipitation (kg/m2sec)
        REAL LWDWN              !Downwelling longwave radiation (W/m2)
	REAL SWDN		!Downwelling shortwave radiation (W/m2)
	REAL SFCPRS          	!Surface Pressure (Pa)
	REAL AT			!Albedo FROM EDAS DATA

C	OUTPUT VARIABLES
	REAL OUTRC		!Canopy Conductance
	REAL OUTRA		!Aerodynamic Conductance

        REAL SOILWM             !Soil wetness variable, max avail volumetric mst
        REAL SOILWW             !Soil wetness variable, actual avail vol. mst

!=== End Variable Definition =============================================


c  RESET was originially used to read data from ETA files
c  since the LDAS Subdriver has its own restart files
c  and doesn't have the capability to read ETA soil moisture
c  data then RESET SHOULD ALWAYS BE ZERO
c    note that if ETA reading ability is installed
c   then reset is only equal to one for the very first call
c   and then should be set to zero again!!!!
          RESET=0

C  THE FOLLOWING SECTIONS BREAKS DOWN THE THREE PARAMETER ARRAYS
C   INTO THEIR INDIVIDUAL VARIABLE NAMES FOR LATER USAGE
C
C  STATIC VEGETATION PARAMETERS

           VGRF11=MOS%VEGP(1)
	     VGRF12=MOS%VEGP(2)
	     VGTR11=MOS%VEGP(3)
	     VGTR12=MOS%VEGP(4)
	     VGZ2  =MOS%VEGP(5)
	     VGROTD=MOS%VEGP(6)
	     VGRDRS=MOS%VEGP(7)
	     VGROCA=MOS%VEGP(8)
	     VKC   =MOS%VEGP(9)
	     VGPH1X=MOS%VEGP(10)
	     VGPH2X=MOS%VEGP(11)
	     VGRPLX=MOS%VEGP(12)
	     VGCHIL=MOS%VEGP(13)
	     VGZMEW=MOS%VEGP(14)
	     VGRST1=MOS%VEGP(15)
	     VGRST2=MOS%VEGP(16)
	     VGRST3=MOS%VEGP(17)
	     VGDFAC=MOS%VEGP(18)
	     VGTLL =MOS%VEGP(19)
	     VGTU  =MOS%VEGP(20)
	     VGTCF1=MOS%VEGP(21)
	     VGTCF2=MOS%VEGP(22)
	     VGTCF3=MOS%VEGP(23)
	     SNWMID=MOS%VEGP(24)
	     
C  VARIABLE VEGETATION PARAMETERS
           GREEN =MOS%VEGIP(1)
	     LAI   =MOS%VEGIP(2)
	     VGZ0  =MOS%VEGIP(3)
	     VGROTL=MOS%VEGIP(4)
           VGRDC =MOS%VEGIP(5)
	     VGDD  =MOS%VEGIP(6)

C  STATIC SOIL PARAMETERS
           VGBEEX   =  MOS%SOILP(1)
	     VGPSAX   =  MOS%SOILP(2)
	     VGCSAX   =  MOS%SOILP(3)
	     VGZDEX(1)=  MOS%SOILP(4)
	     VGZDEX(2)=  MOS%SOILP(5)
	     VGZDEX(3)=  MOS%SOILP(6)
	     VGSLOX   =  MOS%SOILP(7)
	     VGWMAX(1)=  MOS%SOILP(8)
	     VGWMAX(2)=  MOS%SOILP(9)
	     VGWMAX(3)=  MOS%SOILP(10)      

C  THE FOLLOWING BREAKS DOWN THE FORCING VARIABLES
           TMP2M = TILE%FORCING(1)
	     HUMID = TILE%FORCING(2)
	     SWDN  = TILE%FORCING(3)
	     LWDWN = TILE%FORCING(4) 
	     XV2WIND=(TILE%FORCING(5))*(TILE%FORCING(5))
	     YV2WIND=(TILE%FORCING(6))*(TILE%FORCING(6))
	     WND   = SQRT( XV2WIND + YV2WIND )
	     SFCPRS= TILE%FORCING(7)
	     TPCP  = TILE%FORCING(8)
	     CPCP  = TILE%FORCING(9)
	     AT    = TILE%FORCING(10)
	     if(AT.lt.001.or.AT.gt.100.1) AT=20.0

!=== Prevent Numerical Instability
      if(WND.le.0.01) WND=0.01

!=== Prevent Numerical Instability with HUMID
C	Store original variable for output
	HUMIDORIG=HUMID
       if(HUMID.le.0.00001)  HUMID=0.00001    ! 1.0E-5  or 0.1E-4

C Compute number of seconds into day
      SEC=(LDAS%MN*60)+(LDAS%HR*60*60)  !total sec into the day

C Compute Zenith angle of sun, SUNANG
      CALL ASTRO(LDAS%YR,LDAS%MO,LDAS%DA,SEC,
     1  TILE%LAT,TILE%LON,1,SUNANG,RA)

C  Call the calculation of the albedo for the given tile, Sunang, veg type
C  This is in the fortran program calc_albedo.f
      DUMBSNOW=1.0
      DUMBIRUN=1                  !this should stay equal to one
C      CALL oldSIBALB( avisdr,anirdr, !albedo Visible/Near Ir both direct
C     & avisdf,anirdf,             !albedo Vis/Near IR both Diffuse
C     & LAI,GREEN,SUNANG,          !LAI and Green, Solar Zenith Angle
C     & DUMBSNOW,TILE%VEGT,        !dumby for snow (not used), Vegetation type
C     & DUMBIRUN,MOS%CT)           !dumby run variable, Temp of Canopy (not used?)

      call umd_sibalb     
     & ( avisdr,anirdr,         !albedo Visible/Near Ir both direct
     & avisdf,anirdf,           !albedo Vis/Near IR both Diffuse
     & LAI,GREEN,SUNANG,        !LAI and Green, Solar Zenith Angle
     & DUMBSNOW,TILE%VEGT,      !dumby for snow(not used), Vegetation type
     & DUMBIRUN,MOS%CT,         !dumby run variable, Temp of Canopy
     & sib )                    !sibalb look up table of coefficients
 
      ALBEDO1=avisdr
      ALBEDO2=anirdr
      ALBEDO3=avisdf
      ALBEDO4=anirdf
      ALBED=(AT/100.0) !Convert Forcing Source Albedo from 0 to 100, 
                       !to 0 to 1 scale.

!=== Process Section
      Call M_PMONTH(LDAS%DOY,TILE%LAT,GREEN,LAI,
     &  Z0,SQSCAT,Z2,DZM,RSOIL1, RSOIL2,
     &  SATCAP,RDC,U2FAC,VGRDRS,VGZ2,VGROTD,
     &  VGROCA,VGRDC,VGROTL,VGZ0,VGDD,
     &  VGTR12,VGTR11,VGRF11,VGRF12)
       
      TSNOW=0.
      TRAINC=AMIN1(CPCP,TPCP)
      TRAINL=TPCP-TRAINC

!=== Determine if Precip is snow
      if(TMP2M.lt.(273.15))THEN
       tsnow=tpcp
       trainl=0.0
       trainc=0.0     
      endif

      psur=sfcprs/100.0

!=== Compute ESI and ESA to get the Dew Point, DPT2M
      EAI=HUMID/0.622*SFCPRS/1000.0
      ESA=0.6108*EXP((17.27*(TMP2M-273.15))/
     &               ((237.3+(TMP2M-273.15))))
      DPT2M=TMP2M/(1-(LOG(EAI/ESA)/17.27))
      em=m_esat(dpt2m)             !m_esat is a function
      ea=MOS%qa*psur/epsilon
      sunang=amax1(sunang,0.01)

!=== account for fact that forcing can lead the land
      snwfrc=(MOS%snow/(MOS%snow+snwmid))
      esattc=(snwfrc*m_qsat(MOS%ct,psur,alhs))+
     &       (1.-snwfrc)*m_qsat(MOS%ct,psur,alhe)*psur/epsilon
     
      if(ea.gt.esattc.and.ea.gt.em) ea=amax1(em,esattc)
      if(ea.lt.esattc.and.ea.lt.em) ea=amin1(em,esattc)
	    
!=== Duplicate what happens in the Call of Turb {Sub-process is turb}	    
      data grav/9.81/,rgas/287./,cp/1010./
      data vkrmn/0.41/,epsi/0.611/,b/5./,c/5./,d/5./

!=== Special parameters for this section
      grav=9.81
      rgas=287.0
      cp=1010.0
      vkrmn=0.41
      epsi=0.611
      b=5.0
      c=5.0
      d=5.0
      zthick=50.0     
      rhoair=psur*100./(rgas*MOS%ct)
      qm=em*epsi/psur
      qc=ea*epsi/psur
      dtv=TMP2M*(1.+epsi*qm)-MOS%ct*(1.+epsi*qc)
      ri=grav*zthick*dtv/(TMP2M*WND*WND)
      dridtc=-grav*zthick*(1.+epsi*qc)/(TMP2M*WND*WND)
      dridqc=-grav*zthick*epsi*MOS%ct/(TMP2M*WND*WND)
      dummy=alog(zthick/z0+1.)
      cn=vkrmn*vkrmn/(dummy*dummy)

      if(ri.ge.0.) then
       t=sqrt(1.+d*ri)
       fu=1./(1.+2.*b*ri/t)
       ft=1./(1.+3.*b*ri*t)
       dftdri=-3.*b*ft*ft*( (d*ri)/(2.*t) + t )
      endif
	  
      if(ri.lt.0.) then
       cs=cn*sqrt((zthick/z0)+1.)
       r=3.*b*cs*sqrt(-ri)    !   so that it can be used in the sqrt function
       s=1./(1.+c*r)	
       t=b*ri*s	 
       ft=1.-3.*t        
       dftdri=-1.5*b*s*s*(2.+c*r)
      endif

      dftdqc=dftdri*dridqc
      dftdtc=dftdri*dridtc
      ctt=cn*ft
      hsturb=rhoair*cp*ctt*WND*(MOS%ct-TMP2M)
      eturb=rhoair*ctt*WND*(qc-qm)
      dedqa=rhoair*cn*WND*( dftdqc*(qc-qm) + ft )
      dedtc=rhoair*cn*WND*dftdtc*(qc-qm)
      dhsdqa=rhoair*cp*cn*WND*dftdqc*(MOS%ct-TMP2M)
      dhsdtc=rhoair*cp*cn*WND*( dftdtc*(MOS%ct-TMP2M) + ft )
      dedea=dedqa*epsi/psur
      dhsdea=dhsdqa*epsi/psur
      ra=1/(ctt*WND)
      if(dhsdtc.lt.0.) dhsdtc=0.
      if(dhsdea.lt.0.) dhsdea=0.
      if(dedtc.lt.0.) dedtc=0.
      if(dedea.lt.0.) dedea=0.

!=== Adjustment to eturb and hsturb with dtcanal based on PSAS and BC
      if(ldas%rpsas.eq.1)then
       eturb=eturb+dedtc*mos%dtcanal
       hsturb=hsturb+dhsdtc*mos%dtcanal
      endif
       

      if(MOS%snow .gt. 0.) then
       snwfrc=MOS%snow/(MOS%snow+snwmid)
       avisdr=avisdr*(1.-snwfrc) + 0.85*snwfrc
       anirdr=anirdr*(1.-snwfrc) + 0.50*snwfrc
       avisdf=avisdf*(1.-snwfrc) + 0.85*snwfrc
       anirdf=anirdf*(1.-snwfrc) + 0.50*snwfrc
      endif
      albave=0.25*(avisdr+anirdr+avisdf+anirdf)
      swnet=(1.-albave)*swdn
      par=0.5*swdn
      pdir=0.5
      totalb=albave

      zlai=LAI
      qm = em * epsilon / psur
      MOS%qa = ea * epsilon / psur
      snwfrc = MOS%snow / ( MOS%snow + snwmid)
      qsattc = (snwfrc*m_qsat(MOS%ct,psur,alhs)+
     &         (1.-snwfrc)*m_qsat(MOS%ct,psur,alhe))
      dqsdtc = qsattc * 5418. / ( MOS%ct * MOS%ct )
      dedqa = dedea * psur / epsilon
      dhsdqa = dhsdea * psur / epsilon
      pardir = par * pdir
      pardif = par * ( 1. - pdir )
      cd = 1. / ( WND * ra )

!=== Compute constants for longwave radiation linearization
      lwattc = 5.67e-8 * MOS%ct*MOS%ct*MOS%ct*MOS%ct
      dlwdtc =  4. * lwattc / MOS%ct
      alwrad = lwattc - dlwdtc * MOS%ct
      blwrad = dlwdtc

      CALL MOSTILE(RESET,FLOAT(LDAS%TS),TRAINL,TRAINC,TSNOW,WND,ETURB,
     I  DEDQA,DEDTC,HSTURB, DHSDQA, DHSDTC,TMP2M,QM,CD,SUNANG,PARDIR,
     I  PARDIF,SWNET,LWDWN,PSUR,ZLAI,GREEN,Z2,SQSCAT,RSOIL1,RSOIL2,
     I  RDC,U2FAC,QSATTC,DQSDTC,ALWRAD,BLWRAD,MOS%DTCANAL,
     U  MOS%CT,MOS%SoT,MOS%QA,MOS%SoWET(1),MOS%SoWET(2),MOS%SoWET(3),
     U  MOS%ICS,MOS%SNOW,
     O  EVAP,SHFLUX,RUNOFF,BOMB,EINT,ESOI,EVEG,ESNO,SMELT,HLATN,HLWUP,
     O  GDRAIN,RUNSRF,FWSOIL,STRDG1,STRDG2,STRDG3,STRDG4,WSMAX1,WSMAX2,
     O  WSMAX3,CONDRY,GHFLUX,POR1,POR2,POR3,water1,water2,water3,VGBEEX,
     O  VGPSAX,VGCSAX,VGZDEX,VGSLOX,VGPH1X,VGPH2X,VGRPLX,VGWMAX,SNWMID,
     O  VGCHIL,VGZMEW,VGRST1,VGRST2,VGRST3,VGRDRS,VGZ2,VGROTD,VGDFAC,
     O  VGTLL,VGTU,VGTCF1,VGTCF2,VGTCF3,TN,LDAS,TILE)

!=== Restrict snow depth to maximum of 100000.0 kg/m^2 ==================
!=== Sometimes such great values are assigned over Greenland & N. Canada
!=== Grib output files do not handle values over this threshold well ====
      if (MOS%SNOW > 100000.0) MOS%SNOW = 100000.0

!=== Collect the output variables into MOS%RETURN
       MOS%RETURN(1)=MOS%CT
       MOS%RETURN(2)=MOS%SoT
       MOS%RETURN(3)=MOS%QA
       MOS%RETURN(4)=MOS%ICS
       MOS%RETURN(5)=MOS%SNOW
       MOS%RETURN(6)=MOS%SoWET(1)
       MOS%RETURN(7)=MOS%SoWET(2)
       MOS%RETURN(8)=MOS%SoWET(3)

       MOS%RETURN(9)=(WATER1+WATER2+WATER3)/
     1               (VGWMAX(1)+VGWMAX(2)+VGWMAX(3)) 

       MOS%RETURN(10)=WATER1
       MOS%RETURN(11)=WATER2
       MOS%RETURN(12)=WATER3
       MOS%RETURN(13)=WATER1+WATER2+WATER3
       MOS%RETURN(14)=MOS%SoWET(1)*POR1
       MOS%RETURN(15)=MOS%SoWET(2)*POR2
       MOS%RETURN(16)=MOS%SoWET(3)*POR3
       MOS%RETURN(60)=(MOS%SoWET(1)*POR1)*(VGZDEX(1)/
     &                   (VGZDEX(1)+VGZDEX(2)+VGZDEX(3)))+
     &                (MOS%SoWET(2)*POR2)*(VGZDEX(2)/
     &                   (VGZDEX(1)+VGZDEX(2)+VGZDEX(3)))+
     &                (MOS%SoWET(3)*POR3)*(VGZDEX(3)/
     &                   (VGZDEX(1)+VGZDEX(2)+VGZDEX(3)))

!=== Must determine wilting point in terms of SoWET(1)
      Wiltpoint1=((VGPH1X/VGPSAX) ** (1.0 / (-VGBEEX)))
      Wiltpoint2=((VGPH2X/VGPSAX) ** (1.0 / (-VGBEEX)))
      Wiltpoint3=0.0

      TopLayer = (VGWMAX(1)*( MOS%SoWET(1)  - Wiltpoint2 )) /
     &  (VGWMAX(1)-Wiltpoint2*VGWMAX(1))
      MidLayer = (VGWMAX(2)*( MOS%SoWET(2)  - Wiltpoint2 )) /
     &  (VGWMAX(2)-Wiltpoint2*VGWMAX(2))
      MLTOP=( MOS%SoWET(2)  - Wiltpoint2 )
      MLBOT= (VGWMAX(2)-Wiltpoint2*VGWMAX(2))

      BotLayer = ( MOS%SoWET(3)*VGWMAX(3) - Wiltpoint3 ) /
     &  (VGWMAX(3)-Wiltpoint3)

       SOILWM=0.0
       SOILWW=0.0
       SOILWM=POR1*VGZDEX(1)
       SOILWM=SOILWM+(POR2*VGZDEX(2))
       SOILWM=SOILWM+(POR3*VGZDEX(3))

       SOILWW=(MOS%SoWET(1)*POR1)*VGZDEX(1)
       SOILWW=SOILWW+((MOS%SoWET(2)*POR2)*VGZDEX(2))
       SOILWW=SOILWW+((MOS%SoWET(3)*POR3)*VGZDEX(3))

       MOS%RETURN(17)=SOILWW/SOILWM


       MOS%RETURN(18)=(-1.0)*SWNET
       MOS%RETURN(19)=HLWUP-LWDWN
       MOS%RETURN(20)=SHFLUX
       MOS%RETURN(21)=HLATN
       MOS%RETURN(22)=(-1.0)*GHFLUX
       MOS%RETURN(23)=EVAP*LDAS%TS
       MOS%RETURN(24)=EINT*(-1.0)
       MOS%RETURN(25)=ESOI*(-1.0)
       MOS%RETURN(26)=EVEG*(-1.0)
       IF (TSNOW.EQ.0) THEN
         MOS%RETURN(27)=TPCP*LDAS%TS
       ELSE
         MOS%RETURN(27)=0.0
       ENDIF
       MOS%RETURN(28)=TSNOW*LDAS%TS
       MOS%RETURN(29)=SMELT*LDAS%TS
       MOS%RETURN(30)=fWSOIL
       MOS%RETURN(31)=GDRAIN
       MOS%RETURN(32)=RUNOFF*LDAS%TS
       MOS%RETURN(33)=RUNSRF*LDAS%TS
       MOS%RETURN(34)=ZLAI
       MOS%RETURN(35)=GREEN
       MOS%RETURN(36)=( (TopLayer*VGZDEX(1)) + 
     &	(MidLayer*VGZDEX(2)) + (BotLayer*VGZDEX(3)) ) /
     &  (VGZDEX(1)+VGZDEX(2)+VGZDEX(3))

       MOS%RETURN(37)=SWDN
       MOS%RETURN(38)=LWDWN
       MOS%RETURN(39)=ALBAVE*100.0
       MOS%RETURN(40)=ESNO*(-1.0)

       IF (MOS%RETURN(5).GT.0.0) THEN
         MOS%RETURN(41)=1.0
       ELSE
         MOS%RETURN(41)=0.0
       ENDIF	
c        if(MOS%snow .eq. 0.) then 
c	  =ldas%udef
c	endif
       MOS%RETURN(42)=LDAS%UDEF
       MOS%RETURN(43)=TILE%AC
       MOS%RETURN(44)=TILE%CC
       MOS%RETURN(45)=SMELT*382000.
       MOS%RETURN(46)=MOS%RETURN(32)-MOS%RETURN(33)	
       MOS%RETURN(47)=LDAS%UDEF
       MOS%RETURN(48)=MOS%RETURN(40)+MOS%RETURN(24)
C	Convert mm liq snow equiv to Meters of snow depth using
C	subjective snow to water ratio of 1mm liq to 10mm snow
       MOS%RETURN(49)=MOS%RETURN(5)*10.0/1000.0
c	IF (MOS%RETURN(5).GT.0.0) THEN
c	  MOS%RETURN(50)=1.0
c	ELSE
c	  MOS%RETURN(50)=0.0
c	ENDIF
	MOS%RETURN(50)=snwfrc
C	Calculate soil moisture (kg/m2) in top 1meter layer of soil
C	If soil is less than 1 meter deep, scale the soil moisture
C	up so that it has an effective depth of 1 meter of soil

	DEPTH=0.0
	FLAG=0
	DEPTH=DEPTH+VGZDEX(1)
	IF (DEPTH.GE.1) THEN
	  FACTOR=(1.0-((DEPTH-1.0)/DEPTH))
	  IF (DEPTH.EQ.1.0) THEN
	    MOS%RETURN(51)=MOS%RETURN(10)
            FLAG=1
	  ELSE
 	    MOS%RETURN(51)=MOS%RETURN(10)*FACTOR
            FLAG=1
	  ENDIF
	ENDIF
	IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.1.0) THEN
            IF (DEPTH.EQ.1.0) THEN
              MOS%RETURN(51)=MOS%RETURN(10)+MOS%RETURN(11)
              FLAG=1
            ELSE
              TEMP1=DEPTH-1.0
              FACTOR=(1.0-(TEMP1/VGZDEX(2)))
              TEMP2=MOS%RETURN(11)*FACTOR
              MOS%RETURN(51)=MOS%RETURN(10)+TEMP2
	      FLAG=1
            ENDIF
       	  ENDIF
	ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(3)
          IF (DEPTH.GE.1.0) THEN
            IF (DEPTH.EQ.1.0) THEN
              MOS%RETURN(51)=MOS%RETURN(10)+MOS%RETURN(11)+
     &          MOS%RETURN(12)
              FLAG=1
            ELSE
              TEMP1=DEPTH-1.0
              FACTOR=(1.0-(TEMP1/VGZDEX(3)))
              TEMP2=MOS%RETURN(12)*FACTOR
              MOS%RETURN(51)=MOS%RETURN(10)+MOS%RETURN(11)
     &          + TEMP2
              FLAG=1
            ENDIF
          ELSE 
	    MOS%RETURN(51)=(MOS%RETURN(10)+MOS%RETURN(11)+
     &          MOS%RETURN(12))/DEPTH
c	print *,'old value=',MOS%RETURN(10)+MOS%RETURN(11)+
c     &          MOS%RETURN(12)
c	print *,'new value,depth=',MOS%RETURN(51),depth
	    FLAG=1
          ENDIF
        ENDIF
	MOS%RETURN(52)=MOS%RETURN(10)+MOS%RETURN(11)
       SOILWM=0.0
       SOILWW=0.0
       SOILWM=(POR1-(POR1*WILTPOINT2))*VGZDEX(1)
       SOILWM=SOILWM+((POR2-(POR2*WILTPOINT2))*VGZDEX(2))

       SOILWW=((MOS%SoWET(1)*POR1)-(POR1*WILTPOINT2))*VGZDEX(1)
       SOILWW=SOILWW+(((MOS%SoWET(2)*POR2)-
     &        (POR2*WILTPOINT2))*VGZDEX(2))

       MOS%RETURN(53)=SOILWW/SOILWM
        MOS%RETURN(54)=TMP2M
	MOS%RETURN(55)=HUMIDORIG
	MOS%RETURN(56)=TILE%FORCING(5)
	MOS%RETURN(57)=TILE%FORCING(6)
	MOS%RETURN(58)=SFCPRS/100.0
	MOS%RETURN(59)=CPCP*LDAS%TS
C       Calculate volumetric soil moisture at 0.05 m depth

        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.0.05) THEN
            MOS%RETURN(61)=MOS%RETURN(10)/(VGZDEX(1)*1000)
            FLAG=1
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.0.05) THEN
              MOS%RETURN(61)=MOS%RETURN(11)/(VGZDEX(2)*1000)
              FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
              MOS%RETURN(61)=MOS%RETURN(12)/(VGZDEX(3)*1000)
            FLAG=1
        ENDIF

C       Calculate volumetric soil moisture at 0.25 m depth

        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.0.25) THEN
            MOS%RETURN(62)=MOS%RETURN(10)/(VGZDEX(1)*1000)
            FLAG=1
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.0.25) THEN
              MOS%RETURN(62)=MOS%RETURN(11)/(VGZDEX(2)*1000)
              FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
              MOS%RETURN(62)=MOS%RETURN(12)/(VGZDEX(3)*1000)
            FLAG=1
        ENDIF

C       Calculate volumetric soil moisture at 0.60 m depth

        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.0.6) THEN
            MOS%RETURN(63)=MOS%RETURN(10)/(VGZDEX(1)*1000)
            FLAG=1
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.0.6) THEN
              MOS%RETURN(63)=MOS%RETURN(11)/(VGZDEX(2)*1000)
              FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
              MOS%RETURN(63)=MOS%RETURN(12)/(VGZDEX(3)*1000)
            FLAG=1
        ENDIF

C       Calculate volumetric soil moisture at 0.75 m depth
        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.0.75) THEN
            MOS%RETURN(64)=MOS%RETURN(10)/(VGZDEX(1)*1000)
            FLAG=1
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.0.75) THEN
              MOS%RETURN(64)=MOS%RETURN(11)/(VGZDEX(2)*1000)
              FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
              MOS%RETURN(64)=MOS%RETURN(12)/(VGZDEX(3)*1000)
            FLAG=1
        ENDIF

      ea=MOS%QA*psur/epsilon
      MOS%QA=ea*epsilon/psur
      cond=1./ra
    	     
      RETURN
      END


	     
	     
	     


