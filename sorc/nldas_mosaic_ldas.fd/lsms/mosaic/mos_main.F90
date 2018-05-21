!-------------------------------------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
! Released October 2005
!
! See SOFTWARE DISTRIBUTION POLICY for software distribution policies
!
! The LIS source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See COPYRIGHT.TXT for copyright details.
!
!-------------------------------------------------------------------------
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
!=========================================================================

      SUBROUTINE MOS_MAIN()

      USE lisdrv_module, only : lis, tile,grid,gindex
      USE mos_varder      ! MOS tile variables
      use sibalb_module   ! SiB-based coefficients for albedo calculation.
      use tile_spmdMod

      IMPLICIT NONE

!      type (sibalbdec) sib

!=== Local Variables =====================================================
        integer kkk
        integer index
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
	  
	 REAL CONDRY		!
	 
	 REAL GHFLUX		!Ground heat flux (W/m2) ?
	 
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


!   MOSAIC MONTHLY VEGETATION PARAMETERS
!   HOWEVER, these have already been interpolated
!   So the value passed in is singular for each variable
!				
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



				
!   MOSAIC STATIC PARAMETERS (Passed into mos_middle in MOS%VEGP array)
!  THESE ARE THE STATIC VEGETATION ONES
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

! MOSAIC STATIC SOIL PARAMETERS

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



!   MOSAIC TEMPORALLY INTERPOLATED OR OTHERWISE DERIVED PARAMETERS
!  THESE VALUES ARE GENERALLY BASED ON THE VALUES CONTAINED IN THE
!   VEGPARAM ARRAY OR ON FORCING VARIABLES
					
!        REAL GREEN              !Current greenness fraction of vegetation
!                                !  Used among other things to derive
!                                !  the SQSCAT parameter
!        REAL LAI                !Current leaf area index of vegetation
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
	


! MOSAIC FORCING VARIABLES (IN TILE SPACE)
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

!	OUTPUT VARIABLES
	REAL OUTRC		!Canopy Conductance
	REAL OUTRA		!Aerodynamic Conductance

        REAL SOILWM             !Soil wetness variable, max avail volumetric mst
        REAL SOILWW             !Soil wetness variable, actual avail vol. mst
        integer M
!********************************        
        real tempac,tempcc
!*********************************        

      !call mos_debug()

!=== End Variable Definition =============================================
        do M = 1, di_array(iam)

!  RESET was originially used to read data from ETA files
!  since the LDAS Subdriver has its own restart files
!  and doesn't have the capability to read ETA soil moisture
!  data then RESET SHOULD ALWAYS BE ZERO
!    note that if ETA reading ability is installed
!   then reset is only equal to one for the very first call
!   and then should be set to zero again!!!!
          RESET=0

!  THE FOLLOWING SECTIONS BREAKS DOWN THE THREE PARAMETER ARRAYS
!   INTO THEIR INDIVIDUAL VARIABLE NAMES FOR LATER USAGE
!
!  STATIC VEGETATION PARAMETERS

           VGRF11=MOS(M)%VEGP(1)
	     VGRF12=MOS(M)%VEGP(2)
	     VGTR11=MOS(M)%VEGP(3)
	     VGTR12=MOS(M)%VEGP(4)
	     VGZ2  =MOS(M)%VEGP(5)
	     VGROTD=MOS(M)%VEGP(6)
	     VGRDRS=MOS(M)%VEGP(7)
	     VGROCA=MOS(M)%VEGP(8)
	     VKC   =MOS(M)%VEGP(9)
	     VGPH1X=MOS(M)%VEGP(10)
	     VGPH2X=MOS(M)%VEGP(11)
	     VGRPLX=MOS(M)%VEGP(12)
	     VGCHIL=MOS(M)%VEGP(13)
	     VGZMEW=MOS(M)%VEGP(14)
	     VGRST1=MOS(M)%VEGP(15)
	     VGRST2=MOS(M)%VEGP(16)
	     VGRST3=MOS(M)%VEGP(17)
	     VGDFAC=MOS(M)%VEGP(18)
	     VGTLL =MOS(M)%VEGP(19)
	     VGTU  =MOS(M)%VEGP(20)
	     VGTCF1=MOS(M)%VEGP(21)
	     VGTCF2=MOS(M)%VEGP(22)
	     VGTCF3=MOS(M)%VEGP(23)
	     SNWMID=MOS(M)%VEGP(24)
	     
!  VARIABLE VEGETATION PARAMETERS
           GREEN =MOS(M)%VEGIP(1)
	     LAI   =MOS(M)%VEGIP(2)
	     VGZ0  =MOS(M)%VEGIP(3)
	     VGROTL=MOS(M)%VEGIP(4)
           VGRDC =MOS(M)%VEGIP(5)
	     VGDD  =MOS(M)%VEGIP(6)

!  STATIC SOIL PARAMETERS
           VGBEEX   =  MOS(M)%SOILP(1)
	     VGPSAX   =  MOS(M)%SOILP(2)
	     VGCSAX   =  MOS(M)%SOILP(3)
	     VGZDEX(1)=  MOS(M)%SOILP(4)
	     VGZDEX(2)=  MOS(M)%SOILP(5)
	     VGZDEX(3)=  MOS(M)%SOILP(6)
	     VGSLOX   =  MOS(M)%SOILP(7)
	     VGWMAX(1)=  MOS(M)%SOILP(8)
	     VGWMAX(2)=  MOS(M)%SOILP(9)
	     VGWMAX(3)=  MOS(M)%SOILP(10)      

!  THE FOLLOWING BREAKS DOWN THE FORCING VARIABLES
             
             TMP2M = MOS(M)%FORCING(1)



	     HUMID = MOS(M)%FORCING(2)
	     SWDN  = MOS(M)%FORCING(3)
	     LWDWN = MOS(M)%FORCING(4)
	     XV2WIND=(MOS(M)%FORCING(5))*(MOS(M)%FORCING(5))
	     YV2WIND=(MOS(M)%FORCING(6))*(MOS(M)%FORCING(6))
	     WND   = SQRT( XV2WIND + YV2WIND )
	     SFCPRS= MOS(M)%FORCING(7)
	     TPCP  = MOS(M)%FORCING(8)
	     CPCP  = MOS(M)%FORCING(9)
	     AT    = MOS(M)%FORCING(10)

!             print*, 'here.. ',m,tmp2m,humid,at
             if(AT.lt.001.or.AT.gt.100.1) AT=20.0
!=== Prevent Numerical Instability
      if(WND.le.0.01) WND=0.01

!=== Prevent Numerical Instability with HUMID
!	Store original variable for output
	HUMIDORIG=HUMID
       if(HUMID.le.0.00001)  HUMID=0.00001    ! 1.0E-5  or 0.1E-4

! Compute number of seconds into day
      SEC=(LIS%t%MN*60)+(LIS%t%HR*60*60)  !total sec into the day

! Compute Zenith angle of sun, SUNANG
!      CALL ASTRO(LIS%t%YR,LIS%t%MO,LIS%t%DA,SEC, &
!       TILE(M)%LAT,TILE(M)%LON,1,SUNANG,RA)
!       index = gindex(tile(m)%col,tile(m)%row)
       index = tile(m)%index
      CALL ASTRO(LIS%t%YR,LIS%t%MO,LIS%t%DA,SEC, &
       grid(index)%lat,grid(index)%lon,1,SUNANG,RA)


!      if (m .eq. 100000) then
!          print*, LIS%t%YR,LIS%t%MO,LIS%t%DA,SEC
!	  print*, grid(index)%lat,grid(index)%lon
!	  print*, SUNANG,RA
!      endif				  

!  Call the calculation of the albedo for the given tile, Sunang, veg type
!  This is in the fortran program calc_albedo.f
      DUMBSNOW=1.0
      DUMBIRUN=1                  !this should stay equal to one
!      CALL oldSIBALB( avisdr,anirdr, !albedo Visible/Near Ir both direct
!     & avisdf,anirdf,             !albedo Vis/Near IR both Diffuse
!     & LAI,GREEN,SUNANG,          !LAI and Green, Solar Zenith Angle
!     & DUMBSNOW,TILE%VEGT,        !dumby for snow (not used), Vegetation type
!     & DUMBIRUN,MOS(M)%CT)           !dumby run variable, Temp of Canopy (not used?)

!      do i=1,13
!        print*, i, sib%ALVDR(1,1,i)
!      enddo

!      stop
      
!      write(*,63) m, tile(m)%row,tile(m)%col,TILE(M)%VEGT,avisdr
! 63   format(i6,1x,i4,1x,i4,1x,i2,1x,f8.5)
         call umd_sibalb   &   
              ( avisdr,anirdr,  &       !albedo Visible/Near Ir both direct
              avisdf,anirdf,  &         !albedo Vis/Near IR both Diffuse
              LAI,GREEN,SUNANG, &        !LAI and Green, Solar Zenith Angle
              DUMBSNOW,TILE(M)%VEGT, &      !dumby for snow(not used), Vegetation type
              DUMBIRUN,MOS(M)%CT)           !dumby run variable, Temp of Canopy



!                      !sibalb look up table of coefficients
!       print*, "M = ",M

!       if (m .eq. 100000) then
!        print*, TILE(M)%VEGT,LAI,GREEN,SUNANG
!	print*, avisdr,anirdr,avisdf,anirdf
!       endif
       
					   
!      write(*,63) m, tile(m)%row,tile(m)%col,TILE(M)%VEGT,avisdr
! 64   format(i6,1x,i4,1x,i4,1x,i2,1x,f8.5)
      
!      print*, "hereM3 ",M,MOS(M)%CT
      
      ALBEDO1=avisdr
      ALBEDO2=anirdr
      ALBEDO3=avisdf
      ALBEDO4=anirdf
      ALBED=(AT/100.0) !Convert Forcing Source Albedo from 0 to 100, 
                       !to 0 to 1 scale.

!=== Process Section
      Call M_PMONTH(LIS%t%DOY,MOS(M)%LAT,GREEN,LAI, &
       Z0,SQSCAT,Z2,DZM,RSOIL1, RSOIL2, &
       SATCAP,RDC,U2FAC,VGRDRS,VGZ2,VGROTD, &
       VGROCA,VGRDC,VGROTL,VGZ0,VGDD, &
       VGTR12,VGTR11,VGRF11,VGRF12) 

!       if (m .eq. 100000) then
!       print*, MOS(M)%LAT,MOS(M)%QA
!       endif
      TSNOW=0.

      TRAINC=AMIN1(CPCP,TPCP)
      TRAINL=TPCP-TRAINC

!=== Determine if Precip is snow'

      if(TMP2M.lt.(273.15))THEN
       tsnow=tpcp
       trainl=0.0
       trainc=0.0     
      endif


      psur=sfcprs/100.0
!=== Compute ESI and ESA to get the Dew Point, DPT2M
      EAI=HUMID/0.622*SFCPRS/1000.0
      ESA=0.6108*EXP((17.27*(TMP2M-273.15))/ &
                    ((237.3+(TMP2M-273.15))))
      DPT2M=TMP2M/(1-(LOG(EAI/ESA)/17.27))
      em=m_esat(dpt2m)             !m_esat is a function
      ea=MOS(M)%qa*psur/epsilon
      sunang=amax1(sunang,0.01)


!      if (isnan(MOS(M)%snow)) then
!        MOS(M)%snow = 0.0
!      end if

!=== account for fact that forcing can lead the land
      snwfrc=(MOS(M)%snow/(MOS(M)%snow+snwmid))
      esattc=(snwfrc*m_qsat(MOS(M)%ct,psur,alhs))+ &
            (1.-snwfrc)*m_qsat(MOS(M)%ct,psur,alhe)*psur/epsilon
     
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
      rhoair=psur*100./(rgas*MOS(M)%ct)
      qm=em*epsi/psur
      qc=ea*epsi/psur
      dtv=TMP2M*(1.+epsi*qm)-MOS(M)%ct*(1.+epsi*qc)
      ri=grav*zthick*dtv/(TMP2M*WND*WND)
      dridtc=-grav*zthick*(1.+epsi*qc)/(TMP2M*WND*WND)
      dridqc=-grav*zthick*epsi*MOS(M)%ct/(TMP2M*WND*WND)
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
      hsturb=rhoair*cp*ctt*WND*(MOS(M)%ct-TMP2M)
      eturb=rhoair*ctt*WND*(qc-qm)
      dedqa=rhoair*cn*WND*( dftdqc*(qc-qm) + ft )
      dedtc=rhoair*cn*WND*dftdtc*(qc-qm)
      dhsdqa=rhoair*cp*cn*WND*dftdqc*(MOS(M)%ct-TMP2M)
      dhsdtc=rhoair*cp*cn*WND*( dftdtc*(MOS(M)%ct-TMP2M) + ft )
!      if(m.eq.95094) print *, &
!      'dhsdtc,rhoair,cp,cn,wnd,dftdtc,MOS(M)%ct,TMP2M,ft', &
!      dhsdtc,rhoair,cp,cn,wnd,dftdtc,MOS(M)%ct,TMP2M,ft 
      dedea=dedqa*epsi/psur
      dhsdea=dhsdqa*epsi/psur
      ra=1/(ctt*WND)
      if(dhsdtc.lt.0.) dhsdtc=0.
      if(dhsdea.lt.0.) dhsdea=0.
      if(dedtc.lt.0.) dedtc=0.
      if(dedea.lt.0.) dedea=0.

!=== Adjustment to eturb and hsturb with dtcanal based on PSAS and BC
      if(lis%a%rpsas.eq.1)then
       eturb=eturb+dedtc*MOS(M)%dtcanal
       hsturb=hsturb+dhsdtc*MOS(M)%dtcanal
      endif
       


      if(MOS(M)%snow .gt. 0.) then
       snwfrc=MOS(M)%snow/(MOS(M)%snow+snwmid)
       avisdr=avisdr*(1.-snwfrc) + 0.85*snwfrc
       anirdr=anirdr*(1.-snwfrc) + 0.50*snwfrc
       avisdf=avisdf*(1.-snwfrc) + 0.85*snwfrc
       anirdf=anirdf*(1.-snwfrc) + 0.50*snwfrc
      endif

      albave=0.25*(avisdr+anirdr+avisdf+anirdf)
!      print*, mos(m)%snow, avisdr, anirdr, avisdf, anirdf

!      print*, "ALB check1 ", swdn, MOS(M)%snow, snwfrc
!      print*, "ALB check2 ", avisdr, anirdr, avisdf, anirdf
      swnet=(1.-albave)*swdn
      par=0.5*swdn
      pdir=0.5
      totalb=albave

      zlai=LAI
      qm = em * epsilon / psur
      MOS(M)%qa = ea * epsilon / psur
      snwfrc = MOS(M)%snow / ( MOS(M)%snow + snwmid)
      qsattc = (snwfrc*m_qsat(MOS(M)%ct,psur,alhs)+ &
              (1.-snwfrc)*m_qsat(MOS(M)%ct,psur,alhe))
      dqsdtc = qsattc * 5418. / ( MOS(M)%ct * MOS(M)%ct )
      dedqa = dedea * psur / epsilon
      dhsdqa = dhsdea * psur / epsilon
      pardir = par * pdir
      pardif = par * ( 1. - pdir )
      cd = 1. / ( WND * ra )

!=== Compute constants for longwave radiation linearization
      lwattc = 5.67e-8 * MOS(M)%ct*MOS(M)%ct*MOS(M)%ct*MOS(M)%ct
      dlwdtc =  4. * lwattc / MOS(M)%ct
      alwrad = lwattc - dlwdtc * MOS(M)%ct
      blwrad = dlwdtc

!      if (m .eq. 1) then
!      print*, "hereM1 ", M,MOS(M)%CT
!      endif


      CALL MOSTILE(RESET,FLOAT(LIS%t%TS),TRAINL,TRAINC,TSNOW,WND,ETURB,&
       DEDQA,DEDTC,HSTURB, DHSDQA, DHSDTC,TMP2M,QM,CD,SUNANG,PARDIR, &
       PARDIF,SWNET,LWDWN,PSUR,ZLAI,GREEN,Z2,SQSCAT,RSOIL1,RSOIL2, &
       RDC,U2FAC,QSATTC,DQSDTC,ALWRAD,BLWRAD,MOS(M)%DTCANAL, &
       MOS(M)%CT,MOS(M)%SoT,MOS(M)%QA,MOS(M)%SoWET(1),MOS(M)%SoWET(2), &
       MOS(M)%SoWET(3),MOS(M)%ICS,MOS(M)%SNOW, &
       EVAP,SHFLUX,RUNOFF,BOMB,EINT,ESOI,EVEG,ESNO,SMELT,HLATN,HLWUP, &
       GDRAIN,RUNSRF,FWSOIL,STRDG1,STRDG2,STRDG3,STRDG4,WSMAX1,WSMAX2, &
       WSMAX3,CONDRY,GHFLUX,POR1,POR2,POR3,water1,water2,water3,VGBEEX, &
       VGPSAX,VGCSAX,VGZDEX,VGSLOX,VGPH1X,VGPH2X,VGRPLX,VGWMAX,SNWMID, &
       VGCHIL,VGZMEW,VGRST1,VGRST2,VGRST3,VGRDRS,VGZ2,VGROTD,VGDFAC, &
       VGTLL,VGTU,VGTCF1,VGTCF2,VGTCF3,tempcc,tempac,m)









!=== Restrict snow depth to maximum of 100000.0 kg/m^2 ==================
!=== Sometimes such great values are assigned over Greenland & N. Canada
!=== Grib output files do not handle values over this threshold well ====
      if (MOS(M)%SNOW > 100000.0) MOS(M)%SNOW = 100000.0
      mos(m)%swnet = mos(m)%swnet + swnet
      mos(m)%lwnet = mos(m)%lwnet +hlwup-lwdwn
      mos(m)%qle = mos(m)%qle+ hlatn
      mos(m)%qh=mos(m)%qh+ shflux
      mos(m)%qg=mos(m)%qg+ghflux
      mos(m)%swrad=mos(m)%swrad+mos(m)%forcing(3)
      mos(m)%lwrad=mos(m)%lwrad+mos(m)%forcing(4)



      mos(m)%snowf= mos(m)%snowf + tsnow
      if (tsnow.eq.0) then
         mos(m)%rainf=mos(m)%rainf + tpcp
      else
         mos(m)%rainf=mos(m)%rainf+ 0.0
      endif

      mos(m)%laiout = zlai
      mos(m)%green = green
      mos(m)%acond = tempac
      mos(m)%ccond = tempcc

      mos(m)%sbsno =  mos(m)%sbsno+ (esno*(1.0))
      !mos(m)%sbsno =  mos(m)%sbsno+ (esno*(-1.0))
      mos(m)%snohf =  mos(m)%snohf+(SMELT*382000.)
      mos(m)%evap = mos(m)%evap + evap*LIS%t%TS
      mos(m)%qs=mos(m)%qs + runsrf*LIS%t%TS
      mos(m)%qsb=mos(m)%qsb + (runoff-runsrf)*LIS%t%TS

      mos(m)%qsm=mos(m)%qsm + smelt*LIS%t%TS

      mos(m)%avgsurft = mos(m)%ct
      mos(m)%albedo=albave*100.0
!      mos(m)%swe = mos(m)%swe + mos(m)%snow
      mos(m)%swe = mos(m)%snow


      mos(m)%water1 = water1
      mos(m)%water2 = water2
      mos(m)%water3 = water3

!      mos(m)%soilmoist1 = mos(m)%soilmoist1+  water1
!      mos(m)%soilmoist2 = mos(m)%soilmoist2 + water2
!      mos(m)%soilmoist3 = mos(m)%soilmoist3 + water3

      mos(m)%soilmoist1 = water1
      mos(m)%soilmoist2 = water2
      mos(m)%soilmoist3 = water3

      mos(m)%soilmtot = mos(m)%soilmoist1+mos(m)%soilmoist2+mos(m)%soilmoist3
      mos(m)%soilmr = mos(m)%soilmoist1+mos(m)%soilmoist2
!      Calculate soil moisture (kg/m2) in top 1meter layer of soil
!      If soil is less than 1 meter deep, scale the soil moisture
!      up so that it has an effective depth of 1 meter of soil
!
        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.1) THEN
          FACTOR=(1.0-((DEPTH-1.0)/DEPTH))
          IF (DEPTH.EQ.1.0) THEN
            mos(m)%soilm1=mos(m)%soilmoist1
            FLAG=1
          ELSE
            mos(m)%soilm1=mos(m)%soilmoist1*FACTOR
            FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.1.0) THEN
            IF (DEPTH.EQ.1.0) THEN
              mos(m)%soilm1=mos(m)%soilmoist1+mos(m)%soilmoist2
              FLAG=1
            ELSE
              TEMP1=DEPTH-1.0
              FACTOR=(1.0-(TEMP1/VGZDEX(2)))
              TEMP2=mos(m)%soilmoist2*FACTOR
              mos(m)%soilm1=mos(m)%soilmoist1+TEMP2
              FLAG=1
            ENDIF
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(3)
          IF (DEPTH.GE.1.0) THEN
            IF (DEPTH.EQ.1.0) THEN
              mos(m)%soilm1=mos(m)%soilmoist1+mos(m)%soilmoist2+ &
               mos(m)%soilmoist3
              FLAG=1
            ELSE
              TEMP1=DEPTH-1.0
              FACTOR=(1.0-(TEMP1/VGZDEX(3)))
              TEMP2=mos(m)%soilmoist3*FACTOR
              mos(m)%soilm1=mos(m)%soilmoist1+mos(m)%soilmoist2 &
               + TEMP2
              FLAG=1
            ENDIF
          ELSE
            mos(m)%soilm1=(mos(m)%soilmoist1+mos(m)%soilmoist2+ &
               mos(m)%soilmoist3)/DEPTH
!c       print *,'old value=',mos(m)%soilmoist1+mos(m)%soilmoist2+
!c     &          mos(m)%soilmoist3
!c       print *,'new value,depth=',mos(m)%soilmr,depth
            FLAG=1
          ENDIF
        ENDIF





      Wiltpoint1=((VGPH1X/VGPSAX) ** (1.0 / (-VGBEEX)))
      Wiltpoint2=((VGPH2X/VGPSAX) ** (1.0 / (-VGBEEX)))
      Wiltpoint3=0.0

      soilwm=0.0
      soilww=0.0
      soilwm=por1*vgzdex(1)
      soilwm=soilwm+(por2*vgzdex(2))
      soilwm=soilwm+(por3*vgzdex(3))
      soilww=(mos(m)%sowet(1)*por1)*vgzdex(1)
      soilww=soilww+((mos(m)%sowet(2)*por2)*vgzdex(2))
      soilww=soilww+((mos(m)%sowet(3)*por3)*vgzdex(3))
!      mos(m)%soilwet= mos(m)%soilwet + soilww/soilwm
      mos(m)%soilwet= soilww/soilwm


       SOILWM=0.0
       SOILWW=0.0
       SOILWM=(POR1-(POR1*WILTPOINT2))*VGZDEX(1)
       SOILWM=SOILWM+((POR2-(POR2*WILTPOINT2))*VGZDEX(2))
       SOILWW=((MOS(m)%SoWET(1)*POR1)-(POR1*WILTPOINT2))*VGZDEX(1)
       SOILWW=SOILWW+(((MOS(m)%SoWET(2)*POR2)- &
             (POR2*WILTPOINT2))*VGZDEX(2))
       mos(m)%mstavr=SOILWW/SOILWM



      mos(m)%ecanop = mos(m)%ecanop + eint*(-1.0) 
      mos(m)%tveg = mos(m)%tveg + eveg*(-1.0)
      mos(m)%esoil = mos(m)%esoil + esoi*(-1.0)

!      mos(m)%rootmoist = mos(m)%rootmoist + water1 + water2 
      mos(m)%rootmoist = water1 + water2

      mos(m)%canopint= mos(m)%canopint + mos(m)%ics
      mos(m)%acond = tempac
      if(lis%t%tscount == 0 .or. lis%t%tscount ==1 &
           .or.lis%f%rstflag.eq.1) then
         mos(m)%soilm_prev = mos(m)%soilm_prev + mos(m)%soilmoist1+ & 
              mos(m)%soilmoist2 + mos(m)%soilmoist3
         mos(m)%swe_prev = mos(m)%swe
      endif
       MOS(m)%SNOD=MOS(m)%SNOW*10.0/1000.0
       MOS(m)%SNWFRCOUT=snwfrc

!C       Calculate volumetric soil moisture at 0.05 m depth
        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.0.05) THEN
            MOS(m)%soilv1=MOS(m)%soilmoist1/(VGZDEX(1)*1000)
            FLAG=1
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.0.05) THEN
              MOS(m)%soilv1=MOS(m)%soilmoist2/(VGZDEX(2)*1000)
              FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
              MOS(m)%soilv1=MOS(m)%soilmoist3/(VGZDEX(3)*1000)
            FLAG=1
        ENDIF

!C       Calculate volumetric soil moisture at 0.25 m depth

        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.0.25) THEN
            mos(m)%soilv2=MOS(m)%soilmoist1/(VGZDEX(1)*1000)
            FLAG=1
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.0.25) THEN
              mos(m)%soilv2=MOS(m)%soilmoist2/(VGZDEX(2)*1000)
              FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
              mos(m)%soilv2=MOS(m)%soilmoist3/(VGZDEX(3)*1000)
            FLAG=1
        ENDIF

!C       Calculate volumetric soil moisture at 0.60 m depth

        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.0.6) THEN
            mos(m)%soilv3=MOS(m)%soilmoist1/(VGZDEX(1)*1000)
            FLAG=1
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.0.6) THEN
              mos(m)%soilv3=MOS(m)%soilmoist2/(VGZDEX(2)*1000)
              FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
              mos(m)%soilv3=MOS(m)%soilmoist3/(VGZDEX(3)*1000)
            FLAG=1
        ENDIF

!C       Calculate volumetric soil moisture at 0.75 m depth
        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.0.75) THEN
            mos(m)%soilv4=MOS(m)%soilmoist3/(VGZDEX(1)*1000)
            FLAG=1
        ENDIF
        IF (FLAG.EQ.0) THEN
          DEPTH=DEPTH+VGZDEX(2)
          IF (DEPTH.GE.0.75) THEN
              mos(m)%soilv4=MOS(m)%soilmoist3/(VGZDEX(2)*1000)
              FLAG=1
          ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
              mos(m)%soilv4=MOS(m)%soilmoist3/(VGZDEX(3)*1000)
            FLAG=1
        ENDIF

       MOS(m)%soilv5=(MOS(m)%SoWET(1)*POR1)*(VGZDEX(1)/ &
                        (VGZDEX(1)+VGZDEX(2)+VGZDEX(3)))+ &
                     (MOS(m)%SoWET(2)*POR2)*(VGZDEX(2)/ &
                        (VGZDEX(1)+VGZDEX(2)+VGZDEX(3)))+ &
                     (MOS(m)%SoWET(3)*POR3)*(VGZDEX(3)/ &
                        (VGZDEX(1)+VGZDEX(2)+VGZDEX(3)))

       MOS(m)%soilv6=MOS(m)%SoWET(1)*POR1
       MOS(m)%soilv7=MOS(m)%SoWET(2)*POR2
       MOS(m)%soilv8=MOS(m)%SoWET(3)*POR3

      mos(m)%count = mos(m)%count + 1
   enddo
    	     
   return
 end subroutine mos_main


	     
	     
	     


