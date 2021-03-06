#include <misc.h>

module clm1type

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
!  Module for 1-D (vertical) CLM1 variable specification.
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!  27 Nov 2001: Jon Gottschalck; Added variables for AVHRR LAI
!  28 Jan 2002: Jon Gottschalck; Added variables for initialization
!  15 May 2002: Urszula Jambor; Changed LOGICAL to LOGICAL*1 to match new 
!                GRIB libraries
!  20 Nov 2002: Jon Radakovich; Added term for analysis and BC
!=========================================================================
! $Id: clm1type.F90,v 1.1.1.1 2003/02/06 16:10:47 jgottsch Exp $
!=========================================================================

  use precision
  use clm1_varpar
  implicit none
  public clm11d

  type clm11d

!=== Arguments ===========================================================

! time invariant variables

     integer :: kpatch     ! patch index
     integer :: itypwat    ! water type
     logical*1 :: lakpoi     ! true => lake point
     logical*1 :: baresoil   ! true => bare soil  
     logical*1 :: irrig      ! true => for irrigated soil
     integer :: itypprc    ! precipitation type (from met data) 1= rain 2 =snow
     integer :: isoicol    ! color classes for soil albedos
     real(r8):: latdeg     ! latitude  (degrees)
     real(r8):: londeg     ! longitude (degrees)
     real(r8):: lat        ! latitude  (radians)
     real(r8):: lon        ! longitude (radians)
                           
     real(r8) dtime        ! model time step [second]
     integer  istep        ! number of time step

! Leaf constants (read into 2-D grid module variables)

     real(r8) :: dewmx     ! Maximum allowed dew [mm]

! Roughness lengths (read into 2-D grid module variables)

     real(r8) :: zlnd      ! Roughness length for soil [m]
     real(r8) :: zsno      ! Roughness length for snow [m]
     real(r8) :: csoilc    ! Drag coefficient for soil under canopy [-]

! Hydraulic constants of soil (read into 2-D grid module variables)

     real(r8) :: wtfact    ! Fraction of model area with high water table

! Numerical finite-difference(read into 2-D grid module variables)

     real(r8) :: capr      ! Tuning factor to turn first layer T into surface T
     real(r8) :: cnfac     ! Crank Nicholson factor between 0 and 1
     real(r8) :: smpmin    ! Restriction for min of soil poten. (mm)
     real(r8) :: ssi       ! Irreducible water saturation of snow
     real(r8) :: wimp      ! Water impermeable if porosity < wimp
     real(r8) :: pondmx    ! Ponding depth (mm)

! Vegetation static, dynamic, derived parameters

     real(r8) :: fdry      ! fraction of foliage that is green and dry [-]
     real(r8) :: fwet      ! fraction of foliage covered by water [-]
     real(r8) :: tlai      ! time interpolated leaf area index
     real(r8) :: tsai      ! time interpolated stem area index
     real(r8) :: elai      ! exposed leaf area index
     real(r8) :: esai      ! exposed stem area index
     real(r8) :: minlai    ! minimum leaf area index
     real(r8) :: maxlai    ! maximum leaf area index
     integer  :: laiflag   ! Flag to decide which LAI parameterization to use
     real(r8) :: lai_t1_f  ! Time 1 AVHRR LAI
     real(r8) :: lai_t2_f  ! Time 2 AVHRR LAI
     real(r8) :: sai_t1_f  ! Time 1 AVHRR SAI
     real(r8) :: sai_t2_f  ! Time 2 AVHRR SAI

! Soil physical parameters

     real(r8) :: bsw   (nlevsoi) ! Clapp and Hornberger "b"
     real(r8) :: watsat(nlevsoi) ! volumetric soil water at saturation (porosity)
     real(r8) :: hksat (nlevsoi) ! hydraulic conductivity at saturation (mm H2O /s)
     real(r8) :: sucsat(nlevsoi) ! minimum soil suction (mm)
     real(r8) :: csol  (nlevsoi) ! heat capacity, soil solids (J/m**3/Kelvin)
     real(r8) :: tkmg  (nlevsoi) ! thermal conductivity, soil minerals  [W/m-K]  
     real(r8) :: tkdry (nlevsoi) ! thermal conductivity, dry soil       (W/m/Kelvin)
     real(r8) :: tksatu(nlevsoi) ! thermal conductivity, saturated soil [W/m-K]  
     real(r8) :: rootfr(nlevsoi) ! fraction of roots in each soil layer
     real(r8) :: rootr(nlevsoi)  ! effective fraction of roots in each soil layer

! Forcing

     real(r8) :: forc_u           ! wind speed in eastward direction [m/s]
     real(r8) :: forc_v           ! wind speed in northward direction [m/s]
     real(r8) :: forc_t           ! temperature at agcm reference height [kelvin]
     real(r8) :: forc_q           ! specific humidity at agcm reference height [kg/kg]
     real(r8) :: forc_rain        ! rain rate [mm/s]
     real(r8) :: forc_snow        ! snow rate [mm/s]
     real(r8) :: forc_pbot        ! atmosphere pressure at the surface [pa]
     real(r8) :: forc_rho         ! density air [kg/m3]
     real(r8) :: forc_hgt_u       ! observational height of wind [m]
     real(r8) :: forc_hgt_t       ! observational height of temperature [m]
     real(r8) :: forc_hgt_q       ! observational height of humidity [m]
     real(r8) :: forc_lwrad       ! atmospheric infrared (longwave) radiation [W/m2]
     
! Added forcing variables for initialization
     real(r8) :: forc_swc1        ! Layer 1 (0-10 cm) soil water content
     real(r8) :: forc_swc2        ! Layer 2 (10-200 cm) soil water content (m3/m3)
     real(r8) :: forc_stemp1      ! Layer 1 (0-10 cm) soil temperature (K)
     real(r8) :: forc_stemp2      ! Layer 2 (10-200 cm) soil temperature (K)
     real(r8) :: forc_sdepth      ! Model liquid equivalent snow depth (kg/m2)

! Main variables needed for restart

     integer ::  snl              ! number of snow layers
     integer frac_veg_nosno       ! fraction of veg cover, excluding snow-covered veg (now 0 OR 1) [-]

     real(r8) :: zi(-nlevsno+0:nlevsoi)          !interface level below a "z" level (m)
     real(r8) :: dz(-nlevsno+1:nlevsoi)          !layer depth (m)
     real(r8) :: z (-nlevsno+1:nlevsoi)          !layer thickness (m)
     real(r8) :: t_soisno  (-nlevsno+1:nlevsoi)  !soil temperature (Kelvin)
     real(r8) :: h2osoi_liq(-nlevsno+1:nlevsoi)  !liquid water (kg/m2) (new)
     real(r8) :: h2osoi_ice(-nlevsno+1:nlevsoi)  !ice lens (kg/m2) (new)

     real(r8) :: frac_sno        ! fractional snow cover
     real(r8) :: t_veg           ! leaf temperature [K]
     real(r8) :: h2ocan          ! depth of water on foliage [kg/m2/s]
     real(r8) :: snowage         ! non dimensional snow age [-]
     real(r8) :: h2osno          ! snow mass (kg/m2)
     real(r8) :: h2osno_old      ! snow mass for previous time step (kg/m2)
     real(r8) :: snowdp          ! snow depth (m)
     real(r8) :: t_grnd          ! ground surface temperature [k]

! Fluxes

     real(r8) :: taux                 ! wind stress: E-W [kg/m/s**2]
     real(r8) :: tauy                 ! wind stress: N-S [kg/m/s**2]
     real(r8) :: eflx_lh_tot          ! latent heat flux from canopy height to atmosphere [W/2]
     real(r8) :: eflx_sh_tot          ! sensible heat from canopy height to atmosphere [W/m2]
     real(r8) :: eflx_sh_grnd         ! sensible heat flux from ground [W/m2]
     real(r8) :: eflx_sh_veg          ! sensible heat from leaves [W/m2]
     real(r8) :: qflx_evap_tot        ! evapotranspiration from canopy height to atmosphere [mm/s]
     real(r8) :: qflx_evap_veg        ! evaporation+transpiration from leaves [mm/s]
     real(r8) :: qflx_evap_soi        ! evaporation heat flux from ground [mm/s]
     real(r8) :: qflx_tran_veg        ! transpiration rate [mm/s]
     real(r8) :: eflx_lwrad_out       ! outgoing long-wave radiation from ground+canopy
     real(r8) :: eflx_soil_grnd       ! ground heat flux [W/m2]
     real(r8) :: qflx_surf            ! surface runoff (mm h2o/s)
     real(r8) :: t_ref2m              ! 2 m height air temperature [K]
     real(r8) :: t_rad                ! radiative temperature [K]

! Diagnostic Variables

     real(r8), pointer :: diagsurf(:)    ! Surface diagnostics defined by user
     real(r8), pointer :: diagsoil(:,:)  ! Soil layer diagnostics defined by user
     real(r8), pointer :: diagsnow(:,:)  ! Snow layer diagnostics defined by user

     integer :: surfind   !Number of surface diagnostic variables
     integer :: soilind   !Number of soil layer diagnostic variables
     integer :: snowind   !Number of snow layer diagnostic variables

! hydrology 

     integer  :: imelt      (-nlevsno+1:nlevsoi) ! Flag for melting (=1), freezing (=2), Not=0 
     real(r8) :: frac_iceold(-nlevsno+1:nlevsoi) ! fraction of ice relative to the total water 

     real(r8) :: sfact            ! term for implicit correction to evaporation
     real(r8) :: sfactmax         ! maximim of "sfact"
     real(r8) :: qflx_snow_grnd   ! ice onto ground [kg/(m2 s)]
     real(r8) :: qflx_rain_grnd   ! liquid water onto ground [kg/(m2 s)]
     real(r8) :: qflx_evap_grnd   ! ground surface evaporation rate (mm h2o/s)
     real(r8) :: qflx_dew_grnd    ! ground surface dew formation (mm h2o /s) [+]
     real(r8) :: qflx_sub_snow    ! sublimation rate from snow pack (mm h2o /s) [+]
     real(r8) :: qflx_dew_snow    ! surface dew added to snow pack (mm h2o /s) [+]
     real(r8) :: qflx_snomelt     ! rate of snowmelt [kg/(m2 s)]

! Surface solar radiation 

     real(r8) :: rssun          ! sunlit stomatal resistance (s/m)
     real(r8) :: rssha          ! shaded stomatal resistance (s/m)
     real(r8) :: psnsun         ! sunlit leaf photosynthesis (umol CO2 /m**2/ s) 
     real(r8) :: psnsha         ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
     real(r8) :: laisun         ! sunlit leaf area
     real(r8) :: laisha         ! shaded leaf area
     real(r8) :: sabg           ! solar radiation absorbed by ground (W/m**2)
     real(r8) :: sabv           ! solar radiation absorbed by vegetation (W/m**2)
     real(r8) :: fsa            ! solar radiation absorbed (total) (W/m**2)
     real(r8) :: fsr            ! solar radiation reflected (W/m**2)
     real(r8) :: ndvi           ! Normalized Difference Vegetation Index (diagnostic)

! surfacealbedo 

     real(r8) :: parsun         ! average absorbed PAR for sunlit leaves (W/m**2)
     real(r8) :: parsha         ! average absorbed PAR for shaded leaves (W/m**2)
     real(r8) :: albd(numrad)   ! surface albedo (direct)                     
     real(r8) :: albi(numrad)   ! surface albedo (diffuse)                    
     real(r8) :: albgrd(numrad) ! ground albedo (direct)                      
     real(r8) :: albgri(numrad) ! ground albedo (diffuse)                     
     real(r8) :: fabd(numrad)   ! flux absorbed by veg per unit direct flux   
     real(r8) :: fabi(numrad)   ! flux absorbed by veg per unit diffuse flux  
     real(r8) :: ftdd(numrad)   ! down direct flux below veg per unit dir flx 
     real(r8) :: ftid(numrad)   ! down diffuse flux below veg per unit dir flx
     real(r8) :: ftii(numrad)   ! down diffuse flux below veg per unit dif flx
     real(r8) :: fsun           ! sunlit fraction of canopy                   
     real(r8) :: surfalb        ! instantaneous all-wave surface albedo
     real(r8) :: snoalb         ! instantaneous all_wave snow albedo

!hydrology

     real(r8) :: h2osoi_vol(nlevsoi)     ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
     real(r8) :: eff_porosity(nlevsoi)   ! effective porosity = porosity - vol_ice

     real(r8) :: qflx_infl      ! infiltration (mm H2O /s) 
     real(r8) :: qflx_drain     ! sub-surface runoff (mm H2O /s) 
     real(r8) :: qflx_top_soil  ! net water input into soil from top (mm/s)
     real(r8) :: qflx_prec_intr ! interception of precipitation [mm/s]
     real(r8) :: qflx_prec_grnd ! water onto ground including canopy runoff [kg/(m2 s)]
     real(r8) :: qflx_qirr      ! qflx_surf directed to irrig (mm H2O/s)
     real(r8) :: qflx_qrgwl     ! qflx_surf at glaciers, wetlands, lakes
     real(r8) :: btran          ! transpiration wetness factor (0 to 1) 
     real(r8) :: smpmax         ! wilting point potential in mm (new)

     real(r8) :: eflx_snomelt   ! added to be consistent with lsm hybrid code
     real(r8) :: eflx_impsoil   ! implicit evaporation for soil temperature equation (W/m**2)
     real(r8) :: eflx_lh_vege   ! veg evaporation heat flux (W/m**2) [+ to atm]
     real(r8) :: eflx_lh_vegt   ! veg transpiration heat flux (W/m**2) [+ to atm]
     real(r8) :: eflx_lh_grnd   ! ground evaporation heat flux (W/m**2) [+ to atm]   
     real(r8) :: eflx_lwrad_net ! net infrared (longwave) rad (W/m**2) [+ = to atm]

! water and energy balance check

     real(r8) :: begwb                 !water mass begining of the time step
     real(r8) :: endwb                 !water mass end of the time step
     real(r8) :: errh2o                !water conservation error (mm H2O)
     real(r8) :: errsoi                !soil/lake energy conservation error (W/m**2)
     real(r8) :: errseb                !surface energy conservation error (W/m**2)
     real(r8) :: errsol                !solar radiation conservation error (W/m**2)
     real(r8) :: errlon                !longwave radiation conservation error (W/m**2)
     real(r8) :: acc_errseb            !accumulation of surface energy balance error
     real(r8) :: acc_errh2o            !accumulation of water balance error

!forcing

     real(r8) :: forc_solad(numrad)    !direct beam radiation (vis=forc_sols , nir=forc_soll )
     real(r8) :: forc_solai(numrad)    !diffuse radiation     (vis=forc_solsd, nir=forc_solld)

! temperatures

     real(r8) :: dt_veg                !change in t_veg, last iteration (Kelvin)
     real(r8) :: dt_grnd               !change in t_grnd, last iteration (Kelvin)

! new lsm terms from pft_varcon - to avoid indirect indexing

     real(r8) :: z0m            ! aerodynamic roughness length [m]
     real(r8) :: displa         ! displacement height [m]
     real(r8) :: dleaf          ! leaf dimension [m]
     real(r8) :: xl             ! pft_varcon leaf/stem orientation index
     real(r8) :: rhol(numrad)   ! pft_varcon leaf reflectance  : 1=vis, 2=nir 
     real(r8) :: rhos(numrad)   ! pft_varcon stem reflectance  : 1=vis, 2=nir 
     real(r8) :: taul(numrad)   ! pft_varcon leaf transmittance: 1=vis, 2=nir 
     real(r8) :: taus(numrad)   ! pft_varcon stem transmittance: 1=vis, 2=nir 
     real(r8) :: qe25           ! quantum efficiency at 25c (umol co2 / umol photon)
     real(r8) :: ko25           ! o2 michaelis-menten constant at 25c (pa)
     real(r8) :: kc25           ! co2 michaelis-menten constant at 25c (pa)
     real(r8) :: vcmx25         ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
     real(r8) :: ako            ! q10 for ko25
     real(r8) :: akc            ! q10 for kc25
     real(r8) :: avcmx          ! q10 for vcmx25
     real(r8) :: bp             ! minimum leaf conductance (umol/m**2/s)
     real(r8) :: mp             ! slope for conductance-to-photosynthesis relationship
     real(r8) :: folnmx         ! foliage nitrogen concentration when f(n)=1 (%)
     real(r8) :: folnvt         ! foliage nitrogen concentration (%)
     real(r8) :: c3psn          ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8) :: htop           ! canopy top height (m)
     real(r8) :: hbot           ! canopy bottom height (m)

! Variables needed for ALMA output

     real(r8) :: diffusion                 !heat diffusion through layer zero interface 
     real(r8) :: h2osoi_liq_old(1:nlevsoi) !liquid water from previous timestep
     real(r8) :: h2ocan_old                !depth of water on foliage from previous timestep
     real(r8) :: acond                     !aerodynamic conductance (m/s)

! Totalizing arrays for output in LDAS 
     real(r8) :: totfsa                    ! solar absorbed by ground + vegetation [W/m2]
     real(r8) :: toteflx_lwrad_net         ! net longwave radiation [W/m2]
     real(r8) :: toteflx_lh_tot            ! latent heat flux from canopy height to atmosphere [W/2]
     real(r8) :: toteflx_sh_tot            ! sensible heat from canopy height to atmosphere [W/m2]
     real(r8) :: toteflx_soil_grnd         ! ground heat flux [W/m2]  
     real(r8) :: totqflx_snomelt           ! rate of snowmelt [kg/(m2 s)]
     real(r8) :: totsolisbd                ! total downward surface shortwave radiation [W/m2]
     real(r8) :: totforc_lwrad             ! atmospheric infrared (longwave) radiation [W/m2]
     real(r8) :: totsnow                   ! accumulation of snow [mm]
     real(r8) :: totrain                   ! accumulation of rain [mm]
     real(r8) :: totqflx_evap              ! total evaporation [mm]    
     real(r8) :: totqflx_surf              ! surface runoff (mm h2o/s) 
     real(r8) :: totqflx_drain             ! subsurface runoff (mm h2o/s) 
     real(r8) :: totqflx_ecanop            ! interception evaporation [W/m2]
     real(r8) :: totqflx_tran_veg          ! transpiration rate [mm/s]  
     real(r8) :: totqflx_evap_grnd         ! evaporation heat flux from ground [mm/s]
     real(r8) :: totqflx_sub_snow          ! sublimation rate from snow pack (mm h2o /s) [+]
     integer  :: count                     ! counter used for time averaging of output

! Variables for PSAS temperature assimilation
     real(r8) :: dtcanal                   ! Change in skin temperature based on analysis

!=== End Variable List ===================================================

  end type clm11d

end module clm1type



