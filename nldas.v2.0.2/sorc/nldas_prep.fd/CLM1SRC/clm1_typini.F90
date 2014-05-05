#include <misc.h>

subroutine clm1_typini (ntiles, clm1)

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
!  initialize clm variables
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!  27 Nov 2001: Jon Gottschalck; Added varaibles for AVHRR LAI
!  28 Aug 2002: Jon Gottschalck; Added forcing variables for surface initialization
!  20 Nov 2002: Jon Radakovich; Added change in skin temperature term for analysis
!=========================================================================
! $Id: clm1_typini.F90,v 1.1.1.1 2003/02/06 16:10:45 jgottsch Exp $
!=========================================================================

  use precision
  use infnan
  use clm1type
  use drv_module      ! 1-D Land Model Driver variables
  use drv_tilemodule  ! Tile-space variables
  implicit none

!=== Arguments ===========================================================  

  integer, intent(in)         :: ntiles     !number of tiles
  type (clm11d), intent(inout) :: clm1(ntiles)

!=== Local Variables =====================================================

  integer ::  k

!=========================================================================

  do k = 1, ntiles

     clm1(k)%itypwat = bigint   ! water type
     clm1(k)%itypprc = bigint   ! precipitation type (from met data) 1= rain 2 =snow
     clm1(k)%isoicol = bigint   ! color classes for soil albedos
!     clm1(k)%latdeg = NaN       ! latitude (degrees)
!     clm1(k)%londeg = NaN       ! longitude (degrees)
                   
     clm1(k)%dtime = NaN        ! model time step [second]
     clm1(k)%istep = bigint     ! number of time step

! Leaf constants (read into 2-D grid module variables)

     clm1(k)%dewmx  = NaN    ! Maximum allowed dew [mm]

! Roughness lengths (read into 2-D grid module variables)

     clm1(k)%zlnd   = NaN    ! Roughness length for soil [m]
     clm1(k)%zsno   = NaN    ! Roughness length for snow [m]
     clm1(k)%csoilc = NaN    ! Drag coefficient for soil under canopy [-]

! Hydraulic constants of soil (read into 2-D grid module variables)

     clm1(k)%wtfact = NaN   ! Fraction of model area with high water table

! Numerical finite-difference(read into 2-D grid module variables)

     clm1(k)%capr   = NaN    ! Tuning factor to turn first layer T into surface T
     clm1(k)%cnfac  = NaN    ! Crank Nicholson factor between 0 and 1
     clm1(k)%smpmin = NaN    ! Restriction for min of soil poten. (mm)
     clm1(k)%ssi    = NaN    ! Irreducible water saturation of snow
     clm1(k)%wimp   = NaN    ! Water impermeable if porosity < wimp
     clm1(k)%pondmx = NaN    ! Ponding depth (mm)

! Vegetation static, dynamic, derived parameters

     clm1(k)%fdry   = NaN   ! fraction of foliage that is green and dry [-]
     clm1(k)%fwet   = NaN   ! fraction of foliage covered by water [-]
     clm1(k)%tlai   = NaN   ! time interpolated leaf area index
     clm1(k)%tsai   = NaN   ! time interpolated stem area index
     clm1(k)%elai   = NaN   ! exposed leaf area index
     clm1(k)%esai   = NaN   ! exposed stem area index
     clm1(k)%minlai = NaN   ! minimum leaf area index
     clm1(k)%maxlai = NaN   ! maximum leaf area index
     clm1(k)%lai_t1_f = NaN ! time 1 AVHRR lai
     clm1(k)%lai_t2_f = NaN ! time 2 AVHRR lai
     clm1(k)%sai_t1_f = NaN ! time 1 AVHRR sai
     clm1(k)%sai_t2_f = NaN ! time 2 AVHRR sai
     clm1(k)%z0m    = NaN   ! aerodynamic roughness length [m]
     clm1(k)%displa = NaN   ! displacement height [m]
     clm1(k)%dleaf  = NaN   ! inverse sqrt of leaf dimension [m**-0.5]

! Soil physical parameters

     clm1(k)%bsw   (:) = NaN   ! Clapp and Hornberger "b"
     clm1(k)%watsat(:) = NaN   ! volumetric soil water at saturation (porosity)
     clm1(k)%hksat (:) = NaN   ! hydraulic conductivity at saturation (mm H2O /s)
     clm1(k)%sucsat(:) = NaN   ! minimum soil suction (mm)
     clm1(k)%csol  (:) = NaN   ! heat capacity, soil solids (J/m**3/Kelvin)
     clm1(k)%tkmg  (:) = NaN   ! thermal conductivity, soil minerals  [W/m-K]  
     clm1(k)%tkdry (:) = NaN   ! thermal conductivity, dry soil       (W/m/Kelvin)
     clm1(k)%tksatu(:) = NaN   ! thermal conductivity, saturated soil [W/m-K]  
     clm1(k)%rootfr(:) = NaN   ! fraction of roots in each soil layer
     clm1(k)%rootr (:) = NaN   ! effective fraction of roots in each soil layer

! Forcing

     clm1(k)%forc_u     = NaN  ! wind speed in eastward direction [m/s]
     clm1(k)%forc_v     = NaN  ! wind speed in northward direction [m/s]
     clm1(k)%forc_t     = NaN  ! temperature at agcm reference height [kelvin]
     clm1(k)%forc_q     = NaN  ! specific humidity at agcm reference height [kg/kg]
     clm1(k)%forc_rain  = NaN  ! rain rate [mm/s]
     clm1(k)%forc_snow  = NaN  ! snow rate [mm/s]
     clm1(k)%forc_pbot  = NaN  ! atmosphere pressure at the surface [pa]
     clm1(k)%forc_rho   = NaN  ! density air [kg/m3]
     clm1(k)%forc_hgt_u = NaN  ! observational height of wind [m]
     clm1(k)%forc_hgt_t = NaN  ! observational height of temperature [m]
     clm1(k)%forc_hgt_q = NaN  ! observational height of humidity [m]
     clm1(k)%forc_lwrad = NaN  ! atmospheric infrared (longwave) radiation [W/m2]

! Forcing for surface initialization

     clm1(k)%forc_swc1   = NaN ! Soil water content (0-10 cm) [m3/m3]
     clm1(k)%forc_swc2   = NaN ! Soil water content (10-100 cm) [m3/m3]
     clm1(k)%forc_stemp1 = NaN ! Soil temperature (0-10 cm) [K]
     clm1(k)%forc_stemp2 = NaN ! Soil temperature (10-100 cm) [K]
     clm1(k)%forc_sdepth = NaN ! Snow depth [m]

! Main variables needed for restart

     clm1(k)%snl = bigint              ! number of snow layers
     clm1(k)%frac_veg_nosno = bigint   ! fraction of veg cover, excluding snow-covered veg (now 0 OR 1) [-]

     clm1(k)%zi(:)         = NaN   ! interface level below a "z" level (m)
     clm1(k)%dz(:)         = NaN   ! layer depth (m)
     clm1(k)%z(:)          = NaN   ! layer thickness (m)
     clm1(k)%t_soisno(:)   = NaN   ! soil + snow layer temperature [K]
     clm1(k)%h2osoi_liq(:) = NaN   ! liquid water (kg/m2)
     clm1(k)%h2osoi_ice(:) = NaN   ! ice lens (kg/m2)

     clm1(k)%frac_sno       = NaN   ! fractional snow cover
     clm1(k)%t_veg          = NaN   ! leaf temperature [K]
     clm1(k)%h2ocan         = NaN   ! depth of water on foliage [kg/m2/s]
     clm1(k)%snowage        = NaN   ! non dimensional snow age [-]
     clm1(k)%h2osno         = NaN   ! snow mass (kg/m2)
     clm1(k)%h2osno_old     = NaN   ! snow mass for previous time step (kg/m2)
     clm1(k)%snowdp         = NaN   ! snow depth (m)
     clm1(k)%t_grnd         = NaN   ! ground surface temperature [k]

! Fluxes

     clm1(k)%taux           = NaN   ! wind stress: E-W [kg/m/s**2]
     clm1(k)%tauy           = NaN   ! wind stress: N-S [kg/m/s**2]
     clm1(k)%eflx_lh_tot    = NaN   ! latent heat flux from canopy height to atmosphere [W/2]
     clm1(k)%eflx_sh_tot    = NaN   ! sensible heat from canopy height to atmosphere [W/m2]
     clm1(k)%eflx_sh_grnd   = NaN   ! sensible heat flux from ground [W/m2]
     clm1(k)%eflx_sh_veg    = NaN   ! sensible heat from leaves [W/m2]
     clm1(k)%qflx_evap_tot  = NaN   ! evapotranspiration from canopy height to atmosphere [mm/s]
     clm1(k)%qflx_evap_veg  = NaN   ! evaporation+transpiration from leaves [mm/s]
     clm1(k)%qflx_evap_soi  = NaN   ! evaporation heat flux from ground [mm/s]
     clm1(k)%qflx_tran_veg  = NaN   ! transpiration rate [mm/s]
     clm1(k)%eflx_lwrad_out = NaN   ! outgoing long-wave radiation from ground+canopy
     clm1(k)%eflx_soil_grnd = NaN   ! ground heat flux [W/m2]
     clm1(k)%qflx_surf      = NaN   ! surface runoff (mm h2o/s)
     clm1(k)%t_ref2m        = NaN   ! 2 m height air temperature [K]
     clm1(k)%t_rad          = NaN   ! radiative temperature [K]

! Diagnostic Variables

     clm1(k)%diagsurf(:)   = NaN  ! Surface diagnostics defined by user
     clm1(k)%diagsoil(:,:) = NaN  ! Soil layer diagnostics defined by user
     clm1(k)%diagsnow(:,:) = NaN  ! Snow layer diagnostics defined by user
! Initialized in drv_readclm1in in LDAS
!     clm1(k)%surfind = bigint  ! Number of surface diagnostic variables
!     clm1(k)%soilind = bigint  ! Number of soil layer diagnostic variables
!     clm1(k)%snowind = bigint  ! Number of snow layer diagnostic variables

! hydrology 

     clm1(k)%imelt(:) = bigint      ! Flag for melting (=1), freezing (=2), Not=0         
     clm1(k)%frac_iceold(:) = NaN   ! Fraction of ice relative to the total water

     clm1(k)%sfact          = NaN  ! term for implicit correction to evaporation
     clm1(k)%sfactmax       = NaN  ! maximim of "sfact"
     clm1(k)%qflx_snow_grnd = NaN  ! ice onto ground [kg/(m2 s)]
     clm1(k)%qflx_rain_grnd = NaN  ! liquid water onto ground [kg/(m2 s)]
     clm1(k)%qflx_evap_grnd = NaN  ! ground surface evaporation rate (mm h2o/s)
     clm1(k)%qflx_dew_grnd  = NaN  ! ground surface dew formation (mm h2o /s) [+]
     clm1(k)%qflx_sub_snow  = NaN  ! sublimation rate from snow pack (mm h2o /s) [+]
     clm1(k)%qflx_dew_snow  = NaN  ! surface dew added to snow pack (mm h2o /s) [+]
     clm1(k)%qflx_snomelt   = NaN  ! rate of snowmelt [kg/(m2 s)]

! added to be consistent with LSM

     clm1(k)%eflx_snomelt = NaN 

     clm1(k)%rhol(:) = NaN  ! pft_varcon leaf reflectance  : 1=vis, 2=nir 
     clm1(k)%rhos(:) = NaN  ! pft_varcon stem reflectance  : 1=vis, 2=nir 
     clm1(k)%taus(:) = NaN  ! pft_varcon stem transmittance: 1=vis, 2=nir 
     clm1(k)%taul(:) = NaN  ! pft_varcon leaf transmittance: 1=vis, 2=nir 
     clm1(k)%xl      = NaN  ! pft_varcon leaf/stem orientation index
     clm1(k)%htop    = NaN  ! canopy top height (m)
     clm1(k)%hbot    = NaN  ! canopy bottom height (m)

! Surface solar radiation 

     clm1(k)%rssun  = NaN  ! sunlit stomatal resistance (s/m)
     clm1(k)%rssha  = NaN  ! shaded stomatal resistance (s/m)
     clm1(k)%psnsun = NaN  ! sunlit leaf photosynthesis (umol CO2 /m**2/ s) 
     clm1(k)%psnsha = NaN  ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
     clm1(k)%laisun = NaN  ! sunlit leaf area
     clm1(k)%laisha = NaN  ! shaded leaf area
     clm1(k)%sabg   = NaN  ! solar radiation absorbed by ground (W/m**2)
     clm1(k)%sabv   = NaN  ! solar radiation absorbed by vegetation (W/m**2)
     clm1(k)%fsa    = NaN  ! solar radiation absorbed (total) (W/m**2)
     clm1(k)%fsr    = NaN  ! solar radiation reflected (W/m**2)
     clm1(k)%ndvi   = NaN  ! Normalized Difference Vegetation Index (diagnostic)

! surfacealbedo 

     clm1(k)%parsun    = NaN  ! average absorbed PAR for sunlit leaves (W/m**2)
     clm1(k)%parsha    = NaN  ! average absorbed PAR for shaded leaves (W/m**2)
     clm1(k)%albd(:)   = NaN  ! surface albedo (direct)                     
     clm1(k)%albi(:)   = NaN  ! surface albedo (diffuse)                    
     clm1(k)%albgrd(:) = NaN  ! ground albedo (direct)                      
     clm1(k)%albgri(:) = NaN  ! ground albedo (diffuse)                     
     clm1(k)%fabd(:)   = NaN  ! flux absorbed by veg per unit direct flux   
     clm1(k)%fabi(:)   = NaN  ! flux absorbed by veg per unit diffuse flux  
     clm1(k)%ftdd(:)   = NaN  ! down direct flux below veg per unit dir flx 
     clm1(k)%ftid(:)   = NaN  ! down diffuse flux below veg per unit dir flx
     clm1(k)%ftii(:)   = NaN  ! down diffuse flux below veg per unit dif flx
     clm1(k)%fsun      = NaN  ! sunlit fraction of canopy                   
     clm1(k)%surfalb   = NaN  ! instantaneous all-wave surface albedo
     clm1(k)%snoalb    = NaN  ! instantaneous all-wave snow albedo

!hydrology

     clm1(k)%h2osoi_vol(:)   = NaN ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
     clm1(k)%eff_porosity(:) = NaN ! effective porosity = porosity - vol_ice

     clm1(k)%qflx_infl       = NaN  ! infiltration (mm H2O /s) 
     clm1(k)%qflx_drain      = NaN  ! sub-surface runoff (mm H2O /s) 
     clm1(k)%qflx_top_soil   = NaN  ! net water input into soil from top (mm/s)
     clm1(k)%qflx_prec_intr  = NaN  ! interception of precipitation [mm/s]
     clm1(k)%qflx_prec_grnd  = NaN  ! water onto ground including canopy runoff [kg/(m2 s)]
     clm1(k)%qflx_qirr       = NaN  ! qflx_surf directed to irrig (mm H2O/s)
     clm1(k)%qflx_qrgwl      = NaN  ! qflx_surf at glaciers, wetlands, lakes
     clm1(k)%btran           = NaN  ! transpiration wetness factor (0 to 1) 
     clm1(k)%smpmax          = NaN  !wilting point potential in mm (new)

     clm1(k)%eflx_impsoil    = NaN  ! implicit evaporation for soil temperature equation (W/m**2)
     clm1(k)%eflx_lh_vege    = NaN  ! veg evaporation heat flux (W/m**2) [+ to atm]
     clm1(k)%eflx_lh_vegt    = NaN  ! veg transpiration heat flux (W/m**2) [+ to atm]
     clm1(k)%eflx_lh_grnd    = NaN  ! ground evaporation heat flux (W/m**2) [+ to atm]   
     clm1(k)%eflx_lwrad_net  = NaN  ! net infrared (longwave) rad (W/m**2) [+ = to atm]

! water and energy balance check

     clm1(k)%begwb       = NaN  ! water mass begining of the time step
     clm1(k)%endwb       = NaN  ! water mass end of the time step
     clm1(k)%errh2o      = NaN  ! water conservation error (mm H2O)
     clm1(k)%errsoi      = NaN  ! soil/lake energy conservation error (W/m**2)
     clm1(k)%errseb      = NaN  ! surface energy conservation error (W/m**2)
     clm1(k)%errsol      = NaN  ! solar radiation conservation error (W/m**2)
     clm1(k)%errlon      = NaN  ! longwave radiation conservation error (W/m**2)
     clm1(k)%acc_errseb  = NaN  ! accumulation of surface energy balance error
     clm1(k)%acc_errh2o  = NaN  ! accumulation of water balance error

!forcing

     clm1(k)%forc_solad(:) = NaN   ! direct beam radiation (vis=forc_sols , nir=forc_soll )
     clm1(k)%forc_solai(:) = NaN   ! diffuse radiation     (vis=forc_solsd, nir=forc_solld)

! temperatures

     clm1(k)%dt_veg  = NaN         ! change in t_veg, last iteration (Kelvin)
     clm1(k)%dt_grnd = NaN         ! change in t_grnd, last iteration (Kelvin)

! new lsm terms from pft_varcon - to avoid indirect indexing

     clm1(k)%z0m          = NaN    ! aerodynamic roughness length [m]
     clm1(k)%displa       = NaN    ! displacement height [m]
     clm1(k)%dleaf        = NaN    ! leaf dimension [m]
     clm1(k)%qe25         = NaN    ! quantum efficiency at 25c (umol co2 / umol photon)
     clm1(k)%ko25         = NaN    ! o2 michaelis-menten constant at 25c (pa)
     clm1(k)%kc25         = NaN    ! co2 michaelis-menten constant at 25c (pa)
     clm1(k)%vcmx25       = NaN    ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
     clm1(k)%ako          = NaN    ! q10 for ko25
     clm1(k)%akc          = NaN    ! q10 for kc25
     clm1(k)%avcmx        = NaN    ! q10 for vcmx25
     clm1(k)%bp           = NaN    ! minimum leaf conductance (umol/m**2/s)
     clm1(k)%mp           = NaN    ! slope for conductance-to-photosynthesis relationship
     clm1(k)%folnmx       = NaN    ! foliage nitrogen concentration when f(n)=1 (%)
     clm1(k)%folnvt       = NaN    ! foliage nitrogen concentration (%)
     clm1(k)%c3psn        = NaN    ! photosynthetic pathway: 0. = c4, 1. = c3

! alma output

     clm1(k)%diffusion      = NaN  ! heat diffusion through layer zero interface 
     clm1(k)%h2osoi_liq_old = NaN  ! liquid water from previous timestep
     clm1(k)%h2ocan_old     = NaN  ! depth of water on foliage from previous timestep
     clm1(k)%acond          = 0.  ! aerodynamic conductance (m/s)

! LDAS standard output totalizing arrays
     clm1(k)%totfsa             = NaN ! solar absorbed by ground + vegetation [W/m2]
     clm1(k)%toteflx_lwrad_net  = NaN ! net longwave radiation [W/m2]
     clm1(k)%toteflx_lh_tot     = NaN ! latent heat flux from canopy height to atmosphere [W/2]
     clm1(k)%toteflx_sh_tot     = NaN ! sensible heat from canopy height to atmosphere [W/m2]
     clm1(k)%toteflx_soil_grnd  = NaN ! ground heat flux [W/m2]  
     clm1(k)%totqflx_snomelt    = NaN ! rate of snowmelt [kg/(m2 s)]
     clm1(k)%totsolisbd         = NaN ! clm1(k)%total downward surface shortwave radiation [W/m2]
     clm1(k)%totforc_lwrad      = NaN ! atmospheric infrared (longwave) radiation [W/m2]
     clm1(k)%totsnow            = NaN ! accumulation of snow [mm]
     clm1(k)%totrain            = NaN ! accumulation of rain [mm]
     clm1(k)%totqflx_evap       = NaN ! clm1(k)%total evaporation [mm]    
     clm1(k)%totqflx_surf       = NaN ! surface runoff (mm h2o/s) 
     clm1(k)%totqflx_drain      = NaN ! subsurface runoff (mm h2o/s) 
     clm1(k)%totqflx_ecanop     = NaN ! interception evaporation [W/m2]
     clm1(k)%totqflx_tran_veg   = NaN ! transpiration rate [mm/s]  
     clm1(k)%totqflx_evap_grnd  = NaN ! evaporation heat flux from ground [mm/s]
     clm1(k)%totqflx_sub_snow   = NaN ! sublimation rate from snow pack (mm h2o /s) [+]
     clm1(k)%dtcanal             = NaN ! change in skin temperature
     clm1(k)%count              = bigint

  end do

end subroutine clm1_typini



















