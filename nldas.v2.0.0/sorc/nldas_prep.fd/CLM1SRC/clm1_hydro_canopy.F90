#include <misc.h>

subroutine clm1_hydro_canopy (clm1)

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
!  Calculation of 
!  (1) water storage of intercepted precipitation
!  (2) direct throughfall and canopy drainage of precipitation
!  (3) the fraction of foliage covered by water and the fraction
!      of foliage that is dry and transpiring. 
!  (4) snow layer initialization if the snow accumulation exceeds 10 mm.
!
! Note:  The evaporation loss is taken off after the calculation of leaf 
! temperature in the subroutine clm_leaftem.f90 not in this subroutine.
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!=========================================================================

  use precision
  use clm1type
  use clm1_varcon, only : tfrz, istice, istwet, istsoil
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1        !CLM 1-D Module

!=== Local Variables =====================================================

  real(r8)  &
       prcp,              & ! precipitation rate [mm/s]
       h2ocanmx,          & ! maximum allowed water on canopy [mm]
       fpi,               & ! coefficient of interception
       vegt,              & ! frac_veg_nosno*lsai
       xrun,              & ! excess water that exceeds the leaf capacity [mm/s]
       qflx_candrip,      & ! rate of canopy runoff and snow falling off canopy [mm/s]
       qflx_through,      & ! direct throughfall [mm/s]
       dz_snowf,          & ! layer thickness rate change due to precipitation [mm/s]
       flfall,            & ! fraction of liquid water within falling precip.
       bifall               ! bulk density of newly fallen dry snow [kg/m3]

  integer  &
       i,                 & ! do looping index
       newnode              ! signification when new snow node is set, (1=yes, 0=non)

  real(r8)    :: &
       dewmxi               ! inverse of maximum allowed dew [1/mm]

!=== End Variable List ===================================================

! ========================================================================
! [1] Canopy interception and precipitation onto ground surface
! ========================================================================

! 1.1 Add precipitation to leaf water 

  if (clm1%itypwat==istsoil .OR. clm1%itypwat==istwet) then       ! soil or wetland point

     qflx_candrip = 0.                 ! rate of canopy runoff
     qflx_through = 0.                 ! precipitation direct through canopy
     clm1%qflx_prec_intr = 0.           ! intercepted precipitation  

     prcp = clm1%forc_rain + clm1%forc_snow  ! total precipitation

     if (clm1%frac_veg_nosno == 1 .AND. prcp > 0.) then

! The leaf water capacities for solid and liquid are different, 
! generally double for snow, but these are of somewhat less significance
! for the water budget because of lower evap. rate at lower temperature.
! Hence, it is reasonable to assume that vegetation storage of solid water 
! is the same as liquid water.

        h2ocanmx = clm1%dewmx * clm1%frac_veg_nosno * (clm1%elai + clm1%esai)

! Direct throughfall

        fpi = 1. - exp(-0.5*(clm1%elai + clm1%esai))
        qflx_through  = prcp*(1.-fpi)*clm1%frac_veg_nosno

! Water storage of intercepted precipitation and dew

        clm1%qflx_prec_intr = prcp*fpi*clm1%frac_veg_nosno
        clm1%h2ocan = max(0., clm1%h2ocan + clm1%dtime*clm1%qflx_prec_intr)

! Initialize rate of canopy runoff and snow falling off canopy

        qflx_candrip = 0.0

! Excess water that exceeds the leaf capacity

        xrun = (clm1%h2ocan - h2ocanmx)/clm1%dtime

! Test on maximum dew on leaf

        if (xrun > 0.) then
           qflx_candrip = xrun
           clm1%h2ocan = h2ocanmx
        endif

     endif

  else if (clm1%itypwat == istice) then  !land ice

     clm1%qflx_prec_intr = 0.
     clm1%h2ocan = 0.
     qflx_candrip = 0.
     qflx_through = 0.  

  endif

! 1.2 Precipitation onto ground (kg/(m2 s))

  if (clm1%frac_veg_nosno == 0) then
     clm1%qflx_prec_grnd = clm1%forc_rain + clm1%forc_snow
  else
     clm1%qflx_prec_grnd = qflx_through + qflx_candrip  
  endif
  
! 1.3 The percentage of liquid water by mass, which is arbitrarily set to 
!     vary linearly with air temp, from 0% at 273.16 to 40% max at 275.16.

  if (clm1%itypprc <= 1) then
     flfall = 1.                              ! fraction of liquid water within falling precip.
     clm1%qflx_snow_grnd = 0.                  ! ice onto ground (mm/s)
     clm1%qflx_rain_grnd = clm1%qflx_prec_grnd  ! liquid water onto ground (mm/s)
     dz_snowf = 0.                            ! rate of snowfall, snow depth/s (m/s)
  else
     if (clm1%forc_t <= tfrz) then
        flfall = 0.
     else if (clm1%forc_t <= tfrz+2.) then
        flfall = -54.632 + 0.2*clm1%forc_t
     else
        flfall = 0.4
     endif
     
     ! Use Alta relationship, Anderson(1976); LaChapelle(1961), 
     ! U.S.Department of Agriculture Forest Service, Project F, 
     ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.
     
     if (clm1%forc_t > tfrz + 2.) then
        bifall =189.
     else if (clm1%forc_t > tfrz - 15.) then
        bifall=50. + 1.7*(clm1%forc_t - tfrz + 15.)**1.5
     else
        bifall=50.
     endif
     
     clm1%qflx_snow_grnd = clm1%qflx_prec_grnd*(1.-flfall)                 
     clm1%qflx_rain_grnd = clm1%qflx_prec_grnd*flfall
     dz_snowf = clm1%qflx_snow_grnd/bifall                
     clm1%snowdp = clm1%snowdp + dz_snowf*clm1%dtime         
     clm1%h2osno = clm1%h2osno + clm1%qflx_snow_grnd*clm1%dtime      ! snow water equivalent (mm)
     
     if (clm1%itypwat==istwet .AND. clm1%t_grnd>=tfrz) then
        clm1%h2osno=0. 
        clm1%snowdp=0. 
        clm1%snowage=0.
     endif
  endif

! ========================================================================
! [2] Determine the fraction of foliage covered by water and the 
!     fraction of foliage that is dry and transpiring.
! ========================================================================

! fwet is the fraction of all vegetation surfaces which are wet 
! including stem area which contribute to evaporation
! fdry is the fraction of elai (***now in LSM***) which is dry because only leaves
! can transpire.  Adjusted for stem area which does not transpire.

  if (clm1%frac_veg_nosno == 1) then
     if (clm1%h2ocan > 0.) then
        vegt     = clm1%frac_veg_nosno*(clm1%elai + clm1%esai)
        dewmxi   = 1.0/clm1%dewmx
        clm1%fwet = ((dewmxi/vegt)*clm1%h2ocan)**.666666666666
        clm1%fwet = min (clm1%fwet,1.0)     ! Check for maximum limit of fwet
     else
        clm1%fwet = 0.
     endif
     clm1%fdry = (1.-clm1%fwet)*clm1%elai/(clm1%elai+clm1%esai)
#if (defined PERGRO)
     clm1%fwet = 0.
     clm1%fdry = clm1%elai/(clm1%elai+clm1%esai)
#endif
  else
     clm1%fwet = 0.
     clm1%fdry = 0.
  endif

! ========================================================================
! [3] When the snow accumulation exceeds 10 mm, initialize snow layer
! ========================================================================

! Currently, the water temperature for the precipitation is simply set 
! as the surface air temperature

  newnode = 0    ! signification when snow node will be initialized
  if (clm1%snl == 0 .AND. clm1%qflx_snow_grnd > 0.0 .AND. clm1%snowdp >= 0.01) then  
     newnode = 1
     clm1%snl = -1
     clm1%dz(0) = clm1%snowdp                       ! meter
     clm1%z(0) = -0.5*clm1%dz(0)
     clm1%zi(-1) = -clm1%dz(0)
     clm1%snowage = 0.                             ! snow age
     clm1%t_soisno (0) = min(tfrz, clm1%forc_t)     ! K
     clm1%h2osoi_ice(0) = clm1%h2osno               ! kg/m2
     clm1%h2osoi_liq(0) = 0.                       ! kg/m2
     clm1%frac_iceold(0) = 1.
  endif

! The change of ice partial density of surface node due to precipitation
! only ice part of snowfall is added here, the liquid part will be added later

  if (clm1%snl < 0 .AND. newnode == 0) then
     clm1%h2osoi_ice(clm1%snl+1) = clm1%h2osoi_ice(clm1%snl+1)+clm1%dtime*clm1%qflx_snow_grnd
     clm1%dz(clm1%snl+1) = clm1%dz(clm1%snl+1)+dz_snowf*clm1%dtime
  endif

end subroutine clm1_hydro_canopy
