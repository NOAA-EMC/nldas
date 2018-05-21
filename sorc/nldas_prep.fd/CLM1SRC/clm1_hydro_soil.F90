#include <misc.h>

subroutine clm1_hydro_soil (clm1) 

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
!  This is the main subroutine used to execute the calculation of water
!  processes over soil
!
!  (1) Water flow within soil (see clm1_soilwater.f90)
!
!  (2) Runoff 
!      The original code was provide by Robert E. Dickinson based on 
!      following clues:  exponential decrease of Ksat, a water table 
!      level determination level including highland and lowland levels 
!      and fractional area of wetland (water table above the surface). 
!      Runoff is parameterized from the lowlands in terms of precip 
!      incident on wet areas and a base flow, where these are estimated 
!      using ideas from TOPMODEL.
!
!  The original scheme was modified by Z.-L. Yang and G.-Y. Niu,
!  *  using a new method to determine water table depth and
!     the fractional wet area (fcov)
!  *  computing runoff (surface and subsurface) from this
!     fraction and the remaining fraction (i.e. 1-fcov)
!  *  for the 1-fcov part, using BATS1e method to compute
!     surface and subsurface runoff.
!
!   The original code on soil moisture and runoff were provided by 
!   R. E. Dickinson in July 1996.
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  12 November 1999:  Z.-L. Yang and G.-Y. Niu
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!  20 November 2002:  Jon Radakovich; Modified so calculations do 
!                     not depend on the number of soil layers
!=========================================================================
! $Id: clm1_hydro_soil.F90,v 1.1.1.1 2003/02/06 16:10:46 jgottsch Exp $
!=========================================================================

  use precision
  use clm1type
  use clm1_varcon, only : denice, denh2o, hvap, tfrz
  use clm1_varpar, only : nlevsoi
  implicit none

! When multiplied by suction*bsw/(water*4.71047e4*t_grnd) gives derivative of
! evaporation with respect to water in top soil layer
! (= d evap/d rel hum  *d rel hum / d soil suction * d sol suction / d water
! and using d evap / d rel humidity = rho * qsat * aerodynamic resistance,
! d rel hum / d soil suction = (rel hum)/(4.71047e4*t_grnd)
! and d soil suction / d soil water = bsw * (soil suction)/water)

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1	 !CLM 1-D Module

!=== Local Variables =====================================================

  integer i                      ! loop counter

  real(r8) hk(1:nlevsoi),       & ! hydraulic conductivity (mm h2o/s)
       dhkdw(1:nlevsoi),        & ! d(hk)/d(vol_liq)
       dwat(1:nlevsoi),         & ! change in soil water
       fcov,                    & ! fractional area with water table at surface
       s(1:nlevsoi),            & ! wetness of soil (including ice)
       s_node,                  & ! vol_liq/porosity
       sdamp,                   & ! extrapolates soiwat dependence of evaporation
       sdampmax,                & ! maximum of "sdamp"
       smp_node,                & ! matrix potential
       vol_liq(1:nlevsoi),      & ! partial volume of liquid water in layer
       vol_ice(1:nlevsoi),      & ! partial volume of ice lens in layer
       wat1,                    & ! water in top soil layer
       xs,                      & ! excess soil water above saturation
       zwice,                   & ! the sum of ice mass of soil (kg/m2)
       zwt,                     & ! water table depth
       zmm (1:nlevsoi),         & ! layer depth (mm)
       dzmm(1:nlevsoi),         & ! layer thickness (mm)
       watmin,                  & ! minimum soil moisture
       hksum                      ! summation of hydraulic cond for layers 6->9

! For Z.-L. Yang & G.-Y. Niu's modification

  real(r8) zmean               ! The surface soil layers contributing to runoff
  real(r8) wmean               ! The averaged soil wetness in surface soil layers
  real(r8) fz                  ! coefficient for water table depth
  real(r8) zsat                ! hydraulic conductivity weighted soil thickness
  real(r8) wsat                ! hydraulic conductivity weighted soil wetness
  real(r8) qflx_drain_wet      ! subsurface runoff from "wet" part (mm h2o/s)
  real(r8) qflx_drain_dry      ! subsurface runoff from "dry" part (mm h2o/s)
  real(r8) dzksum              ! hydraulic conductivity weighted soil thickness
  integer  jwtable             ! soil layer index where water table lies
  
  real(r8) wtold, wtnew        ! temporary variables for irrigation calc

!=== End Variable List ===================================================

!=========================================================================
! [1] Surface runoff
!=========================================================================

! Porosity of soil, partial volume of ice and liquid

  zwice = 0.
  do i = 1,nlevsoi 
     zwice = zwice + clm1%h2osoi_ice(i)
     vol_ice(i) = min(clm1%watsat(i), clm1%h2osoi_ice(i)/(clm1%dz(i)*denice))
     clm1%eff_porosity(i) = clm1%watsat(i)-vol_ice(i)
     vol_liq(i) = min(clm1%eff_porosity(i), clm1%h2osoi_liq(i)/(clm1%dz(i)*denh2o))
  enddo

! Calculate wetness of soil

  do i = 1,nlevsoi
     s(i) = min(1.,(vol_ice(i)+vol_liq(i))/clm1%watsat(i))
  end do

! Determine water table 

  wmean = 0.                                                  
  fz    = 1.0                                                
  do i  = 1, nlevsoi                                        
     wmean = wmean + s(i)*clm1%dz(i)                          
  enddo
  zwt = fz * (clm1%zi(nlevsoi) - wmean)                   
  
! Saturation fraction

  fcov = clm1%wtfact*min(1.,exp(-zwt))
  
! Currently no overland flow parameterization in code is considered
! qflx_surf = 0.   Zong-Liang Yang & G.-Y. Niu                         
!*modified surface runoff according to the concept of TOPMODEL   

  wmean = 0.                                               
  zmean = 0.                                              
  do i = 1, 3                                          
     zmean = zmean + clm1%dz(i)                          
     wmean = wmean + s(i) * clm1%dz(i)                 
  enddo
  wmean = wmean / zmean                           
  clm1%qflx_surf =  max(0., fcov*clm1%qflx_top_soil) + &
                   max(0., (1.-fcov)*min(1.,wmean**4)*clm1%qflx_top_soil)        

! Infiltration into surface soil layer (minus the evaporation)

  if (clm1%snl+1 >= 1) then
     clm1%qflx_infl = clm1%qflx_top_soil - clm1%qflx_surf - clm1%qflx_evap_grnd
  else
     clm1%qflx_infl = clm1%qflx_top_soil - clm1%qflx_surf
  endif
  
!=========================================================================
! [2] Set up r, a, b, and c vectors for tridiagonal solution and renew bl
!=========================================================================

! Following length units are all in millimeters

  do i = 1,nlevsoi 
     zmm(i) = clm1%z(i)*1.e3
     dzmm(i) = clm1%dz(i)*1.e3
  enddo

  sdamp = 0.
  call clm1_soilwater (vol_liq, clm1%eff_porosity, clm1%qflx_infl, sdamp, &
                      dwat   , hk              , dhkdw        , clm1)

! Renew the mass of liquid water

  do i= 1,nlevsoi 
     clm1%h2osoi_liq(i) = clm1%h2osoi_liq(i) + dwat(i)*dzmm(i)
  enddo

!=========================================================================
! [3] Streamflow and total runoff
!=========================================================================

! The amount of streamflow is assumed maintained by flow from the 
! lowland water table with different levels contributing according to 
! their thickness and saturated hydraulic conductivity, i.e. a given 
! layer below the water table interface loses water at a rate per unit 
! depth given by qflx_drain*hk/(sum over all layers below this water table 
! of hk*dz). Because this is a slow smooth process, and not strongly 
! coupled to water in any one layer, it should remain stable for 
! explicit time differencing. Hence, for simplicity it is removed
! explicitly prior to the main soil water calculation.
! Another assumption: no subsurface runoff for ice mixed soil 
! Zong-Liang Yang & G.-Y. Niu                                         

  clm1%qflx_drain = 0.                      ! subsurface runoff
  qflx_drain_wet = 0.                      ! subsurface runoff        
  qflx_drain_dry = 0.                      ! subsurface runoff        
  
! -----------------------------------------------------------------
! Find the layer index corresponding to the water table depth
! -----------------------------------------------------------------

   do i = 1,nlevsoi
      if (zwt <= clm1%zi(i)) goto 10
   end do
10 jwtable = i

  if (jwtable > nlevsoi-1) goto 20

  hksum = 0.
  do i = jwtable,nlevsoi-1                  
     hksum = hksum + hk(i)
  enddo
  if (zwice <= 0. .AND. hksum > 0.) then
     zsat = 0.                                                    
     wsat = 0.                                                    
     dzksum = 0.                                                  
     do i = jwtable,nlevsoi-1                  
        zsat = zsat + clm1%dz(i)*hk(i)                               
        wsat = wsat + s(i)*clm1%dz(i)*hk(i)                         
        dzksum  = dzksum   + hk(i)*clm1%dz(i)                       
     enddo
     wsat = wsat / zsat                                         
     
     qflx_drain_dry = (1.-fcov)*4.e-2* wsat ** (2.*clm1%bsw(1)+3.)  ! mm/s
     qflx_drain_wet = fcov * 1.e-5 * exp(-zwt)                     ! mm/s
     clm1%qflx_drain = qflx_drain_dry + qflx_drain_wet
     
     do i = jwtable, nlevsoi-1                 
        clm1%h2osoi_liq(i) = clm1%h2osoi_liq(i) &
             - clm1%dtime*clm1%qflx_drain*clm1%dz(i)*hk(i)/dzksum                                  
     enddo
  endif
  
20 continue

! --------------------------------------------------------------------
! Limit h2osoi_liq to be greater than or equal to watmin. 
! Get water needed to bring h2osoi_liq equal watmin from lower layer. 
! --------------------------------------------------------------------

  watmin = 0.0
  do i = 1, nlevsoi-1
     if (clm1%h2osoi_liq(i) < 0.) then
        xs = watmin-clm1%h2osoi_liq(i)
     else
        xs = 0.
     end if
     clm1%h2osoi_liq(i  ) = clm1%h2osoi_liq(i  ) + xs
     clm1%h2osoi_liq(i+1) = clm1%h2osoi_liq(i+1) - xs
  end do
  i = nlevsoi
  if (clm1%h2osoi_liq(i) < watmin) then
     xs = watmin-clm1%h2osoi_liq(i)
  else
     xs = 0.
  end if
  clm1%h2osoi_liq(i) = clm1%h2osoi_liq(i) + xs
  clm1%qflx_drain = clm1%qflx_drain - xs/clm1%dtime

! Determine water in excess of saturation

  xs = max(0., clm1%h2osoi_liq(1)-(clm1%pondmx+clm1%eff_porosity(1)*dzmm(1)))
  if (xs > 0.) clm1%h2osoi_liq(1) = clm1%pondmx+clm1%eff_porosity(1)*dzmm(1)
  
  do i = 2,nlevsoi 
     xs = xs + max(clm1%h2osoi_liq(i)-clm1%eff_porosity(i)*dzmm(i), 0.)     ! [mm]
     clm1%h2osoi_liq(i) = min(clm1%eff_porosity(i)*dzmm(i), clm1%h2osoi_liq(i))
  enddo

! Sub-surface runoff and drainage 

  clm1%qflx_drain = clm1%qflx_drain + xs/clm1%dtime  &
       + hk(nlevsoi) + dhkdw(nlevsoi)*dwat(nlevsoi) ! [mm/s]
  
! set imbalance for glacier, lake and wetland to 0. 

  clm1%qflx_qrgwl = 0.   !only set for lakes, wetlands and glaciers 
  
! for now set implicit evaporation to zero

  clm1%eflx_impsoil = 0.

! Renew the ice and liquid mass due to condensation
  
  if (clm1%snl+1 >= 1) then
     clm1%h2osoi_liq(1) = clm1%h2osoi_liq(1) + clm1%qflx_dew_grnd*clm1%dtime
     clm1%h2osoi_ice(1) = clm1%h2osoi_ice(1) + (clm1%qflx_dew_snow-clm1%qflx_sub_snow)*clm1%dtime
  endif

end subroutine clm1_hydro_soil




