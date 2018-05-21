#include <misc.h>

subroutine clm1_leaftem (z0mv,       z0hv,       z0qv,           &
                        thm,        th,         thv,            & 
                        tg,         qg,         dqgdT,          &
                        htvp,       sfacx,      dqgmax,         &
                        emv,        emg,        dlrad,          &
                        ulrad,      cgrnds,     cgrndl,         &
                        cgrnd,      clm1)

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
!  This subroutine:
!  1. Calculates the leaf temperature: 
!     Use the Newton-Raphson iteration to solve for the foliage 
!     temperature that balances the surface energy budget:
!
!     f(t_veg) = Net radiation - Sensible - Latent = 0
!     f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
!
!  Note:
!  (1) In solving for t_veg, t_grnd is given from the previous timestep.
!  (2) The partial derivatives of aerodynamical resistances, which cannot 
!      be determined analytically, are ignored for d(H)/dT and d(LE)/dT
!  (3) The weighted stomatal resistance of sunlit and shaded foliage is used 
!  (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
!                                                           => Ec + Eg = Ea
!  (5) Energy loss is due to: numerical truncation of energy budget equation
!      (*); and "ecidif" (see the code) which is dropped into the sensible 
!      heat 
!  (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n) and 
!      del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference of 
!      water flux from the leaf between the iteration step (n+1) and (n) 
!      less than 0.1 W/m2; or the iterative steps over 40.
!
!  2. Calculates the leaf fluxes, transpiration, photosynthesis and 
!     updates the dew accumulation due to evaporation.
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!  20 November 2002:  Jon Radakovich; Added term for assimilation and bias correction
!=========================================================================
! $Id: clm1_leaftem.F90,v 1.1.1.1 2003/02/06 16:10:45 jgottsch Exp $
!=========================================================================

! Declare Modules and data structures

  use precision
  use clm1type
  use clm1_varcon, only : sb, cpair, hvap, vkc, grav
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1	 !CLM 1-D Module

  real(r8), intent(in) :: &
       htvp              ! latent heat of evaporation (/sublimation) [J/kg]

! vegetation parameters
  real(r8), intent(in) ::    &
       z0mv,              & ! roughness length, momentum [m]
       z0hv,              & ! roughness length, sensible heat [m]
       z0qv                 ! roughness length, latent heat [m]

! input variables
  real(r8), intent(in) ::    &
       thm,               & ! intermediate variable (forc_t+0.0098*forc_hgt_t)
       th,                & ! potential temperature (kelvin)
       thv                  ! virtual potential temperature (kelvin)

  real(r8), intent(in) ::    &
       tg,                & ! ground surface temperature [K]
       qg,                & ! specific humidity at ground surface [kg/kg]
       dqgdT,             & ! temperature derivative of "qg"
       sfacx,             & ! coefficient for "sfact"
       dqgmax,            & ! max of d(qg)/d(theta)
       emv,               & ! ground emissivity
       emg                  ! vegetation emissivity

  real(r8), intent(inout) :: &
       cgrnd,             & ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
       cgrndl,            & ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
       cgrnds               ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]

  real(r8), intent(out) ::   &
       dlrad,             & ! downward longwave radiation blow the canopy [W/m2]
       ulrad                ! upward longwave radiation above the canopy [W/m2]

!=== Local Variables =====================================================

  real(r8) zldis,         & ! reference height "minus" zero displacement heght [m]
       zii,               & ! convective boundary layer height [m]
       zeta,              & ! dimensionless height used in Monin-Obukhov theory
       beta,              & ! coefficient of conective velocity [-]
       wc,                & ! convective velocity [m/s]
       dth,               & ! diff of virtual temp. between ref. height and surface 
       dthv,              & ! diff of vir. poten. temp. between ref. height and surface
       dqh,               & ! diff of humidity between ref. height and surface
       obu,               & ! Monin-Obukhov length (m)
       um,                & ! wind speed including the stablity effect [m/s]
       ur,                & ! wind speed at reference height [m/s]
       uaf,               & ! velocity of air within foliage [m/s]
       temp1,             & ! relation for potential temperature profile
       temp2,             & ! relation for specific humidity profile
       ustar,             & ! friction velocity [m/s]
       tstar,             & ! temperature scaling parameter
       qstar,             & ! moisture scaling parameter
       thvstar,           & ! virtual potential temperature scaling parameter
       taf,               & ! air temperature within canopy space [K]
       qaf                  ! humidity of canopy air [kg/kg]

  real(r8) rpp,           & ! fraction of potential evaporation from leaf [-]
       rppdry,            & ! fraction of potential evaporation through transp [-]
       cf,                & ! heat transfer coefficient from leaves [-]
       rb,                & ! leaf boundary layer resistance [s/m]
       ram(2),            & ! aerodynamical resistance [s/m]
       rah(2),            & ! thermal resistance [s/m]
       raw(2),            & ! moisture resistance [s/m]
       wta,               & ! heat conduactance for air [m/s]
       wtg,               & ! heat conduactance for ground [m/s]
       wtl,               & ! heat conduactance for leaf [m/s]
       wta0,              & ! normalized heat conduactance for air [-]
       wtl0,              & ! normalized heat conduactance for leaf [-]
       wtg0,              & ! normalized heat conduactance for ground [-]
       wtal,              & ! normalized heat conductance for air and leaf [-]
       wtgl,              & ! normalized heat conductance for leaf and ground [-]
       wtga,              & ! normalized heat cond. for air and ground  [-]
       wtaq,              & ! latent heat conduactance for air [m/s]
       wtlq,              & ! latent heat conduactance for leaf [m/s]
       wtgq,              & ! latent heat conduactance for ground [m/s]
       wtaq0,             & ! normalized latent heat conduactance for air [-]
       wtlq0,             & ! normalized latent heat conduactance for leaf [-]
       wtgq0,             & ! normalized heat conduactance for ground [-]
       wtalq,             & ! normalized latent heat cond. for air and leaf [-]
       wtglq,             & ! normalized latent heat cond. for leaf and ground [-]
       wtgaq,             & ! normalized latent heat cond. for air and ground [-]
       el,                & ! vapor pressure on leaf surface [pa]
       deldT,             & ! derivative of "el" on "t_veg" [pa/K]
       qsatl,             & ! leaf specific humidity [kg/kg]
       qsatldT,           & ! derivative of "qsatl" on "t_veg"
       air,bir,cir,       & ! atmos. radiation temporay set
       dc1,dc2              ! derivative of energy flux [W/m2/K]

  real(r8) delt,          & ! temporary
       delq                 ! temporary

  integer                 & !
       itlef,             & ! counter for leaf temperature iteration [-]
       itmax,             & ! maximum number of iteration [-]
       itmin                ! minimum number of iteration [-]

  real(r8) del,           & ! absolute change in leaf temp in current iteration [K]
       del2,              & ! change in leaf temperature in previous iteration [K]
       dele,              & ! change in latent heat flux from leaf [K]
       delmax,            & ! maximum change in  leaf temperature [K]
       dels,              & ! change in leaf temperature in current iteration [K]
       det,               & ! maximum leaf temp. change in two consecutive iter [K]
       dlemin,            & ! max limit for energy flux convergence [w/m2]
       dtmin,             & ! max limit for temperature convergence [K]
       efeb                 ! latent heat flux from leaf (previous iter) [mm/s]

  real(r8) efeold           ! latent heat flux from leaf (previous iter) [mm/s]

  real(r8) efpot,         & ! potential latent energy flux [kg/m2/s]
       efe,               & ! water flux from leaf [mm/s]
       efsh                 ! sensible heat from leaf [mm/s]

  integer nmozsgn           ! number of times moz changes sign

  real(r8) obuold,           & ! monin-obukhov length from previous iteration
       tlbef,                & ! leaf temperature from previous iteration [K]
       ecidif,               & ! excess energies [W/m2]
       err                     ! balance error
  real(r8) erre                ! balance error

! Constant atmospheric co2 and o2
  real(r8) po2                 ! partial pressure  o2 (mol/mol)
  real(r8) pco2                ! partial pressure co2 (mol/mol)

  data po2,pco2 /0.209,355.e-06/

  real(r8) co2                 ! atmospheric co2 concentration (pa)
  real(r8) o2                  ! atmospheric o2 concentration (pa)

  real(r8) svpts               ! saturation vapor pressure at t_veg (pa)
  real(r8) eah                 ! canopy air vapor pressure (pa)
  real(r8) foln                ! foliage nitrogen (%)

  real(r8) gdir                ! relative projected leaf+stem area in solar direction
  real(r8) ext                 ! optical depth of direct beam per unit leaf area
  real(r8) wl                  ! fraction of lai+sai that is lai

  real(r8) :: mpe = 1.e-6      ! prevents overflow error if division by zero

!=== End Variable List ===================================================

! Initialization

  del   = 0.0  ! change in leaf temperature from previous iteration
  itlef = 0    ! counter for leaf temperature iteration
  efeb  = 0.0  ! latent head flux from leaf for previous iteration

  wtlq = 0.0
  wtlq0 = 0.0
  wtgq0 = 0.0
  wtalq = 0.0
  wtgaq = 0.0
  wtglq = 0.0
  wtaq = 0.0
  wtgq = 0.0
  wtaq0 = 0.0
  wtlq0 = 0.0
  wtgq0 = 0.0
  wtalq = 0.0
  wtgaq = 0.0
  wtglq = 0.0

! Assign iteration parameters

  delmax = 1.0  ! maximum change in  leaf temperature
  itmax  = 40   ! maximum number of iteration
  itmin  = 2    ! minimum number of iteration
  dtmin  = 0.01 ! max limit for temperature convergence
  dlemin = 0.1  ! max limit for energy flux convergence

! Net absorbed longwave radiation by canopy and ground
! =air+bir*t_veg**4+cir*t_grnd**4

  air =   clm1%frac_veg_nosno * emv * (1.+(1.-emv)*(1.-emg)) * clm1%forc_lwrad
  bir = - clm1%frac_veg_nosno * (2.-emv*(1.-emg)) * emv * sb
  cir =   clm1%frac_veg_nosno * emv*emg*sb

! Saturated vapor pressure and humidity and their derivation

  call clm1_qsadv (clm1%t_veg, clm1%forc_pbot, el, deldT, qsatl, qsatldT)

! For use Bonan's stomatal resistance scheme
! atmospheric co2 and o2 are currently constants

  co2 = pco2*clm1%forc_pbot
  o2  = po2*clm1%forc_pbot

! Initialize flux profile

  nmozsgn = 0
  obuold = 0.
  zii=1000.         ! m  (pbl height)
  beta=1.           ! -  (in computing W_*)

  taf = (tg + thm)/2.
  qaf = (clm1%forc_q+qg)/2.

  ur = max(1.0,sqrt(clm1%forc_u*clm1%forc_u+clm1%forc_v*clm1%forc_v))    ! limit must set to 1.0, otherwise
  dth = thm-taf
  dqh = clm1%forc_q-qaf
  dthv = dth*(1.+0.61*clm1%forc_q)+0.61*th*dqh
  zldis = clm1%forc_hgt_u - clm1%displa

  call clm1_obuini(ur,thv,dthv,zldis,z0mv,um,obu)

!                 +----------------------------+
!>----------------| BEGIN stability iteration: |-----------------------<
!                 +----------------------------+

  ITERATION : do while (itlef <= itmax) 

     tlbef = clm1%t_veg
     del2 = del

! Evaluate stability-dependent variables using moz from prior iteration

     call clm1_obult (clm1%displa, z0mv, z0hv, z0qv, obu, um, ustar, temp1, temp2, clm1)

! Aerodynamic resistance

     ram(1)=1./(ustar*ustar/um)
     rah(1)=1./(temp1*ustar) 
     raw(1)=1./(temp2*ustar) 

! Bulk boundary layer resistance of leaves

     uaf = um*sqrt( 1./(ram(1)*um) )
     cf = 0.01/(sqrt(uaf)*sqrt(clm1%dleaf))
     rb = 1./(cf*uaf)

! Aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.
! if no vegetation, rah(2)=0 because zpd+z0h = z0hg.
! (Dickinson et al., 1993, pp.54)

     ram(2) = 0.               ! not used
     rah(2) = 1./(clm1%csoilc*uaf)
     raw(2) = rah(2) 

! Stomatal resistances for sunlit and shaded fractions of canopy.
! should do each iteration to account for differences in eah, tv.

     svpts = el                        ! pa
     eah = clm1%forc_pbot * qaf / 0.622 ! pa
     foln = clm1%folnvt

     call clm1_stomata(mpe      , clm1%parsun, svpts     , eah       ,    &
                      thm      , o2        , co2       ,                &
                      clm1%btran, rb        , clm1%rssun , clm1%psnsun,    &
                      clm1%qe25 , clm1%kc25  , clm1%ko25  , clm1%vcmx25,    &
                      clm1%akc  , clm1%ako   , clm1%avcmx , clm1%bp    ,    &
                      clm1%mp   , foln      , clm1%folnmx, clm1%c3psn , clm1)
                                                      
     call clm1_stomata(mpe      , clm1%parsha, svpts     , eah       ,    &
                      thm      , o2        , co2       ,                &
                      clm1%btran, rb        , clm1%rssha , clm1%psnsha,    &
                      clm1%qe25 , clm1%kc25  , clm1%ko25  , clm1%vcmx25,    &
                      clm1%akc  , clm1%ako   , clm1%avcmx , clm1%bp    ,    &
                      clm1%mp   , foln      , clm1%folnmx, clm1%c3psn , clm1)

! Heat conductance for air, leaf and ground  

     call clm1_condch(rah(1),rb,rah(2),wta,wtl,wtg,wta0,wtl0,wtg0, &
          wtal,wtga,wtgl,clm1)

! Fraction of potential evaporation from leaf

     if (clm1%fdry .gt. 0.0) then
        rppdry  = clm1%fdry*rb*(clm1%laisun/(rb+clm1%rssun) + clm1%laisha/(rb+clm1%rssha))/clm1%elai
     else
        rppdry = 0.0
     endif
     efpot = clm1%forc_rho*wtl*(qsatl-qaf)

     if (efpot > 0.) then
       if (clm1%btran > 1.e-10) then

        clm1%qflx_tran_veg = efpot*rppdry
        rpp = rppdry + clm1%fwet

! No transpiration if btran below 1.e-10

       else
        rpp = clm1%fwet
        clm1%qflx_tran_veg = 0.
       endif

! Check total evapotranspiration from leaves

       rpp = min(rpp, (clm1%qflx_tran_veg+clm1%h2ocan/clm1%dtime)/efpot)

     else

! No transpiration if potential evaporation less than zero
       rpp = 1.
       clm1%qflx_tran_veg = 0.

     endif

! Update conductances for changes in rpp 
! Latent heat conductances for ground and leaf.
! Air has same conductance for both sensible and latent heat.

     call clm1_condcq(raw(1),rb,raw(2),rpp,wtaq,wtlq,wtgq,wtaq0, &
          wtlq0,wtgq0,wtalq,wtgaq,wtglq,clm1) 

! The partial derivatives of aerodynamical resistance are ignored 
! which cannot be determined analytically. 

     dc1 = clm1%forc_rho*cpair*wtl
     dc2 = hvap*clm1%forc_rho*wtlq

     efsh = dc1*(wtga*clm1%t_veg-wtg0*tg-wta0*thm)
     efe = dc2*(wtgaq*qsatl-wtgq0*qg-wtaq0*clm1%forc_q)

! Evaporation flux from foliage

     erre = 0.
     if (efe*efeb < 0.0) then
        efeold = efe
        efe  = 0.1*efeold
        erre = efe - efeold
     endif
     clm1%dt_veg = (clm1%sabv + air + bir*clm1%t_veg**4 + cir*tg**4 - efsh - efe) &
                  / (- 4.*bir*clm1%t_veg**3 +dc1*wtga +dc2*wtgaq*qsatldT) + &
                  clm1%frac_veg_nosno*(clm1%dtcanal) ! Added for PSAS Assimilation

     clm1%t_veg = tlbef + clm1%dt_veg
     dels = clm1%t_veg-tlbef
     del  = abs(dels)
     err = 0.
     if (del > delmax) then
        clm1%dt_veg = delmax*dels/del
        clm1%t_veg = tlbef + clm1%dt_veg
        err = clm1%sabv + air + bir*tlbef**3*(tlbef + 4.*clm1%dt_veg) &
             + cir*tg**4 - (efsh + dc1*wtga*clm1%dt_veg)          &
             - (efe + dc2*wtgaq*qsatldT*clm1%dt_veg)
     endif

! Fluxes from leaves to canopy space
! "efe" was limited as its sign changes frequently.  This limit may
! result in an imbalance in "hvap*qflx_evap_veg" and "efe + dc2*wtgaq*qsatldT*dt_veg" 

     efpot = clm1%forc_rho*wtl*(wtgaq*(qsatl+qsatldT*clm1%dt_veg) &
          -wtgq0*qg-wtaq0*clm1%forc_q)
     clm1%qflx_evap_veg = rpp*efpot

! Calculation of evaporative potentials (efpot) and
! interception losses; flux in kg m**-2 s-1.  ecidif 
! holds the excess energy if all intercepted water is evaporated
! during the timestep.  This energy is later added to the
! sensible heat flux.

     ecidif = 0.
     if (efpot > 0. .AND. clm1%btran > 1.e-10) then
        clm1%qflx_tran_veg = efpot*rppdry
     else
        clm1%qflx_tran_veg = 0.
     endif
     ecidif = max(0., clm1%qflx_evap_veg-clm1%qflx_tran_veg-clm1%h2ocan/clm1%dtime)
     clm1%qflx_evap_veg = min(clm1%qflx_evap_veg,clm1%qflx_tran_veg+clm1%h2ocan/clm1%dtime)

! The energy loss due to above two limits is added to 
! the sensible heat flux.

     clm1%eflx_sh_veg = efsh + dc1*wtga*clm1%dt_veg + err + erre +hvap*ecidif

! Recalculate leaf saturated vapor pressure (eg) for updated leaf temperature
! and adjust specific humidity (qsatl) proportionately 

     call clm1_qsadv(clm1%t_veg,clm1%forc_pbot,el,deldT,qsatl,qsatldT)

! Update vegetation/ground surface temperature, canopy air temperature, 
! canopy vapor pressure, aerodynamic temperature, and
! Monin-Obukhov stability parameter moz for next iteration. 

     taf = wtg0*tg + wta0*thm + wtl0*clm1%t_veg
     qaf = wtlq0*qsatl+wtgq0*qg+clm1%forc_q*wtaq0

! Update Monin-Obukhov length and wind speed including the stability effect

     dth = thm-taf       
     dqh = clm1%forc_q-qaf

     tstar=temp1*dth
     qstar=temp2*dqh

     dthv=dth*(1.+0.61*clm1%forc_q)+0.61*th*dqh

     thvstar=tstar*(1.+0.61*clm1%forc_q) + 0.61*th*qstar
     zeta=zldis*vkc*grav*thvstar/(ustar**2*thv)

     if (zeta >= 0.) then     !stable
        zeta = min(2.,max(zeta,0.01))
     else                     !unstable
        zeta = max(-100.,min(zeta,-0.01))
     endif
     obu = zldis/zeta

     if (dthv >= 0.) then
        um=max(ur,0.1)
     else
        wc=beta*(-grav*ustar*thvstar*zii/thv)**0.333
        um=sqrt(ur*ur+wc*wc)
     endif

     if (obuold*obu < 0.) nmozsgn = nmozsgn+1
     if (nmozsgn >= 4) then 
        obu = zldis/(-0.01)
     endif

     obuold = obu

! Test for convergence

     itlef = itlef+1
     if (itlef > itmin) then
        dele = abs(efe-efeb)
        efeb = efe
        det  = max(del,del2)
        if (det < dtmin .AND. dele < dlemin) exit 
     endif

! Repeat iteration

  enddo ITERATION

!                   +--------------------------+
!>------------------| END stability iteration: |-----------------------<
!                   +--------------------------+
! Balance check

  err = clm1%sabv + air + bir*tlbef**3*(tlbef + 4.*clm1%dt_veg) &
       + cir*tg**4 - clm1%eflx_sh_veg - hvap*clm1%qflx_evap_veg
  if (abs(err) > 0.1) then
     write(79,*) 'energy balance in canopy X',err
  endif

! Fluxes from ground to canopy space 

  delt  = wtal*tg-wtl0*clm1%t_veg-wta0*thm
  delq  = wtalq*qg-wtlq0*qsatl-wtaq0*clm1%forc_q
  clm1%taux  = clm1%taux - clm1%frac_veg_nosno*clm1%forc_rho*clm1%forc_u/ram(1)
  clm1%tauy  = clm1%tauy - clm1%frac_veg_nosno*clm1%forc_rho*clm1%forc_v/ram(1)
  clm1%eflx_sh_grnd = clm1%eflx_sh_grnd + cpair*clm1%forc_rho*wtg*delt
  clm1%qflx_evap_soi = clm1%qflx_evap_soi +   clm1%forc_rho*wtgq*delq

! 2 m height air temperature

  clm1%t_ref2m   = clm1%t_ref2m + clm1%frac_veg_nosno*(taf + temp1*dth * &
       1./vkc *log((2.+z0hv)/z0hv))
  if (delq > 0.) then
     clm1%sfact = clm1%sfact + clm1%forc_rho*(wtgq*wtalq)*sfacx
     clm1%sfactmax = clm1%sfactmax + clm1%forc_rho*(wtgq*wtalq)*dqgmax
  endif

! Downward longwave radiation below the canopy    

  dlrad = clm1%frac_veg_nosno*(1.-emv)*emg*clm1%forc_lwrad &
       + clm1%frac_veg_nosno * emv*emg * sb * &
       tlbef**3*(tlbef + 4.*clm1%dt_veg)

! Upward longwave radiation above the canopy    

  ulrad = clm1%frac_veg_nosno* ( (1.-emg)*(1.-emv)*(1.-emv)*clm1%forc_lwrad &
       + emv*(1.+(1.-emg)*(1.-emv))*sb * tlbef**3 &
       *(tlbef + 4.*clm1%dt_veg) + emg *(1.-emv) *sb * tg**4)

! Derivative of soil energy flux with respect to soil temperature (cgrnd) 

  cgrnds = cgrnds + cpair*clm1%forc_rho*wtg*wtal
  cgrndl = cgrndl + clm1%forc_rho*wtgq*wtalq*dqgdT
  cgrnd  = cgrnds + cgrndl*htvp

! Update dew accumulation (kg/m2) 

  clm1%h2ocan = max(0.,clm1%h2ocan + (clm1%qflx_tran_veg-clm1%qflx_evap_veg)*clm1%dtime)

end subroutine clm1_leaftem

