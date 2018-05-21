#include <misc.h>

subroutine clm1_thermal (clm1)

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
!  This is the main subroutine to execute the calculation of thermal 
!  processes and surface fluxes.
!  (1) Leaf temperature
!      Foliage energy conservation is given by the foliage energy budget 
!      equation:
!                     Rnet - Hf - LEf = 0 
!      The equation is solved by Newton-Raphson iteration, in which this 
!      iteration includes the calculation of the photosynthesis and 
!      stomatal resistance, and the integration of turbulent flux profiles. 
!      The sensible and latent heat transfer between foliage and atmosphere 
!      and ground is linked by the equations:  
!                     Ha = Hf + Hg and Ea = Ef + Eg
!
!  (2) Snow and soil temperatures
!      o The volumetric heat capacity is calculated as a linear combination 
!        in terms of the volumetric fraction of the constituent phases. 
!      o The thermal conductivity of soil is computed from 
!        the algorithm of Johansen (as reported by Farouki 1981), and the 
!        conductivity of snow is from the formulation used in
!        SNTHERM (Jordan 1991).
!      o Boundary conditions:  
!        F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
!      o Soil / snow temperature is predicted from heat conduction 
!        in 10 soil layers and up to 5 snow layers. 
!        The thermal conductivities at the interfaces between two 
!        neighboring layers (j, j+1) are derived from an assumption that 
!        the flux across the interface is equal to that from the node j 
!        to the interface and the flux from the interface to the node j+1. 
!        The equation is solved using the Crank-Nicholson method and 
!        results in a tridiagonal system equation.
!
!  (3) Phase change (see clm1_meltfreeze.F90)
!
!  FLOW DIAGRAM FOR clm1_thermal.F90
!
!  thermal ===> clm1_qsadv
!               clm1_obuini
!               clm1_obult
!               clm1_leaftem  
!                  ===> clm1_qsadv    
!                       clm1_obuini  
!                       clm1_obult   
!                       clm1_stomata 
!                       clm1_condch  
!                       clm1_condcq  
!               clm1_thermalk          
!               clm1_tridia            
!               clm1_meltfreeze        
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!  20 November 2002:  Jon Radakovich; Added term for assimilation and BC
!=========================================================================
! $Id: clm1_thermal.F90,v 1.1.1.1 2003/02/06 16:10:46 jgottsch Exp $
!=========================================================================

  use precision
  use clm1type
  use clm1_varcon, only : denh2o, denice, roverg, hvap, hsub, &
                        rair, cpair, grav, vkc, tfrz, sb, istice, istwet 
  use clm1_varpar, only : nlevsoi
  implicit none

!=== Arguments  =====================================================

  type (clm11d), intent(inout)  :: clm1	 !CLM 1-D Module

!=== Local Variables =====================================================

  integer i,j

  real(r8)  &
       at(clm1%snl+1 : nlevsoi),    & ! "a" vector for tridiagonal matrix
       bt(clm1%snl+1 : nlevsoi),    & ! "b" vector for tridiagonal matrix
       ct(clm1%snl+1 : nlevsoi),    & ! "c" vector for tridiagonal matrix
       rt(clm1%snl+1 : nlevsoi),    & ! "r" vector for tridiagonal solution
       tg,                         & ! ground surface temperature [K]
       cv(clm1%snl+1 : nlevsoi),    & ! heat capacity [J/(m2 K)]
       tk(clm1%snl+1 : nlevsoi),    & ! thermal conductivity [W/(m K)]
       tssbef(clm1%snl+1 : nlevsoi),& ! soil/snow temperature before update
       qred,                       & ! soil surface relative humidity
       z0mg,                       & ! roughness length over ground, momentum [m]
       z0hg,                       & ! roughness length over ground, sensible heat [m]
       z0qg,                       & ! roughness length over ground, latent heat [m]
       z0mv,                       & ! roughness length over vegetation, momentum [m]
       z0hv,                       & ! roughness length over vegetation, sensible heat [m]
       z0qv                          ! roughness length over vegetation, latent heat [m]

  real(r8)  htvp,                  & ! latent heat of vapor of water (or sublimation) [j/kg]
       fact(clm1%snl+1 : nlevsoi),  & ! used in computing tridiagonal matrix
       fn  (clm1%snl+1 : nlevsoi),  & ! heat diffusion through the layer interface [W/m2]
       fn1 (clm1%snl+1 : nlevsoi),  & ! heat diffusion through the layer interface [W/m2]
       dzm,                        & ! used in computing tridiagonal matrix
       dzp                           ! used in computing tridiagonal matrix

  integer &
       niters,                     & ! maximum number of iterations for surface temperature
       iter,                       & ! iteration index
       nmozsgn                       ! number of times moz changes sign

  real(r8)  beta,                     & ! coefficient of conective velocity [-]
       zii,                        & ! convective boundary height [m]
       zldis,                      & ! reference height "minus" zero displacement heght [m]
       ur,                         & ! wind speed at reference height [m/s]
       th,                         & ! potential temperature (kelvin)
       thm,                        & ! intermediate variable (forc_t+0.0098*forc_hgt_t) 
       thv,                        & ! virtual potential temperature (kelvin)
       dth,                        & ! diff of virtual temp. between ref. height and surface
       dqh,                        & ! diff of humidity between ref. height and surface
       dthv,                       & ! diff of vir. poten. temp. between ref. height and surface
       thvstar,                    & ! virtual potential temperature scaling parameter
       obu,                        & ! monin-obukhov length (m)
       zeta,                       & ! dimensionless height used in Monin-Obukhov theory
       wc,                         & ! convective velocity [m/s]
       um,                         & ! wind speed including the stability effect [m/s]
       temp1,                      & ! relation for potential temperature profile
       temp2,                      & ! relation for specific humidity profile
       ustar,                      & ! friction velocity [m/s]
       tstar,                      & ! temperature scaling parameter
       qstar,                      & ! moisture scaling parameter
       ram,                        & ! aerodynamical resistance [s/m]
       rah,                        & ! thermal resistance [s/m]
       raw,                        & ! moisture resistance [s/m]
       raih,                       & ! temporary variable [kg/m2/s]
       raiw,                       & ! temporary variable [kg/m2/s]
       emg,                        & ! ground emissivity (0.97 for snow, glaciers and water surface; 0.96 for soil and wetland)
       emv,                        & ! vegetation emissivity
       avmuir                        ! ir inverse optical depth per unit leaf area

  real(r8)                         & 
       cgrnd,                      & ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
       cgrndl,                     & ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
       cgrnds,                     & ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
       hs,                         & ! net energy flux into the surface (w/m2)
       dhsdt,                      & ! d(hs)/dT
       eg,                         & ! water vapor pressure at temperature T [pa]
       qsatg,                      & ! saturated humidity [kg/kg]
       degdT,                      & ! d(eg)/dT
       qsatgdT,                    & ! d(qsatg)/dT
       fac,                        & ! soil wetness of surface layer
       psit,                       & ! negative potential of soil
       hr,                         & ! relative humidity
       dqgmax,                     & ! maximum of d(qg)/d(theta)
       sfacx,                      & ! coefficient for "sfact"
       qg,                         & ! ground specific humidity [kg/kg]
       dqgdT,                      & ! d(qg)/dT
       wice0(clm1%snl+1 : nlevsoi), & ! ice mass from previous time-step
       wliq0(clm1%snl+1 : nlevsoi), & ! liquid mass from previous time-step
       wx,                         & ! patitial volume of ice and water of surface layer
       egsmax,                     & ! max. evaporation which soil can provide at one time step
       egidif,                     & ! the excess of evaporation over "egsmax"
       brr(clm1%snl+1 : nlevsoi),   & ! temporay set 
       xmf,                        & ! total latent heat of phase change of ground water
       dlrad,                      & ! downward longwave radiation blow the canopy [W/m2]
       ulrad,                      & ! upward longwave radiation above the canopy [W/m2]
       tinc,                       & ! temperature difference of two time step
       obuold                        ! monin-obukhov length from previous iteration

  real(r8) temp                      !temporary variable                                      
  real(r8) cf                        !s m**2/umol -> s/m

  real(r8) vol_ice(1:nlevsoi)        ! partial volume of ice lens in layer
  real(r8) eff_porosity(1:nlevsoi)   ! effective porosity in layer
  real(r8) vol_liq(1:nlevsoi)        ! partial volume of liquid water in layer
  real(r8) s_node                    ! vol_liq/eff_porosity
  real(r8) smp_node                  ! matrix potential
  real(r8) rresis(1:nlevsoi)         ! soil water contribution to root resistance


!=== End Variable List ===================================================

!=========================================================================
! [1] Initial set 
!=========================================================================

! Fluxes 

  clm1%taux     = 0.
  clm1%tauy     = 0.
  clm1%eflx_sh_tot    = 0.  
  clm1%qflx_evap_tot    = 0.  
  clm1%eflx_lh_tot   = 0.  
  clm1%eflx_sh_veg    = 0.  
  clm1%qflx_evap_veg    = 0.  
  clm1%qflx_tran_veg      = 0.  
  clm1%eflx_sh_grnd    = 0.
  clm1%qflx_evap_soi    = 0.  
  dlrad    = 0.
  ulrad    = 0.
  cgrnds   = 0.
  cgrndl   = 0.
  cgrnd    = 0.
  clm1%sfact    = 0.
  clm1%sfactmax = 0.
  clm1%t_ref2m     = 0.

! Effective porosity of soil, partial volume of ice and liquid

  do i = 1,nlevsoi
     vol_ice(i) = min(clm1%watsat(i), clm1%h2osoi_ice(i)/(clm1%dz(i)*denice))
     eff_porosity(i) = clm1%watsat(i)-vol_ice(i)
     vol_liq(i) = min(eff_porosity(i), clm1%h2osoi_liq(i)/(clm1%dz(i)*denh2o))
  enddo

!Temperature and water mass from previous time step

  tg = clm1%t_soisno(clm1%snl+1)
  do i = clm1%snl+1, nlevsoi
     tssbef(i) = clm1%t_soisno(i)
     wice0(i) = clm1%h2osoi_ice(i)
     wliq0(i) = clm1%h2osoi_liq(i)
  enddo


!=========================================================================
! [2] Specific humidity and its derivative at ground surface
!=========================================================================

  qred = 1.
  if (clm1%itypwat/=istwet .AND. clm1%itypwat/=istice) then ! NOT wetland and ice land
     wx   = (clm1%h2osoi_liq(1)/denh2o+clm1%h2osoi_ice(1)/denice)/clm1%dz(1)
     fac  = min(1., wx/clm1%watsat(1))
     fac  = max( fac, 0.01 )
     psit = -clm1%sucsat(1) * fac ** (- clm1%bsw(1))
     psit = max(clm1%smpmin, psit)
     hr   = exp(psit/roverg/tg)
     qred = (1.-clm1%frac_sno)*hr + clm1%frac_sno
  endif

  call clm1_qsadv(tg,clm1%forc_pbot,eg,degdT,qsatg,qsatgdT)

  qg = qred*qsatg  
  dqgdT = qred*qsatgdT

  sfacx = 0.
  dqgmax = 0.
  if (clm1%itypwat/=istwet .AND. clm1%itypwat/=istice) then ! NOT wetland and ice land
     sfacx = (1.-clm1%frac_sno)*hr*qsatg*clm1%bsw(1)/(roverg*tg)
     dqgmax = (1.-qred)/clm1%watsat(1) * qsatg
  endif

  if (qsatg > clm1%forc_q .AND. clm1%forc_q > qred*qsatg) then
     qg = clm1%forc_q
     dqgdT = 0.
     sfacx = 0.
     dqgmax = 0.
  endif


!=========================================================================
! [3] Leaf and ground surface temperature and fluxes
!=========================================================================

! 3.1 Propositional variables

! Emissivity

  if (clm1%h2osno>0. .OR.clm1%itypwat==istice) then
     emg = 0.97
  else
     emg = 0.96
  endif
  avmuir=1.
  emv=1.-exp(-(clm1%elai+clm1%esai)/avmuir)

! Latent heat, we arbitrarily assume that the sublimation occurs 
! only as h2osoi_liq = 0

  htvp = hvap
  if (clm1%h2osoi_liq(clm1%snl+1) <= 0. .AND. clm1%h2osoi_ice(clm1%snl+1) > 0.) htvp = hsub

#if (defined PERGRO)
  htvp = hvap
#endif

! Roughness length

  if (clm1%frac_sno > 0.) then
     z0mg = clm1%zsno
     z0hg = z0mg            ! initial set
     z0qg = z0mg            ! initial set
  else
     z0mg = clm1%zlnd
     z0hg = z0mg
     z0qg = z0mg
  endif

  z0mv = clm1%z0m
  z0hv = z0mv
  z0qv = z0mv

! Potential temperature at the reference height

  beta=1.        ! -  (in computing W_*)
  zii = 1000.    ! m  (pbl height)
  thm = clm1%forc_t + 0.0098*clm1%forc_hgt_t              
  th = clm1%forc_t*(100000./clm1%forc_pbot)**(rair/cpair)  ! potential T  (forc_t*(forc_psrf/forc_pbot)**(rair/cp))
  thv = th*(1.+0.61*clm1%forc_q)                          ! virtual potential T
  ur = max(1.0,sqrt(clm1%forc_u*clm1%forc_u+clm1%forc_v*clm1%forc_v)) ! limit must set to 1.0, otherwise,

! 3.2 BARE PART
! Ground fluxes and temperatures
! NOTE: in the current scheme clm1%frac_veg_nosno is EITHER 1 or 0

! Compute sensible and latent fluxes and their derivatives with repect 
! to ground temperature using ground temperatures from previous time step.

  if (clm1%frac_veg_nosno == 0) then  

! Initialization variables

     nmozsgn = 0
     obuold = 0.
     dth   = thm-tg
     dqh   = clm1%forc_q-qg
     dthv  = dth*(1.+0.61*clm1%forc_q)+0.61*th*dqh
     zldis = clm1%forc_hgt_u-0.
     call clm1_obuini(ur,thv,dthv,zldis,z0mg,um,obu)

! Evaluated stability-dependent variables using moz from prior iteration

     niters=3
     do iter = 1, niters         ! begin stability iteration
        call clm1_obult(0.,z0mg,z0hg,z0qg,obu,um,ustar,temp1,temp2,clm1)
        tstar = temp1*dth
        qstar = temp2*dqh
        z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
        z0qg = z0hg

        thvstar=tstar*(1.+0.61*clm1%forc_q) + 0.61*th*qstar
        zeta=zldis*vkc*grav*thvstar/(ustar**2*thv)
        if (zeta >= 0.) then     !stable
           zeta = min(2.,max(zeta,0.01))
        else                     !unstable
           zeta = max(-100.,min(zeta,-0.01))
        endif

        obu = zldis/zeta

        if (dthv >= 0.) then
           um = max(ur,0.1)
        else
           wc = beta*(-grav*ustar*thvstar*zii/thv)**0.333
           um = sqrt(ur*ur+wc*wc)
        endif

        if (obuold*obu < 0.) nmozsgn = nmozsgn+1
        if (nmozsgn >= 4) EXIT

        obuold = obu
     enddo                       ! end stability iteration

! Get derivative of fluxes with repect to ground temperature

     clm1%acond = ustar*ustar/um ! Add-in for ALMA output

     ram    = 1./(ustar*ustar/um)
     rah    = 1./(temp1*ustar) 
     raw    = 1./(temp2*ustar) 
     raih   = (1-clm1%frac_veg_nosno)*clm1%forc_rho*cpair/rah
     raiw   = (1-clm1%frac_veg_nosno)*clm1%forc_rho/raw          
     cgrnds = raih
     cgrndl = raiw*dqgdT
     cgrnd  = cgrnds + htvp*cgrndl
     clm1%sfact  = raiw*sfacx
     if (dqh >= 0.) clm1%sfact = 0.
     clm1%sfactmax = raiw*dqgmax

! Surface fluxes of momentum, sensible and latent heat
! using ground temperatures from previous time step

     clm1%taux   = -(1-clm1%frac_veg_nosno)*clm1%forc_rho*clm1%forc_u/ram        
     clm1%tauy   = -(1-clm1%frac_veg_nosno)*clm1%forc_rho*clm1%forc_v/ram
     clm1%eflx_sh_grnd  = -raih*dth
     clm1%qflx_evap_soi  = -raiw*dqh 
     clm1%eflx_sh_tot  = clm1%eflx_sh_grnd
     clm1%qflx_evap_tot  = clm1%qflx_evap_soi

! 2 m height air temperature

     clm1%t_ref2m=(1-clm1%frac_veg_nosno)*(tg+temp1*dth * 1./vkc *log((2.+z0hg)/z0hg))

! Equate canopy temperature to air over bareland.
! Needed as frac_veg_nosno=0 carried over to next time step

     clm1%t_veg = clm1%forc_t

     clm1%btran = 0.     !needed for history file for bare soil
     clm1%rootr(:) = 0.

     cf = clm1%forc_pbot/(8.314*thm)*1.e06 
     clm1%rssun = 1./clm1%bp * cf
     clm1%rssha = 1./clm1%bp * cf

! 3.3 VEGETATED PART
! Calculate canopy temperature, latent and sensible fluxes from the canopy,
! and leaf water change by evapotranspiration 

  else    

! Potential and root resistance factors

     clm1%btran = 1.e-10
     do i = 1,nlevsoi
        if (clm1%t_soisno(i) > tfrz) then
           s_node = max(vol_liq(i)/eff_porosity(i),0.01)
           smp_node = max(clm1%smpmax, -clm1%sucsat(i)*s_node**(-clm1%bsw(i)))
           rresis(i) = (1.-smp_node/clm1%smpmax)/(1.+clm1%sucsat(i)/clm1%smpmax)
           clm1%rootr(i) = clm1%rootfr(i)*rresis(i)
           clm1%btran = clm1%btran + clm1%rootr(i)
        else
           clm1%rootr(i) = 0.
        endif
     enddo

! Normalize root resistances to get layer contribution to ET

     do i = 1,nlevsoi
        clm1%rootr(i)  = clm1%rootr(i)/clm1%btran
     enddo

     call clm1_leaftem(z0mv,z0hv,z0qv,thm,th,thv,tg,qg,dqgdT,htvp,sfacx,     &
          dqgmax,emv,emg,dlrad,ulrad,cgrnds,cgrndl,cgrnd,clm1)

  endif

!=========================================================================
! [4] Ground temperature
!=======================================================================

! 4.1 Thermal conductivity and Heat capacity

  call clm1_thermalk(tk,cv,clm1)

  j       = clm1%snl+1
  fact(j) = clm1%dtime / cv(j) &
       * clm1%dz(j) / (0.5*(clm1%z(j)-clm1%zi(j-1)+clm1%capr*(clm1%z(j+1)-clm1%zi(j-1))))

  do j = clm1%snl+1 + 1, nlevsoi
     fact(j) = clm1%dtime/cv(j)
  enddo

! 4.2 Net ground heat flux into the surface and its temperature derivative
  hs    = clm1%sabg + dlrad &
       + (1-clm1%frac_veg_nosno)*emg*clm1%forc_lwrad - emg*sb*tg**4 &
       - (clm1%eflx_sh_grnd+clm1%qflx_evap_soi*htvp) &
       + (1-clm1%frac_veg_nosno)*(clm1%dtcanal*(clm1%dtime/fact(clm1%snl+1)))/clm1%dtime    ! Add-in for PSAS assimilation     
  dhsdT = - cgrnd - 4.*emg * sb * tg**3


  do j = clm1%snl+1, nlevsoi - 1
     fn(j) = tk(j)*(clm1%t_soisno(j+1)-clm1%t_soisno(j))/(clm1%z(j+1)-clm1%z(j))
  enddo
  fn(nlevsoi) = 0.

! Add in for ALMA output
  if(clm1%snl < 0)then
     clm1%diffusion=fn(0)
  else
     clm1%diffusion=0.
  endif

! 4.3 Set up vector r and vectors a, b, c that define tridiagonal matrix

  j     = clm1%snl+1
  dzp   = clm1%z(j+1)-clm1%z(j)
  at(j) = 0.
  bt(j) = 1+(1.-clm1%cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
  ct(j) =  -(1.-clm1%cnfac)*fact(j)*tk(j)/dzp
  rt(j) = clm1%t_soisno(j) +  fact(j)*( hs - dhsdT*clm1%t_soisno(j) + clm1%cnfac*fn(j) )

  do j    = clm1%snl+1 + 1, nlevsoi - 1
     dzm   = (clm1%z(j)-clm1%z(j-1))
     dzp   = (clm1%z(j+1)-clm1%z(j))

     at(j) =   - (1.-clm1%cnfac)*fact(j)* tk(j-1)/dzm
     bt(j) = 1.+ (1.-clm1%cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
     ct(j) =   - (1.-clm1%cnfac)*fact(j)* tk(j)/dzp

     rt(j) = clm1%t_soisno(j) + clm1%cnfac*fact(j)*( fn(j) - fn(j-1) )
  enddo

  j     =  nlevsoi
  dzm   = (clm1%z(j)-clm1%z(j-1))
  at(j) =   - (1.-clm1%cnfac)*fact(j)*tk(j-1)/dzm
  bt(j) = 1.+ (1.-clm1%cnfac)*fact(j)*tk(j-1)/dzm
  ct(j) = 0.
  rt(j) = clm1%t_soisno(j) - clm1%cnfac*fact(j)*fn(j-1)


! 4.4 Solve for t_soisno

  i = size(at)
  call clm1_tridia (i ,at ,bt ,ct ,rt ,clm1%t_soisno(clm1%snl+1:nlevsoi))

!=========================================================================
! [5] Melting or Freezing 
!=========================================================================
  do j = clm1%snl+1, nlevsoi - 1
     fn1(j) = tk(j)*(clm1%t_soisno(j+1)-clm1%t_soisno(j))/(clm1%z(j+1)-clm1%z(j))
  enddo
  fn1(nlevsoi) = 0.

  j = clm1%snl+1
  brr(j) = clm1%cnfac*fn(j) + (1.-clm1%cnfac)*fn1(j)

  do j = clm1%snl+1 + 1, nlevsoi
     brr(j) = clm1%cnfac*(fn(j)-fn(j-1)) + (1.-clm1%cnfac)*(fn1(j)-fn1(j-1))
  enddo

  call clm1_meltfreeze (fact(clm1%snl+1), brr(clm1%snl+1), hs, dhsdT, &
       tssbef(clm1%snl+1),xmf, clm1)

  tg = clm1%t_soisno(clm1%snl+1)


!=========================================================================
! [6] Correct fluxes to present soil temperature
!========================================================================= 

  tinc = clm1%t_soisno(clm1%snl+1) - tssbef(clm1%snl+1)
  clm1%eflx_sh_grnd =  clm1%eflx_sh_grnd + tinc*cgrnds 
  clm1%qflx_evap_soi =  clm1%qflx_evap_soi + tinc*cgrndl

! Calculation of evaporative potential; flux in kg m**-2 s-1.  
! egidif holds the excess energy if all water is evaporated
! during the timestep.  This energy is later added to the
! sensible heat flux.

  egsmax = (clm1%h2osoi_ice(clm1%snl+1)+clm1%h2osoi_liq(clm1%snl+1)) / clm1%dtime

  egidif = max( 0., clm1%qflx_evap_soi - egsmax )
  clm1%qflx_evap_soi = min ( clm1%qflx_evap_soi, egsmax )
  clm1%eflx_sh_grnd = clm1%eflx_sh_grnd + htvp*egidif

! Ground heat flux

  clm1%eflx_soil_grnd = clm1%sabg + dlrad + (1-clm1%frac_veg_nosno)*emg*clm1%forc_lwrad &
       - emg*sb*tssbef(clm1%snl+1)**3*(tssbef(clm1%snl+1) + 4.*tinc) &
       - (clm1%eflx_sh_grnd+clm1%qflx_evap_soi*htvp)

  clm1%eflx_sh_tot = clm1%eflx_sh_veg + clm1%eflx_sh_grnd
  clm1%qflx_evap_tot = clm1%qflx_evap_veg + clm1%qflx_evap_soi
  clm1%eflx_lh_tot= hvap*clm1%qflx_evap_veg + htvp*clm1%qflx_evap_soi   ! W/m2 (accouting for sublimation)

  clm1%qflx_evap_grnd = 0.
  clm1%qflx_sub_snow = 0.
  clm1%qflx_dew_snow = 0.
  clm1%qflx_dew_grnd = 0.

  if (clm1%qflx_evap_soi >= 0.) then
     ! Do not allow for sublimation in melting (melting ==> evap. ==> sublimation)
     clm1%qflx_evap_grnd = min(clm1%h2osoi_liq(clm1%snl+1)/clm1%dtime, clm1%qflx_evap_soi)
     clm1%qflx_sub_snow = clm1%qflx_evap_soi - clm1%qflx_evap_grnd
  else
     if (tg < tfrz) then
        clm1%qflx_dew_snow = abs(clm1%qflx_evap_soi)
     else
        clm1%qflx_dew_grnd = abs(clm1%qflx_evap_soi)
     endif
  endif

! Outgoing long-wave radiation from canopy + ground

  clm1%eflx_lwrad_out = ulrad &
       + (1-clm1%frac_veg_nosno)*(1.-emg)*clm1%forc_lwrad &
       + (1-clm1%frac_veg_nosno)*emg*sb * tssbef(clm1%snl+1)**4 &
       ! For conservation we put the increase of ground longwave to outgoing
  + 4.*emg*sb*tssbef(clm1%snl+1)**3*tinc

! Radiative temperature

  clm1%t_rad = (clm1%eflx_lwrad_out/sb)**0.25

!=========================================================================
![7] Soil Energy balance check
!=========================================================================

  clm1%errsoi = 0. 
  do j = clm1%snl+1, nlevsoi
     clm1%errsoi = clm1%errsoi - (clm1%t_soisno(j)-tssbef(j))/fact(j) 
  enddo
  clm1%errsoi = clm1%errsoi + clm1%eflx_soil_grnd - xmf

!=========================================================================
![8] Variables needed by history tap
!=========================================================================

 clm1%dt_grnd        = tinc
 clm1%eflx_lh_vege   = (clm1%qflx_evap_veg - clm1%qflx_tran_veg) * hvap
 clm1%eflx_lh_vegt   = clm1%qflx_tran_veg * hvap       
 clm1%eflx_lh_grnd   = clm1%qflx_evap_soi * htvp
 clm1%eflx_lwrad_net = clm1%eflx_lwrad_out -  clm1%forc_lwrad  

end subroutine clm1_thermal











