#include <misc.h>

subroutine clm1_twoStream (clm1, ib , ic , coszen, vai, &
                          rho, tau, fab, fre   , ftd, &
                          fti, gdir)
       
!----------------------------------------------------------------------- 
! 
! Purpose: 
! two-stream fluxes for canopy radiative transfer
! 
! Method: 
! Use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse 
! flux given an underlying surface with known albedo.
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: clm1_twostream.F90,v 1.1.1.1 2003/02/06 16:10:47 jgottsch Exp $
!-----------------------------------------------------------------------

  use precision
  use clm1type
  use clm1_varpar, only : numrad
  use clm1_varcon, only : omegas, tfrz, betads, betais
  implicit none

! ------------------------ arguments ------------------------------
  type (clm11d), intent(inout) :: clm1   !CLM 1-D Module
  integer , intent(in)  :: ib          !waveband number 
  integer , intent(in)  :: ic          !0=unit incoming direct; 1=unit incoming diffuse
  real(r8), intent(in)  :: coszen      !cosine solar zenith angle for next time step
  real(r8), intent(in)  :: vai         !elai+esai
  real(r8), intent(in)  :: rho(numrad) !leaf/stem refl weighted by fraction LAI and SAI
  real(r8), intent(in)  :: tau(numrad) !leaf/stem tran weighted by fraction LAI and SAI
  real(r8), intent(out) :: fab(numrad) !flux abs by veg layer (per unit incoming flux)   
  real(r8), intent(out) :: fre(numrad) !flux refl above veg layer (per unit incoming flx)
  real(r8), intent(out) :: ftd(numrad) !down dir flux below veg layer (per unit in flux) 
  real(r8), intent(out) :: fti(numrad) !down dif flux below veg layer (per unit in flux) 
  real(r8), intent(out) :: gdir        !aver projected leaf/stem area in solar direction
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  integer j               !array index
  integer i               !array index 
  real(r8) cosz        !0.001 <= coszen <= 1.000
  real(r8) asu         !single scattering albedo
  real(r8) chil        ! -0.4 <= xl <= 0.6
  real(r8) tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  real(r8) p1,p2,p3,p4,s1,s2,u1,u2,u3
  real(r8) b,c,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
  real(r8) phi1,phi2,sigma
  real(r8) ftds,ftis,fres
  real(r8) twostext !optical depth of direct beam per unit leaf area 
  real(r8) avmu     !average diffuse optical depth
  real(r8) omega    !fraction of intercepted radiation that is scattered
  real(r8) omegal   !omega for leaves
  real(r8) betai    !upscatter parameter for diffuse radiation 
  real(r8) betail   !betai for leaves
  real(r8) betad    !upscatter parameter for direct beam radiation 
  real(r8) betadl   !betad for leaves
  real(r8) dewmxi,vegt
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Calculate two-stream parameters omega, betad, betai, avmu, gdir, twostext.
! Omega, betad, betai are adjusted for snow. Values for omega*betad 
! and omega*betai are calculated and then divided by the new omega
! because the product omega*betai, omega*betad is used in solution. 
! Also, the transmittances and reflectances (tau, rho) are linear 
! weights of leaf and stem values.
! -----------------------------------------------------------------
  cosz = max(0.001, coszen)
  chil = min( max(clm1%xl, -0.4), 0.6 )
!  if (abs(chil) <= 0.01) chil = 0.01
  phi1 = 0.5 - 0.633*chil - 0.330*chil*chil
  phi2 = 0.877 * (1.-2.*phi1)
  gdir = phi1 + phi2*cosz
  twostext = gdir/cosz
! JDR Added bug fix from Dai
  if(phi1 == 0.) then
   avmu = 1.0/phi2
  else if(phi2 == 0.) then
   avmu = 0.5/phi1
  else
   avmu = ( 1. - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
  endif
  omegal = rho(ib) + tau(ib)
  tmp0 = gdir + phi2*cosz
  tmp1 = phi1*cosz
  if(tmp0 == 0..or.tmp1 == 0.)then
   asu = 0.5*gdir * 0.5
  else
   asu = 0.5*omegal*gdir/tmp0 * ( 1. - tmp1/tmp0 * log((tmp1+tmp0)/tmp1) )
  endif
  betadl = (1.+avmu*twostext)/(omegal*avmu*twostext)*asu
  betail = 0.5 * ((rho(ib)+tau(ib)) + (rho(ib)-tau(ib)) * ((1.+chil)/2.)**2) / omegal

! Adjust omega, betad, and betai for intercepted snow
  if (clm1%t_veg > tfrz) then                             !no snow
     tmp0 = omegal           
     tmp1 = betadl 
     tmp2 = betail  
  else if(clm1%t_veg <= tfrz) then
     if (clm1%h2ocan > 0.) then
      vegt     = clm1%frac_veg_nosno*(clm1%elai + clm1%esai)
      dewmxi   = 1.0/clm1%dewmx
      clm1%fwet = ((dewmxi/vegt)*clm1%h2ocan)**.666666666666
      clm1%fwet = min (clm1%fwet,1.0)     ! Check for maximum limit of fwet
     else
      clm1%fwet = 0.
     endif

     tmp0 =   (1.-clm1%fwet)*omegal        + clm1%fwet*omegas(ib)          
     tmp1 = ( (1.-clm1%fwet)*omegal*betadl + clm1%fwet*omegas(ib)*betads ) / tmp0
     tmp2 = ( (1.-clm1%fwet)*omegal*betail + clm1%fwet*omegas(ib)*betais ) / tmp0
  end if
  omega = tmp0           
  betad = tmp1 
  betai = tmp2  
! -----------------------------------------------------------------
! Absorbed, reflected, transmitted fluxes per unit incoming radiation
! -----------------------------------------------------------------
  b = 1. - omega + omega*betai
  c = omega*betai
  tmp0 = avmu*twostext
  d = tmp0 * omega*betad
  f = tmp0 * omega*(1.-betad)
  tmp1 = b*b - c*c
  sigma = tmp0*tmp0 - tmp1
!JDR Bug fix
  if(abs(sigma) <= 1.e-6)then
   omega = 0.98*omega
   betad = 0.98*betad
   c = 0.98*c
   b = 1. - omega + c
   tmp1 = b*b - c*c
   sigma = tmp0*tmp0 - tmp1
  endif
  h = sqrt(tmp1) / avmu
  p1 = b + avmu*h
  p2 = b - avmu*h
  p3 = b + tmp0
  p4 = b - tmp0
  s1 = exp(-h*vai)
  s2 = exp(-twostext*vai)
  if (ic == 0) then
     u1 = b - c/clm1%albgrd(ib)
     u2 = b - c*clm1%albgrd(ib)
     u3 = f + c*clm1%albgrd(ib)
  else
     u1 = b - c/clm1%albgri(ib)
     u2 = b - c*clm1%albgri(ib)
     u3 = f + c*clm1%albgri(ib)
  end if
  tmp2 = u1 - avmu*h
  tmp3 = u1 + avmu*h
  d1 = p1*tmp2/s1 - p2*tmp3*s1
  tmp4 = u2 + avmu*h
  tmp5 = u2 - avmu*h
  d2 = tmp4/s1 - tmp5*s1
  h1 = -d*p4 - c*f
  tmp6 = d - h1*p3/sigma
  tmp7 = ( d - c - h1/sigma*(u1+tmp0) ) * s2
  h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
  h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
  h4 = -f*p3 - c*d
  tmp8 = h4/sigma
  tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
  h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
  h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
  h7 = (c*tmp2) / (d1*s1)
  h8 = (-c*tmp3*s1) / d1
  h9 = tmp4 / (d2*s1)
  h10 = (-tmp5*s1) / d2

! Downward direct and diffuse fluxes below vegetation

  if (ic == 0) then
     ftds = s2
     ftis = h4*s2/sigma + h5*s1 + h6/s1
  else
     ftds = 0.
     ftis = h9*s1 + h10/s1
  end if
  ftd(ib) = ftds
  fti(ib) = ftis

! Flux reflected by vegetation

  if (ic == 0) then
     fres = h1/sigma + h2 + h3
  else
     fres = h7 + h8
  end if
  fre(ib) = fres
  
! Flux absorbed by vegetation

  fab(ib) = 1. - fre(ib) - (1.-clm1%albgrd(ib))*ftd(ib) - (1.-clm1%albgri(ib))*fti(ib)
  return
end subroutine clm1_twostream
