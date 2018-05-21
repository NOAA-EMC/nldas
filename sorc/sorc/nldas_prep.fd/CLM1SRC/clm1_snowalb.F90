#include <misc.h>

subroutine clm1_snowalb (clm1, coszen, nband, ind, alb)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Determine snow albedos
! 
!-----------------------------------------------------------------------
! $Id: clm1_snowalb.F90,v 1.1.1.1 2003/02/06 16:10:44 jgottsch Exp $
!-----------------------------------------------------------------------

  use precision
  use clm1type
  implicit none

! ------------------------ arguments ------------------------------ 
  type (clm11d), intent(inout) :: clm1       !CLM 1-D Module
  real(r8)    , intent(in) :: coszen       !cosine solar zenith angle for next time step
  integer     , intent(in) :: nband        !number of solar radiation waveband classes
  integer     , intent(in) :: ind          !0=direct beam, 1=diffuse radiation
  real(r8)    , intent(out):: alb(numrad)  !snow albedo by waveband
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------

  integer  :: ib          !waveband class

  real(r8) :: snal0 = 0.95 !vis albedo of new snow for sza<60
  real(r8) :: snal1 = 0.65 !nir albedo of new snow for sza<60
  real(r8) :: conn  = 0.5  !constant for visible snow alb calculation [-]
  real(r8) :: cons  = 0.2  !constant (=0.2) for nir snow albedo calculation [-]
  real(r8) :: sl    = 2.0  !factor that helps control alb zenith dependence [-]

  real(r8) :: age          !factor to reduce visible snow alb due to snow age [-]
  real(r8) :: albs         !temporary vis snow albedo
  real(r8) :: albl         !temporary nir snow albedo
  real(r8) :: cff          !snow alb correction factor for zenith angle > 60 [-]
  real(r8) :: czf          !solar zenith correction for new snow albedo [-]

! -----------------------------------------------------------------

! zero albedos

  do ib = 1, nband
     alb(ib) = 0._r8
  end do

! =========================================================================
! CLM Albedo for snow cover.
! Snow albedo depends on snow-age, zenith angle, and thickness of snow age
! gives reduction of visible radiation
! =========================================================================

! Correction for snow age

  age = 1.-1./(1.+clm1%snowage)
  albs = snal0*(1.-cons*age)
  albl = snal1*(1.-conn*age)

  if (ind == 0) then

! Czf corrects albedo of new snow for solar zenith

    cff    = ((1.+1./sl)/(1.+max(0.001,coszen)*2.*sl )- 1./sl)
    cff    = max(cff,0.)
    czf    = 0.4*cff*(1.-albs)
    albs = albs+czf
    czf    = 0.4*cff*(1.-albl)
    albl = albl+czf

  endif

  alb(1) = albs
  alb(2) = albl

  return
end subroutine clm1_snowalb


