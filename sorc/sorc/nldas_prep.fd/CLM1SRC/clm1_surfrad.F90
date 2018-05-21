#include <misc.h>

subroutine clm1_surfrad (clm1)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Solar fluxes absorbed by vegetation and ground surface
! 
! Method: 
! Note possible problem when land is on different grid than atmosphere.
!
! Land may have sun above the horizon (coszen > 0) but atmosphere may
! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
! because all fluxes (absorbed, reflected, transmitted) are multiplied
! by the incoming flux and all will equal zero.
!
! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0) but
! land may have sun below horizon. This is okay because fabd, fabi,
! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words, all
! the radiation is reflected. NDVI should equal zero in this case.
! However, the way the code is currently implemented this is only true
! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
!
! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: clm1_surfrad.F90,v 1.1.1.1 2003/02/06 16:10:47 jgottsch Exp $
!-----------------------------------------------------------------------

  use precision
  use clm1type
  implicit none

! ------------------------ arguments-------------------------------
  type (clm11d), intent(inout) :: clm1    !CLM 1-D Module
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  integer ib              !waveband number (1=vis, 2=nir)
  real(r8) abs            !absorbed solar radiation (W/m**2) 
  real(r8) rnir           !reflected solar radiation [nir] (W/m**2)
  real(r8) rvis           !reflected solar radiation [vis] (W/m**2)
  real(r8) laifra         !leaf area fraction of canopy
  real(r8) trd            !transmitted solar radiation: direct (W/m**2)
  real(r8) tri            !transmitted solar radiation: diffuse (W/m**2)
  real(r8) cad(numrad)    !direct beam absorbed by canopy (W/m**2)
  real(r8) cai(numrad)    !diffuse radiation absorbed by canopy (W/m**2)
  integer  nband          !number of solar radiation waveband classes              
  real(r8) fsha           !shaded fraction of canopy
  real(r8) vai            !total leaf area index + stem area index, one sided
  real(r8) mpe            !prevents overflow for division by zero                  
! -----------------------------------------------------------------

  mpe   = 1.e-06
  nband = numrad

  fsha = 1.-clm1%fsun
  clm1%laisun = clm1%elai*clm1%fsun
  clm1%laisha = clm1%elai*fsha
  vai = clm1%elai+ clm1%esai

! zero summed solar fluxes

  clm1%sabg = 0.
  clm1%sabv = 0.
  clm1%fsa  = 0.

! loop over nband wavebands

  do ib = 1, nband

! absorbed by canopy

     cad(ib)  = clm1%forc_solad(ib)*clm1%fabd(ib)
     cai(ib)  = clm1%forc_solai(ib)*clm1%fabi(ib)
     clm1%sabv = clm1%sabv + cad(ib) + cai(ib)
     clm1%fsa  = clm1%fsa  + cad(ib) + cai(ib)

! transmitted = solar fluxes incident on ground

     trd = clm1%forc_solad(ib)*clm1%ftdd(ib)
     tri = clm1%forc_solad(ib)*clm1%ftid(ib) + clm1%forc_solai(ib)*clm1%ftii(ib)

! solar radiation absorbed by ground surface

     abs = trd*(1.-clm1%albgrd(ib)) + tri*(1.-clm1%albgri(ib)) 
     clm1%sabg = clm1%sabg + abs
     clm1%fsa  = clm1%fsa  + abs

  end do

! partion visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves

  laifra = clm1%elai / max(vai,mpe)
  if (clm1%fsun > 0.) then
     clm1%parsun = (cad(1)+cai(1)) * laifra
     clm1%parsha = 0._r8
  else
     clm1%parsun = 0._r8 
     clm1%parsha = 0._r8
  endif

! NDVI and reflected solar radiation

  rvis = clm1%albd(1)*clm1%forc_solad(1) + clm1%albi(1)*clm1%forc_solai(1) 
  rnir = clm1%albd(2)*clm1%forc_solad(2) + clm1%albi(2)*clm1%forc_solai(2)
  clm1%fsr = rvis + rnir
  clm1%ndvi = (rnir-rvis) / max(rnir+rvis,mpe)

  return
end subroutine clm1_surfrad

