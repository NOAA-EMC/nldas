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
#include "misc.h"

subroutine SurfaceAlbedo (clm, caldayp1, eccen, obliqr, lambm0, mvelpp)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                        M  available land surface process model.
!  M --COMMON LAND MODEL--  C
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List:
!
!-----------------------------------------------------------------------
! Purpose: 
! Surface albedo and two-stream fluxes
! 
! Method: 
! Surface albedos. Also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation. 
! Also sunlit fraction of the canopy. 
!
! The calling sequence is:
!   -> SurfaceAlbedo:     albedos for next time step
!        -> SnowAge:      snow age
!        -> SnowAlbedo:   snow albedos: direct beam
!        -> SnowAlbedo:   snow albedos: diffuse
!        -> SoilAlbedo:   soil/lake/glacier/wetland albedos
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dir)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dif)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (nir dir)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (nir dif)
! 
! Author: 
! Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: SurfaceAlbedo.F90,v 1.7 2004/11/24 22:56:45 jim Exp $
!-----------------------------------------------------------------------

  use precision
  use clmtype
  use shr_orb_mod
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm  !CLM 1-D Module

  real(r8), intent(in) :: caldayp1 ! calendar day at Greenwich (1.00, ..., 365.99)
  real(r8), intent(in) :: eccen    ! Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr   ! Earth's obliquity in radians
  real(r8), intent(in) :: lambm0   ! Mean longitude of perihelion at the vernal equinox
                                   ! (radians)
  real(r8), intent(in) :: mvelpp   ! Earth's moving vernal equinox long. of perihelion
                                   ! + pi (radians)

!----Local Variables----------------------------------------------------

  integer  :: ib              ! band index
  integer  :: ic              ! direct beam: ic=0; diffuse: ic=1
  integer  :: nband = numrad  ! number of solar radiation wave bands
  real(r8) :: wl              ! fraction of LAI+SAI that is LAI
  real(r8) :: ws              ! fraction of LAI+SAI that is SAI
  real(r8) :: mpe = 1.e-06    ! prevents overflow for division by zero 
  real(r8) :: vai             ! elai+esai
  real(r8) :: rho(numrad)     ! leaf/stem refl weighted by fraction LAI and SAI
  real(r8) :: tau(numrad)     ! leaf/stem tran weighted by fraction LAI and SAI
  real(r8) :: ftdi(numrad)    ! down direct flux below veg per unit dif flux = 0
  real(r8) :: albsnd(numrad)  ! snow albedo (direct)
  real(r8) :: albsni(numrad)  ! snow albedo (diffuse)
  real(r8) :: gdir            ! aver projected leaf/stem area in solar direction
  real(r8) :: ext             ! optical depth direct beam per unit LAI+SAI
  real(r8) :: delta           ! solar declination angle in radians
  real(r8) :: eccf            ! earth orbit eccentricity factor
  real(r8) :: coszen          ! cosine solar zenith angle for next time step
  real(r8) :: albd(2)
  real(r8) :: albi(2)

!----End Variable List--------------------------------------------------

!
! Solar declination  for next time step
!

  call shr_orb_decl (caldayp1, eccen, mvelpp, lambm0, obliqr, &
                     delta   , eccf )

!
! Cosine solar zenith angle for next time step
!

  coszen = shr_orb_cosz(caldayp1, clm%lat, clm%lon, delta)

!
! Initialize output because solar radiation only done if coszen > 0
!

  do ib = 1, nband
     albd(ib)   = 1.
     albi(ib)   = 1.
     clm%albgrd(ib) = 0._r8
     clm%albgri(ib) = 0._r8
     clm%fabd(ib)   = 0._r8
     clm%fabi(ib)   = 0._r8
     clm%ftdd(ib)   = 0._r8
     clm%ftid(ib)   = 0._r8
     clm%ftii(ib)   = 0._r8
     if (ib==1) clm%fsun = 0.
  end do

!===LDAS modification
!Original CLM2 code for solar zenith angle less than 0 is replaced with that in CLM1
!
! Return if coszen is not positive
!
!  if (coszen <= 0._r8) RETURN

!===LDAS modification
! IF COSZEN IS NOT POSITIVE - NEW for CLM_OFFLINE
! NOTE: All NEW changes for CLM offline assumes that the incoming solar
! radiation is 70% direct, 30% diffuse, and 50% visible, 50% near-infrared
! (as is assumed in atmdrvMod.F90).

  do ib = 1, nband
    albsnd(ib)     = 0._r8
    albsni(ib)     = 0._r8
  enddo
  if (coszen <= 0._r8) then
   clm%surfalb = 35./100.*(albd(1)+albd(2)) &
		 +15./100.*(albi(1)+albi(2))
   clm%snoalb = 35./100.*(albsnd(1)+albsnd(2)) &
                 +15./100.*(albsni(1)+albsni(2))
   RETURN
  endif

!
! weight reflectance/transmittance by lai and sai
!

  do ib = 1, nband
     vai = clm%elai + clm%esai
     wl = clm%elai / max( vai,mpe )
     ws = clm%esai / max( vai,mpe )
     rho(ib) = max( clm%rhol(ib)*wl + clm%rhos(ib)*ws, mpe )
     tau(ib) = max( clm%taul(ib)*wl + clm%taus(ib)*ws, mpe )
  end do

!
! Snow albedos: only if h2osno > 0
!

  if ( clm%h2osno > 0._r8 ) then
     ic=0; call SnowAlbedo (clm, coszen, nband, ic, albsnd)
     ic=1; call SnowAlbedo (clm, coszen, nband, ic, albsni)  
  else
     albsnd(:) = 0._r8
     albsni(:) = 0._r8
  endif

!===LDAS modification
!===NEW for CLM offline

!  clm%snoalb = 35./100.*(albsnd(1)+albsnd(2)) + 15./100.*(albsni(1)+albsni(2))

!
! Ground surface albedos
!
  call SoilAlbedo (clm, coszen, nband, albsnd, albsni)      
  if (vai /= 0.) then  ! vegetated patch
!   if(clm%itypveg .ne.12) then 
!
! Loop over nband wavebands to calculate surface albedos and solar 
! fluxes for vegetated patch for unit incoming direct 
! (ic=0) and diffuse flux (ic=1)
!
     do ib = 1, nband
        ic = 0
        call TwoStream (clm,      ib,  ic,       coszen,   vai,      &
                        rho,      tau, clm%fabd, albd, clm%ftdd, &
                        clm%ftid, gdir )
        ic = 1
        call TwoStream (clm,      ib,  ic,       coszen,   vai,  &
                        rho,      tau, clm%fabi, albi, ftdi, &
                        clm%ftii, gdir )
     end do
!     print*, 'aft twostream.. ',clm%kpatch, clm%fabd,clm%fabi
!
! Sunlit fraction of canopy. Set fsun = 0 if fsun < 0.01.
!
     
     ext = gdir/coszen * sqrt(1.-rho(1)-tau(1))
     clm%fsun = (1.-exp(-ext*vai)) / max(ext*vai,mpe)
     ext = clm%fsun                                       !temporary fsun
     if (ext < 0.01) then 
        wl = 0._r8                                        !temporary fsun
     else
        wl = ext                                          !temporary fsun
     end if
     clm%fsun = wl

  else     ! non-vegetated patch

     do ib = 1,numrad
        clm%fabd(ib) = 0.
        clm%fabi(ib) = 0.
        clm%ftdd(ib) = 1.
        clm%ftid(ib) = 0.
        clm%ftii(ib) = 1.
        albd(ib) = clm%albgrd(ib)
        albi(ib) = clm%albgri(ib)
        clm%fsun     = 0.
     end do

  endif

!===LDAS modification
!===NEW for CLM offline:
  
    clm%surfalb = 35./100.*(albd(1)+albd(2)) +15./100.*(albi(1)+albi(2))
!    print*, 'end twostream.. ',clm%kpatch, clm%fabd,clm%fabi
  return
end subroutine SurfaceAlbedo




