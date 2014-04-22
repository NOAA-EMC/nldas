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

subroutine SoilAlbedo (clm, coszen, nband, albsnd, albsni)

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
! Determine ground surface albedo, accounting for snow
! 
! Method: 
! 
! Author: 
! Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: SoilAlbedo.F90,v 1.7 2004/11/24 22:56:39 jim Exp $
!-----------------------------------------------------------------------

  use precision
  use clmtype
!  use clm_varpar, only : numrad
  use clm_varcon, only : albsat, albdry, alblak, albice, tfrz, istice, istsoil
!  use spmdMod
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm            !CLM 1-D Module

  integer, intent(in) :: nband           ! number of solar radiation waveband classes
  real(r8), intent(in) :: coszen         ! cosine solar zenith angle for next time step
  real(r8), intent(in) :: albsnd(numrad) ! snow albedo (direct)
  real(r8), intent(in) :: albsni(numrad) ! snow albedo (diffuse)

!----Local Variables----------------------------------------------------

  integer  ib      !waveband number (1=vis, 2=nir)
  real(r8) inc     !soil water correction factor for soil albedo
  real(r8) albsod  !soil albedo (direct)
  real(r8) albsoi  !soil albedo (diffuse)
!----End Variable List--------------------------------------------------

  do ib = 1, nband
     if (clm%itypwat == istsoil)  then               !soil
        inc    = max(0.11-0.40*clm%h2osoi_vol(1), 0._r4)
        albsod = min(albsat(clm%isoicol,ib)+inc, albdry(clm%isoicol,ib))
        albsoi = albsod
     else if (clm%itypwat == istice)  then           !land ice
        albsod = albice(ib)
        albsoi = albsod
     else if (clm%t_grnd > tfrz) then                !unfrozen lake, wetland
        albsod = 0.05/(max(0.001_r4,coszen) + 0.15)
        albsoi = albsod
     else                                            !frozen lake, wetland
        albsod = alblak(ib)
        albsoi = albsod
     end if
     clm%albgrd(ib) = albsod*(1.-clm%frac_sno) + albsnd(ib)*clm%frac_sno
     clm%albgri(ib) = albsoi*(1.-clm%frac_sno) + albsni(ib)*clm%frac_sno
  end do
  return
end subroutine SoilAlbedo
