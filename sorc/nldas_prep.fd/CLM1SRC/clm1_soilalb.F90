#include <misc.h>

subroutine clm1_soilalb (clm1, coszen, nband, albsnd, albsni)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Determine ground surface albedo, accounting for snow
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: clm1_soilalb.F90,v 1.1.1.1 2003/02/06 16:10:46 jgottsch Exp $
!-----------------------------------------------------------------------

  use precision
  use clm1type
  use clm1_varpar, only : numrad
  use clm1_varcon, only : albsat, albdry, alblak, albice, tfrz, istice, istsoil
  implicit none

! ------------------------- arguments ----------------------------
  type (clm11d), intent(inout) :: clm1            !CLM 1-D Module
  real(r8)    , intent(in)    :: coszen         !cosine solar zenith angle for next time step
  integer     , intent(in)    :: nband          !number of solar radiation waveband classes
  real(r8)    , intent(in)    :: albsnd(numrad) !snow albedo (direct)
  real(r8)    , intent(in)    :: albsni(numrad) !snow albedo (diffuse)
! -----------------------------------------------------------------

! ------------------------- local variables -----------------------
  integer  ib      !waveband number (1=vis, 2=nir)
  real(r8) inc     !soil water correction factor for soil albedo
  real(r8) albsod  !soil albedo (direct)
  real(r8) albsoi  !soil albedo (diffuse)
! -----------------------------------------------------------------

  do ib = 1, nband
     if (clm1%itypwat == istsoil)  then               !soil
        inc    = max(0.11-0.40*clm1%h2osoi_vol(1), 0._r8)
        albsod = min(albsat(clm1%isoicol,ib)+inc, albdry(clm1%isoicol,ib))
        albsoi = albsod
     else if (clm1%itypwat == istice)  then           !land ice
        albsod = albice(ib)
        albsoi = albsod
     else if (clm1%t_grnd > tfrz) then                !unfrozen lake, wetland
        albsod = 0.05/(max(0.001,coszen) + 0.15)
        albsoi = albsod
     else                                            !frozen lake, wetland
        albsod = alblak(ib)
        albsoi = albsod
     end if

     clm1%albgrd(ib) = albsod*(1.-clm1%frac_sno) + albsnd(ib)*clm1%frac_sno
     clm1%albgri(ib) = albsoi*(1.-clm1%frac_sno) + albsni(ib)*clm1%frac_sno

  end do

  return
end subroutine clm1_soilalb


