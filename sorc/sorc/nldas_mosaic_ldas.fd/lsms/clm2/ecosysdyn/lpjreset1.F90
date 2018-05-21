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

subroutine lpjreset1 (caldayp1, eccen, obliqr, lambm0, mvelpp)
 
  use precision
  use clmtype
  use clm_varder
  use clm_varmap, only : begpatch,endpatch
  implicit none

! ----------------------------------------------------------------------
! purpose           : to reset variables related to lpj
! date first created: November 2000 - lsm version 2 
! by whom           : Sam Levis
! date last revised : 
! by whom           : 
! ----------------------------------------------------------------------
! $Id: lpjreset1.F90,v 1.6 2004/11/24 22:57:05 jim Exp $
! ----------------------------------------------------------------------

! -------------------------- arguments ---------------------------------
  real(r4), intent(in) :: caldayp1 !calendar day at Greenwich (1.00, ..., 365.99) for nstep+1
  real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox long. of perihelion + pi (radians)
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
  integer  :: k                    !indices 
! ----------------------------------------------------------------------

! reset a few variables here at the very end of the year
#if 0 
  do k = begpatch,endpatch
     clm(k)%annpsn     = 0.
     clm(k)%annpsnpot  = 0.
     clm(k)%bm_inc     = 0.
     clm(k)%afmicr     = 0.
     clm(k)%firelength = 0.
     clm(k)%agddtw     = 0.
     clm(k)%agdd       = 0.
     clm(k)%t10min     = 1.0e+36
     clm(k)%t_mo_min   = 1.0e+36
  end do

! call EcosystemDyn because need info for first timestep of next year

  do k = begpatch,endpatch

     call EcosystemDyn (clm(k), .false., .true.)

     call SurfaceAlbedo (clm(k), caldayp1, eccen, obliqr, lambm0, mvelpp)

  end do

  call iniTimeConstDGVM ()
#endif
  return
end subroutine lpjreset1
