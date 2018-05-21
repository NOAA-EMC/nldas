#include <misc.h>

subroutine drv_getforce (ldas,drv,ldas_grid,tile,clm1)

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
!  Access meteorological data - this current version reads 1D forcing
!  and distributes it to the clm domain (spatially constant).  This routine
!  must be modified to allow for spatially variable forcing, or coupling to
!  a GCM.
!
!  The user may likely want to modify this subroutine significantly,
!  to include such things as space/time intrpolation of forcing to the
!  CLM grid, reading of spatially variable binary data, etc.
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December  1999: Paul Houser and Jon Radakovich; F90 Revision
!  28 January   2002: Jon Gottschalck; Added variables used in surface initialization
!  03 June      2003: Jon Gottschalck, Modified initialization section for GEOS
!=========================================================================
! $Id: drv_getforce.F90,v 1.2 2003/06/03 21:22:29 jgottsch Exp $
!=========================================================================

  use precision
  use ldas_module         ! LDAS non-model-specific 1-D variables
  use drv_module          ! 1-D Land Model Driver variables
  use grid_module         ! LDAS non-model-specific grid variables
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! 1-D CLM variables
  use clm1_varcon, only : tfrz, tcrit
  implicit none

!=== Arguments ===========================================================

  type (ldasdec)                   :: ldas
  type (drvdec) ,intent(inout)     :: drv              
  type (griddec),intent(inout)     :: ldas_grid(drv%nc,drv%nr)
  type (clm_tiledec),intent(inout) :: tile(drv%nch)
  type (clm11d)  ,intent(inout)     :: clm1 (drv%nch)

!=== Local Variables =====================================================

  real(r8) solar(drv%nch)     ! incident solar radiation [w/m2]
  real(r8) prcp(drv%nch)      ! precipitation [mm/s]
  integer t          ! Tile looping variable

!=== End Variable List ===================================================

!=== Increment Time Step Counter

  clm1%istep=clm1%istep+1 

  do t=1,drv%nch
     clm1(t)%forc_t   =  ldas_grid(tile(t)%col,tile(t)%row)%forcing(1)
     clm1(t)%forc_q   =  ldas_grid(tile(t)%col,tile(t)%row)%forcing(2)
     solar(t)        =  ldas_grid(tile(t)%col,tile(t)%row)%forcing(3)
     clm1(t)%forc_solad(1) =  solar(t)*35./100.
     clm1(t)%forc_solad(2) =  solar(t)*35./100.
     clm1(t)%forc_solai(1) =  solar(t)*15./100.
     clm1(t)%forc_solai(2) =  solar(t)*15./100.
     clm1(t)%forc_lwrad   = ldas_grid(tile(t)%col,tile(t)%row)%forcing(4)
     clm1(t)%forc_u       = ldas_grid(tile(t)%col,tile(t)%row)%forcing(5)
     clm1(t)%forc_v       = ldas_grid(tile(t)%col,tile(t)%row)%forcing(6)
     clm1(t)%forc_pbot    = ldas_grid(tile(t)%col,tile(t)%row)%forcing(7)
     prcp(t)             = ldas_grid(tile(t)%col,tile(t)%row)%forcing(8)
     clm1(t)%forc_rho     = clm1(t)%forc_pbot/(clm1(t)%forc_t*2.8704e2)

     if (drv%startcode .eq. 4 .and. ldas%tscount .eq. 0) then
      select case (ldas%force)
      case(1)
      clm1(t)%forc_swc1   = ldas_grid(tile(t)%col,tile(t)%row)%forcing(11)
      clm1(t)%forc_swc2   = ldas_grid(tile(t)%col,tile(t)%row)%forcing(12)
      clm1(t)%forc_stemp1 = ldas_grid(tile(t)%col,tile(t)%row)%forcing(13)
      clm1(t)%forc_stemp2 = ldas_grid(tile(t)%col,tile(t)%row)%forcing(14)
      clm1(t)%forc_sdepth = ldas_grid(tile(t)%col,tile(t)%row)%forcing(15)
      case(2)
      clm1(t)%forc_swc1   = ldas_grid(tile(t)%col,tile(t)%row)%forcing(13)
      clm1(t)%forc_sdepth = ldas_grid(tile(t)%col,tile(t)%row)%forcing(12) 
      case default
       print*, "ONLY FORCING OPTIONS AVIALBLE FOR THIS METHOD ARE GDAS, GEOS"
       STOP
      end select      
     endif

! Next set upper limit of air temperature for snowfall at 275.65K.
! This cut-off was selected based on Fig. 1, Plate 3-1, of Snow
! Hydrology (1956).

     if (prcp(t) > 0.) then
        if(clm1(t)%forc_t > (tfrz + tcrit))then
           clm1(t)%itypprc = 1
           clm1(t)%forc_rain = prcp(t)

           clm1(t)%forc_snow = 0.
        else
           clm1(t)%itypprc = 2
           clm1(t)%forc_rain = 0.
           clm1(t)%forc_snow = prcp(t)
        endif
     else
        clm1(t)%itypprc = 0
        clm1(t)%forc_rain = 0.
        clm1(t)%forc_snow = 0
     endif
   enddo

end subroutine drv_getforce
