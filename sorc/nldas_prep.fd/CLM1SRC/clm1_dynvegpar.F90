#include <misc.h>

subroutine clm1_dynvegpar (clm1)

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
!  Vegetation dynamic parameters and snow cover fraction as subgrid vectors 
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December  1999: Paul Houser and Jon Radakovich; F90 Revision 
!  15 November  2000: Mariana Vertenstein
!  01 October   2002: Jon Gottschalck; Modifed to account for MODIS LAI data 
!  20 November  2002: Jon Radakovich; Modified calculation of seasb for
!                     5 soil layers.
!=========================================================================
! $Id: clm1_dynvegpar.F90,v 1.1.1.1 2003/02/06 16:10:46 jgottsch Exp $
!=========================================================================

  use precision
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! CLM tile variables
  use clm1_varcon, only :istice
  implicit none

!=== Arguments ===========================================================

  type (clm11d)   :: clm1

!=== Local Variables =====================================================

  real(r8) seasb   !temperature dependence of vegetation cover [-]
  real(r8) fb      !fraction of canopy layer covered by snow
  real(r8) ol      !thickness of canopy layer covered by snow (m)

!=== End Variable List ===================================================

! Note: temporarily set, they can be given by measurement, or dynamic ecosystem model
! only nonzero if NOT glacier/ice, water or bare soil

  if ((.not. clm1%lakpoi) .AND. (.not. clm1%baresoil) .AND. (clm1%itypwat/=istice)) then
    IF (clm1%LAIFLAG .GE. 2) THEN
       clm1%tlai  = clm1%tlai
       clm1%tsai  = clm1%tsai
     ELSE
       seasb = max(0., 1. - 0.0016*max(298.-clm1%t_soisno(5), .0)**2)
       clm1%tlai  = clm1%maxlai + (clm1%minlai-clm1%maxlai)*(1.-seasb)
       clm1%tsai  = clm1%tsai
     ENDIF
  else
     clm1%tlai  = 0.
     clm1%tsai  = 0.
  endif

! Adjust lai and sai for burying by snow. if exposed lai and sai are less than 0.05,
! set equal to zero to prevent numerical problems associated with very small lai,sai

  ol = min( max(clm1%snowdp-clm1%hbot,0._r8), clm1%htop-clm1%hbot)
  fb = 1. - ol / max(1.e-06, clm1%htop-clm1%hbot)
  clm1%elai = clm1%tlai*fb
  clm1%esai = clm1%tsai*fb

  if (clm1%elai < 0.05) clm1%elai = 0._r8
  if (clm1%esai < 0.05) clm1%esai = 0._r8

! Fraction of vegetation free of snow

  if ((clm1%elai + clm1%esai) >= 0.05) then
     clm1%frac_veg_nosno = 1
  else
     clm1%frac_veg_nosno = 0
  endif
  
! Fraction of soil covered by snow

  clm1%frac_sno = clm1%snowdp/(10.*clm1%zlnd + clm1%snowdp)  

end subroutine clm1_dynvegpar
