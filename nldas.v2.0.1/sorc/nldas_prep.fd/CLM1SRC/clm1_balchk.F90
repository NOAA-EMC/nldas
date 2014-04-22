#include <misc.h>

subroutine clm1_balchk (clm1) 

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
!  Water and energy balance check  
!
!  This subroutine accumulates the numerical truncation errors of the water
!  and energy balance calculation. It is helpful to see the performance of 
!  the process of integration.
!
!  The error for energy balance: 
!  error = abs(Net radiation - the change of internal energy - Sensible heat
!              - Latent heat) 
!  The error should be less than 0.02 W/m2 in each time integration interval;
!
!  The error for water balance:
!  error = abs(precipitation - change of water storage - evaporation - runoff)
!  The error should be less than 0.001 mm in  each time integration interval.
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!  10 November 2000: Mariana Vertenstein
!  03 June 2003: Jon Gottschalck; Modified water balance error cutoff to 0.20 from 0.10
!=========================================================================
! $Id: clm1_balchk.F90,v 1.2 2003/06/03 21:23:01 jgottsch Exp $
!=========================================================================

  use precision
  use clm1type
  use clm1_varpar, only : nlevsoi
  use clm1_varcon, only : istsoil, tfrz
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1	 !CLM 1-D Module

!=== Local Variables =====================================================

  integer j                    ! do loop index
  logical :: constop = .false. ! true => stop if energy balance err too great

!=== End Variable List ===================================================

!------------------------------------------------------------------------
! Water balance 
!------------------------------------------------------------------------

  clm1%errh2o = clm1%endwb - clm1%begwb - &
           ( clm1%forc_rain  + clm1%forc_snow - clm1%qflx_evap_tot - clm1%qflx_surf &
           - clm1%qflx_qrgwl + clm1%qflx_qirr - clm1%qflx_drain ) * clm1%dtime

  if (abs(clm1%errh2o) > .20) then
     write(79,200)'water balance error',clm1%istep,clm1%kpatch,clm1%errh2o
     write(79,*)'clm1 model is stopping'
     call endrun
  endif

!------------------------------------------------------------------------
! Solar radiation energy balance
!------------------------------------------------------------------------

  clm1%errsol = clm1%fsa + clm1%fsr - &
              (clm1%forc_solad(1) + clm1%forc_solad(2) + clm1%forc_solai(1) + clm1%forc_solai(2))

  if (abs(clm1%errsol) > .10 ) then
     write(79,100)'solar radiation balance error',clm1%istep,clm1%kpatch,clm1%errsol
     write(79,*)'clm1 model is stopping'
     call endrun
  endif

!------------------------------------------------------------------------
! Longwave radiation energy balance
!------------------------------------------------------------------------

  clm1%errlon = clm1%eflx_lwrad_out - clm1%eflx_lwrad_net - clm1%forc_lwrad

  if (abs(clm1%errlon) > .10 ) then
     write(79,100)'longwave energy balance error',clm1%istep,clm1%kpatch,clm1%errlon
     write(79,*)'clm1 model is stopping'
     call endrun
  endif

!------------------------------------------------------------------------
! Surface energy balance
!------------------------------------------------------------------------

  clm1%errseb = clm1%sabv + clm1%sabg  &
             + clm1%forc_lwrad - clm1%eflx_lwrad_out &
             - clm1%eflx_sh_tot &
             - clm1%eflx_lh_tot &
             - clm1%eflx_soil_grnd 
                
  if (abs(clm1%errseb) > .10 ) then
     write(79,100)'surface flux energy balance error',clm1%istep,clm1%kpatch,clm1%errseb
     write(79,*)'clm1 model is stopping'
     call endrun
  endif

!------------------------------------------------------------------------
! Accumulation of water and surface energy balance error
!------------------------------------------------------------------------

  clm1%acc_errh2o = clm1%acc_errh2o + clm1%errh2o
  clm1%acc_errseb = clm1%acc_errseb + clm1%errseb

100 format (1x,a14,' istep =',i10,' point =',i6,' imbalance =',f8.2,' W/m2') 
200 format (1x,a14,' istep =',i10,' point =',i6,' imbalance =',f8.2,' mm') 

end subroutine clm1_balchk











