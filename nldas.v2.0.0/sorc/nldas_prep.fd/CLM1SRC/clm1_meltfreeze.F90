#include <misc.h>

subroutine clm1_meltfreeze (fact,     brr,       hs,     dhsdT,  &   
                           tssbef,   xmf,       clm1) 

!=========================================================================
!
!  clm1clm1CLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
!  Calculation of the phase change within snow and soil layers:
!
!  (1) Check the conditions for which the phase change may take place, 
!      i.e., the layer temperature is great than the freezing point 
!      and the ice mass is not equal to zero (i.e. melting), 
!      or the layer temperature is less than the freezing point 
!      and the liquid water mass is not equal to zero (i.e. freezing).
!  (2) Assess the rate of phase change from the energy excess (or deficit) 
!      after setting the layer temperature to freezing point.
!  (3) Re-adjust the ice and liquid mass, and the layer temperature
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!=========================================================================
! $Id: clm1_meltfreeze.F90,v 1.1.1.1 2003/02/06 16:10:46 jgottsch Exp $
!=========================================================================

! Declare Modules and data structures

  use precision
  use clm1type
  use clm1_varcon, only : tfrz, hfus
  use clm1_varpar, only : nlevsoi
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1   ! CLM 1-D Module

  real(r8), intent(in) ::     &
       tssbef(clm1%snl+1 : nlevsoi),   & ! temperature at previous time step [K]
       brr   (clm1%snl+1 : nlevsoi),   & ! 
       fact  (clm1%snl+1 : nlevsoi),   & ! temporary variables
       hs,                            & ! net ground heat flux into the surface
       dhsdT                            ! temperature derivative of "hs"

  real(r8), intent(out) ::    &
       xmf                              ! total latent heat of phase change

!=== Local Variables =====================================================

  integer j
  real(r8)  hm(clm1%snl+1 : nlevsoi),     & ! energy residual [W/m2]
       xm(clm1%snl+1 : nlevsoi),       & ! melting or freezing within a time step [kg/m2]
       heatr,                         & ! energy residual or loss after melting or freezing
       temp1                            ! temporary variables [kg/m2]

  real(r8), dimension(clm1%snl+1 : nlevsoi) :: wmass0, wice0, wliq0
  real(r8)  propor,tinc  

!=== End Variable List ===================================================

! Initial 

  clm1%qflx_snomelt = 0.
  xmf = 0.
  do j = clm1%snl+1, nlevsoi
     clm1%imelt(j) = 0
     hm(j) = 0.
     xm(j) = 0.
     wice0(j) = clm1%h2osoi_ice(j)
     wliq0(j) = clm1%h2osoi_liq(j)
     wmass0(j) = clm1%h2osoi_ice(j) + clm1%h2osoi_liq(j)
  enddo

! Melting identification
! If ice exists above melt point, melt some to liquid.

  do j = clm1%snl+1, nlevsoi
     if (clm1%h2osoi_ice(j) > 0. .AND. clm1%t_soisno(j) > tfrz) then
        clm1%imelt(j) = 1
        clm1%t_soisno(j) = tfrz
     endif

     ! Freezing identification
     ! If liquid exists below melt point, freeze some to ice.

     if (clm1%h2osoi_liq(j) > 0. .AND. clm1%t_soisno(j) < tfrz) then
        clm1%imelt(j) = 2
        clm1%t_soisno(j) = tfrz
     endif
  enddo

! If snow exists, but its thickness is less than the critical value (0.01 m)

  if (clm1%snl+1 == 1 .AND. clm1%h2osno > 0.) then
     if (clm1%t_soisno(1) > tfrz) then
        clm1%imelt(1) = 1
        clm1%t_soisno(1) = tfrz
     endif
  endif

! Calculate the energy surplus and loss for melting and freezing

  do j = clm1%snl+1, nlevsoi
     if (clm1%imelt(j) > 0) then
        tinc = clm1%t_soisno(j)-tssbef(j)
        if (j > clm1%snl+1) then
           hm(j) = brr(j) - tinc/fact(j) 
        else
           hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j) 
        endif
     endif
  enddo

! These two errors were checked carefully.  They result from the the computed error
! of "Tridiagonal-Matrix" in subroutine "thermal".

  do j = clm1%snl+1, nlevsoi
     if (clm1%imelt(j) == 1 .AND. hm(j) < 0.) then
        hm(j) = 0.
        clm1%imelt(j) = 0
     endif

     if (clm1%imelt(j) == 2 .AND. hm(j) > 0.) then
        hm(j) = 0.
        clm1%imelt(j) = 0
     endif
  enddo

! The rate of melting and freezing

  do j = clm1%snl+1, nlevsoi

     if (clm1%imelt(j) > 0 .and. abs(hm(j)) > .0) then
        xm(j) = hm(j)*clm1%dtime/hfus                        ! kg/m2

        ! If snow exists, but its thickness is less than the critical value (1 cm)
        ! Note: more work is needed to determine how to tune the snow depth for this case

        if (j == 1) then
           if ((clm1%snl+1 == 1) .AND. (clm1%h2osno > 0.) .AND. (xm(j) > 0.)) then
              temp1 = clm1%h2osno                                        ! kg/m2
              clm1%h2osno = max(0.,temp1-xm(j))
              propor = clm1%h2osno/temp1
              clm1%snowdp = propor * clm1%snowdp
              heatr = hm(j) - hfus*(temp1-clm1%h2osno)/clm1%dtime         ! W/m2
              if (heatr > 0.) then
                 xm(j) = heatr*clm1%dtime/hfus                           ! kg/m2
                 hm(j) = heatr                                          ! W/m2
              else
                 xm(j) = 0.
                 hm(j) = 0.
              endif
              clm1%qflx_snomelt = max(0.,(temp1-clm1%h2osno))/clm1%dtime   ! kg/(m2 s)
              xmf = hfus*clm1%qflx_snomelt
           endif
        endif

        heatr = 0.
        if (xm(j) > 0.) then
           clm1%h2osoi_ice(j) = max(0., wice0(j)-xm(j))
           heatr = hm(j) - hfus*(wice0(j)-clm1%h2osoi_ice(j))/clm1%dtime
        else if (xm(j) < 0.) then
           clm1%h2osoi_ice(j) = min(wmass0(j), wice0(j)-xm(j))
           heatr = hm(j) - hfus*(wice0(j)-clm1%h2osoi_ice(j))/clm1%dtime  
        endif

        clm1%h2osoi_liq(j) = max(0.,wmass0(j)-clm1%h2osoi_ice(j))

        if (abs(heatr) > 0.) then
           if (j > clm1%snl+1) then
              clm1%t_soisno(j) = clm1%t_soisno(j) + fact(j)*heatr
           else
              clm1%t_soisno(j) = clm1%t_soisno(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
           endif
           if (clm1%h2osoi_liq(j)*clm1%h2osoi_ice(j)>0.) clm1%t_soisno(j) = tfrz
        endif

        xmf = xmf + hfus * (wice0(j)-clm1%h2osoi_ice(j))/clm1%dtime

        if (clm1%imelt(j) == 1 .AND. j < 1) then
           clm1%qflx_snomelt = clm1%qflx_snomelt + max(0.,(wice0(j)-clm1%h2osoi_ice(j)))/clm1%dtime  
        endif

     endif

  enddo

! needed for output to lsm history file

  clm1%eflx_snomelt = clm1%qflx_snomelt * hfus  

end subroutine clm1_meltfreeze
