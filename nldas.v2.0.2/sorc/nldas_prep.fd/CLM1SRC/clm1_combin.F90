#include <misc.h>

subroutine clm1_combin (clm1) 

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
!  This subroutine checks for elements which are below the prescribed 
!  minimum for thickness or mass.  If the snow element thickness or mass 
!  is less than a prescribed minimum, then it is combined with a 
!  neighboring element.  The subroutine clm_combo.f90 then executes the 
!  combination of mass and energy.
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!=========================================================================
! $Id: clm1_combin.F90,v 1.1.1.1 2003/02/06 16:10:45 jgottsch Exp $
!=========================================================================

! Declare Modules and data structures

  use precision
  use clm1type
  implicit none

!=== Arguments ===========================================================

  type (clm11d), intent(inout) :: clm1      !CLM 1-D Module

!=== Local Variables =====================================================

  real(r8)                    &
       dzmin(5),              & ! minimum of snow layer 1 (top) to msn0 (bottom)
       zwice,                 & ! total ice mass in snow
       zwliq                    ! total liquid water in snow

  integer                     & !
       i,                     & ! number of do looping
       j,                     & ! node index
       k,                     & ! number of do looping
       l,                     & ! node index
       msn_old,               & ! number of snow layer 1 (top) to msn0 (bottom)
       mssi,                  & ! node index
       neibor                   ! adjacent node selected for combination

  data dzmin /0.010, 0.015, 0.025, 0.055, 0.115/

!=== End Variable List ===================================================

! Check the mass of ice lens of snow, when the total is less than a small value,
! combine it with the underlying neighbor.

  msn_old = clm1%snl
  do j = msn_old+1, 0
     if(clm1%h2osoi_ice(j) <= .1)then
        clm1%h2osoi_liq(j+1) = clm1%h2osoi_liq(j+1) + clm1%h2osoi_liq(j)
        clm1%h2osoi_ice(j+1) = clm1%h2osoi_ice(j+1) + clm1%h2osoi_ice(j)

! shift all elements above this down one.
        if(j > clm1%snl+1 .AND. clm1%snl < -1)then
           do i =  j, clm1%snl+2, -1
              clm1%t_soisno(i) = clm1%t_soisno(i-1)
              clm1%h2osoi_liq(i) = clm1%h2osoi_liq(i-1)
              clm1%h2osoi_ice(i) = clm1%h2osoi_ice(i-1)
              clm1%dz(i) = clm1%dz(i-1)
           enddo
        endif
        clm1%snl = clm1%snl + 1
!*      write(79,*) 'one snow layer is gone'
     endif
  enddo

  if(clm1%snl == 0)then
     clm1%h2osno = 0.
     clm1%snowdp = 0.
!*     write(79,*) 'all snow has gone'
     return
  else
     clm1%h2osno = 0.
     clm1%snowdp = 0.
     zwice = 0.
     zwliq = 0.
     do j = clm1%snl + 1, 0
        clm1%h2osno = clm1%h2osno + clm1%h2osoi_ice(j) + clm1%h2osoi_liq(j)
        clm1%snowdp = clm1%snowdp + clm1%dz(j)
        zwice = zwice + clm1%h2osoi_ice(j)
        zwliq = zwliq + clm1%h2osoi_liq(j)
     enddo
  endif

! Check the snow depth
  if(clm1%snowdp < 0.01)then       !!! all snow gone 
     clm1%snl = 0
     clm1%h2osno = zwice
     if(clm1%h2osno <= 0.) clm1%snowdp = 0.

! The liquid water assumed ponding on soil surface.
     clm1%h2osoi_liq(1) = clm1%h2osoi_liq(1) + zwliq
!**    write(79,'(17h all snow is gone)')
     return
  else                        !!! snow layers combined

! two or more layers 
     if(clm1%snl < -1)then
        msn_old = clm1%snl
        mssi = 1
        do i = msn_old+1, 0

! If top node is removed, combine with bottom neighbor.
           if(clm1%dz(i) < dzmin(mssi))then
              if(i == clm1%snl+1)then
                 neibor = i + 1

! If the bottom neighbor is not snow, combine with the top neighbor.
              else if(i == 0)then
                 neibor = i - 1

! If none of the above special cases apply, combine with the thinnest neighbor
              else
                 neibor = i + 1
                 if((clm1%dz(i-1)+clm1%dz(i)) < (clm1%dz(i+1)+clm1%dz(i))) neibor = i-1
              endif

! Node l and j are combined and stored as node j.
              if(neibor > i)then
                 j = neibor
                 l = i
              else
                 j = i
                 l = neibor
              endif

              call clm1_combo ( clm1%dz(j), clm1%h2osoi_liq(j), clm1%h2osoi_ice(j), clm1%t_soisno(j),&
                   clm1%dz(l), clm1%h2osoi_liq(l), clm1%h2osoi_ice(l), clm1%t_soisno(l) )

! Now shift all elements above this down one.
              if(j-1 > clm1%snl+1) then
                 do k = j-1, clm1%snl+2, -1
                    clm1%t_soisno(k) = clm1%t_soisno(k-1)
                    clm1%h2osoi_ice(k) = clm1%h2osoi_ice(k-1)
                    clm1%h2osoi_liq(k) = clm1%h2osoi_liq(k-1)
                    clm1%dz(k) = clm1%dz(k-1)
                 enddo
              endif

              clm1%snl = clm1%snl + 1

!**    write(79,'(7h Nodes ,i4,4h and,i4,14h combined into,i4)') l,j,j

              if(clm1%snl >= -1) EXIT

! The layer thickness is greater than the prescribed minimum value
           else
              mssi = mssi + 1 
           endif
        enddo

     endif

! Reset the node depth and the depth of layer interface
     do k = 0, clm1%snl+1, -1
        clm1%z(k) = clm1%zi(k) - 0.5*clm1%dz(k)
        clm1%zi(k-1) = clm1%zi(k) - clm1%dz(k)
     enddo

  endif                       !!! snow layers combined 

end subroutine clm1_combin

