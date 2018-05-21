!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! clm_bias.F90: 
!
! DESCRIPTION:
!  This subroutine corrects the analysis increment with the forecast bias.  
!
! REVISION HISTORY:
!  26 July 2001: Jon Radakovich; Initial Code
!
!=========================================================================
subroutine clm1_bias(ldas,grid,tile,clm1)
  use ldas_module      ! ldas non-model-specific 1-d variables
  use grid_module      ! ldas non-model-specific grid variables
  use tile_module      ! ldas non-model-specific tile variables
  use clm1type          ! 1-D CLM variables
  use precision

  implicit none
  type (ldasdec)     :: ldas
  type (griddec)     :: grid(ldas%nc,ldas%nr)
  type (tiledec)     :: tile(ldas%nch)
  type (clm11d)       :: clm1(ldas%nch)

!=== Local Variables =====================================================
  integer t,c,r        ! Loop counters
  real(r8) time            ! Model time of day in seconds
  real(r8) gamma           ! Tunable weight parameter for bias correction
  real(r8) omega1,omega2   ! Frequencies of harmonics
  real(r8) alpha           ! Amplitude of harmonic
  real(r8) fbiast(ldas%nch)       ! Forecast bias in tile space
  real(r8), parameter :: pi = 3.14159265

!=== End Variables =====================================================

  time=(float(ldas%hr))*3600.+  &
       (float(ldas%mn))*60.  +  &
       float(ldas%ss)

  gamma=0.2
  omega1=2.*pi/86400.
  omega2=2.*pi/43200.

  do r=1,ldas%nr
   do c=1,ldas%nc
    alpha=gamma*grid(c,r)%clmdelt
    if(ldas%ribc.eq.1)then
     grid(c,r)%clmbetak(1)=grid(c,r)%clmbetak(1) - alpha
     grid(c,r)%clmbetak(2)=0.
     grid(c,r)%clmbetak(3)=0.
     grid(c,r)%clmbetak(4)=0.
     grid(c,r)%clmbetak(5)=0.
    else if(ldas%rdbc.eq.1)then
     grid(c,r)%clmbetak(1)=grid(c,r)%clmbetak(1) - alpha
     grid(c,r)%clmbetak(2)=grid(c,r)%clmbetak(2) -    &
                           alpha*cos(omega1*time)
     grid(c,r)%clmbetak(3)=grid(c,r)%clmbetak(3) -    &
                           alpha*sin(omega1*time)
     grid(c,r)%clmbetak(4)=0.
     grid(c,r)%clmbetak(5)=0.        
    else if(ldas%rsdbc.eq.1)then
     grid(c,r)%clmbetak(1)=(grid(c,r)%clmbetak(1) - alpha)
     grid(c,r)%clmbetak(2)=(grid(c,r)%clmbetak(2) -   &
                           alpha*cos(omega1*time))
     grid(c,r)%clmbetak(3)=(grid(c,r)%clmbetak(3) -   &
                           alpha*sin(omega1*time))
     grid(c,r)%clmbetak(4)=(grid(c,r)%clmbetak(4) -   &
                           alpha*cos(omega2*time))
     grid(c,r)%clmbetak(5)=(grid(c,r)%clmbetak(5) -   &
                           alpha*sin(omega2*time))
    endif
   enddo
  enddo

  do r=1,ldas%nr
   do c=1,ldas%nc
    grid(c,r)%fbias=(grid(c,r)%clmbetak(1)+                     &
                     grid(c,r)%clmbetak(2)*cos(omega1*time)+    &
                     grid(c,r)%clmbetak(3)*sin(omega1*time)+    &
                     grid(c,r)%clmbetak(4)*cos(omega2*time)+    &
                     grid(c,r)%clmbetak(5)*sin(omega2*time))/86400.

    if(grid(c,r)%imask.eq.0)grid(c,r)%fbias=ldas%udef
   enddo
  enddo
      

!=== Transfer forecast bias to tile space
      call g2tr(fbiast,grid%fbias,ldas%nc,ldas%nr,ldas%nch,     &
                tile%col,tile%row)

!=== Add forecast bias term to analysis increment dtc 
      do t=1,ldas%nch
       clm1(t)%dtcanal=clm1(t)%dtcanal-float(ldas%ts)*fbiast(t)
      enddo

end subroutine clm1_bias
