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
! mos_bias.f: 
!
! DESCRIPTION:
!  This subroutine corrects the analysis increment with the forecast bias.  
!
! REVISION HISTORY:
!  7  March    2001: Jon Radakovich; Initial Code
!  20 November 2002: Jon Radakovich; Updated with fbias variable from grid module
!
!=========================================================================
      subroutine mos_bias(ldas,grid,tile,mos)

      use ldas_module      ! ldas non-model-specific 1-d variables
      use grid_module      ! ldas non-model-specific grid variables
      use tile_module      ! ldas non-model-specific tile variables
      use mos_module       ! mosaic tile variables
      implicit none
      type (ldasdec) :: ldas 
      type (griddec) :: grid(ldas%nc,ldas%nr)           
      type (tiledec) :: tile(ldas%nch)
      type (mosdec)  :: mos(ldas%nch)

!=== Local Variables =====================================================
      integer t,c,r        ! Loop counters
      real time            ! Model time of day in seconds
      real gamma           ! Tunable weight parameter for bias correction
      real omega1,omega2   ! Frequencies of harmonics
      real alpha           ! Amplitude of harmonic
      real fbiast(ldas%nch)       ! Forecast bias in tile space
      real pi
      parameter (pi = 3.14159265)

!=== End Variables =====================================================
      time=(float(ldas%hr))*3600.+
     &     (float(ldas%mn))*60.  +
     &     float(ldas%ss)

      gamma=0.2
      omega1=2.*pi/86400.
      omega2=2.*pi/43200.

      do r=1,ldas%nr
       do c=1,ldas%nc
        alpha=gamma*grid(c,r)%mosdelt
        if(ldas%ribc.eq.1)then
         grid(c,r)%mosbetak(1)=grid(c,r)%mosbetak(1) - alpha
         grid(c,r)%mosbetak(2)=0.
         grid(c,r)%mosbetak(3)=0.
         grid(c,r)%mosbetak(4)=0.
         grid(c,r)%mosbetak(5)=0.
        else if(ldas%rdbc.eq.1)then
         grid(c,r)%mosbetak(1)=grid(c,r)%mosbetak(1) - alpha
         grid(c,r)%mosbetak(2)=grid(c,r)%mosbetak(2) - 
     &                         alpha*cos(omega1*time)
         grid(c,r)%mosbetak(3)=grid(c,r)%mosbetak(3) - 
     &                         alpha*sin(omega1*time)
         grid(c,r)%mosbetak(4)=0.
         grid(c,r)%mosbetak(5)=0.
        else if(ldas%rsdbc.eq.1)then
         grid(c,r)%mosbetak(1)=grid(c,r)%mosbetak(1) - alpha
         grid(c,r)%mosbetak(2)=grid(c,r)%mosbetak(2) - 
     &                         alpha*cos(omega1*time)
         grid(c,r)%mosbetak(3)=grid(c,r)%mosbetak(3) - 
     &                         alpha*sin(omega1*time)
         grid(c,r)%mosbetak(4)=grid(c,r)%mosbetak(4) - 
     &                         alpha*cos(omega2*time)
         grid(c,r)%mosbetak(5)=grid(c,r)%mosbetak(5) - 
     &                         alpha*sin(omega2*time)
        endif
       enddo
      enddo

      do r=1,ldas%nr
       do c=1,ldas%nc
        grid(c,r)%fbias=(grid(c,r)%mosbetak(1)+
     &              grid(c,r)%mosbetak(2)*cos(omega1*time)+
     &              grid(c,r)%mosbetak(3)*sin(omega1*time)+
     &              grid(c,r)%mosbetak(4)*cos(omega2*time)+
     &              grid(c,r)%mosbetak(5)*sin(omega2*time))/86400.
        if(grid(c,r)%imask.eq.0)grid(c,r)%fbias=ldas%udef
       enddo
      enddo

C=== Transfer forecast bias to tile space
      call g2tr(fbiast,grid%fbias,ldas%nc,ldas%nr,ldas%nch,     
     &          tile%col,tile%row)


C=== Add forecast bias term to analysis increment dtc 
      do t=1,ldas%nch
       mos(t)%dtcanal=mos(t)%dtcanal-float(ldas%ts)*fbiast(t)
      enddo

      end subroutine mos_bias
