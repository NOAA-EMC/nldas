#include <misc.h>

subroutine drv_clmini (ldas, drv, grid, tile, clm1)

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
!  Setup initial CLM information in tile space.
!  
!  The observational initial values are the first choice.  If there are no
!  observational initial data, arbitrarily set the initial fields, and then
!  spin up the first model year to equilibrium state, and take the 
!  equilibrium state variables as the initial values.
!
!   The arbitrary initial data are created based on the following rules:
!    (1) Foliage temperature is initially set to lowest the atmospheric 
!        model air-temperature.
!    (2) Canopy water storage is set to zero.
!    (3) Soil temperatures are initialized as in bucket type
!        parameterizations using the lowest atmospheric model
!        air-temperature and a climatological deep-ground temperature.
!    (4) Soil moistures are initialized to a percentage of field capacity, 
!        and the percent of liquid water and ice lens are determined by the 
!        layer temperatures.
!    (5) If the depth of snow is known, then subdivide the snow pack 
!        up to five layers based on the following rules: 
!         From top layer to bottom layer
!         minimum thickness:     0.010, 0.015, 0.025, 0.055, 0.115 (m),
!         and maximum thickness: 0.02, 0.05, 0.11, 0.23, and  >0.23 m,
!    (6) The snow layer temperature is set to the surface air temperature.
!        If the air temperature is greater than freezing point, then it is 
!        set to 273.16 K.  If no information on snow is available, then 
!        the snow mass for areas of permanent land ice is initially set to 
!        50000 kg m-2. For other areas, all snow related variables 
!        are set to 0.
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December  1999: Paul Houser and Jon Radakovich; F90 Revision 
!   3 March     2000: Jon Radakovich; Revision for diagnostic output
!  28 January   2002: Jon Gottschalck; Added multiple sections to use model
!                     forcing to initialize land surface
!  28 March     2002: Brian Cosgrove; Added code to make use of Dag 
!                     Lohmann's initial conditions for NLDAS runs.
!                     His soil wetness (0=wilt,1=porosity) are transformed
!                     into a wetness factor (0=no mst,1=poros)
!  20 November  2002: Jon Radakovich; Modified for 5 soil and lake layers
!                     Modified calculation of xksat so it is based on an
!                     e-folding depth in order to reduce spin-up time.
!  03 June      2003: Jon Gottschalck, Modified initialization section for GEOS
!=========================================================================
! $Id: drv_clmini.F90,v 1.2 2003/06/03 21:23:36 jgottsch Exp $
!=========================================================================

  use precision
  use infnan
  use ldas_module         ! LDAS run variables
  use drv_module          ! 1-D Land Model Driver variables
  use drv_gridmodule      ! Grid-space variables
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! CLM tile variables 
  use clm1_varcon, only : denice, denh2o, istwet, istice, istsoil, istslak, istdlak
  use clm1_varpar, only : nlevsoi, nlevsno
  implicit none

!=== Arguments ===========================================================

  type (ldasdec)     :: ldas
  type (drvdec)      :: drv              
  type (clm_griddec) :: grid(drv%nc,drv%nr)   
  type (clm_tiledec) :: tile(drv%nch)
  type (clm11d)       :: clm1(drv%nch)

!=== Local Variables =====================================================

  integer  i, j, L, t, m         !loop indices
  integer  noaaic            ! 0=Use IC's from card file, 1=Use NOAA IC's from
                             ! Dag's data files
  real(r8) bd                !bulk density of dry soil material [kg/m^3]
  real(r8) dmvol             !fractional volume of dry soil material
  real(r8) tkm               !mineral conductivity
  real(r8) zlak(1:nlevlak)   !temporary z
  real(r8) dzlak(1:nlevlak)  !temporary dz
  real(r8) pi                !3.14159...
  real(r8) xksat             !maximum hydraulic conductivity of soil [mm/s]
  REAL(R8) WT1(drv%nch,1:nlevsoi),WT2(drv%nch,1:nlevsoi)
  REAL SOILMOIST(LDAS%NC,LDAS%NR)
  REAL SOILTEMP1(LDAS%NC,LDAS%NR),SOILTEMP2(LDAS%NC,LDAS%NR)
  real totaldepth(ldas%nch)
  real rootfr
  real swetwilt(ldas%nch)
  real avgwatsat(ldas%nch),avgrootsat(ldas%nch),temp
!	real aaa,bbb,ccc,soilwr(ldas%nch),mstav(464,224)
!	integer x,y
  PARAMETER (NOAAIC=0)

      ! Weights to be given for the GDAS soil layers
      ! when using the model initialization

!=== End Variable List ===================================================

 if (ldas%startcode == 3.and.noaaic.eq.1) then

               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'USING DAGS IC DATA, NOT CARD!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'USING DAGS IC DATA, NOT CARD!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

!       READ IN INITIAL SOIL TEMP AND SOIL WETNESS FROM NOAA's initial
!       conditions from Dag
        open (unit=22,file=    &
       './BCS/N0.125/1996093012.moist_scale.bin',  &
       form='unformatted')
        read(22) ((soilmoist(i,j),i=1,464),j=1,224)
        close (22)
        open (unit=22,file=  &
       './BCS/N0.125/1996093012.temp1.bin',  &
       form='unformatted')
        read(22)((soiltemp1(i,j),i=1,464),j=1,224)
        close(22)
        open (unit=22,file=  &
       './BCS/N0.125/1996093012.temp2.bin',  &
       form='unformatted')
        read(22)((soiltemp2(i,j),i=1,464),j=1,224)
        close(22)
        do t=1,drv%nch
        tile(t)%it1=SOILTEMP1(TILE(T)%COL,TILE(T)%ROW)
        tile(t)%it2=SOILTEMP2(TILE(T)%COL,TILE(T)%ROW)
        enddo

 	endif ! the noaaic=1 and startcode=3	

! ========================================================================
! Initialize Time Parameters - additional init done in drv_tick 
! ========================================================================
  do t=1,drv%nch
  clm1(t)%kpatch = t  

  clm1(t)%istep = 0
  clm1(t)%dtime = drv%ts

! ========================================================================
! TIME CONSTANT [1]
!  Define layer structure (here "0" refers to soil surface and 
! "nlevsoi" refers to the bottom of model soil)
! ========================================================================

  if (nlevlak /= nlevsoi) then
     write(6,*)'number of soil levels and number of lake levels must be the same'
     write(6,*)'nlevsoi= ',nlevsoi,' nlevlak= ',nlevlak
     stop
  endif

  if (clm1(t)%itypwat == istdlak) then        !assume all lakes are deep lakes

     if (nlevlak > 10) then
        write(6,*) ' must set new values for lake levels > 10'
        stop
     endif

     dzlak(1) = 1.               
     dzlak(2) = 2.               
     dzlak(3) = 3.               
     dzlak(4) = 4.               
     dzlak(5) = 5.
!     dzlak(6) = 7.
!     dzlak(7) = 7.
!     dzlak(8) = 7.
!     dzlak(9) = 7.
!     dzlak(10)= 7.

     zlak(1) =  0.5
     zlak(2) =  1.5
     zlak(3) =  4.5
     zlak(4) =  8.0
     zlak(5) = 12.5
!     zlak(6) = 18.5
!     zlak(7) = 25.5
!     zlak(8) = 32.5
!     zlak(9) = 39.5
!     zlak(10)= 46.5

     do j = 1,nlevlak
        clm1(t)%z(j) = zlak(j)
        clm1(t)%dz(j) = dzlak(j)
     end do
     
  else if (clm1(t)%itypwat == istslak) then   !shallow lake (not used)
     
     clm1(t)%dz(1:nlevlak) = NaN
     clm1(t)%z(1:nlevlak)  = NaN

  else                                    !soil, ice, wetland
     
     do j = 1, nlevsoi
        clm1(t)%z(j) = tile(t)%scalez*(exp(0.5*(j-0.5))-1.)     !node depths
     enddo

     
     clm1(t)%dz(1)  = 0.5*(clm1(t)%z(1)+clm1(t)%z(2))         !thickness b/n two interfaces
     do j = 2,nlevsoi-1
        clm1(t)%dz(j)= 0.5*(clm1(t)%z(j+1)-clm1(t)%z(j-1)) 
     enddo
     clm1(t)%dz(nlevsoi)= clm1(t)%z(nlevsoi)-clm1(t)%z(nlevsoi-1)


     clm1(t)%zi(0)   = 0.                             !interface depths 
     do j = 1, nlevsoi-1
        clm1(t)%zi(j)= 0.5*(clm1(t)%z(j)+clm1(t)%z(j+1))     
     enddo
     clm1(t)%zi(nlevsoi) = clm1(t)%z(nlevsoi) + 0.5*clm1(t)%dz(nlevsoi) 

  endif

!=== Setting weights for model initialization for soil model temp. and moisture
!=== GDAS has 2 levels so depending on where clm1 soil layers are located
!=== weight then accordingly (linear weights)

     do j = 1, nlevsoi
       if (j == 1) then
	 wt1(t,j) = (clm1(t)%zi(j) - 0.00) / clm1(t)%zi(j)
	 wt2(t,j) = (clm1(t)%zi(j) - 0.10) / clm1(t)%zi(j)
       else
	 if ((clm1(t)%zi(j) - clm1(t)%zi(j-1)) == 0.0) then
          wt1(t,j) = (0.10 - clm1(t)%zi(j-1)) / ((clm1(t)%zi(j) + 0.0001) &
&                     - clm1(t)%zi(j-1))
	  wt2(t,j) = (clm1(t)%zi(j) - 0.10)   / ((clm1(t)%zi(j) + 0.0001) &
&                     - clm1(t)%zi(j-1))
        else
	  wt1(t,j) = (0.10 - clm1(t)%zi(j-1)) / (clm1(t)%zi(j) - &
&                 clm1(t)%zi(j-1))
          wt2(t,j) = (clm1(t)%zi(j) - 0.10)   / (clm1(t)%zi(j) - &
&                 clm1(t)%zi(j-1))
        endif
      endif

      if (wt1(t,j) .lt. 0.0) wt1(t,j) = 0.0
      if (wt2(t,j) .lt. 0.0) wt2(t,j) = 0.0
      if (wt1(t,j) .gt. 1.0) wt1(t,j) = 1.0
      if (wt2(t,j) .gt. 1.0) wt2(t,j) = 1.0

    enddo

! ========================================================================
! TIME CONSTANT [2]
! Initialize root fraction (computing from surface, d is depth in meter):
! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with beta & d_obs
! given in Zeng et al. (1998).
! ========================================================================

  do j = 1, nlevsoi-1
     clm1(t)%rootfr(j) = .5*( exp(-tile(t)%roota*clm1(t)%zi(j-1))  &
                        + exp(-tile(t)%rootb*clm1(t)%zi(j-1))  &
                        - exp(-tile(t)%roota*clm1(t)%zi(j  ))  &
                        - exp(-tile(t)%rootb*clm1(t)%zi(j  )) )
  enddo
  clm1(t)%rootfr(nlevsoi)=.5*( exp(-tile(t)%roota*clm1(t)%zi(nlevsoi-1))&
                         + exp(-tile(t)%rootb*clm1(t)%zi(nlevsoi-1)))

! reset depth variables assigned by user in clm1in file 

  do l=1,nlevsoi
     if (grid(tile(t)%col,tile(t)%row)%rootfr /= drv%udef) &
          clm1(t)%rootfr(l)=grid(tile(t)%col,tile(t)%row)%rootfr    
  enddo

! ========================================================================
! TIME CONSTANT [3]
! Initialize soil thermal and hydraulic properties
! ========================================================================

! Define the vertical profile of soil thermal and hydraulic properties

  if (clm1(t)%itypwat == istsoil) then  ! soil

     do j = 1, nlevsoi
        clm1(t)%bsw(j)    = 2.91 + 0.159*tile(t)%clay(j)*100.0
        clm1(t)%watsat(j) = 0.489 - 0.00126*tile(t)%sand(j)*100.0
!        xksat         = 0.0070556 *( 10.**(-0.884+0.0153*tile(t)%sand(j)*100.0)) ! mm/s
         xksat         = 0.0070556 *( 10.**(-0.884+0.0153*tile(t)%sand(j)*100.0)) & ! mm/s
                                   *exp(1.0/tile(t)%hkdepth)                
        clm1(t)%hksat(j)  = xksat * exp(-clm1(t)%zi(j)/tile(t)%hkdepth)
        clm1(t)%sucsat(j) = 10. * ( 10.**(1.88-0.0131*tile(t)%sand(j)*100.0) )
        tkm           = (8.80*tile(t)%sand(j)*100.0+2.92*tile(t)%clay(j)*100.0) /  &
                        (tile(t)%sand(j)*100.0+tile(t)%clay(j)*100.0)          ! W/(m K)
        bd            = (1.-clm1(t)%watsat(j))*2.7e3
        clm1(t)%tkmg(j)   = tkm ** (1.- clm1(t)%watsat(j))
        clm1(t)%tksatu(j) = clm1(t)%tkmg(j)*0.57**clm1(t)%watsat(j)
        clm1(t)%tkdry(j)  = (0.135*bd + 64.7) / (2.7e3 - 0.947*bd)
        clm1(t)%csol(j)   = (2.128*tile(t)%sand(j)*100.0+2.385*tile(t)%clay(j)*100.0)/ &
                        (tile(t)%sand(j)*100.0+tile(t)%clay(j)*100.0)*1.e6     ! J/(m3 K)
     enddo

  else                                ! ice/glacier, lakes, wetlands

     do j = 1, nlevsoi
        clm1(t)%bsw(j)    = drv%udef
        clm1(t)%watsat(j) = 1.
        clm1(t)%hksat(j)  = drv%udef
        clm1(t)%sucsat(j) = drv%udef
        clm1(t)%tkmg(j)   = drv%udef
        clm1(t)%tksatu(j) = drv%udef
        clm1(t)%tkdry(j)  = drv%udef
        clm1(t)%csol(j)   = drv%udef
     enddo

  endif
 if (ldas%startcode == 3.and.noaaic.eq.1) then
!	Convert Dag's soil wetness (0=wilt pt, 1=porosity) to 
!       wetness (0 is totally dry, 1 is porosity) for use in initializing model
	do m=1,nlevsoi
         swetwilt(t)=0.
         avgwatsat(t)=0.
         swetwilt(t)=swetwilt(t) + clm1(t)%dz(m)*(clm1(t)%watsat(m)*   &
         ((-1)*clm1(t)%smpmax/clm1(t)%sucsat(m))**(-1/clm1(t)%bsw(m)))
	 swetwilt(t)=swetwilt(t)/clm1(t)%dz(m)

         tile(t)%iwet(m)= ((soilmoist(TILE(T)%COL,TILE(T)%ROW)*  &
         clm1(t)%watsat(m))  &
         -(soilmoist(TILE(T)%COL,TILE(T)%ROW)*clm1(t)%watsat(m)* &
          swetwilt(t))+  &
         (clm1(t)%watsat(m)*swetwilt(t)))/clm1(t)%watsat(m)
	enddo
  endif !endif the startcode=3 and noaaic=1	


! ========================================================================
! TIME VARIANT [1]
! Temperatures and snow cover fraction are initialized in clm1
! according to atmospheric temperatures and snow cover.
! ========================================================================

! set water and temperatures to constant values: all points
!=== Used in model initialization option

 if (ldas%startcode == 4) then
  clm1(t)%h2ocan  = 0.
  clm1(t)%snowage = 0.
  clm1(t)%h2osno  = clm1(t)%forc_sdepth
  clm1(t)%snowdp  = clm1(t)%forc_sdepth / 250.0    !the arbitary snow density = 250 kg/m3
  clm1(t)%t_veg   = clm1(t)%forc_t
  clm1(t)%t_grnd  = clm1(t)%forc_t
 else
  clm1(t)%h2ocan  = 0.
  clm1(t)%snowage = 0. 
  clm1(t)%h2osno  = drv%h2osno_ini	  
  clm1(t)%snowdp  = drv%h2osno_ini/250.  !the arbitary snow density = 250 kg/m3
   if (ldas%startcode == 3.and.noaaic.eq.1) then
    clm1(t)%t_veg   = tile(t)%it1
    clm1(t)%t_grnd  = tile(t)%it1
   else
    clm1(t)%t_veg   = drv%t_ini
    clm1(t)%t_grnd  = drv%t_ini
   endif
 endif

! For lake points only:

  if (clm1(t)%lakpoi) then  !lake 
     if (clm1(t)%t_grnd <= 273.16) then
        clm1(t)%h2osno = 0.
        clm1(t)%snowdp = 0.
     endif
  endif

  clm1(t)%acc_errseb = 0.
  clm1(t)%acc_errh2o = 0.

! ========================================================================
! TIME VARIANT [2]
! Snow layer number, depth and thickiness 
! ========================================================================

  if (.not. clm1(t)%lakpoi) then  !not lake
     if (clm1(t)%snowdp < 0.01)then
        clm1(t)%snl = 0
        clm1(t)%dz(-nlevsno+1:0) = 0.
        clm1(t)%z (-nlevsno+1:0) = 0.
        clm1(t)%zi(-nlevsno+0:0) = 0.
     else
        if ((clm1(t)%snowdp >= 0.01) .AND. (clm1(t)%snowdp <= 0.03))then
           clm1(t)%snl = -1
           clm1(t)%dz(0)  = clm1(t)%snowdp
        else if ((clm1(t)%snowdp > 0.03) .AND. (clm1(t)%snowdp <= 0.04))then
           clm1(t)%snl = -2
           clm1(t)%dz(-1) = clm1(t)%snowdp/2.
           clm1(t)%dz( 0) = clm1(t)%dz(-1)
        else if ((clm1(t)%snowdp > 0.04) .AND. (clm1(t)%snowdp <= 0.07))then
           clm1(t)%snl = -2
           clm1(t)%dz(-1) = 0.02
           clm1(t)%dz( 0) = clm1(t)%snowdp - clm1(t)%dz(-1)
        else if ((clm1(t)%snowdp > 0.07) .AND. (clm1(t)%snowdp <= 0.12))then
           clm1(t)%snl = -3
           clm1(t)%dz(-2) = 0.02
           clm1(t)%dz(-1) = (clm1(t)%snowdp - 0.02)/2.
           clm1(t)%dz( 0) = clm1(t)%dz(-1)
        else if ((clm1(t)%snowdp > 0.12) .AND. (clm1(t)%snowdp <= 0.18))then
           clm1(t)%snl = -3
           clm1(t)%dz(-2) = 0.02
           clm1(t)%dz(-1) = 0.05
           clm1(t)%dz( 0) = clm1(t)%snowdp - clm1(t)%dz(-2) - clm1(t)%dz(-1)
        else if ((clm1(t)%snowdp > 0.18) .AND. (clm1(t)%snowdp <= 0.29))then
           clm1(t)%snl = -4
           clm1(t)%dz(-3) = 0.02
           clm1(t)%dz(-2) = 0.05
           clm1(t)%dz(-1) = (clm1(t)%snowdp - clm1(t)%dz(-3) - clm1(t)%dz(-2))/2.
           clm1(t)%dz( 0) = clm1(t)%dz(-1)
        else if ((clm1(t)%snowdp > 0.29) .AND. (clm1(t)%snowdp <= 0.41))then
           clm1(t)%snl = -4
           clm1(t)%dz(-3) = 0.02
           clm1(t)%dz(-2) = 0.05
           clm1(t)%dz(-1) = 0.11
           clm1(t)%dz( 0) = clm1(t)%snowdp - clm1(t)%dz(-3) - clm1(t)%dz(-2) - clm1(t)%dz(-1)
        else if ((clm1(t)%snowdp > 0.41) .AND. (clm1(t)%snowdp <= 0.64))then
           clm1(t)%snl = -5
           clm1(t)%dz(-4) = 0.02
           clm1(t)%dz(-3) = 0.05
           clm1(t)%dz(-2) = 0.11
           clm1(t)%dz(-1) = (clm1(t)%snowdp - clm1(t)%dz(-4) - clm1(t)%dz(-3) - clm1(t)%dz(-2))/2.
           clm1(t)%dz( 0) = clm1(t)%dz(-1)
        else if (clm1(t)%snowdp > 0.64)then
           clm1(t)%snl = -5
           clm1(t)%dz(-4) = 0.02
           clm1(t)%dz(-3) = 0.05
           clm1(t)%dz(-2) = 0.11
           clm1(t)%dz(-1) = 0.23
           clm1(t)%dz( 0)=clm1(t)%snowdp-clm1(t)%dz(-4)-clm1(t)%dz(-3)-clm1(t)%dz(-2)-clm1(t)%dz(-1)
        endif
        do i = 0, clm1(t)%snl+1, -1
           clm1(t)%z(i) = clm1(t)%zi(i) - 0.5*clm1(t)%dz(i)
           clm1(t)%zi(i-1) = clm1(t)%zi(i) - clm1(t)%dz(i)
        enddo
     endif
  else   ! lake points
     clm1(t)%snl = 0
     clm1(t)%dz(-nlevsno+1:0) = 0.
     clm1(t)%z (-nlevsno+1:0) = 0.
     clm1(t)%zi(-nlevsno+0:0) = 0.
  endif

! ========================================================================
! TIME VARIANT [3]
! Snow/soil temperature (Modified to allow for model initialization [startcode=4])
! ========================================================================

  if (.not. clm1(t)%lakpoi) then  !not lake
     if (clm1(t)%snl < 0) then
        do i = clm1(t)%snl+1, 0
	  if (ldas%startcode == 4) then
	   if (clm1(t)%forc_t  < 273.16) then
	    clm1(t)%t_soisno(i) = clm1(t)%forc_t
	   else
	    clm1(t)%t_soisno(i) = 273.16 - 1.
	   endif
	  else
           if (ldas%startcode == 3.and.noaaic.eq.1) then
             if (tile(t)%it1  < 273.16) then
                clm1(t)%t_soisno(i) = tile(t)%it1
             else
                clm1(t)%t_soisno(i) = 273.16 - 1.
             endif
           else
             if (drv%t_ini  < 273.16) then
                clm1(t)%t_soisno(i) = drv%t_ini
             else
                clm1(t)%t_soisno(i) = 273.16 - 1.
             endif
           endif
	 endif
        enddo
     endif
     do i = 1, nlevsoi
       if (ldas%startcode == 4) then
         if (clm1(t)%itypwat == istice) then
	  clm1(t)%t_soisno(i) = clm1(t)%forc_t
	 else if (clm1(t)%itypwat == istwet) then
	  clm1(t)%t_soisno(i) = clm1(t)%forc_t
	 else
          select case (ldas%force)
          case(1)
             clm1(t)%t_soisno(i) = clm1(t)%forc_t
!!!          clm1(t)%t_soisno(i) = wt1(t,i) * clm1(t)%forc_stemp1 + &
!!!&                              wt2(t,i) * clm1(t)%forc_stemp2
          case(2)
          clm1(t)%t_soisno(i) = clm1(t)%forc_t
          case default
           print*, "Only can use GDAS, GEOS for initialization"
           stop
          end select
         endif
       else
       if (ldas%startcode == 3.and.noaaic.eq.1) then
        if (clm1(t)%itypwat == istice) then
           if (i.eq.1.or.i.eq.2.or.i.eq.3) then
             clm1(t)%t_soisno(i) = tile(t)%it1
           else
             clm1(t)%t_soisno(i) = tile(t)%it2
           endif
        else if (clm1(t)%itypwat == istwet) then
           if (i.eq.1.or.i.eq.2.or.i.eq.3) then
             clm1(t)%t_soisno(i) = tile(t)%it1
           else
             clm1(t)%t_soisno(i) = tile(t)%it2
           endif
        else
           if (i.eq.1.or.i.eq.2.or.i.eq.3) then
             clm1(t)%t_soisno(i) = tile(t)%it1
           else
             clm1(t)%t_soisno(i) = tile(t)%it2
           endif
        endif
       else
        if (clm1(t)%itypwat == istice) then
           clm1(t)%t_soisno(i) = drv%t_ini
        else if (clm1(t)%itypwat == istwet) then
           clm1(t)%t_soisno(i) = drv%t_ini
        else
           clm1(t)%t_soisno(i) = drv%t_ini
        endif
       endif 
      endif	
    enddo
  else
     do i = 1, nlevlak
       if (ldas%startcode == 3.and.noaaic.eq.1) then
        if (i.eq.1.or.i.eq.2.or.i.eq.3) then
          clm1(t)%t_soisno(i) = tile(t)%it1
        else
          clm1(t)%t_soisno(i) = tile(t)%it2
        endif
	if (ldas%startcode == 4) clm1(t)%t_soisno(i) = clm1(t)%forc_t
       else
        clm1(t)%t_soisno(i) = drv%t_ini
        if (ldas%startcode == 4) clm1(t)%t_soisno(i) = clm1(t)%forc_t
       endif
     enddo
  endif

! ========================================================================
! TIME VARIANT [4]
! Snow/soil ice and liquid mass (Modified to allow for model initialization [startcode=4])
! ========================================================================


  if (.not. clm1(t)%lakpoi) then  !not lake
     if (clm1(t)%snl < 0)then
        do i = clm1(t)%snl+1, 0
           clm1(t)%h2osoi_ice(i) = clm1(t)%dz(i)*250.
           clm1(t)%h2osoi_liq(i) = 0.
        enddo
     endif
     do i = 1, nlevsoi
        if (clm1(t)%t_soisno(i) <= 273.16) then
	 if (ldas%startcode == 4) then
           select case (ldas%force)
           case(1)
           clm1(t)%h2osoi_ice(i) = clm1(t)%dz(i) * &
&          (wt1(t,i)*clm1(t)%forc_swc1/clm1(t)%watsat(i)  &
&          +wt2(t,i)*clm1(t)%forc_swc2/clm1(t)%watsat(i)) &
&          *clm1(t)%watsat(i)*denice
           case(2)
           clm1(t)%h2osoi_ice(i) = clm1(t)%dz(i)*clm1(t)%forc_swc1*clm1(t)%watsat(i)*denice
           case default
            print*, "Only can use GDAS, GEOS initialization"
            stop
           end select
	 else
          if (ldas%startcode == 3.and.noaaic.eq.1) then
          clm1(t)%h2osoi_ice(i) = clm1(t)%dz(i)* tile(t)%iwet(i)*clm1(t)%watsat(i)*denice
          else
           clm1(t)%h2osoi_ice(i) = clm1(t)%dz(i)* drv%sw_ini*clm1(t)%watsat(i)*denice
          endif
	 endif	   
           clm1(t)%h2osoi_liq(i) = 0.
           if (clm1(t)%itypwat==istwet .or. clm1(t)%itypwat==istice) clm1(t)%h2osoi_ice(i)=clm1(t)%dz(i)*denice
        else
          clm1(t)%h2osoi_ice(i) = 0.
	  if (ldas%startcode == 4) then
          select case (ldas%force)
          case(1)
            clm1(t)%h2osoi_liq(i) = clm1(t)%dz(i)* &
&           (wt1(t,i)*clm1(t)%forc_swc1/clm1(t)%watsat(i)  &
&           +wt2(t,i)*clm1(t)%forc_swc2/clm1(t)%watsat(i)) &
&           *clm1(t)%watsat(i)*denh2o
          case(2)
            clm1(t)%h2osoi_liq(i) = clm1(t)%dz(i)*clm1(t)%forc_swc1*clm1(t)%watsat(i)*denh2o
          case default
           print*, "Only can use GDAS, GEOS initialization"
           stop
          end select
          else
          if (ldas%startcode == 3.and.noaaic.eq.1) then
           clm1(t)%h2osoi_liq(i) = clm1(t)%dz(i)*tile(t)%iwet(i)*clm1(t)%watsat(i)*denh2o

          else
            clm1(t)%h2osoi_liq(i) = clm1(t)%dz(i)* drv%sw_ini*clm1(t)%watsat(i)*denh2o
	  endif
	  endif
          if (clm1(t)%itypwat==istwet .or. clm1(t)%itypwat==istice) clm1(t)%h2osoi_liq(i)=clm1(t)%dz(i)*denh2o
        endif
     enddo
  else    !not used for lake
     do i = -nlevsno+1, nlevlak
        clm1(t)%h2osoi_liq(i) = NaN
        clm1(t)%h2osoi_ice(i) = NaN
     end do
  endif



! ========================================================================
! TIME VARIANT [5]
! need to set h2osoi_vol (needed by clm_soilalb) -  this is also needed
! upon restart since the albedo calculation is called before h2osoi_vol
! is computed
! ========================================================================

  do l = 1,nlevsoi
     if (.not. clm1(t)%lakpoi) then
        clm1(t)%h2osoi_vol(l) = clm1(t)%h2osoi_liq(l)/(clm1(t)%dz(l)*denh2o) &
                          + clm1(t)%h2osoi_ice(l)/(clm1(t)%dz(l)*denice)
     else
        clm1(t)%h2osoi_vol(l) = 1.0
     endif
  end do


! ========================================================================
! TIME VARIANT [6]
! Ecosystem dynamics: phenology, vegetation, soil carbon and snow cover 
! fraction
! ========================================================================

  call clm1_dynvegpar (clm1(t))

! ========================================================================
! TIME VARIANT [7]
! Initialize DIAG arrays 
! ========================================================================

  do i = 1, drv%surfind
     clm1(t)%diagsurf(i) = 0.
  enddo

  do i = 1, drv%soilind
     do j = 1, nlevsoi
        clm1(t)%diagsoil(i,j) = 0.
     enddo
  enddo

  do i = 1, drv%snowind
     do j = -nlevsno+1,0
        clm1(t)%diagsnow(i,j) = 0.
     enddo
  enddo

  enddo

!        do T=1,LDAS%NCH
!	MSTAV(TILE(T)%COL,TILE(T)%ROW) =0.0
!	enddo

!	do T=1,LDAS%NCH
!        swetwilt(t)=clm1(t)%dz(m)*(clm1(t)%watsat(m)* &
!         ((-1)*clm1(t)%smpmax/clm1(t)%sucsat(m))**(-1/clm1(t)%bsw(m)))
!        swetwilt(t)=swetwilt(t)/clm1(t)%dz(m)
!!     temp is volumetric
!	temp=soilmoist(TILE(T)%COL,TILE(T)%ROW)* &
!        (clm1(t)%watsat(m)-clm1(t)%watsat(m)*swetwilt(t))+ &
!        (clm1(t)%watsat(m)*swetwilt(t))
!
!	aaa=((clm1(t)%h2osoi_liq(1)+clm1(t)%h2osoi_ice(1))/(clm1(t)%dz(1)*  &
!       1000.0*clm1(t)%watsat(1)))*clm1(t)%watsat(1)
!	bbb=(clm1(t)%watsat(1))*swetwilt(t)
!	ccc=clm1(t)%watsat(1)-(clm1(t)%watsat(1)*swetwilt(t))
!       soilwr(t)=(aaa-bbb)/ccc
!
!        MSTAV(TILE(T)%COL,TILE(T)%ROW) =  & 
!       MSTAV(TILE(T)%COL,TILE(T)%ROW)+soilwr(t)*   &
!         TILE(T)%FGRD
!	ENDDO
!        OPEN(UNIT=68,FILE='test2.bin',FORM='UNFORMATTED')
!        write(68) ((MSTAV(x,y),X=1,464),Y=1,224)
!	close(68)

  return
end subroutine drv_clmini

