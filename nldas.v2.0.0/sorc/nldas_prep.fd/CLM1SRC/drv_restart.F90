#include <misc.h>

subroutine drv_restart (rw, drv, ldas, grid, tile, clm1)

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
!  This program reads and writes restart files for CLM.  This
!   includes all relevant water/energy storages, tile information,
!   and time information.  It also rectifies changes in the tile space.  
!
! REVISION HISTORY:
!  22  Oct 1999: Jon Radakovich and Paul Houser; Initial code
!  28  Jan 2002: Jon Gottschalck; Added model initialization section (startcode =4)
!  09  Apr 2002: Brian Cosgrove; Changed code to match new output directory structure
!  20  Nov 2002: Jon Radakovich; Updated for PSAS temperature assimilation and BC so
!                the forecast bias is included in the restart file.  Added conditionals
!                based on startcode=5 (Using spun-up restart so model time comes
!                from card file) and startcode=6 (Restarting a bias correction run).
!  14  Jan 2003: Urszula Jambor; Added deallocation statements near end of routine
!                and changed pointer variables to allocatable.
!=========================================================================
! RESTART FILE FORMAT(fortran sequential binary):
!  yr,mo,da,hr,mn,ss,vclass,nch !Restart time,Veg class,no.tiles, no.soil lay 
!  tile(nch)%col        !Grid Col of Tile   
!  tile(nch)%row        !Grid Row of Tile
!  tile(nch)%fgrd       !Fraction of Grid covered by tile
!  tile(nch)%vegt       !Vegetation Type of Tile
!  clm1(nch)%states      !Model States in Tile Space
!=========================================================================

  use precision
  use drv_module          ! 1-D Land Model Driver variables
  use ldas_module         ! LDAS non-model-specific 1-D variables
  use grid_module         ! LDAS non-model-specific Grid variables
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! 1-D CLM variables
  use clm1_varpar, only : nlevsoi, nlevsno
  use clm1_varcon, only : denh2o, denice
  implicit none

!=== Arguments ===========================================================  

  integer, intent(in) :: rw   ! 1=read restart, 2=write restart
  type (drvdec)       :: drv              
  type (ldasdec)      :: ldas
  type (griddec)      :: grid(ldas%nc,ldas%nr)
  type (clm_tiledec)  :: tile(drv%nch)
  type (clm11d)        :: clm1 (drv%nch)

!=== Local Variables =====================================================

  integer :: c,r,t,l,n,i       ! Loop counters
  integer :: found         ! Counting variable

!=== Temporary tile space transfer files (different than in DRV_module)

  integer :: yr,mo,da,hr,mn,ss        ! Time variables
  integer :: vclass,nc,nr,nch
  integer, allocatable  :: col(:)         ! Column
  integer, allocatable  :: row(:)         ! Row
  integer, allocatable  :: vegt(:)        ! Tile veg type
  real(r8), allocatable :: fgrd(:)        ! Grid Fraction of Tile
  real(r8), allocatable :: t_grnd(:)      ! CLM Soil Surface Temperature [K]
  real(r8), allocatable :: t_veg(:)       ! CLM Leaf Temperature [K]
  real(r8), allocatable :: h2osno(:)      ! CLM Snow Cover, Water Equivalent [mm]
  real(r8), allocatable :: snowage(:)     ! CLM Non-dimensional snow age [-] 
  real(r8), allocatable :: snowdp(:)      ! CLM Snow Depth [m] 
  real(r8), allocatable :: h2ocan(:)      ! CLM Depth of Water on Foliage [mm]

  real(r8), allocatable :: frac_sno(:)            ! CLM Fractional Snow Cover [-]
  real(r8), allocatable :: elai(:)                ! CLM Leaf Area Index
  real(r8), allocatable :: esai(:)                ! CLM Stem Area Index

  integer, allocatable  :: snl(:)       ! CLM Actual number of snow layers
  integer           :: istep        ! number of time step
  real(r8), allocatable :: xerr(:)      ! accumulation of water balance error
  real(r8), allocatable :: zerr(:)      ! accumulation of energy balnce error

  real(r8), allocatable :: dz(:,:)           ! CLM Layer Depth [m]
  real(r8), allocatable :: z(:,:)            ! CLM Layer Thickness [m]
  real(r8), allocatable :: zi(:,:)           ! CLM Interface Level Below a "z" Level [m]
  real(r8), allocatable :: t_soisno(:,:)     ! CLM Soil + Snow Layer Temperature [K]
  real(r8), allocatable :: h2osoi_liq(:,:)   ! CLM Average Soil Water Content [kg/m2]
  real(r8), allocatable :: h2osoi_ice(:,:)   ! CLM Average Ice Content [kg/m2]

  real(r8), allocatable :: tmptileot(:) ! Temporary Transfer Array   
  real(r8), allocatable :: tmptileow(:) ! Temporary Transfer Array 
  real(r8), allocatable :: tmptileoi(:) ! Temporary Transfer Array
  real(r8), allocatable :: tmptileor(:) ! Temporary Transfer Array
  real(r8), allocatable :: tmptileoa(:) ! Temporary Transfer Array
  real(r8) :: tmptilent(drv%nch)    ! Temporary Transfer Array   
  real(r8) :: tmptilenw(drv%nch)    ! Temporary Transfer Array
  real(r8) :: tmptileni(drv%nch)    ! Temporary Transfer Array
  real(r8) :: tmptilenr(drv%nch)    ! Temporary Transfer Array
  real(r8) :: tmptilena(drv%nch)    ! Temporary Transfer Array

  real(r8) :: g_t_grnd(drv%nc,drv%nr)         ! CLM Soil Surface Temperature [K]
  real(r8) :: g_t_veg(drv%nc,drv%nr)          ! CLM Leaf Temperature [K] 
  real(r8) :: g_h2osno(drv%nc,drv%nr)         ! CLM Snow Cover, Water Equivalent [mm] 
  real(r8) :: g_snowage(drv%nc,drv%nr)        ! CLM Non-dimensional snow age [-] 
  real(r8) :: g_snowdp(drv%nc,drv%nr)         ! CLM Snow Depth [m] 
  real(r8) :: g_h2ocan(drv%nc,drv%nr)         ! CLM Depth of Water on Foliage [mm]
  real(r8) :: g_frac_sno(drv%nc,drv%nr)       ! CLM Fractional Snow Cover [-]
  real(r8) :: g_elai(drv%nc,drv%nr)           ! CLM Leaf + Stem Area Index
  real(r8) :: g_esai(drv%nc,drv%nr)           ! CLM Leaf + Stem Area Index

  real(r8) :: g_dz(drv%nc,drv%nr,-nlevsno+1:nlevsoi)    ! CLM Layer Depth [m]
  real(r8) :: g_z (drv%nc,drv%nr,-nlevsno+1:nlevsoi)    ! CLM Layer Thickness [m]
  real(r8) :: g_zi(drv%nc,drv%nr,-nlevsno:nlevsoi)      ! CLM Interface Level Below a "z" Level [m]
  real(r8) :: g_t_soisno  (drv%nc,drv%nr,nlevsoi)       ! CLM Soil + Snow Layer Temperature [K]
  real(r8) :: g_h2osoi_liq(drv%nc,drv%nr,nlevsoi)       ! CLM Average Soil Water Content [kg/m2]
  real(r8) :: g_h2osoi_ice(drv%nc,drv%nr,-nlevsno+1:nlevsoi) ! CLM Average Ice Content [kg/m2]
  real(r8) :: g_fbias(drv%nc,drv%nr)                    ! CLM skin temperature bias

  character*80 :: filen,mkfyrmo
  character*1  :: fname(80),fbase(80),fsubs(80),fmkdir(80)
  character*1  :: ftime(10),fyrmodir(80)



!=== End Variable Definition =============================================

!=== Read Active Archive File ============================================

  if((rw.eq.1.and.drv%startcode.eq.1).or.  &
     (rw.eq.1.and.drv%startcode.eq.5.and.  &
      ldas%tscount.eq.0).or.  &
     (rw.eq.1.and.drv%startcode.eq.6))then

     open(40,file=drv%rstf,form='unformatted')

     read(40) yr,mo,da,hr,mn,ss,vclass,nc,nr,nch  !Time, veg class, no. tiles

     allocate (col(nch),row(nch),fgrd(nch),vegt(nch))
     allocate (t_grnd(nch),t_veg(nch),h2osno(nch),snowage(nch),         &
               snowdp(nch),h2ocan(nch),frac_sno(nch))
     allocate (elai(nch), esai(nch), snl(nch),xerr(nch),zerr(nch))
     allocate (dz(nch,-nlevsno+1:nlevsoi),        &
               z(nch,-nlevsno+1:nlevsoi),         &
               zi(nch,-nlevsno:nlevsoi),          &           
               t_soisno(nch,-nlevsno+1:nlevsoi),  &
               tmptileot(nch),                    &
               h2osoi_liq(nch,-nlevsno+1:nlevsoi),&
               tmptileow(nch),                    &
               h2osoi_ice(nch,-nlevsno+1:nlevsoi),&
               tmptileoi(nch),tmptileoa(nch))

     read(40) col                  !Grid Col of Tile   
     read(40) row                  !Grid Row of Tile
     read(40) fgrd                 !Fraction of Grid covered by tile
     read(40) vegt                 !Vegetation Type of Tile
     read(40) t_grnd               !CLM Soil Surface Temperature [K] 
     read(40) t_veg                !CLM Leaf Temperature [K] 
     read(40) h2osno               !CLM Snow Cover, Water Equivalent [mm] 
     read(40) snowage              !CLM Non-dimensional snow age [-] 
     read(40) snowdp               !CLM Snow Depth [m]
     read(40) h2ocan               !CLM Depth of Water on Foliage [mm]
     read(40) frac_sno             !CLM Fractional Snow Cover [-]
     read(40) elai                 !CLM Leaf Area Index
     read(40) esai                 !CLM Stem Area Index
     read(40) snl                  !CLM Actual number of snow layers
     read(40) xerr                 !CLM Accumulation of water balance error
     read(40) zerr                 !CLM Accumulation of energy balnce error
     read(40) istep                !CLM Number of time step
     if(ldas%rbias.eq.1.and.ldas%startcode.eq.6)then
      read(40) g_fbias             !CLM forecast bias
      read(40) grid%clmdelt        !CLM analysis increment
      do r=1,ldas%nr
       do c=1,ldas%nc
        grid(c,r)%clmbetak(1)=g_fbias(c,r)
       enddo
      enddo
     endif

     do l = -nlevsno+1,nlevsoi
        read(40) tmptileoa  !CLM Layer Depth [m]
        do t = 1,drv%nch
           dz(t,l) = tmptileoa(t) 
        enddo
     enddo
     do l = -nlevsno+1,nlevsoi
        read(40) tmptileoa  !CLM Layer Thickness [m]
        do t = 1,drv%nch
           z(t,l) = tmptileoa(t) 
        enddo
     enddo
     do l = -nlevsno,nlevsoi
        read(40) tmptileoa  !CLM Interface Level Below a "z" Level [m]
        do t = 1,drv%nch
           zi(t,l) = tmptileoa(t) 
        enddo
     enddo

     do l = -nlevsno+1,nlevsoi
        read(40) tmptileot  !CLM Soil + Snow Layer Temperature [K]
        do t = 1,drv%nch
           t_soisno(t,l) = tmptileot(t) 
        enddo
     enddo
     do l = -nlevsno+1,nlevsoi
        read(40) tmptileow  !Average Soil Water Content [kg/m2]
        do t = 1,drv%nch
           h2osoi_liq(t,l) = tmptileow(t)
        enddo
     enddo
     do l = -nlevsno+1,nlevsoi
        read(40) tmptileoi  !CLM Average Ice Content [kg/m2]
        do t = 1,drv%nch
           h2osoi_ice(t,l) = tmptileoi(t)
        enddo
     enddo

     close(40)
     write(*,*)'CLM Restart File Read: ',drv%rstf
     write(79,*)'CLM Restart File Read: ',drv%rstf
!=== Establish Model Restart Time  

     if(drv%startcode.eq.1.or.drv%startcode.eq.6)then
        drv%yr = yr
        drv%mo = mo 
        drv%da = da
        drv%hr = hr
        drv%mn = mn
        drv%ss = ss
        call drv_date2time(drv%time,drv%doy,drv%day,drv%gmt,yr,mo,da,hr,mn,ss) 
        drv%ctime = drv%time
        write(*,*)'CLM Restart File Time Used: ',drv%rstf
	write(79,*)'CLM Restart File Time Used: ',drv%rstf
	
        ldas%yr=drv%yr
        ldas%mo=drv%mo
        ldas%da=drv%da
        ldas%hr=drv%hr
        ldas%mn=drv%mn
        ldas%ss=drv%ss
        ldas%time=drv%time
     endif

!=== Using spunup IC file, then don't use timestamp in restart    
     if(rw.eq.1.and.drv%startcode.eq.5)then 
        drv%yr = drv%syr
        drv%mo = drv%smo 
        drv%da = drv%sda
        drv%hr = drv%shr
        drv%mn = drv%smn
        drv%ss = drv%sss
        call drv_date2time(drv%time,drv%doy,drv%day,drv%gmt, &
             drv%yr,drv%mo,drv%da,drv%hr,drv%mn,drv%ss) 
        write(*,*)'Using ldas.crd start time for spun-up IC ',drv%time
        write(79,*)'Using ldas.crd start time for spun-up IC ',drv%time
        ldas%yr=drv%yr
        ldas%mo=drv%mo
        ldas%da=drv%da
        ldas%hr=drv%hr
        ldas%mn=drv%mn
        ldas%ss=drv%ss
        ldas%time=drv%time
     endif

!=== Rectify Restart Tile Space to DRV Tile Space =====================

     if(drv%clm_ic.eq.1)then

! Check for Vegetation Class Conflict 

        if(vclass.ne.drv%vclass)then
           write(*,*)drv%rstf,' Vegetation class conflict - clm1 HALTED'
	   write(79,*)drv%rstf,' Vegetation class conflict - clm1 HALTED'
           stop
        endif

! Check for Grid Space Conflict 

        if(nc.ne.drv%nc.or.nr.ne.drv%nr)then
           write(*,*)drv%rstf,'Grid space mismatch - clm1 HALTED'
           write(79,*)drv%rstf,'Grid space mismatch - clm1 HALTED'
           stop
        endif

! Transfer Restart tile space to DRV tile space

        if(nch.ne.drv%nch)then
           write(*,*)'Restart Tile Space Mismatch-Transfer in Progress'
           write(79,*)'Restart Tile Space Mismatch-Transfer in Progress'

!  Start by finding grid averages

           call drv_t2gr(t_grnd         ,g_t_grnd         ,drv%nc,drv%nr,nch,fgrd,col,row)
           call drv_t2gr(t_veg          ,g_t_veg          ,drv%nc,drv%nr,nch,fgrd,col,row)
           call drv_t2gr(h2osno         ,g_h2osno         ,drv%nc,drv%nr,nch,fgrd,col,row)
           call drv_t2gr(snowage        ,g_snowage        ,drv%nc,drv%nr,nch,fgrd,col,row)
           call drv_t2gr(snowdp         ,g_snowdp         ,drv%nc,drv%nr,nch,fgrd,col,row)
           call drv_t2gr(h2ocan         ,g_h2ocan         ,drv%nc,drv%nr,nch,fgrd,col,row)
           call drv_t2gr(frac_sno       ,g_frac_sno       ,drv%nc,drv%nr,nch,fgrd,col,row)
           call drv_t2gr(elai           ,g_elai           ,drv%nc,drv%nr,nch,fgrd,col,row)
           call drv_t2gr(esai           ,g_esai           ,drv%nc,drv%nr,nch,fgrd,col,row)

           do l = -nlevsno+1,nlevsoi
              call drv_t2gr(dz(:,l),g_dz(:,:,l),drv%nc,drv%nr,nch,fgrd,col,row)
           enddo
           do l = -nlevsno+1,nlevsoi
              call drv_t2gr(z(:,l),g_z(:,:,l),drv%nc,drv%nr,nch,fgrd,col,row)
           enddo
           do l = -nlevsno,nlevsoi
              call drv_t2gr(zi(:,l),g_zi(:,:,l),drv%nc,drv%nr,nch,fgrd,col,row)
           enddo

           do l = -nlevsno+1,nlevsoi
              call drv_t2gr(t_soisno(:,l),g_t_soisno(:,:,l),drv%nc,drv%nr,nch,fgrd,col,row)
           enddo

           do l = -nlevsno+1,nlevsoi
              call drv_t2gr(h2osoi_liq(:,l),g_h2osoi_liq(:,:,l),drv%nc,drv%nr,nch,fgrd,col,row)
           enddo

           do l = -nlevsno+1,nlevsoi
              call drv_t2gr(h2osoi_ice(:,l),g_h2osoi_ice(:,:,l),drv%nc,drv%nr,nch,fgrd,col,row)
           enddo

! Perform state transfer

           c = 0
           do 555 t = 1,drv%nch 

              if(amod(float(t),10000.0).eq.0.0)write(*,23)'  Transferred ', &
                   100.0*float(t)/float(drv%nch),' Percent of Tiles'

23            format(a14,f5.2,a17)
              found = 0
              do n = 1,nch
                 if ( tile(t)%vegt.eq.vegt(n) .and.   &
                      tile(t)%col .eq. col(n) .and.   &
                      tile(t)%row .eq. row(n) )then
                    clm1(t)%t_grnd = t_grnd(n)
                    clm1(t)%t_veg = t_veg(n)
                    clm1(t)%h2osno = h2osno(n)
                    clm1(t)%snowage = snowage(n)
                    clm1(t)%snowdp = snowdp(n)
                    clm1(t)%h2ocan = h2ocan(n)
                    clm1(t)%frac_sno = frac_sno(n)
                    clm1(t)%elai = elai(n)
                    clm1(t)%esai = esai(n)

                    do l = -nlevsno+1,nlevsoi
                       clm1(t)%dz(l) = dz(n,l)
                    enddo
                    do l = -nlevsno+1,nlevsoi
                       clm1(t)%z(l) = z(n,l)
                    enddo
                    do l = -nlevsno,nlevsoi
                       clm1(t)%zi(l) = zi(n,l)
                    enddo
                    do l = -nlevsno+1,nlevsoi
                       clm1(t)%t_soisno(l) = t_soisno(n,l)
                    enddo
                    do l = -nlevsno+1,nlevsoi
                       clm1(t)%h2osoi_liq(l) = h2osoi_liq(n,l)
                    enddo
                    do l = -nlevsno+1,nlevsoi
                       clm1(t)%h2osoi_ice(l) = h2osoi_ice(n,l)
                    enddo

                    found = 1
                    goto 555 
                 endif
              enddo

              if(found.eq.0)then        
                 clm1(t)%t_grnd = g_t_grnd(tile(t)%col,tile(t)%row)
                 clm1(t)%t_veg = g_t_veg(tile(t)%col,tile(t)%row)
                 clm1(t)%h2osno = g_h2osno(tile(t)%col,tile(t)%row)
                 clm1(t)%snowage = g_snowage(tile(t)%col,tile(t)%row)
                 clm1(t)%snowdp = g_snowdp(tile(t)%col,tile(t)%row)
                 clm1(t)%h2ocan = g_h2ocan(tile(t)%col,tile(t)%row)
                 clm1(t)%frac_sno = g_frac_sno(tile(t)%col,tile(t)%row)
                 clm1(t)%elai = g_elai(tile(t)%col,tile(t)%row)
                 clm1(t)%esai = g_esai(tile(t)%col,tile(t)%row)
                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%dz(l) = g_dz(tile(t)%col,tile(t)%row,l)
                 enddo
                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%z(l) = g_z(tile(t)%col,tile(t)%row,l)
                 enddo
                 do l = -nlevsno,nlevsoi
                    clm1(t)%zi(l) = g_zi(tile(t)%col,tile(t)%row,l)
                 enddo
                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%t_soisno(l) = g_t_soisno(tile(t)%col,tile(t)%row,l)
                 enddo
                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%h2osoi_liq(l) = g_h2osoi_liq(tile(t)%col,tile(t)%row,l)
                 enddo
                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%h2osoi_ice(l) = g_h2osoi_ice(tile(t)%col,tile(t)%row,l)
                 enddo
                 c = 0
              endif
555           continue
              write(*,*)'Tile Space Transfer Complete'
              write(*,*)'clm1 Restart NCH:',nch,'Current NCH:',drv%nch
              write(*,*) c, ' Tiles not found in old clm1 restart'        
              write(*,*)
	      write(79,*)'Tile Space Transfer Complete'
              write(79,*)'clm1 Restart NCH:',nch,'Current NCH:',drv%nch
              write(79,*) c, ' Tiles not found in old clm1 restart'        
              write(79,*)

           else  !The number of tiles is a match

              clm1%istep=istep

              do t = 1,drv%nch
                 clm1(t)%t_grnd = t_grnd(t)
                 clm1(t)%t_veg = t_veg(t)
                 clm1(t)%h2osno = h2osno(t)
                 clm1(t)%snowage = snowage(t)
                 clm1(t)%snowdp = snowdp(t)
                 clm1(t)%h2ocan = h2ocan(t)
                 clm1(t)%frac_sno = frac_sno(t)
                 clm1(t)%elai = elai(t)
                 clm1(t)%esai = esai(t)
                 clm1(t)%snl=snl(t)
                 clm1(t)%acc_errh2o=xerr(t)
                 clm1(t)%acc_errseb=zerr(t)

                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%dz(l) = dz(t,l)
                 enddo
                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%z(l) = z(t,l)
                 enddo
                 do l = -nlevsno,nlevsoi
                    clm1(t)%zi(l) = zi(t,l)
                 enddo

                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%t_soisno(l) = t_soisno(t,l)
                 enddo
                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%h2osoi_liq(l) = h2osoi_liq(t,l)
                 enddo
                 do l = -nlevsno+1,nlevsoi
                    clm1(t)%h2osoi_ice(l) = h2osoi_ice(t,l)
                 enddo
              enddo
           endif

        endif

! Determine h2osoi_vol(1) - needed for soil albedo calculation

        do t = 1,drv%nch
           clm1(t)%h2osoi_vol(1) = clm1(t)%h2osoi_liq(1)/(clm1(t)%dz(1)*denh2o) &
                                + clm1(t)%h2osoi_ice(1)/(clm1(t)%dz(1)*denice)
        end do

     endif !RW option 1

! === Set starttime to clm1in Stime when STARTCODE = 2 

     if(rw.eq.1.and.drv%startcode.eq.3)then 
        drv%yr = drv%syr
        drv%mo = drv%smo 
        drv%da = drv%sda
        drv%hr = drv%shr
        drv%mn = drv%smn
        drv%ss = drv%sss
        call drv_date2time(drv%time,drv%doy,drv%day,drv%gmt, &
             drv%yr,drv%mo,drv%da,drv%hr,drv%mn,drv%ss) 
        write(*,*)'Using ldas.crd start time for cold IC ',drv%time
	write(79,*)'Using ldas.crd start time for cold IC ',drv%time
        ldas%yr=drv%yr
        ldas%mo=drv%mo
        ldas%da=drv%da
        ldas%hr=drv%hr
        ldas%mn=drv%mn
        ldas%ss=drv%ss
        ldas%time=drv%time
     endif
     
!=== Set time to card file when doing the initialization.
     if(rw.eq.1.and.drv%startcode.eq.4)then
       drv%yr = drv%syr
       drv%mo = drv%smo
       drv%da = drv%sda
       drv%hr = drv%shr
       drv%mn = drv%smn
       drv%ss = drv%sss
       call drv_date2time(drv%time,drv%doy,drv%day,drv%gmt, &
            drv%yr,drv%mo,drv%da,drv%hr,drv%mn,drv%ss)
       write(*,*)'Using ldas.crd start time for cold IC ',drv%time
       write(79,*)'Using ldas.crd start time for cold IC ',drv%time

       ldas%yr=drv%yr
       ldas%mo=drv%mo
       ldas%da=drv%da
       ldas%hr=drv%hr
       ldas%mn=drv%mn
       ldas%ss=drv%ss
       ldas%time=drv%time
   endif

     if(rw.eq.1)then
        write(*,*)'clm1 Start Time: ',drv%yr,drv%mo,drv%da,drv%hr,drv%mn,drv%ss
        write(*,*)
	write(79,*)'clm1 Start Time: ',drv%yr,drv%mo,drv%da,drv%hr,drv%mn,drv%ss
        write(79,*)
     endif

!=== Restart Writing (2 file are written - active and archive)

     if(rw.eq.2)then

        open(40,file=drv%rstf,form='unformatted') !Active archive restart

        write(40) drv%yr,drv%mo,drv%da,drv%hr,drv%mn,drv%ss,	&
             drv%vclass,drv%nc,drv%nr,drv%nch  !Veg class, no tiles       
        write(40) tile%col                  !Grid Col of Tile   
        write(40) tile%row                  !Grid Row of Tile
        write(40) tile%fgrd                 !Fraction of Grid covered by tile
        write(40) tile%vegt                 !Vegetation Type of Tile
        write(40) clm1%t_grnd                !clm1 Soil Surface Temperature [K]
        write(40) clm1%t_veg                 !CLM Leaf Temperature [K]
        write(40) clm1%h2osno                !CLM Snow Cover, Water Equivalent [mm]
        write(40) clm1%snowage               !CLM Non-dimensional snow age [-]
        write(40) clm1%snowdp                !CLM Snow Depth [m]
        write(40) clm1%h2ocan                !CLM Depth of Water on Foliage [mm]
        write(40) clm1%frac_sno              !CLM Fractional Snow Cover [-]
        write(40) clm1%elai                  !CLM Leaf Area Index
        write(40) clm1%esai                  !CLM Stem Area Index
        write(40) clm1%snl                   !CLM Actual number of snow layers
        write(40) clm1%acc_errh2o            !CLM Accumulation of water balance error
        write(40) clm1%acc_errseb            !CLM Accumulation of energy balance error
        write(40) clm1%istep
        if(ldas%rbias.eq.1)then
         write(40) grid%clmbetak(1)         !CLM forecast bias
	 write(40) grid%clmdelt             !CLM analysis increment
        endif

        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptilent(t) = clm1(t)%dz(l)  
           enddo
           write(40) tmptilent     !CLM Layer Depth [m]
        enddo
        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptilent(t) = clm1(t)%z(l)  
           enddo
           write(40) tmptilent     !CLM Layer Thickness [m]
        enddo
        do l = -nlevsno,nlevsoi
           do t = 1,drv%nch
              tmptilent(t) = clm1(t)%zi(l)  
           enddo
           write(40) tmptilent     !CLM Interface Level Below a "z" Level [m]
        enddo

        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptilent(t) = clm1(t)%t_soisno(l)  
           enddo
           write(40) tmptilent     !CLM Soil + Snow Layer Temperature [K] 
        enddo
        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptilenw(t) = clm1(t)%h2osoi_liq(l)
           enddo
           write(40) tmptilenw     !CLM Average Soil Water Content [kg/m2]
        enddo
        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptileni(t) = clm1(t)%h2osoi_ice(l)
           enddo
           write(40) tmptileni     !CLM Average Ice Content [kg/m2] 
        enddo

        close(40)

        write(*,*)'clm1 Active Restart Written: ',drv%rstf
        write(79,*)'clm1 Active Restart Written: ',drv%rstf

!=== Now write the archived restart file
 91    FORMAT(a4,i3,a5,i4,a1,i4,i2,i2,a7,i3,a1)
 92    format(80a1)
 93    format(a80)
 94    format(i4,i2,i2,i2)
 95    format(10a1)
 96    format(a40)
 98    format(a4,i3,a5,i4,a1,i4,i2,i2)
100    format(a9)
101    format(a7)
       open(90,file='temp',form='formatted',access='direct',recl=80)
       write(90,94,rec=1)drv%yr,drv%mo,drv%da,drv%hr
       read(90,95,rec=1)ftime
       do i=1,10
        if(ftime(i).eq.(' '))ftime(i)='0'
       enddo
       write(90,91,rec=1)'/EXP',ldas%expcode,'/CLM/',drv%yr,  &
       '/',drv%yr,drv%mo,drv%da,'/LDAS.E',ldas%expcode,'.'
       read(90,92,rec=1) (fname(i),i=1,36)
       do i=1,29
        if(fname(i).eq.(' '))fname(i)='0'
       enddo

       write(90,100,rec=1)'mkdir -p '
       read(90,92,rec=1)(fmkdir(i),i=1,9)
       write(90,98,rec=1)'/EXP',ldas%expcode,'/CLM/',   &
        drv%yr,'/',drv%yr,drv%mo,drv%da
       read(90,92,rec=1) (fyrmodir(i),i=1,25)
       do i=1,25
        if(fyrmodir(i).eq.(' '))fyrmodir(i)='0'
       enddo

       write(90,101,rec=1)'.CLMrst'
       read(90,92,rec=1) (fsubs(i),i=1,7)

       write(90,96,rec=1) ldas%odir                       
       read(90,92,rec=1) (fbase(i),i=1,80)
       c=0
       do i=1,80
        if(fbase(i).eq.(' ').and.c.eq.0)c=i-1
       enddo
       write(90,92,rec=1)(fbase(i),i=1,c), (fname(i),i=1,36),  &
                         (ftime(i),i=1,10),(fsubs(i),i=1,7 ) 
       read(90,93,rec=1)filen

 
       write(90,92,rec=1)(fmkdir(i),i=1,9),(fbase(i),i=1,c),   &
                         (fyrmodir(i),i=1,25)
       read(90,93,rec=1)mkfyrmo

!== Archive File Name Generation Complete
!== Make the directories for the clm1 restart file
       if(ldas%gmt.le.ldas%writeintc1)then             
        call system(mkfyrmo)
       endif
       open(40,file=filen,status='unknown',form='unformatted')
        write(40) drv%yr,drv%mo,drv%da,drv%hr,drv%mn,drv%ss,	&
             drv%vclass,drv%nc,drv%nr,drv%nch  !Veg class, no tiles       
        write(40) tile%col                  !Grid Col of Tile   
        write(40) tile%row                  !Grid Row of Tile
        write(40) tile%fgrd                 !Fraction of Grid covered by tile
        write(40) tile%vegt                 !Vegetation Type of Tile
        write(40) clm1%t_grnd                !CLM Soil Surface Temperature [K]
        write(40) clm1%t_veg                 !CLM Leaf Temperature [K]
        write(40) clm1%h2osno                !CLM Snow Cover, Water Equivalent [mm]
        write(40) clm1%snowage               !CLM Non-dimensional snow age [-]
        write(40) clm1%snowdp                !CLM Snow Depth [m]
        write(40) clm1%h2ocan                !CLM Depth of Water on Foliage [mm]
        write(40) clm1%frac_sno              !CLM Fractional Snow Cover [-]
        write(40) clm1%elai                  !CLM Leaf Area Index
        write(40) clm1%esai                  !CLM Stem Area Index
        write(40) clm1%snl                   !CLM Actual number of snow layers
        write(40) clm1%acc_errh2o            !CLM Accumulation of water balance error
        write(40) clm1%acc_errseb            !CLM Accumulation of energy balance error
        write(40) clm1%istep
        if(ldas%rbias.eq.1)then
         write(40) grid%clmbetak(1)         !CLM forecast bias
	 write(40) grid%clmdelt             !CLM analysis increment
        endif

        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptilent(t) = clm1(t)%dz(l)  
           enddo
           write(40) tmptilent     !CLM Layer Depth [m]
        enddo
        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptilent(t) = clm1(t)%z(l)  
           enddo
           write(40) tmptilent     !CLM Layer Thickness [m]
        enddo
        do l = -nlevsno,nlevsoi
           do t = 1,drv%nch
              tmptilent(t) = clm1(t)%zi(l)  
           enddo
           write(40) tmptilent     !CLM Interface Level Below a "z" Level [m]
        enddo

        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptilent(t) = clm1(t)%t_soisno(l)  
           enddo
           write(40) tmptilent     !CLM Soil + Snow Layer Temperature [K] 
        enddo
        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptilenw(t) = clm1(t)%h2osoi_liq(l)
           enddo
           write(40) tmptilenw     !CLM Average Soil Water Content [kg/m2]
        enddo
        do l = -nlevsno+1,nlevsoi
           do t = 1,drv%nch
              tmptileni(t) = clm1(t)%h2osoi_ice(l)
           enddo
           write(40) tmptileni     !CLM Average Ice Content [kg/m2] 
        enddo

        close(40)

        write(*,*)'clm1 Archive Restart Written: ',filen
        write(79,*)'clm1 Archive Restart Written: ',filen

     endif

!== Deallocate arrays if necessary
     if (allocated(col))        deallocate(col)
     if (allocated(row))        deallocate(row)
     if (allocated(vegt))       deallocate(vegt)
     if (allocated(fgrd))       deallocate(fgrd)
     if (allocated(t_grnd))     deallocate(t_grnd)
     if (allocated(t_veg))      deallocate(t_veg)
     if (allocated(h2osno))     deallocate(h2osno)
     if (allocated(snowage))    deallocate(snowage)
     if (allocated(snowdp))     deallocate(snowdp)
     if (allocated(h2ocan))     deallocate(h2ocan)
     if (allocated(frac_sno))   deallocate(frac_sno)
     if (allocated(elai))       deallocate(elai)
     if (allocated(esai))       deallocate(esai)
     if (allocated(snl))        deallocate(snl)
     if (allocated(xerr))       deallocate(xerr)
     if (allocated(zerr))       deallocate(zerr)
     if (allocated(dz))         deallocate(dz)
     if (allocated(z))          deallocate(z)
     if (allocated(zi))         deallocate(zi)
     if (allocated(t_soisno))   deallocate(t_soisno)
     if (allocated(h2osoi_liq)) deallocate(h2osoi_liq)
     if (allocated(h2osoi_ice)) deallocate(h2osoi_ice)
     if (allocated(tmptileot))  deallocate(tmptileot)
     if (allocated(tmptileow))  deallocate(tmptileow)
     if (allocated(tmptileoi))  deallocate(tmptileoi)
     if (allocated(tmptileor))  deallocate(tmptileor)
     if (allocated(tmptileoa))  deallocate(tmptileoa)

     return
   end subroutine drv_restart


