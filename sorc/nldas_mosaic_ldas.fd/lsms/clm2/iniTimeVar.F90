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
!BOP
! !ROUTINE: iniTimeVar.F90
! 
! !DESCRIPTION: 
! 
! Initialize the following time varying variables:
!
!   water      : h2osno, h2ocan, h2osoi$_-$liq, h2osoi$_-$ice, h2osoi$_-$vol\\
!   snow       : snowdp, snowage, snl, dz, z, zi \\
!   temperature: t$_-$soisno, t$_-$veg, t$_-$grnd\\
!
! Note - h2osoi$_-$vol is needed by clm$_-$soilalb -this is not needed on 
! restart since it is computed before the soil albedo computation is 
! called
! 
! Note -  remaining variables are initialized by calls to ecosystem 
! dynamics and albedo subroutines. 
!
! Method: 
! Initial data is saved to instantaneous initial data files for each of
! the [maxpatch] subgrid patches for each of the [numland] land points. 
! If a subgrid patch is not active (e.g., 3 patches rather 
! than [maxpatch]), the inactive subgrid patches have data values for the 
! first subgrid patch. 
! This way, as long as the land mask DOES NOT change among runs 
! (i.e., [numland] is the same), an initial data file can be used in 
! numerous experiments even if the surface types (and hence [numpatch]) 
! differ
! 
!
! Author: Mariana Vertenstein
! 
! !INTERFACE:
subroutine iniTimeVar (readini, eccen, obliqr, lambm0 , mvelpp, lis, tile)
! !USES:
#if ( defined USE_NETCDF )
  use netcdf
#endif
  use lisdrv_module, only : gindex
  use lis_module        
  use tile_module       
  use precision
  use clm_varder
  use clm_varctl, only : clmdrv
  use clm_varmap  , only : numpatch
  use clm_varcon  , only : bdsno, istice, istwet, istsoil, denice, denh2o, tfrz, spval, doalb
  use inicFileMod , only : type_inidat, inicrd, histrd 
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use time_manager, only : get_nstep, get_curr_calday    
#if (defined SPMD)
  use mpishorthand, only : mpicom, mpichar
#endif
!EOP
  implicit none

  type(lisdec)      :: lis
  type(tiledec)      :: tile(lis%d%nch)

!  real :: soiltemp(lis%d%lnc, lis%d%lnr)
!  real :: soiltemp1(lis%d%glbnch)
  real, allocatable, dimension(:,:) :: soiltemp
  real, allocatable, dimension(:)   :: soiltemp1
  integer :: c,r

! ------------------------ arguments ---------------------------------
  logical , intent(in) :: readini  !true if read in initial data set
  real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox long. of perihelion + pi (radians)
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer  :: i,j,k              !loop indices
  integer  :: ier                       !MPI return code
!  real(r8) :: pi                        !3.14159...
!  logical  :: doalb                     !true => albedo time step
  character(len=16) :: initype          !type of initial dataset
  real(r8) :: calday                    !calendar day
  integer :: tmp_nstep
  integer :: ncid, status
  integer :: stid
  integer :: cindex,rindex
! --------------------------------------------------------------------

! ----------------------------------------------------------------------
! Initialize water and temperature based on:
! o readini = true : read initial data set -- requires netCDF codes
! o readini = false: arbitrary initialization
! ----------------------------------------------------------------------
!BOC
  if (readini) then
     if ( masterproc ) write (6,*) 'Reading initial data '
     call type_inidat(initype)
     if (trim(initype) == 'INICFILE') then
        call inicrd () 
     else if (trim(initype) == 'HISTFILE') then
        call histrd ()
     else
        call shr_sys_abort('initial data type is limited to INIC or HIST file only')
     endif
     do k = 1,numpatch
        do j = 1,nlevsoi
           clm(k)%h2osoi_vol(j) = clm(k)%h2osoi_liq(j)/(clm(k)%dz(j)*denh2o) &
                                + clm(k)%h2osoi_ice(j)/(clm(k)%dz(j)*denice)
        end do
     end do
  else
     if ( masterproc ) write (6,*) 'Setting initial data to non-spun up values'

! ========================================================================
! Set snow water 
! ========================================================================

! NOTE: h2ocan, h2osno, snowdp and snowage has valid values everywhere

     do k = 1,numpatch
        if (lis%o%startcode == 4) then
           clm(k)%h2ocan  = 0.
           clm(k)%snowage = 0.
        else
           clm(k)%h2ocan = 0.
           if (clm(k)%itypwat == istice) then
              clm(k)%h2osno = 1000.
           else
              clm(k)%h2osno = clmdrv%clm2_iscv
           endif
           clm(k)%snowdp  = clm(k)%h2osno/bdsno
           clm(k)%snowage = 0.
        endif
     end do
     
! ========================================================================
! Set snow layer number, depth and thickiness 
! ========================================================================
     call snowdp2lev ()

! ========================================================================
! Set snow/soil temperature
! ========================================================================
        
! NOTE: 
! t$_-$soisno only has valid values over non-lake
! t$_-$lake   only has valid values over lake
! t$_-$grnd has valid values over all land
! t$_-$veg  has valid values over all land
#if ( defined USE_NETCDF )
      if ( lis%o%startcode == 3 ) then
         allocate(soiltemp(lis%d%lnc, lis%d%lnr))
         allocate(soiltemp1(lis%d%glbnch))

         soiltemp = -9999.0
   
         call lis_log_msg('MSG: iniTimeVar -- Reading initial soil temp: ' &
                           //trim(lis%p%soiltemp_init))
         status = nf90_open(path=trim(lis%p%soiltemp_init), mode=nf90_nowrite, &
                            ncid=ncid)
         status = nf90_inq_varid(ncid, "SoilTemp_init", stid)
         status = nf90_get_var(ncid, stid, soiltemp1)
         status = nf90_close(ncid)
   
         do c=1,lis%d%lnc
            do r=1,lis%d%lnr
               rindex = 150-r+1
               cindex = c
               if(gindex(cindex,rindex).ne.-1) then 
                  soiltemp(cindex,rindex) = soiltemp1(gindex(cindex,rindex))
               endif
            enddo
         enddo
         deallocate(soiltemp1)
      endif
#endif

     do k =1,numpatch
        clm(k)%t_soisno(-nlevsno+1:nlevsoi) = 0
        if (lis%o%startcode == 4) then
           clm(k)%t_veg   = clm(k)%forc_t
        else         
           clm(k)%t_veg = clmdrv%clm2_it
        endif
        if (.not. clm(k)%lakpoi) then  !not lake
           clm(k)%t_soisno(-nlevsno+1:0) = spval
           if (clm(k)%snl < 0) then    !snow layer temperatures
              do i = clm(k)%snl+1, 0     
                 if (lis%o%startcode == 4) then
                    if (clm(k)%forc_t  < 273.15) then
                       clm(k)%t_soisno(i) = clm(k)%forc_t
                    else                  
                       clm(k)%t_soisno(i) = 273.15 - 1.
                    endif
                 else      
                    if (clmdrv%clm2_it  < 273.15) then
                       clm(k)%t_soisno(i) = clmdrv%clm2_it
                    else
                       clm(k)%t_soisno(i) = 273.15 - 1.
                    endif
                 endif
              enddo
           endif

           do i = 1, nlevsoi
              if (lis%o%startcode == 4) then    
                 if (clm(k)%itypwat == istice) then
                    clm(k)%t_soisno(i) = clm(k)%forc_t
                 else if (clm(k)%itypwat ==istwet) then
                    clm(k)%t_soisno(i) = clm(k)%forc_t
                 else
                    clm(k)%t_soisno(i) = clm(k)%forc_t
                 endif
              else    
                 if (clm(k)%itypwat == istice) then
                    clm(k)%t_soisno(i) = clmdrv%clm2_it
                 else if (clm(k)%itypwat == istwet) then
                    clm(k)%t_soisno(i) = clmdrv%clm2_it
                 else
                    if ( lis%o%startcode == 3 ) then
                       if(soiltemp(tile(k)%col,tile(i)%row).ne.-9999.0) then 
                          clm(k)%t_soisno(i) = soiltemp(tile(k)%col,tile(k)%row)
                       else 
                          clm(k)%t_soisno(i) = clmdrv%clm2_it
                       endif
                    else
                       clm(k)%t_soisno(i) = clmdrv%clm2_it
                    endif
                 endif
              endif
           enddo
           clm(k)%t_grnd = clm(k)%t_soisno(clm(k)%snl+1)
        else                           !lake
           if (lis%o%startcode == 4) then
              clm(k)%t_grnd = clm(k)%forc_t
           else
              clm(k)%t_grnd = clmdrv%clm2_it
           endif
        endif
     end do

     if ( lis%o%startcode == 3 ) then
        deallocate(soiltemp)
     endif

! ========================================================================
! Set snow/soil ice and liquid mass
! ========================================================================
        
! volumetric water is set first and liquid content and ice lens are
! then obtained
! NOTE: h2osoi$_-$vol, h2osoi$_-$liq and h2osoi$_-$ice only have valid values 
! over soil
     do k = 1,numpatch
        clm(k)%h2osoi_vol(         1:nlevsoi) = spval
        clm(k)%h2osoi_liq(-nlevsno+1:nlevsoi) = spval
        clm(k)%h2osoi_ice(-nlevsno+1:nlevsoi) = spval

        if (.not. clm(k)%lakpoi) then  !not lake
           if (clm(k)%snl < 0) then    !snow 
              do i = clm(k)%snl+1, 0
                 clm(k)%h2osoi_ice(i) = clm(k)%dz(i)*250.
                 clm(k)%h2osoi_liq(i) = 0.
              enddo
           endif
           do i = 1, nlevsoi           !soil layers
              if (clm(k)%t_soisno(i) <= tfrz) then       
                 if (lis%o%startcode == 4) then
              else
                 clm(k)%h2osoi_ice(i) = clm(k)%dz(i)* &
                      clmdrv%clm2_ism*clm(k)%watsat(i)*denice
              endif
              clm(k)%h2osoi_liq(i) = 0.
              if (clm(k)%itypwat==istwet .or. clm(k)%itypwat==istice) & 
              clm(k)%h2osoi_ice(i)=clm(k)%dz(i)*denice
           else
              if (lis%o%startcode == 4) then
                 print*, 'Not supposed to be called..'
               else
                  clm(k)%h2osoi_liq(i) = clm(k)%dz(i)* &
                       clmdrv%clm2_ism*clm(k)%watsat(i)*denh2o
               endif
               clm(k)%h2osoi_ice(i) = 0.
               if (clm(k)%itypwat==istwet .or. clm(k)%itypwat==istice) &
               clm(k)%h2osoi_liq(i)=clm(k)%dz(i)*denh2o
            endif
         enddo
         
         do i = 1,nlevsoi
            if (clm(k)%itypwat == istsoil) then
               clm(k)%h2osoi_vol(i) = 0.3_r8
               clm(k)%h2osoi_vol(i) = clm(k)%h2osoi_liq(i)/&
                    (clm(k)%dz(i)*denh2o) &
                    + clm(k)%h2osoi_ice(i)/(clm(k)%dz(i)*denice)
            else
               clm(k)%h2osoi_vol(i) = 1.0_r8
            endif
              clm(k)%h2osoi_vol(i) = min(clm(k)%h2osoi_vol(i),clm(k)%watsat(i))
           end do
        endif
     end do
  end if  ! end of arbitrary initialization if-block

! ========================================================================
! Remaining variables are initialized by calls to ecosystem dynamics and
! albedo subroutines. 
! Note: elai, esai, frac_veg_nosno are computed in Ecosysdyn and needed
! by Fwet and SurfaceAlbedo
! Note: fwet is needed in routine clm_twostream (called by clm_surfalb)
! ========================================================================

  if ( masterproc ) then
     calday = get_curr_calday(lis%t)
  endif
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
  call MPI_BCAST(calday,1,MPI_REAL,0,MPI_COMM_WORLD,ier)
#endif
  doalb = .true.

    print*, 'DBG: iniTimeVar -- calling clm2lairead',' (',iam,')'
    call clm2lairead

#if (defined DGVM)
  call iniTimeConstDGVM()
#endif

  if ( masterproc ) then
     tmp_nstep = get_nstep(lis%t)
  endif
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
  call MPI_BCAST(tmp_nstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
#endif
  do k = 1,numpatch
     clm(k)%nstep = tmp_nstep
     call EcosystemDyn (clm(k), doalb, .false.)
     clm(k)%frac_sno = clm(k)%snowdp/(0.1 + clm(k)%snowdp)  
     clm(k)%frac_veg_nosno = clm(k)%frac_veg_nosno_alb
     if ( lis%d%landcover == 1 ) then     !for UMD
        if( clm(k)%itypveg == 12 ) then
           clm(k)%frac_veg_nosno = 0
        endif
     elseif ( lis%d%landcover == 2 ) then
        if ( clm(k)%itypveg == 11 .or. &
             clm(k)%itypveg == 15 .or. &
             clm(k)%itypveg == 16 ) then  !for IGBP
           clm(k)%frac_veg_nosno = 0
        endif
     endif
     call Fwet(clm(k))
     call SurfaceAlbedo (clm(k), calday, eccen, obliqr, lambm0, mvelpp)
  end do
  return
!EOP
end subroutine iniTimeVar
