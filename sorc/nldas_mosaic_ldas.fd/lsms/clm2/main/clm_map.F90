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

subroutine clm_map (gindex, nc,nr) 

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Build subgrid patch <-> grid mapping indices and weights
! 
! Method: 
! Build mapping indices and weights: [lsmlon] x [lsmlat] grid <->
! [numland] vector of land points <-> [numpatch] vector of subgrid patches. 
! Allow for variable longitudinal resolution: [numlon] <= [lsmlon]
!
! Example: convert of vector of [numpatch] subgrid patches [t_vec(numpatch)]
! to a 2-d lon x lat array [t_xy(lsmlon,lsmlat)]
!
! Algorithm: for 1 <= i <= numland
!
! t_xy( ixy(i), jxy(i) ) =  t_vec( patch(i,       1) ) * landvec%wtxy(i,       1) + 
!                                              ...                         +
!                           t_vec( patch(i,maxpatch) ) * landvec%wtxy(i,maxpatch) 
!
! If the grid cell has less than [maxpatch] subgrid patches, the value of 
! [patch] used to index the subgrid vector is the value for the first valid
! patch in that grid cell. The weight is zero.
!
! Example: convert a 2-d lon x lat array [t_xy(lsmlon,lsmlat)] to a
! vector of [numpatch] subgrid patches [t_vec(numpatch)]
!
! Algorithm: for 1 <= k <= numpatch
!
! t_vec(k) = t_xy( ixy(land(k)), jxy(land(k)) )
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
#include "misc.h"
  use precision
  use lisdrv_module, only : lis
  use clm_varpar , only : maxpatch
!  use clm_varsur , only : landmask
  use clm_varmap , only : numland, numpatch, begpatch, endpatch, &
                          begland, endland, landvec, patchvec, mapvar_ini
!  use clm_varder , only : clm_varder_init
  use histFileMod, only : histvar_ini
  use mvegFileMod, only : monthveg_ini   
#if (defined ACCUM)
  use accumulMod , only : accumvar_ini
#endif
  use spmdMod    , only : masterproc, iam, npes
#if (defined SPMD)

!                          proc_landi, proc_landf, proc_landpts, &
!                          proc_patchi, proc_patchf, proc_patchpts, spmd_init_patch  
  use mpishorthand
  use tile_spmdMod
#endif
  implicit none

! ------------------------ arguments----------------------------------
  integer :: nc, nr
!  real(r8), intent(in) :: wtxy(nc,nr,maxpatch) !subgrid patch weights
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer  :: i,j,l,m,n                 !indices
!  integer  :: li,lf                       !land indices
  integer  :: nland                       !number of land points
  integer  :: npatch                      !number of patch points 
  real(r8) :: sumwts                      !sum of wtxy over patches
!  real(r8) :: procwt                      !weight of processor patches
!  real(r8), allocatable :: land_wt(:)     !weight of patches for land point
!  integer , allocatable :: land_patchs(:) !numger of patches for land point
!  integer :: ier,ierr,p
  integer :: gindex(nc,nr)
! --------------------------------------------------------------------

#if (defined SPMD)

! ----------------------------------------------------------------------
! Initialize spmd arrays for number of land/patch points per processor
! ----------------------------------------------------------------------

!  call spmd_init_patch

#endif

! --------------------------------------------------------------------
! Find number of land grid cells [numland] and total number of patches
! allowing for multiple subgrid patches in a grid cell.
! --------------------------------------------------------------------
!<test masterproc>
!  if(masterproc) then
!</test masterproc>
  npatch = 0
  nland = 0
#if 0
  do j = 1, nr
     do i = 1, nc
        if (gindex(i,j).ne. -1) then         
           nland = nland+1                               !land point number
           do m = 1, maxpatch                
              if (wtxy(i,j,m) > 0.)  npatch = npatch+1   !subgrid patch number
           end do
           !if (j .eq. 169 .and. i .eq. 1431) then
           !   print*,'DBG:', i,j,nland,npatch,(wtxy(i,j,m),m=1,17)
           !   print*,'DBG: STOPPING'
           !   stop 666
           !endif
        end if
     end do
  end do
#endif
#if ( defined OPENDAP ) && ( defined SPMD )
  numland = di_array(iam)
  numpatch = di_array(iam) 
#else
  numland = lis%d%glbnch
  numpatch = lis%d%glbnch
#endif
!<test masterproc>
!  endif
!</test masterproc>
!broadcast numpatch..
!#if (defined SPMD)
!  call MPI_BCAST(numpatch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!#endif
  if (masterproc) then
     write (6,*)' Surface Grid Characteristics'
     write (6,*)'   longitude points         = ',nc
     write (6,*)'   latitude points          = ',nr
!     write (6,*)'   minimum longitude points = ',minval(numlon)
!     write (6,*)'   maximum longitude points = ',maxval(numlon)
     write (6,*)'   land points on grid                    = ',numland
     write (6,*)'   total points including subgrid patches = ',numpatch
     write (6,*)
  endif
! --------------------------------------------------------------------
! Find beginning and ending patch index for processor
! --------------------------------------------------------------------
#if (defined SPMD)
  begland = 1
  endland = di_array(iam)
  begpatch = 1
  endpatch = di_array(iam)
  print*, 'DBG: clm_map -- begpatch,endpatch',begpatch,endpatch,' (',iam,')'
#else
  begland   = 1
  endland   = numland
  begpatch  = 1
  endpatch  = numpatch
  print*, 'DBG: clm_map -- begpatch,endpatch',begpatch,endpatch,' (',iam,')'
#endif
! --------------------------------------------------------------------
! Allocate dynamic memory
! --------------------------------------------------------------------
!  if(masterproc) then
     call mapvar_ini
!  endif
!  call clm_varder_init
!  if(masterproc) then
!     call histvar_ini
!  endif
  print*, 'past inis..',iam
#if (defined ACCUM)
  call accumvar_ini
#endif
  call monthveg_ini


! --------------------------------------------------------------------
! Build land vector and patch vector mapping components
! --------------------------------------------------------------------

! Determine land vector and patch vector mapping components
!<test masterproc>
! if(masterproc) then
!</test masterproc>
#if 0 
  landvec%wtxy(:,:) = 0._r8
  patchvec%wtxy(:)  = 0._r8
  npatch = 0
  nland  = 0
  do j = 1, nr
     do i = 1, nc
        if (gindex(i,j) .ne. -1) then                 
           nland = nland+1                           
           landvec%ixy(nland) = i                       !land longitude index
           landvec%jxy(nland) = j                       !land latitude index
           do m = 1, maxpatch                           
              if (wtxy(i,j,m) > 0.) then                
                 npatch = npatch+1                      
                 landvec%patch(nland,m) = npatch        !land subgrid patch number
                 landvec%wtxy(nland,m)  = wtxy(i,j,m)   !land subgrid weights
                 patchvec%ixy(npatch)   = i             !patch longitude index
                 patchvec%jxy(npatch)   = j             !patch latitude index
                 patchvec%mxy(npatch)   = m             !patch subgrid index of lnd point
                 patchvec%wtxy(npatch)  = wtxy(i,j,m)   !patch weight
                 patchvec%land(npatch)  = nland         !patch land index
              end if
           end do
        end if
     end do
  end do


! Initialize land vector patches with zero weights to first patch with 
! non-zero weight

  do l = 1, numland
     n = 0
     do m = maxpatch, 1, -1
        if (landvec%wtxy(l,m) > 0.) n = m
     end do
     if (n == 0) then
        write (6,*) 'CLM_MAP error: n = 0'
        call endrun
     end if
     do m = 1, maxpatch
        if (landvec%wtxy(l,m) == 0.) then
           landvec%patch(l,m) = landvec%patch(l,n)
        endif
     end do
  end do

! Error check: make sure weights sum to one for each land cell

  do l = 1, numland
     i = landvec%ixy(l)
     j = landvec%jxy(l)
     sumwts = 0
     do m = 1, maxpatch
        sumwts = sumwts + landvec%wtxy(l,m)
     end do
     if (abs(sumwts-1.) > 1.e-05) then
        write (6,*) 'CLM_MAP error: weights do not sum to 1'
        write (6,*) 'lon = ',i,' lat = ',j,' :sum = ',sumwts
        call endrun
     end if
  end do
!<test masterproc>
!endif
!</test masterproc>
#endif
  print*, 'done clm_map', iam
  return
end subroutine clm_map
