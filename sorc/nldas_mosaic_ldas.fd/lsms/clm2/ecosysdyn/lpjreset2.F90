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

subroutine lpjreset2 ()
 
  use precision
  use clmtype
  use clm_varder
  use clm_varctl, only : archive_dir, mss_wpass, mss_irt
  use clm_varpar, only : maxpatch_pft, maxpatch
  use clm_varmap, only : landvec, begland, endland, numland, &
	                 patchvec, begpatch, endpatch, numpatch 
!  use clm_varsur, only : pctlak, pctwet, pctgla, pcturb, landmask
  use clm_varcon, only : cpliq, cpice, istsoil !<- for diagnostic
  use histFileDGVM, only : ncid, spval, dgvm_fn, beg3d, len3d, beg4d, len4d, &
                           fpcgrid_id, itypveg_id, lmind_id, rmind_id, &
                           smind_id, hmind_id, nind_id, begwater_id, endwater_id, &
                           begenergy_id, endenergy_id
  use fileutils, only : set_filename, putfil
#if (defined SPMD)
  use spmdMod, only : masterproc, npes, compute_mpigs_land, compute_mpigs_patch
  use mpishorthand, only : mpir8, mpiint, mpicom
#else
  use spmdMod, only : masterproc
#endif
  use time_manager, only : get_nstep
  implicit none

! ----------------------------------------------------------------------
! purpose           : to reset variables related to lpj
! date first created: November 2000 - lsm version 2 
! by whom           : Sam Levis
! date last revised : 
! by whom           : 
! ----------------------------------------------------------------------
! $Id: lpjreset2.F90,v 1.6 2004/11/24 22:57:06 jim Exp $
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
  integer  :: i,j,k,l,m,ii,ji               !indices
  real(r8) :: sumscl_land(numland)
  real(r8) :: sumwt_land(numland)
  real(r8) :: rbuf_land(maxpatch_pft,numland)
  integer  :: ibuf_land(maxpatch_pft,numland)
  real(r8) :: rbuf1d(numpatch)		
  real(r8) :: rbuf2d(maxpatch,numland)
!  real(r8) :: rbuf3d_xy(lsmlon,lsmlat,maxpatch_pft)
!  integer  :: ibuf3d_xy(lsmlon,lsmlat,maxpatch_pft)
  real(r8) :: begwater(numpatch)          !diagnostic output
  real(r8) :: endwater(numpatch)          !diagnostic output
  real(r8) :: begenergy(numpatch)         !diagnostic output
  real(r8) :: endenergy(numpatch)         !diagnostic output 
!  real(r8) :: begwater_xy(lsmlon,lsmlat)  !diagnostic output
!  real(r8) :: endwater_xy(lsmlon,lsmlat)  !diagnostic output 
!  real(r8) :: begenergy_xy(lsmlon,lsmlat) !diagnostic output
!  real(r8) :: endenergy_xy(lsmlon,lsmlat) !diagnostic output
  real(r8) :: cv(-nlevsno+1:nlevsoi)      !temporary for diagnostics
  character(len=256) :: rem_dir           !remote (archive) directory
  character(len=256) :: rem_fn            !remote (archive) filename
#if (defined SPMD)
  integer  :: numrecvv(0:npes-1)    !vector of items to be received  
  integer  :: displsv(0:npes-1)     !displacement vector
  integer  :: numsend               !number of items to be sent
  integer  :: ier                   !MPI error status
#endif
! ----------------------------------------------------------------------

!------------------------------------------------------------------
! diagnostic for energy and water balance before lpj
!------------------------------------------------------------------
  print*,' lpjreset not supposed to be called..'
#if 0 
  do k = begpatch, endpatch
     if (clm(k)%itypwat == istsoil) then
        begenergy(k) = 0.0
        begwater(k)=clm(k)%h2ocan+clm(k)%h2osno
        do j = 1, nlevsoi
           begwater(k) = begwater(k) + clm(k)%h2osoi_ice(j) + clm(k)%h2osoi_liq(j)
           cv(j) = clm(k)%csol(j)*(1-clm(k)%watsat(j))*clm(k)%dz(j) +   &
                (clm(k)%h2osoi_ice(j)*cpice + clm(k)%h2osoi_liq(j)*cpliq)
           begenergy(k) = begenergy(k) + cv(j)*clm(k)%t_soisno(j)
        enddo
        if (clm(k)%snl+1 == 1 .AND. clm(k)%h2osno > 0.) then
           cv(1) = cv(1) + cpice*clm(k)%h2osno
           begenergy(k) = begenergy(k) + cv(1)*clm(k)%t_soisno(1)
        end if
        if (clm(k)%snl+1 < 1)then
           do j = clm(k)%snl+1, 0
              cv(j) = cpliq*clm(k)%h2osoi_liq(j) + cpice*clm(k)%h2osoi_ice(j)
              begenergy(k) = begenergy(k) + cv(j)*clm(k)%t_soisno(j)
           enddo
        end if
     else
        begwater(k) = spval
        begenergy(k) = spval
     end if
  enddo

!------------------------------------------------------------------
! Reset patchvec%wtxy and landvec%wtxy and do error check
!------------------------------------------------------------------

! In CLM2 with satellite data, the number of veg patches is determined once
! and is less than maxpatch_pft (4) in some cells.
! In LSM with LPJ, the number of veg patches could be dynamic. Until we
! implement it as such, we will make all grid cells have 10 veg patches.


  sumscl_land(:) = 0.
  sumwt_land(:)  = 0.
  do l = begland,endland
     i = landvec%ixy(l)                
     j = landvec%jxy(l)                
     sumscl_land(l) = pcturb(i,j)+pctlak(i,j)+pctwet(i,j)+pctgla(i,j)
     if (sumscl_land(l) < 100.) then
        do m = 1, maxpatch_pft
           k = landvec%patch(l,m)
           patchvec%wtxy(k) = clm(k)%fpcgrid * (100.-sumscl_land(l))/100.
           landvec%wtxy(l,m) = patchvec%wtxy(k)
           sumwt_land(l) = sumwt_land(l) + patchvec%wtxy(k)
        end do
     end if
     sumwt_land(l) = sumwt_land(l) + sumscl_land(l)/100.
     if (abs(sumwt_land(l) - 1.0) > 1.0e-6) then
        write(6,*) 'Error in lpjreset2: sumwt .ne. 1 at', i, j
        write(6,*) 'sumwt =', sumwt_land(l)
        call endrun
     end if
  end do

#if (defined SPMD)

! determine landvec%wtxy for all land points on master processor
! this is needed for global sums 

    do l = begland,endland
       do m = 1,maxpatch
          rbuf2d(m,l) = landvec%wtxy(l,m)
       end do
    end do
    call compute_mpigs_land(maxpatch, numsend, numrecvv, displsv)
    call mpi_gatherv (rbuf2d(1,begland), numsend , mpir8, &
                      rbuf2d           , numrecvv, displsv, mpir8, 0, mpicom, ier)
    if (masterproc) then
       do l = 1,numland
          do m = 1,maxpatch
             landvec%wtxy(l,m) = rbuf2d(m,l)
          end do
       end do
    endif

! determine patchvec%wtxy for all patch points on master processor
! this is needed for calls to v2xy

    do k = begpatch,endpatch
       rbuf1d(k) = patchvec%wtxy(k)
    end do
    call compute_mpigs_patch(1, numsend, numrecvv, displsv)
    call mpi_gatherv (rbuf1d(begpatch), numsend , mpir8, &
                      rbuf1d          , numrecvv, displsv, mpir8, 0, mpicom, ier)
    if (masterproc) then
       do k = 1,numpatch
          patchvec%wtxy(k) = rbuf1d(k)
       end do
    endif

#endif

!------------------------------------------------------------------
! Write out more variables to DGVM history output
!------------------------------------------------------------------

! Write out pctxy and pftxy

  rbuf_land(:,:) = 0.
  ibuf_land(:,:) = 0
  do l = begland,endland
     if (sumscl_land(l) < 100.) then
        do m = 1, maxpatch_pft
           k = landvec%patch(l,m)
           rbuf_land(m,l) = clm(k)%fpcgrid * 100.0
           ibuf_land(m,l) = clm(k)%itypveg
        end do
     end if
  end do
#if (defined SPMD)
  call compute_mpigs_land(maxpatch_pft, numsend, numrecvv, displsv)
  call mpi_gatherv (rbuf_land(1,begland), numsend , mpir8, &
                    rbuf_land           , numrecvv, displsv, mpir8, 0, mpicom, ier)
  call compute_mpigs_land(maxpatch_pft, numsend, numrecvv, displsv)
  call mpi_gatherv (ibuf_land(1,begland), numsend , mpiint, &
                    ibuf_land           , numrecvv, displsv, mpiint, 0, mpicom, ier)
#endif           
  if (masterproc) then
     rbuf3d_xy(:,:,:)= spval
     ibuf3d_xy(:,:,:)= 9999
     do l = 1,numland
        i = landvec%ixy(l)                
        j = landvec%jxy(l)                
        do m=1,maxpatch_pft
           rbuf3d_xy(i,j,m) = rbuf_land(m,l)
           ibuf3d_xy(i,j,m) = ibuf_land(m,l)
        end do
     end do
     call wrap_put_vara_int (ncid, itypveg_id, beg4d, len4d, ibuf3d_xy)
     call wrap_put_vara_realx (ncid, fpcgrid_id, beg4d, len4d, rbuf3d_xy)
endif
  
! Write out lm_ind

  rbuf_land(:,:) = 0.
  do l = begland,endland
     if (sumscl_land(l) < 100.) then
        do m = 1, maxpatch_pft
           k = landvec%patch(l,m)
           rbuf_land(m,l) = clm(k)%lm_ind
        end do
     end if
  end do
#if (defined SPMD)
  call compute_mpigs_land(maxpatch_pft, numsend, numrecvv, displsv)
  call mpi_gatherv (rbuf_land(1,begland), numsend , mpir8, &
                    rbuf_land           , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     rbuf3d_xy(:,:,:) = spval
     do l = 1,numland
        i = landvec%ixy(l)                
        j = landvec%jxy(l)                
        do m=1,maxpatch_pft
           rbuf3d_xy(i,j,m) = rbuf_land(m,l)
        end do
     end do
     call wrap_put_vara_realx (ncid, lmind_id, beg4d, len4d, rbuf3d_xy)
  endif

! Write out rm_ind

  rbuf_land(:,:) = 0.
  do l = begland,endland
     if (sumscl_land(l) < 100.) then
        do m = 1, maxpatch_pft
           k = landvec%patch(l,m)
           rbuf_land(m,l) = clm(k)%rm_ind
        end do
     end if
  end do
#if (defined SPMD)
  call compute_mpigs_land(maxpatch_pft, numsend, numrecvv, displsv)
  call mpi_gatherv (rbuf_land(1,begland), numsend , mpir8, &
                    rbuf_land           , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     rbuf3d_xy(:,:,:) = spval
     do l = 1,numland
        i = landvec%ixy(l)                
        j = landvec%jxy(l)                
        do m=1,maxpatch_pft
           rbuf3d_xy(i,j,m) = rbuf_land(m,l)
        end do
     end do
     call wrap_put_vara_realx (ncid, rmind_id, beg4d, len4d, rbuf3d_xy)
  endif

! Write out sm_ind

  rbuf_land(:,:) = 0.
  do l = begland,endland
     if (sumscl_land(l) < 100.) then
        do m = 1, maxpatch_pft
           k = landvec%patch(l,m)
           rbuf_land(m,l) = clm(k)%sm_ind
        end do
     end if
  end do
#if (defined SPMD)
  call compute_mpigs_land(maxpatch_pft, numsend, numrecvv, displsv)
  call mpi_gatherv (rbuf_land(1,begland), numsend , mpir8, &
                    rbuf_land           , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     rbuf3d_xy(:,:,:) = spval
     do l = 1,numland
        i = landvec%ixy(l)                
        j = landvec%jxy(l)                
        do m=1,maxpatch_pft
           rbuf3d_xy(i,j,m) = rbuf_land(m,l)
        end do
     end do
     call wrap_put_vara_realx (ncid, smind_id, beg4d, len4d, rbuf3d_xy)
  endif

! Write out hm_ind

  rbuf_land(:,:) = 0.
  do l = begland,endland
     if (sumscl_land(l) < 100.) then
        do m = 1, maxpatch_pft
           k = landvec%patch(l,m)
           rbuf_land(m,l) = clm(k)%hm_ind
        end do
     end if
  end do
#if (defined SPMD)
  call compute_mpigs_land(maxpatch_pft, numsend, numrecvv, displsv)
  call mpi_gatherv (rbuf_land(1,begland), numsend , mpir8, &
                    rbuf_land           , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     rbuf3d_xy(:,:,:) = spval
     do l = 1,numland
        i = landvec%ixy(l)                
        j = landvec%jxy(l)                
        do m=1,maxpatch_pft
           rbuf3d_xy(i,j,m) = rbuf_land(m,l)
        end do
     end do
     call wrap_put_vara_realx (ncid, hmind_id, beg4d, len4d, rbuf3d_xy)
  endif

! Write out nind

  rbuf_land(:,:) = 0.
  do l = begland,endland
     if (sumscl_land(l) < 100.) then
        do m = 1, maxpatch_pft
           k = landvec%patch(l,m)
           rbuf_land(m,l) = clm(k)%nind
        end do
     end if
  end do
#if (defined SPMD)
  call compute_mpigs_land(maxpatch_pft, numsend, numrecvv, displsv)
  call mpi_gatherv (rbuf_land(1,begland), numsend , mpir8, &
                    rbuf_land           , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     rbuf3d_xy(:,:,:) = spval
     do l = 1,numland
        i = landvec%ixy(l)                
        j = landvec%jxy(l)                
        do m=1,maxpatch_pft
           rbuf3d_xy(i,j,m) = rbuf_land(m,l)
        end do
     end do
     call wrap_put_vara_realx (ncid, nind_id, beg4d, len4d, rbuf3d_xy)
  endif

! Write out diagnostic for energy and water balance after lpj

  do k = begpatch, endpatch
     if (clm(k)%itypwat == istsoil) then
        endenergy(k) = 0.0
        endwater(k)=clm(k)%h2ocan+clm(k)%h2osno
        do j = 1, nlevsoi
           endwater(k) = endwater(k) + clm(k)%h2osoi_ice(j) + clm(k)%h2osoi_liq(j)
           cv(j) = clm(k)%csol(j)*(1-clm(k)%watsat(j))*clm(k)%dz(j) +   &
                (clm(k)%h2osoi_ice(j)*cpice + clm(k)%h2osoi_liq(j)*cpliq)
           endenergy(k) = endenergy(k) + cv(j)*clm(k)%t_soisno(j)
        enddo
        if (clm(k)%snl+1 == 1 .AND. clm(k)%h2osno > 0.) then
           cv(1) = cv(1) + cpice*clm(k)%h2osno
           endenergy(k) = endenergy(k) + cv(1)*clm(k)%t_soisno(1)
        end if
        if (clm(k)%snl+1 < 1)then
           do j = clm(k)%snl+1, 0
              cv(j) = cpliq*clm(k)%h2osoi_liq(j) + cpice*clm(k)%h2osoi_ice(j)
              endenergy(k) = endenergy(k) + cv(j)*clm(k)%t_soisno(j)
           enddo
        end if
     else
        endwater(k) = spval
        endenergy(k) = spval
     end if
  enddo
#if (defined SPMD)
  call compute_mpigs_patch(1, numsend, numrecvv, displsv)
  call mpi_gatherv (begwater(begpatch), numsend , mpir8, &
                    begwater          , numrecvv, displsv, mpir8, 0, mpicom, ier)
  call mpi_gatherv (endwater(begpatch), numsend , mpir8, &
                    endwater          , numrecvv, displsv, mpir8, 0, mpicom, ier)
  call mpi_gatherv (begenergy(begpatch), numsend , mpir8, &
                    begenergy          , numrecvv, displsv, mpir8, 0, mpicom, ier)
  call mpi_gatherv (endenergy(begpatch), numsend , mpir8, &
                    endenergy          , numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
  if (masterproc) then
     begwater_xy(:,:) = 0.0
     endwater_xy(:,:) = 0.0
     begenergy_xy(:,:) = 0.0
     endenergy_xy(:,:) = 0.0

     call v2xy (begwater, 0._r4, begwater_xy)
     call v2xy (begenergy, 0._r4, begenergy_xy)
     call v2xy (endwater, 0._r4, endwater_xy)
     call v2xy (endenergy, 0._r4, endenergy_xy)

     where (landmask(:,:) /= 1)
        begwater_xy(:,:) = spval
        endwater_xy(:,:) = spval
        begenergy_xy(:,:) = spval
        endenergy_xy(:,:) = spval
     endwhere

     call wrap_put_vara_realx (ncid, begwater_id, beg3d, len3d, begwater_xy)
     call wrap_put_vara_realx (ncid, endwater_id, beg3d, len3d, endwater_xy)
     call wrap_put_vara_realx (ncid, begenergy_id, beg3d, len3d, begenergy_xy)
     call wrap_put_vara_realx (ncid, endenergy_id, beg3d, len3d, endenergy_xy)
  endif
     
! close and archive netcdf DGVM history file

  if (masterproc) then
     call wrap_close(ncid)
     write(6,*)'(LPJRESET2): Finished writing clm2 DGVM history dataset ',&
          trim(dgvm_fn), 'at nstep = ',get_nstep()
     if (mss_irt > 0) then
        rem_dir = trim(archive_dir) // '/hist/'
        rem_fn = set_filename(rem_dir, dgvm_fn)
        call putfil (dgvm_fn, rem_fn, mss_wpass, mss_irt, .true.)
     endif
  endif
#endif
  return
end subroutine lpjreset2
