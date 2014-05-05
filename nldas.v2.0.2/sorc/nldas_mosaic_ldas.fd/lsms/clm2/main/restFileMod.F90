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

module restFileMod

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read/Write CLM restart files
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: restFileMod.F90,v 1.16 2005/02/01 18:53:53 jim Exp $
!-----------------------------------------------------------------------

  use lisdrv_module, only : lis
  use precision
  use clm_varder
  use clm_varpar, only : nlevsoi, numrad, maxhist, &
                         maxpatch_pft, maxpatch
  use clm_varmap, only : patchvec, &
                         landvec
  use clm_varctl, only : nsrest, rpntdir, rpntfil, nrevsn, locfnh, & 
                         archive_dir, mss_irt, mss_wpass, csm_doflxave, caseid
  use fileutils , only : opnfil, putfil, getfil, getavu, relavu, set_filename 
  use histFileMod , only : slfld, mlsoifld, ntim, nbeghis, nhist, ncgetid, ncid
  use histHandlerMod, only : mcdate_i, mcsec_i, mdcur_i, mscur_i       
#if (defined ACCUM)
  use accumulMod, only : naccflds, accval
#endif
#if (defined DGVM)
  use shr_const_mod, only: SHR_CONST_CDAY
#endif
#if (defined SPMD)
  use spmdMod, only : masterproc, npes, compute_mpigs_patch, compute_mpigs_land
  use mpishorthand, only : mpir8, mpilog, mpiint, mpicom 
#else
  use spmdMod, only : masterproc
#endif
#if (defined RTM)
  use RtmMod      , only : volr, ncount_rtm, totrunin_ave, prec_ave, evap_ave, &
                           qchan2, qchocn2, ocnrof_vec, prec_global, evap_global, &
                           runlnd_global, runrtm_global, volrtm_global, ocnrtm_global,&
                           ncount_global, yrold
#endif
#if (defined COUP_CSM)
  use clm_csmMod
#endif
  implicit none

! Methods

  public  :: restrd
  public  :: restwrt
  private :: write_rest_pfile
!  private :: set_restart_filename ! now using create_restart_filename

! Generic procedures

  PRIVATE :: readin, wrtout 

  INTERFACE readin
     MODULE procedure readinsc_log
     MODULE procedure readinsc_int
     MODULE procedure readinsc_real
     MODULE procedure readin1d_int
     MODULE procedure readin2d_int
     MODULE procedure readin3d_int
     MODULE procedure readin1d_real
     MODULE procedure readin2d_real
     MODULE procedure readin3d_real
  END INTERFACE

  INTERFACE wrtout
     MODULE procedure wrtoutsc_log
     MODULE procedure wrtoutsc_int
     MODULE procedure wrtoutsc_real
     MODULE procedure wrtout1d_int
     MODULE procedure wrtout2d_int
     MODULE procedure wrtout3d_int
     MODULE procedure wrtout1d_real
     MODULE procedure wrtout2d_real
     MODULE procedure wrtout3d_real

  END INTERFACE

  integer ,allocatable :: ibuf1d(:)     !temporary integer buffer
  real(r8),allocatable ::  buf1d(:)     !temporary buffer
  integer ,allocatable :: ibuf2d(:,:)   !temporary integer buffer
  real(r8),allocatable ::  buf2d(:,:)   !temporary buffer
  integer ,allocatable :: ibuf3d(:,:,:) !temporary integer buffer
  real(r8),allocatable ::  buf3d(:,:,:) !temporary buffer

  integer, private, parameter :: rest_id = 4

!=======================================================================
CONTAINS
!=======================================================================

  subroutine restrd()
  
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read CLM restart file. Open history file if needed
! 
! Method: 
! This code reads the clm restart file. If the current history file(s) are
! not full, file(s) are opened so that subsequent time samples are added
! until the file is full. A new history file is used on a branch run. 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use clm_varcon , only : denice, denh2o
#if (defined COUP_CAM)
    use time_manager, only : get_step_size, get_nstep 
#else
    use time_manager, only : get_step_size, get_nstep, &
                             timemgr_read_restart, timemgr_restart
#endif  
#if (defined COUP_CSM)
    use controlMod  , only : csm_dtime       !dtime from input namelist
#endif
!    include 'netcdf.inc'

! ------------------------ local variables ------------------------
    integer :: i,j,k,l,m,n                   !indices
    integer :: nio                           !Fortran unit number
    integer :: nstep                         !time index
    integer :: num                           !number of fields (temporary)
    character(len=256) :: fnamer             !full name of restart file
    character(len=256) :: fnameh(maxhist)    !full name of history file
    character(len=256) :: locfn              !local file name 
    character(len= 16) :: casename           !case name read in from restart
    logical :: flxave_res                    !flux averaging flag read from restart file
    integer ierr                              !temporaries 
#if (defined DGVM)
    real(r8) :: sumscl(lsmlon,lsmlat)
    real(r8) :: sumwt(lsmlon,lsmlat)
    real(r8) :: r1
    real(r8) :: r2
#endif
#if (defined SPMD)
    integer :: numsendv(0:npes-1)            !vector of items to be sent
    integer :: numsend                       !number of items to be sent
    integer :: numrecv                       !number of items to be received
    integer :: numrecvv(0:npes-1)            !vector of items to be received  
    integer :: displsv(0:npes-1)             !displacement vector
#endif
    real(r8):: dtime                         !step size (seconds)      
#if (defined COUP_CAM)
    integer :: clm_nstep                     !nstep from restart file
#endif
    integer :: testd
! -----------------------------------------------------------------
    if (masterproc) then
       write(6,*) 'Attempting to read restart data .....'
       write(6,'(72a1)') ("-",i=1,60)
    endif

! -----------------------------------------------------------------
! Get appropriate restart and history file names if continuation run. 
!    o restart: the restart pointer file contains the full mass store or 
!      local names of the current restart file clmr_xxxx and 
!      history files clmh_xxxx, clmha_xxxx, etc.
!    o branch: the variable [nrevsn] is the full mass store or local name
!      of the appropriate restart file clmr_xxxx. do not need to have 
!      history file names because will open new history files starting 
!      with clmh_0001, clmha_0001, etc
! -----------------------------------------------------------------

    if (masterproc) then

! Determine restart filename

!       if (nsrest == 1) then     !restart
!          nio = getavu()
!          !locfn = trim(rpntdir) //'/'// trim(rpntfil)
!          locfn = trim(rpntdir) // trim(rpntfil)
!          call opnfil (locfn, nio, 'f')
!          read (nio,'(a80)') fnamer
!          call relavu (nio)
!       else                      !branch run
!          fnamer = nrevsn
!       end if

! Open and get restart file

       locfn = trim(rpntdir) // trim(rpntfil)
       nio = getavu()
!       call getfil (fnamer, locfn, 0)
       call opnfil (locfn, nio, 'u')

    endif

! -----------------------------------------------------------------
! Read restart file
! -----------------------------------------------------------------

! Check that restart file id matches that for this code version
    
    call readin (nio, n) 
    if (masterproc) then
       if (n == rest_id) then
          write(6,*)'RESTRD: using restart file id ',rest_id
       else
          write(6,*)'RESTRD: input restart file id ',n, &
               ' does not match required restart id ',rest_id
          call endrun
       endif
    endif

#if (defined COUP_CAM)

! read in time step - make check to see that clm restart time step 
! is consistent with cam restart time step
! note that in cam mode - the time manager restart variables have
! already been read in by calls to routines timemgr_read_restart and
! timemgr_restart from the cam restart files

    call readin (nio, clm_nstep)

    if ((clm_nstep+1) /= get_nstep()) then
       write(6,*)'(RESTRD): incompatibility in clm and cam restart dates'
       write(6,*)'  restart step from cam = ',get_nstep()
       write(6,*)'  restart step from clm = ',clm_nstep + 1
       call endrun
    endif

#else

! restart the time manager

    if (masterproc) then
       call timemgr_read_restart(nio)
       !call timemgr_restart()
    endif

#endif

#if (defined COUP_CSM)
    if (csm_dtime /= get_step_size(lis%t)) then
       write(6,*)'(RESTRD): error '
       write(6,*)'namelist dtime on restart does not match input dtime'
       call endrun
    endif
#endif

! set dtime derived type component (for an initial run, this is set
! in routine iniTimeConst.F90) 

    clm(1:lis%d%glbnch)%dtime = get_step_size(lis%t)

! read in clm necessary 1d fields

    allocate ( buf1d(1:lis%d%glbnch))
    allocate (ibuf1d(1:lis%d%glbnch))

    call readin (nio, ibuf1d)
    clm(1:lis%d%glbnch)%snl = ibuf1d(1:lis%d%glbnch)

    call readin (nio, ibuf1d)
    clm(1:lis%d%glbnch)%frac_veg_nosno_alb = ibuf1d(1:lis%d%glbnch)

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%h2osno = buf1d(1:lis%d%glbnch)  

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%h2ocan = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%snowdp = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%snowage = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%frac_sno = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%t_veg = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%t_grnd = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%fwet = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%tlai = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%tsai = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%elai = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%esai = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d) 
    clm(1:lis%d%glbnch)%fsun = buf1d(1:lis%d%glbnch) 

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%htop = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%hbot = buf1d(1:lis%d%glbnch)

#if (defined DGVM)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%wf = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%t_mo_min = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%annpsn = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%annpsnpot = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%fmicr = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%bm_inc = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%afmicr = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%t10min = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%tmomin20 = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%agdd20 = buf1d(1:lis%d%glbnch)

    call readin (nio, ibuf1d)
    clm(1:lis%d%glbnch)%itypveg = ibuf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%fpcgrid = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%lai_ind = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%crownarea = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%dphen = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%leafon = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%leafof = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%firelength = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%litterag = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%litterbg = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%cpool_fast = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%cpool_slow = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%k_fast_ave = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%k_slow_ave = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%litter_decom_ave = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    do k = 1, lis%d%glbnch
       if (buf1d(k) == 1.0) then
          clm(k)%present = .true.
       else
          clm(k)%present = .false.
       end if
    end do

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%nind = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%lm_ind = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%sm_ind = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%hm_ind = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    clm(1:lis%d%glbnch)%rm_ind = buf1d(1:lis%d%glbnch)

    call readin (nio, buf1d)
    do k = 1,lis%d%glbnch
       l = patchvec%land(k)
       m = patchvec%mxy(k)
       landvec%wtxy(l,m) = buf1d(k)
       patchvec%wtxy(k) = buf1d(k)
    end do

#endif

    deallocate (ibuf1d)
    deallocate ( buf1d)

#if (defined DGVM) && (defined SPMD)

! determine landvec%wtxy for all land points on master processor
! this is needed for global sums 

    allocate (buf2d(maxpatch,lis%d%glbnch))
    do l = 1,lis%d%glbnch
       do m = 1,maxpatch
          buf2d(m,l) = landvec%wtxy(l,m)
       end do
    end do
    call compute_mpigs_land(maxpatch, numsend, numrecvv, displsv)
    call mpi_gatherv (buf2d(1,1), numsend , mpir8, &
                      buf2d           , numrecvv, displsv, mpir8, 0, mpicom, ierr)
    if (masterproc) then
       do l = 1,lis%d%glbnch
          do m = 1,maxpatch
             landvec%wtxy(l,m) = buf2d(m,l)
          end do
       end do
    endif
    deallocate (buf2d)

! determine patchvec%wtxy for all patch points on master processor
! this is needed for calls to v2xy

    allocate (buf1d(lis%d%glbnch))
    do k = 1,lis%d%glbnch
       buf1d(k) = patchvec%wtxy(k)
    end do
    call compute_mpigs_patch(1, numsend, numrecvv, displsv)
    call mpi_gatherv (buf1d(1), numsend , mpir8, &
                      buf1d          , numrecvv, displsv, mpir8, 0, mpicom, ierr)
    if (masterproc) then
       do k = 1,lis%d%glbnch
          patchvec%wtxy(k) = buf1d(k)
       end do
    endif
    deallocate (buf1d)

#endif

! read in multi level snow only fields

    allocate (buf2d(-nlevsno+1:0,1:lis%d%glbnch))
    call readin (nio, buf2d, nlevsno)
    do j = -nlevsno+1,0
       clm(1:lis%d%glbnch)%dz(j) = buf2d(j,1:lis%d%glbnch)
    end do

    call readin (nio, buf2d, nlevsno)
    do j = -nlevsno+1,0
       clm(1:lis%d%glbnch)%z(j) = buf2d(j,1:lis%d%glbnch)
    end do
    deallocate(buf2d)

    allocate (buf2d(-nlevsno:0,1:lis%d%glbnch))
    call readin (nio, buf2d, nlevsno+1)
    do j = -nlevsno,0
       clm(1:lis%d%glbnch)%zi(j) = buf2d(j,1:lis%d%glbnch)
    end do
    deallocate(buf2d)

! read in multi level snow-soil fields

    allocate (buf2d(-nlevsno+1:nlevsoi,1:lis%d%glbnch))
    call readin (nio, buf2d, nlevsno+nlevsoi)
    do j = -nlevsno+1,nlevsoi
       clm(1:lis%d%glbnch)%t_soisno(j) = buf2d(j,1:lis%d%glbnch)
    end do

    call readin (nio, buf2d, nlevsoi+nlevsno)
    do j = -nlevsno+1,nlevsoi
       clm(1:lis%d%glbnch)%h2osoi_liq(j) = buf2d(j,1:lis%d%glbnch)
    end do

    call readin (nio, buf2d, nlevsoi+nlevsno)
    do j = -nlevsno+1,nlevsoi
       clm(1:lis%d%glbnch)%h2osoi_ice(j) = buf2d(j,1:lis%d%glbnch)
    end do
    deallocate(buf2d)

!    allocate (buf2d(1:nlevlak,1:lis%d%glbnch))
!    call readin (nio, buf2d, nlevlak)
!    do j = 1,nlevlak
!       clm(1:lis%d%glbnch)%t_lake(j) = buf2d(j,1:lis%d%glbnch)
!    end do
!    deallocate(buf2d)

! Determine volumetric soil water

     do k = 1,lis%d%glbnch
        do j = 1,nlevsoi
           clm(k)%h2osoi_vol(j) = clm(k)%h2osoi_liq(j)/(clm(k)%dz(j)*denh2o) &
                                + clm(k)%h2osoi_ice(j)/(clm(k)%dz(j)*denice)
        end do
     end do

! read in multi level surface albedo and radiation data

    allocate (buf2d(numrad,1:lis%d%glbnch))
!    call readin (nio, buf2d  , numrad)
!    do j = 1,numrad
!       clm(1:lis%d%glbnch)%albd(j) = buf2d(j,1:lis%d%glbnch)
!    end do

!    call readin (nio, buf2d  , numrad)
!    do j = 1,numrad
!       clm(1:lis%d%glbnch)%albi(j) = buf2d(j,1:lis%d%glbnch)
!    end do

    call readin (nio, buf2d, numrad)
    do j = 1,numrad
       clm(1:lis%d%glbnch)%albgrd(j) = buf2d(j,1:lis%d%glbnch)
    end do
    call readin (nio, buf2d, numrad)

    do j = 1,numrad
       clm(1:lis%d%glbnch)%albgri(j) = buf2d(j,1:lis%d%glbnch)
    end do
    call readin (nio, buf2d  , numrad)

    do j = 1,numrad
       clm(1:lis%d%glbnch)%fabd(j) = buf2d(j,1:lis%d%glbnch)
    end do
    call readin (nio, buf2d  , numrad)
    do j = 1,numrad
       clm(1:lis%d%glbnch)%fabi(j) = buf2d(j,1:lis%d%glbnch)
    end do
    call readin (nio, buf2d  , numrad)

    do j = 1,numrad
       clm(1:lis%d%glbnch)%ftdd(j) = buf2d(j,1:lis%d%glbnch)
    end do
    call readin (nio, buf2d  , numrad)

    do j = 1,numrad
       clm(1:lis%d%glbnch)%ftid(j) = buf2d(j,1:lis%d%glbnch)
    end do
    call readin (nio, buf2d  , numrad)

    do j = 1,numrad
       clm(1:lis%d%glbnch)%ftii(j) = buf2d(j,1:lis%d%glbnch)
    end do
    deallocate(buf2d)

#if (defined ACCUM)
! read data for accumulated fields

    call readin (nio, naccflds)
    allocate(buf2d(naccflds,1:lis%d%glbnch))
    call readin (nio, buf2d, naccflds)
    do j = 1, naccflds
       do i = 1, lis%d%glbnch
         accval(i,j) = buf2d(j,i)
       end do
    end do
    deallocate(buf2d)
#endif

! read data for river routing model, if appropriate

#if (defined RTM)

    if (masterproc) then
       read (nio) volr
       read (nio) ncount_rtm
       read (nio) ncount_global, yrold
       read (nio) prec_global,evap_global,runlnd_global,runrtm_global,volrtm_global,ocnrtm_global
       read (nio) (totrunin_ave(j),j=1,lis%d%glbnch)
       read (nio) (prec_ave(j)    ,j=1,lis%d%glbnch)
       read (nio) (evap_ave(j)    ,j=1,lis%d%glbnch)
       read (nio) (qchan2(j)      ,j=1,lis%d%glbnch)
       read (nio) (qchocn2(j)     ,j=1,lis%d%glbnch)
#if (defined COUP_CSM)
       read (nio) ocnrof_vec 
#endif
    endif
!#if (defined SPMD)
!    call compute_mpigs_patch(1, numrecv, numsendv, displsv)
!    call mpi_scatterv (totrunin_ave          , numsendv, displsv, mpir8, &
!                       totrunin_ave(1), numrecv , mpir8, 0, mpicom, ierr)
!    call mpi_scatterv (prec_ave              , numsendv, displsv, mpir8, &
!                       prec_ave(1)    , numrecv , mpir8, 0, mpicom, ierr)
!    call mpi_scatterv (evap_ave              , numsendv, displsv, mpir8, &
!                       evap_ave(1)    , numrecv , mpir8, 0, mpicom, ierr)
!    call mpi_scatterv (qchocn2               , numsendv, displsv, mpir8, &
!                       qchocn2(1)     , numrecv , mpir8, 0, mpicom, ierr)
!    call mpi_scatterv (qchan2                , numsendv, displsv, mpir8, &
!                       qchan2(1)      , numrecv , mpir8, 0, mpicom, ierr)
!#endif

#endif

! read data for flux coupled case, if appropriate

#if (defined COUP_CSM)
    call readin (nio, flxave_res)
    call readin (nio, dosend)
    if (masterproc) then
       if ((flxave_res .and. .not.csm_doflxave).or.(.not.flxave_res .and. csm_doflxave)) then
          write(6,*)'(RESTRD): flxave value from namelist ',csm_doflxave, &
               ' must be the same as flxave value from restart dataset ',flxave_res
          call endrun
       endif
       if (flxave_res .and. .not. dosend) then
          write(6,*)'(RESTRD): assume that current flux coupled model with flux ', &
               'averaging must stop on a time step where dosend (doalb) is true'
          call endrun
       end if
    endif
#endif

! read case name 

    if (masterproc) then
       read (nio) casename
    endif
    
! -----------------------------------------------------------------
! If branch run - check that case name is different and
! return now (no history file data reads for branch run)
! -----------------------------------------------------------------

    if (nsrest == 3) then
       if (masterproc) then
          if (casename==caseid) then
             write(6,*) 'RESTRD ERROR: Must change case name on branch run'
             write(6,*) 'Prev case name ',trim(casename)
             write(6,*)' Current case name ',trim(caseid)
             call endrun
          end if
          call relavu (nio)
          write(6,'(72a1)') ("-",i=1,60)
          write(6,*) 'Successfully read restart data for branch run'
          write(6,*)
       endif
       RETURN
    end if

! -----------------------------------------------------------------
! If restart run - read history file related data 
! -----------------------------------------------------------------

! to reduce file size, only read history accumulators and counters if not
! end of history interval
#if 0 
    if (nsrest == 1) then
       do m = 1, nhist
          call readin (nio, ntim(m))
          call readin (nio, mcdate_i(m))
          call readin (nio, mcsec_i(m))
          call readin (nio, mdcur_i(m))
          call readin (nio, mscur_i(m))
          call readin (nio, nbeghis(m))
          call readin (nio, slfld%num(m))
          call readin (nio, mlsoifld%num(m))
          if (nbeghis(m) /= 1) then
             num = slfld%num(m)
             if (num > 0) then
                allocate(ibuf2d(num,1:lis%d%glbnch))
                call readin (nio, ibuf2d, num)
                do j = 1, num
                do i = 1, lis%d%glbnch
                   slfld%count(i,j,m) = ibuf2d(j,i)
                end do
                end do
                deallocate(ibuf2d)

                allocate( buf2d(num,1:lis%d%glbnch))
                call readin (nio, buf2d, num)
                do j = 1, num
                do i = 1, lis%d%glbnch
                   slfld%value(i,j,m) = buf2d(j,i)
                end do
                end do
                deallocate(buf2d)
             endif
             num = mlsoifld%num(m)
             if (num > 0) then
                allocate(ibuf3d(nlevsoi,num,1:lis%d%glbnch))
                call readin (nio, ibuf3d, nlevsoi, num)
                do j = 1, num
                do n = 1, nlevsoi
                do i = 1, lis%d%glbnch
                   mlsoifld%count(i,n,j,m) = ibuf3d(n,j,i)
                end do
                end do
                end do
                deallocate(ibuf3d)

                allocate( buf3d(nlevsoi,num,1:lis%d%glbnch))
                call readin (nio, buf3d, nlevsoi, num)
                do j = 1, num
                do n = 1, nlevsoi
                do i = 1, lis%d%glbnch
                   mlsoifld%value(i,n,j,m) = buf3d(n,j,i)
                end do
                end do
                end do
                deallocate(buf3d)
             endif
          end if 
       end do

! Read names of history files. If history file is not full, 
! open netCDF history file and set flag to obtain time 
! dependent netCDF variable id's 
! Note id's must be read in from routine histwrt since
! routine histini is called after the routine restrd

       if (masterproc) then
          do m = 1, nhist
             read (nio) fnameh(m)
          end do
          do m = 1, nhist
             if (ntim(m) /= 0) then
                ncgetid(m) = .true.
                call getfil (fnameh(m), locfnh(m), 0)
                call wrap_open (locfnh(m), nf_write, ncid(m))
             end if
          end do
       end if
#endif
! Close unit

       if (masterproc) then
          call relavu (nio)
          write(6,'(72a1)') ("-",i=1,60)
          write(6,*) 'Successfully read restart data for restart run'
          write(6,*)
       endif
          

    return
  end subroutine restrd

!=======================================================================

  subroutine restwrt ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Write CLM restart files
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: restFileMod.F90,v 1.16 2005/02/01 18:53:53 jim Exp $
!-----------------------------------------------------------------------


    use time_manager, only : get_nstep, timemgr_write_restart, is_last_step
    use lis_openfileMod, only : create_restart_filename, create_output_directory
    use tile_spmdMod
!    include 'netcdf.inc'

    implicit none
! ------------------------ local variables ------------------------
    integer :: i,j,k,m,n           !indices
    integer :: nio                 !fortran unit number
    integer :: num                 !number of fields (temporary)
    logical :: lremove             !true => remove file after archive
    character(len=256) :: rem_fn   !remote (archive) filename
    character(len=256) :: rem_dir  !remote (archive) directory
    character(len=256) :: loc_fn   !local restart filename
    character(len=256) :: filename !generic filename
    integer :: ierr                 !error code
#if (defined SPMD)                
    integer :: numrecvv(0:npes-1)  !vector of items to be received  
    integer :: displsv(0:npes-1)   !displacement vector
    integer :: numsend             !number of items to be sent
#endif                            
    integer :: alloc_size, local_size
! -----------------------------------------------------------------

! All processes may be calling this routine.
! The master process must allocate a temporary array that is large
! enough to hold all the land points.  The slave processes should only
! allocate a temporary array the same size as the number of land points
! they computed.
   if ( masterproc ) then
      alloc_size = lis%d%glbnch  ! all land points
   else
      alloc_size = di_array(iam) ! number of land points slave process processed
   endif
#if ( defined OPENDAP )
   local_size = di_array(iam) ! number of land points each process (master or
                              ! slave) processed
#else
   local_size = lis%d%glbnch  ! number of land points master process processed
#endif



! -----------------------------------------------------------------
! Open main restart file. write data. close. dispose to mass store
! -----------------------------------------------------------------
    if (masterproc) then
       write(6,*)
       write(6,'(72a1)') ("-",i=1,60)
!       write(6,*) 'nstep = ',get_nstep()

!       loc_fn = set_restart_filename()
       call create_output_directory()
       call create_restart_filename(loc_fn, '.CLM2rst')
       nio = getavu()
       call opnfil (loc_fn, nio, 'u')
    endif

! -----------------------------------------------------------------
! Write data
! -----------------------------------------------------------------

! main restart data
    call wrtout (nio, rest_id) 
#if (defined COUP_CAM)

    call wrtout (nio, get_nstep())

#else

    if (masterproc) call timemgr_write_restart(nio)

#endif
! write out required 1d fields

    allocate ( buf1d(alloc_size))
    allocate (ibuf1d(alloc_size))

    ibuf1d(1:local_size)= clm(1:local_size)%snl             
    call wrtout (nio,ibuf1d)

    ibuf1d(1:local_size)= clm(1:local_size)%frac_veg_nosno_alb  
    call wrtout (nio,ibuf1d)

    buf1d(1:local_size) = clm(1:local_size)%h2osno          
    call wrtout (nio, buf1d) 

    buf1d(1:local_size) = clm(1:local_size)%h2ocan          
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%snowdp          
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%snowage         
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%frac_sno        
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%t_veg           
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%t_grnd          
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%fwet            
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%tlai            
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%tsai            
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%elai            
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%esai            
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%fsun            
    call wrtout (nio, buf1d)  

    buf1d(1:local_size) = clm(1:local_size)%htop
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%hbot
    call wrtout (nio, buf1d)

#if (defined DGVM)

    buf1d(1:local_size) = clm(1:local_size)%wf
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%t_mo_min
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%annpsn
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%annpsnpot
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%fmicr
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%bm_inc
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%afmicr
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%t10min
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%tmomin20
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%agdd20
    call wrtout (nio, buf1d)

    ibuf1d(1:local_size) = clm(1:local_size)%itypveg
    call wrtout (nio, ibuf1d)

    buf1d(1:local_size) = clm(1:local_size)%fpcgrid
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%lai_ind
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%crownarea
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%dphen
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%leafon
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%leafof
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%firelength
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%litterag
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%litterbg
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%cpool_fast
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%cpool_slow
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%k_fast_ave
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%k_slow_ave
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%litter_decom_ave
    call wrtout (nio, buf1d)

    do k = 1, local_size
       if (clm(k)%present) then
          buf1d(k) = 1.0
       else
          buf1d(k) = 0.0
       end if
    end do
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%nind
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%lm_ind
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%sm_ind
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%hm_ind
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = clm(1:local_size)%rm_ind
    call wrtout (nio, buf1d)

    buf1d(1:local_size) = patchvec%wtxy(1:local_size)
    call wrtout (nio, buf1d)

#endif

    deallocate (ibuf1d)
    deallocate ( buf1d)

! write out multi level snow only fields

    allocate (buf2d(-nlevsno+1:0,1:local_size))
    do j = -nlevsno+1,0
       buf2d(j,1:local_size) = clm(1:local_size)%dz(j)  
    end do
    call wrtout (nio, buf2d, nlevsno)

    do j = -nlevsno+1,0
       buf2d(j,1:local_size) = clm(1:local_size)%z(j)  
    end do
    call wrtout (nio, buf2d, nlevsno)
    deallocate(buf2d)

    allocate (buf2d(-nlevsno:0,1:local_size))
    do j = -nlevsno,0
       buf2d(j,1:local_size) = clm(1:local_size)%zi(j)  
    end do
    call wrtout (nio, buf2d, nlevsno+1)
    deallocate(buf2d)

! write out multi level snow-soil fields

    allocate (buf2d(-nlevsno+1:nlevsoi,1:local_size))
    do j = -nlevsno+1,nlevsoi
       buf2d(j,1:local_size) = clm(1:local_size)%t_soisno(j)  
    end do
    call wrtout (nio, buf2d, nlevsoi+nlevsno)
    do j = -nlevsno+1,nlevsoi
       buf2d(j,1:local_size) = clm(1:local_size)%h2osoi_liq(j)  
    end do
    call wrtout (nio, buf2d, nlevsoi+nlevsno)
    do j = -nlevsno+1,nlevsoi
       buf2d(j,1:local_size) = clm(1:local_size)%h2osoi_ice(j)  
    end do
    call wrtout (nio, buf2d, nlevsoi+nlevsno)
    deallocate(buf2d)

!    allocate (buf2d(1:nlevlak,1:local_size))
!    do j = 1,nlevlak
!       buf2d(j,1:local_size) = clm(1:local_size)%t_lake(j)  
!    end do
!    call wrtout (nio, buf2d, nlevlak)
!    deallocate(buf2d)

! write out multi level albedo and surface radiation related fields

    allocate (buf2d(numrad,1:local_size))
!    do j = 1,numrad
!       buf2d(j,1:local_size) = clm(1:local_size)%albd(j)  
!    end do
!    call wrtout (nio, buf2d  , numrad)

!    do j = 1,numrad
!       buf2d(j,1:local_size) = clm(1:local_size)%albi(j)  
!    end do
!    call wrtout (nio, buf2d  , numrad)

    do j = 1,numrad
       buf2d(j,1:local_size) = clm(1:local_size)%albgrd(j)  
    end do
    call wrtout (nio, buf2d, numrad)

    do j=1,local_size
       if(clm(j)%kpatch == 38) print*, numrad, clm(j)%albgrd
    enddo

    do j = 1,numrad
       buf2d(j,1:local_size) = clm(1:local_size)%albgri(j)  
    end do
    call wrtout (nio, buf2d, numrad)

    do j = 1,numrad
       buf2d(j,1:local_size) = clm(1:local_size)%fabd(j)  
    end do
    call wrtout (nio, buf2d  , numrad)

    do j = 1,numrad
       buf2d(j,1:local_size) = clm(1:local_size)%fabi(j)  
    end do
    call wrtout (nio, buf2d  , numrad)

    do j = 1,numrad
       buf2d(j,1:local_size) = clm(1:local_size)%ftdd(j)  
    end do
    call wrtout (nio, buf2d  , numrad)

    do j = 1,numrad
       buf2d(j,1:local_size) = clm(1:local_size)%ftid(j)  
    end do
    call wrtout (nio, buf2d  , numrad)

    do j = 1,numrad
       buf2d(j,1:local_size) = clm(1:local_size)%ftii(j)  
    end do
    call wrtout (nio, buf2d  , numrad)
    deallocate(buf2d)

#if (defined ACCUM)
! write data for accumulation variables

    call wrtout (nio, naccflds)
    allocate(buf2d(naccflds,1:local_size))
    do j = 1, naccflds
       do i = 1, local_size
          buf2d(j,i) = accval(i,j)
       end do
    end do
    call wrtout (nio, buf2d, naccflds)
    deallocate(buf2d)
#endif

#if (defined RTM)

! write river routing data (if applicable)

!#if (defined SPMD)
!    call compute_mpigs_patch(1, numsend, numrecvv, displsv)
!    call mpi_gatherv (totrunin_ave(1), numsend , mpir8, &
!                      totrunin_ave          , numrecvv, displsv, mpir8, 0, mpicom, ierr)
!    call mpi_gatherv (prec_ave(1)    , numsend , mpir8, &
!                      prec_ave              , numrecvv, displsv, mpir8, 0, mpicom, ierr)
!    call mpi_gatherv (evap_ave(1)    , numsend , mpir8, &
!                      evap_ave              , numrecvv, displsv, mpir8, 0, mpicom, ierr)
!    call mpi_gatherv (qchocn2(1)     , numsend , mpir8, &
!                      qchocn2               , numrecvv, displsv, mpir8, 0, mpicom, ierr)
!    call mpi_gatherv (qchan2(1)      , numsend , mpir8, &
!                      qchan2                , numrecvv, displsv, mpir8, 0, mpicom, ierr)
!#endif
    if (masterproc) then
       write(nio) volr
       write(nio) ncount_rtm
       write(nio) ncount_global, yrold
       write(nio) prec_global,evap_global,runlnd_global,runrtm_global,volrtm_global,ocnrtm_global
       write(nio) (totrunin_ave(j),j=1,local_size)
       write(nio) (prec_ave(j)    ,j=1,local_size)
       write(nio) (evap_ave(j)    ,j=1,local_size)
       write(nio) (qchan2(j)      ,j=1,local_size)
       write(nio) (qchocn2(j)     ,j=1,local_size)
#if (defined COUP_CSM)
       write(nio) ocnrof_vec 
#endif
    endif

#endif

! write coupled model related data (if applicable)

#if (defined COUP_CSM)
    call wrtout (nio, csm_doflxave)
    call wrtout (nio, dosend)
#endif

! write case name

    if (masterproc) then
       write(nio) caseid
    endif

! -----------------------------------------------------------------
! Write history file related data
! -----------------------------------------------------------------
!
! Only needed for history files to restart on any time step. 
! NOTE - to reduce file size, only write history accumulators and
! history counters if not end of history interval.
! The following is NOT read in for branch runs
#if 0 
    do m = 1, nhist
       call wrtout (nio, ntim(m))
       call wrtout (nio, mcdate_i(m))
       call wrtout (nio, mcsec_i(m))
       call wrtout (nio, mdcur_i(m))
       call wrtout (nio, mscur_i(m))
       call wrtout (nio, nbeghis(m))
       call wrtout (nio, slfld%num(m))
       call wrtout (nio, mlsoifld%num(m))
       if (nbeghis(m) /= 1) then
          num = slfld%num(m)
          if (num > 0) then
             allocate(ibuf2d(num,1:local_size))
             do j = 1, num
             do i = 1, local_size
                ibuf2d(j,i) = slfld%count(i,j,m)
             end do
             end do
             call wrtout (nio, ibuf2d, num)
             deallocate(ibuf2d)

             allocate(buf2d(num,1:local_size))
             do j = 1, num
             do i = 1, local_size
                buf2d(j,i) = slfld%value(i,j,m)
             end do
             end do
             call wrtout (nio, buf2d, num)
             deallocate(buf2d)
          endif
          num = mlsoifld%num(m)
          if (num > 0) then
             allocate(ibuf3d(nlevsoi,num,1:local_size))
             do j = 1, num
             do n = 1, nlevsoi
             do i = 1, local_size
                ibuf3d(n,j,i) = mlsoifld%count(i,n,j,m)
             end do
             end do
             end do
             call wrtout (nio, ibuf3d, nlevsoi, num)
             deallocate(ibuf3d)

             allocate(buf3d(nlevsoi,num,1:local_size))
             do j = 1, num
             do n = 1, nlevsoi
             do i = 1, local_size
                buf3d(n,j,i) = mlsoifld%value(i,n,j,m) 
             end do
             end do
             end do
             call wrtout (nio, buf3d, nlevsoi, num)
             deallocate(buf3d)
          endif
       end if
    end do

! write name of current history file 

    if (masterproc) then
       do m = 1, nhist
          if (mss_irt == 0) then
             filename = locfnh(m)
          else
             rem_dir = trim(archive_dir) // '/hist/'
             filename = set_filename(rem_dir, locfnh(m))
          endif
          write(nio) filename
       end do
    endif
#endif
! -----------------------------------------------------------------
! Close and dispose restart file to mass store
! -----------------------------------------------------------------

    if (masterproc) then
!       lremove = .true.
       lremove = .false.
!#if (defined OFFLINE) || (defined COUP_CAM)
!       if (is_last_step()) lremove = .false.
!#elif (defined COUP_CSM)
!       if (csmstop_next) lremove = .false.
!#endif
       call relavu (nio)
       write(6,*) 'Successfully wrote local restart file ',trim(loc_fn)
       if (mss_irt > 0) then
          rem_dir = trim(archive_dir) // '/rest/'
          rem_fn = set_filename(rem_dir, loc_fn)
          call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, lremove)
       endif
    endif

#if 0
! -----------------------------------------------------------------
! Write restart pointer file
! -----------------------------------------------------------------

    if (masterproc) then
       call write_rest_pfile(loc_fn)
       write(6,'(72a1)') ("-",i=1,60)
       write(6,*)
    endif
#endif
    return
  end subroutine restwrt

!=======================================================================

  subroutine write_rest_pfile (loc_fn)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Open restart pointer file. Write names of current restart and 
! history files. If using mass store, these are the mass store
! names except if mss_irt=0 (no mass store files written). Close. 
! 
! Method: 
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

! ------------------------ arguments -----------------------------------
    character(len=*), intent(in) :: loc_fn !local restart filename
!-----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
    integer :: m                           !index
    integer :: nio                         !Fortran unit number
    character(len=256) :: filename         !local file name
    character(len=256) :: rem_dir          !remote directory
!-----------------------------------------------------------------------

    nio = getavu()
    !filename= trim(rpntdir) //'/'// trim(rpntfil)
    filename= trim(rpntdir) // trim(rpntfil)
    call opnfil (filename, nio, 'f')

    ! write name of restart file to pointer file 
    if (mss_irt == 0) then
       write(nio,'(a)') loc_fn
    else
       rem_dir = trim(archive_dir) // '/rest/'
       write(nio,'(a)') set_filename(rem_dir, loc_fn)
    endif

    ! add comments to pointer file of all files that are needed for restart
    write(nio,*)'The following lines list files needed for restart - do not edit'

    ! only write names of open history files
    do m = 1, nhist
       if (locfnh(m) /= ' ') then
          if (mss_irt == 0) then
             filename = locfnh(m)
          else
             rem_dir = trim(archive_dir) // '/hist/'
             filename = set_filename(rem_dir, locfnh(m))
          endif
          write(nio, '(a)') filename
       end if
    end do

    call relavu (nio)
    write(6,*)'Successfully wrote local restart pointer file'

  end subroutine write_rest_pfile

!=======================================================================

#if 0
! now using create_restart_filename
  character(len=256) function set_restart_filename ()

    use clm_varctl  ,  only : caseid, caseoutdir
    use time_manager,  only : get_curr_date
    use lisdrv_module, only : lis

! ------------------------ local variables ------------------------
    character(len=256) :: cdate       !date char string
    character(len=256) :: cyear       !year char string
    character(len=256) :: cnhdate     !no hour date char string
    character(len=256) :: fdir        !
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec, testsec           !seconds into current day
! -----------------------------------------------------------------

    call get_curr_date (lis%t, yr, mon, day, sec) 
!!    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
!!    set_restart_filename = "./"//trim(caseid)//".clm2.r."//trim(cdate)
     testsec = sec / 3600
     write(cdate,'(i4.4,i2.2,i2.2,i2.2)') yr,mon,day,testsec
     write(cnhdate,'(i4.4,i2.2,i2.2)') yr,mon,day
     write(cyear,'(i4)') yr
     
     fdir = "mkdir -p "//trim(caseoutdir)//"/EXP"//trim(caseid)//"/CLM2/"// &
            trim(cyear)//"/"//trim(cnhdate)
     call system(trim(fdir))
     
     set_restart_filename = trim(caseoutdir)//"/EXP"//trim(caseid)//    &
                            "/CLM2/"//trim(cyear)//"/"//trim(cnhdate)// &
                            "/LIS.E"//trim(caseid)//"."//trim(cdate)// &
                            ".CLM2rst"

  end function set_restart_filename
#endif

!=======================================================================
! BEGIN GENERIC PROCEDURE DEFINITIONS
!=======================================================================

  subroutine readinsc_int(iu, scalar)
!-----------------------------------------------------------------------
! Wrapper routine to read real scalar variable from restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in)  :: iu                 !input unit
    integer, intent(out) :: scalar             !scalar to read in
!---------------------------Local variables-----------------------------
    integer ierr               ! errorcode
!-----------------------------------------------------------------------
    if (masterproc) then
       read(iu, iostat=ierr) scalar
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!#if ( defined SPMD ) 
!    call mpi_bcast(scalar,1,mpiint,0,mpicom,ierr)
!#endif
    return
  end subroutine readinsc_int

!=======================================================================

  subroutine readinsc_log(iu, scalar)
!-----------------------------------------------------------------------
! Wrapper routine to read real scalar variable from restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in)  :: iu                 !input unit
    logical, intent(out) :: scalar             !scalar to read in
!---------------------------Local variables-----------------------------
    integer ierr               ! errorcode
!-----------------------------------------------------------------------
    if (masterproc) then
       read(iu, iostat=ierr) scalar
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!#if ( defined SPMD ) 
!    call mpi_bcast(scalar,1,mpilog,0,mpicom,ierr)
!#endif
    return
  end subroutine readinsc_log

!=======================================================================

  subroutine readinsc_real(iu, scalar)
!-----------------------------------------------------------------------
! Wrapper routine to read real scalar variable from restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in)  :: iu            !input unit
    real(r8), intent(out) :: scalar        !scalar to read in
!---------------------------Local variables-----------------------------
    integer ierr               ! errorcode
!-----------------------------------------------------------------------
    if (masterproc) then
       read(iu, iostat=ierr) scalar
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!#if ( defined SPMD ) 
!    call mpi_bcast(scalar,1,mpir8,0,mpicom,ierr)
!#endif
    return
  end subroutine readinsc_real

!=======================================================================

  subroutine readin1d_int (iu, iarr)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in) :: iu                        !input unit
    integer , intent(out):: iarr(1:lis%d%glbnch)   !read data
!---------------------------Local variables-----------------------------
    integer :: ierr            
#if (defined SPMD)
    integer :: i
    integer :: ibuf(lis%d%glbnch)
    integer :: numsendv(0:npes-1)    !vector of items to be sent
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numrecv               !number of items to be received
!-----------------------------------------------------------------------
    if (masterproc) then
       read (iu,iostat=ierr) (ibuf(i),i=1,lis%d%glbnch)
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!    call compute_mpigs_patch(1, numrecv, numsendv, displsv)
!    call mpi_scatterv (ibuf          , numsendv, displsv, mpiint, &
!                       ibuf(1), numrecv , mpiint , 0, mpicom, ierr)
    iarr(1:lis%d%glbnch) = ibuf(1:lis%d%glbnch)
#else 
    read (iu,iostat=ierr) iarr
    if (ierr /= 0 ) then
       write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine readin1d_int

!=======================================================================

  subroutine readin2d_int (iu, iarr, ndim1)
!-----------------------------------------------------------------------
! Wrapper routine to read integer array from restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in)  :: iu                            !input unit
    integer, intent(in)  :: ndim1                         !dimension
    integer, intent(out) :: iarr(ndim1,1:lis%d%glbnch) !read data
!---------------------------Local variables-----------------------------
    integer :: ierr              
#if (defined SPMD)
    integer :: i,j                 
    integer :: ibuf(ndim1,lis%d%glbnch)
    integer :: numsendv(0:npes-1)    !vector of items to be sent
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numrecv               !number of items to be received
!-----------------------------------------------------------------------
! Read in data master processor, scatter from master processor, copy to arr
    if (masterproc) then
       read (iu,iostat=ierr) ((ibuf(j,i),j=1,ndim1),i=1,lis%d%glbnch)
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!    call compute_mpigs_patch(ndim1, numrecv, numsendv, displsv)
!    call mpi_scatterv (ibuf            , numsendv, displsv, mpiint, &
!                       ibuf(1,1), numrecv ,  mpiint, 0, mpicom, ierr)
    iarr(:,1:lis%d%glbnch) = ibuf(:,1:lis%d%glbnch)
#else 
! Read in array directly
    read (iu,iostat=ierr) iarr
    if (ierr /= 0 ) then
       write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine readin2d_int

!=======================================================================

  subroutine readin3d_int (iu, iarr, ndim1, ndim2)
!-----------------------------------------------------------------------
! Wrapper routine to read integer arry from restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in)  :: iu                                  !input unit
    integer, intent(in)  :: ndim1                               !dimension
    integer, intent(in)  :: ndim2                               !dimension
    integer, intent(out) :: iarr(ndim1,ndim2,1:lis%d%glbnch) !read data
!---------------------------Local variables-----------------------------
    integer :: ierr               
#if (defined SPMD)
    integer :: i,j,k                 
    integer :: ibuf(ndim1,ndim2,lis%d%glbnch)
    integer :: numsendv(0:npes-1)    !vector of items to be sent
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numrecv               !number of items to be received
!-----------------------------------------------------------------------
! Read in data master processor, scatter from master processor, copy to arr
    if (masterproc) then
       read (iu, iostat=ierr)(((ibuf(k,j,i),k=1,ndim1),j=1,ndim2),i=1,lis%d%glbnch)
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!    call compute_mpigs_patch(ndim1*ndim2, numrecv, numsendv, displsv)
!    call mpi_scatterv (ibuf              , numsendv, displsv, mpiint, &
!                       ibuf(1,1,1), numrecv ,  mpiint, 0, mpicom, ierr)
    iarr(:,:,1:lis%d%glbnch) = ibuf(:,:,1:lis%d%glbnch)
#else 
! Read in data directly
    read (iu,iostat=ierr) iarr
    if (ierr /= 0 ) then
       write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine readin3d_int

!=======================================================================

  subroutine readin1d_real (iu, arr)
!-----------------------------------------------------------------------
! Wrapper routine to read real array from restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in) :: iu                        !input unit
    real(r8), intent(out):: arr(1:lis%d%glbnch)    !read data
!---------------------------Local variables-----------------------------
    integer :: ierr               ! errorcode
#if (defined SPMD)
    integer  :: i
    real(r8) :: buf(lis%d%glbnch)
    integer  :: numsendv(0:npes-1)    !vector of items to be sent
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numrecv               !number of items to be received
!-----------------------------------------------------------------------
! Read in data master processor, scatter from master processor, copy to arr
    if (masterproc) then
       read (iu,iostat=ierr) (buf(i),i=1,lis%d%glbnch)
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!    call compute_mpigs_patch(1, numrecv, numsendv, displsv)
!    call mpi_scatterv (buf          , numsendv, displsv, mpir8, &
!                       buf(1), numrecv , mpir8, 0, mpicom, ierr)
    arr(1:lis%d%glbnch) = buf(1:lis%d%glbnch)
#else 
! Read in data directly
    read (iu,iostat=ierr) arr
    if (ierr /= 0 ) then
       write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine readin1d_real

!=======================================================================

  subroutine readin2d_real (iu, arr, ndim1)
!-----------------------------------------------------------------------
! Wrapper routine to read real array from restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in)  :: iu                           !input unit
    integer , intent(in)  :: ndim1                        !dimension
    real(r8), intent(out) :: arr(ndim1,1:lis%d%glbnch) !read data
!---------------------------Local variables-----------------------------
    integer  :: ierr              
#if (defined SPMD)
    integer  :: i,j                 
    real(r8) :: buf(ndim1,lis%d%glbnch)
    integer  :: numsendv(0:npes-1)    !vector of items to be sent
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numrecv               !number of items to be received
!-----------------------------------------------------------------------
! Read in data master processor, scatter from master processor, copy to arr
    if (masterproc) then
       read (iu,iostat=ierr) ((buf(j,i),j=1,ndim1),i=1,lis%d%glbnch)
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!    call compute_mpigs_patch(ndim1, numrecv, numsendv, displsv)
!    call mpi_scatterv (buf            , numsendv, displsv, mpir8, &
!                       buf(1,1), numrecv , mpir8, 0, mpicom, ierr)
    arr(:,1:lis%d%glbnch) = buf(:,1:lis%d%glbnch)
#else 
! Read in array directly
    read (iu,iostat=ierr) arr
    if (ierr /= 0 ) then
       write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine readin2d_real

!=======================================================================

  subroutine readin3d_real (iu, arr, ndim1, ndim2)
!-----------------------------------------------------------------------
! Wrapper routine to read real array fom restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in)  :: iu              !input unit
    integer , intent(in)  :: ndim1           !dimension
    integer , intent(in)  :: ndim2           !dimension
    real(r8), intent(out) :: arr(ndim1,ndim2,1:lis%d%glbnch) !read data
!---------------------------Local variables-----------------------------
    integer :: ierr               
#if (defined SPMD)
    integer  :: i,j,k                 
    real(r8) :: buf(ndim1,ndim2,lis%d%glbnch)
    integer  :: numsendv(0:npes-1)    !vector of items to be sent
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numrecv               !number of items to be received
!-----------------------------------------------------------------------
! Read in data master processor, scatter from master processor, copy to arr
    if (masterproc) then
       read (iu, iostat=ierr)(((buf(k,j,i),k=1,ndim1),j=1,ndim2),i=1,lis%d%glbnch)
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
!    call compute_mpigs_patch(ndim1*ndim2, numrecv, numsendv, displsv)
!    call mpi_scatterv (buf              , numsendv, displsv, mpir8, &
!                       buf(1,1,1), numrecv , mpir8, 0, mpicom, ierr)
    arr(:,:,1:lis%d%glbnch) = buf(:,:,1:lis%d%glbnch)
#else 
! Read in data directly
    read (iu,iostat=ierr) arr
    if (ierr /= 0 ) then
       write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine readin3d_real

!=======================================================================

  subroutine wrtoutsc_int(iu, scalar)
!-----------------------------------------------------------------------
! Wrapper routine to read real scalar variable from restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in) :: iu             !input unit
    integer, intent(in) :: scalar         !scalar to read in
!---------------------------Local variables-----------------------------
    integer ierr                           !errorcode
!-----------------------------------------------------------------------
    if (masterproc) then
       write(iu, iostat=ierr) scalar
       if (ierr /= 0 ) then
          write(6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
    return
  end subroutine wrtoutsc_int

!=======================================================================

  subroutine wrtoutsc_log(iu, scalar)
!-----------------------------------------------------------------------
! Wrapper routine to read real scalar variable from restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in) :: iu             !input unit
    logical, intent(in) :: scalar         !scalar to read in
!---------------------------Local variables-----------------------------
    integer ierr                           !errorcode
!-----------------------------------------------------------------------
    if (masterproc) then
       write(iu, iostat=ierr) scalar
       if (ierr /= 0 ) then
          write(6,*) 'READIN ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
    return
  end subroutine wrtoutsc_log

!=======================================================================

  subroutine wrtoutsc_real(iu, scalar)
!-----------------------------------------------------------------------
! Wrapper routine to read real scalar variable from restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in) :: iu           !input unit
    real(r8), intent(in) :: scalar       !scalar to read in
!---------------------------Local variables-----------------------------
    integer ierr                          !errorcode
!-----------------------------------------------------------------------
    if (masterproc) then
       write(iu, iostat=ierr) scalar
       if (ierr /= 0 ) then
          write(6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
    return
  end subroutine wrtoutsc_real

!=======================================================================

  subroutine wrtout1d_int (iu, iarr)
    use spmdMod
    use tile_spmdMod

    implicit none

!-----------------------------------------------------------------------
! Wrapper routine to write integer array to restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in)  :: iu                     !input unit
    integer, dimension(*), intent(inout)  :: iarr  !output data
!---------------------------Local variables-----------------------------
    integer :: ierr            ! errorcode
#if ( defined OPENDAP )
!-----------------------------------------------------------------------
#if ( defined SPMD )
    if ( npes > 1 ) then
       call MPI_GATHERV(iarr(1:di_array(iam)), di_array(iam),             & 
                        MPI_INTEGER, iarr, di_array, displs, MPI_INTEGER, & 
                        0, MPI_COMM_WORLD, ierr)
    endif

#endif
    if (masterproc) then
       write (iu, iostat=ierr) iarr(1:lis%d%glbnch)
       if (ierr /= 0 ) then
          write(6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
#else 
! Write out array directly
    write (iu,iostat=ierr) iarr(1:lis%d%glbnch)
    if (ierr /= 0 ) then
       write(6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout1d_int

!=======================================================================

  subroutine wrtout2d_int (iu, iarr, ndim1)
    use spmdMod
    use tile_spmdMod

    implicit none
!-----------------------------------------------------------------------
! Wrapper routine to write integer array to restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in) :: iu                            !input unit
    integer, intent(in) :: ndim1                         !dimension
    integer, dimension(ndim1,*), intent(in) :: iarr         !output data
!---------------------------Local variables-----------------------------

    integer :: ierr               ! errorcode
#if ( defined OPENDAP )
    integer, dimension(ndim1,lis%d%glbnch) :: ibuf2d
    integer, allocatable, dimension(:) :: ibuf
    integer :: i,j

!-----------------------------------------------------------------------
#if ( defined SPMD )
    if ( npes > 1 ) then
    
       if ( masterproc ) then
          allocate(ibuf(lis%d%glbnch))
       else
          allocate(ibuf(di_array(iam)))
       endif
       do i = 1, ndim1
          ibuf(1:di_array(iam)) = iarr(i,1:di_array(iam))
          call MPI_GATHERV(ibuf, di_array(iam),           &
                           MPI_INTEGER, ibuf, di_array, displs, MPI_INTEGER, &
                           0, MPI_COMM_WORLD, ierr)
   
          if ( masterproc ) then
             ibuf2d(i,:) = ibuf
          endif
       enddo
    endif
#endif

    if (masterproc) then
       if ( npes > 1 ) then
          write (iu,iostat=ierr) ((ibuf2d(j,i),j=1,ndim1),i=1,lis%d%glbnch)
       else
          write (iu,iostat=ierr) iarr(1:ndim1,1:lis%d%glbnch)
       endif
       if (ierr /= 0 ) then
          write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
    deallocate(ibuf)
#else 
! Write out array directly
    write (iu,iostat=ierr) iarr(1:ndim1,1:lis%d%glbnch)
    if (ierr /= 0 ) then
       write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout2d_int

!=======================================================================

  subroutine wrtout3d_int (iu, iarr, ndim1, ndim2)
    use spmdMod
    use tile_spmdMod

    implicit none
!-----------------------------------------------------------------------
! Wrapper routine to write integer array to restart binary file 
!------------------------------Arguments--------------------------------
    integer, intent(in) :: iu                                  !input unit
    integer, intent(in) :: ndim1                               !dimension
    integer, intent(in) :: ndim2                               !dimension
    integer, dimension(ndim1,ndim2,*), intent(in) :: iarr         !output data
!---------------------------Local variables-----------------------------
    integer :: ierr
#if ( defined OPENDAP )
    integer, dimension(ndim1,ndim2,lis%d%glbnch) :: ibuf3d
    integer, allocatable, dimension(:) :: ibuf
    integer :: i,j,k
!-----------------------------------------------------------------------
#if ( defined SPMD )
    if ( npes > 1 ) then
    
       if ( masterproc ) then
          allocate(ibuf(lis%d%glbnch))
       else
          allocate(ibuf(di_array(iam)))
       endif
       do i = 1, ndim1
          do j = 1, ndim2
             ibuf(1:di_array(iam)) = iarr(i,j,1:di_array(iam))
             call MPI_GATHERV(ibuf, di_array(iam),                          &
                              MPI_INTEGER,ibuf,di_array,displs,MPI_INTEGER, &
                              0, MPI_COMM_WORLD, ierr)
   
             if ( masterproc ) then
                ibuf3d(i,j,:) = ibuf
             endif
          enddo
       enddo
    endif
#endif

    if (masterproc) then
       if ( npes > 1 ) then
          write (iu,iostat=ierr)(((ibuf3d(k,j,i),k=1,ndim1),j=1,ndim2),i=1,lis%d%glbnch)
       else
          write (iu,iostat=ierr) iarr(1:ndim1,1:ndim2,1:lis%d%glbnch)
       endif
       if (ierr /= 0 ) then
          write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
    deallocate(ibuf)
#else 
! Write out array directly
    write (iu,iostat=ierr) iarr(1:ndim1,1:ndim2,1:lis%d%glbnch)
    if (ierr /= 0 ) then
       write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout3d_int

!=======================================================================

  subroutine wrtout1d_real (iu, arr)
    use spmdMod
    use tile_spmdMod

    implicit none
!-----------------------------------------------------------------------
! Wrapper routine to write real array to restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in)  :: iu                     !input unit
    real, dimension(*), intent(inout)  :: arr       !output data
!---------------------------Local variables-----------------------------
    integer :: ierr            
#if ( defined OPENDAP )
!-----------------------------------------------------------------------
#if ( defined SPMD )
    if ( npes > 1 ) then
       call MPI_GATHERV(arr(1:di_array(iam)), di_array(iam),       & 
                        MPI_REAL, arr, di_array, displs, MPI_REAL, & 
                        0, MPI_COMM_WORLD, ierr)
    endif
#endif

    if (masterproc) then
       write (iu, iostat=ierr) arr(1:lis%d%glbnch)
       if (ierr /= 0 ) then
          write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
#else 
! Write out array directly
    write (iu,iostat=ierr) arr(1:lis%d%glbnch)
    if (ierr /= 0 ) then
       write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout1d_real

!=======================================================================

  subroutine wrtout2d_real (iu, arr, ndim1)
    use spmdMod
    use tile_spmdMod

    implicit none
!-----------------------------------------------------------------------
! Wrapper routine to write real array to restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in) :: iu                           !input unit
    integer , intent(in) :: ndim1                        !dimension
    real, dimension(ndim1,*), intent(in) :: arr              !output data
!---------------------------Local variables-----------------------------
    integer :: ierr
#if ( defined OPENDAP )
    real, dimension(ndim1,lis%d%glbnch) :: buf2d
    real, allocatable, dimension(:) :: buf
    integer  :: i,j
!-----------------------------------------------------------------------
#if ( defined SPMD )
    if ( npes > 1 ) then
    
       if ( masterproc ) then
          allocate(buf(lis%d%glbnch))
       else
          allocate(buf(di_array(iam)))
       endif
       do i = 1, ndim1
          buf(1:di_array(iam)) = arr(i,1:di_array(iam))
          call MPI_GATHERV(buf, di_array(iam),           &
                           MPI_REAL, buf, di_array, displs, MPI_REAL, &
                           0, MPI_COMM_WORLD, ierr)
   
          if ( masterproc ) then
             buf2d(i,:) = buf
          endif
       enddo
    endif
#endif
    if (masterproc) then
       if ( npes > 1 ) then
          write (iu,iostat=ierr) ((buf2d(j,i),j=1,ndim1),i=1,lis%d%glbnch)
       else
          write (iu,iostat=ierr) arr(1:ndim1,1:lis%d%glbnch)
       endif
       if (ierr /= 0 ) then
          write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
    deallocate(buf)
#else 
! Write out array directly
    write (iu,iostat=ierr) arr(1:ndim1,1:lis%d%glbnch)
    if (ierr /= 0 ) then
       write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout2d_real

!=======================================================================

  subroutine wrtout3d_real (iu, arr, ndim1, ndim2)
    use spmdMod
    use tile_spmdMod

    implicit none
!-----------------------------------------------------------------------
! Wrapper routine to write real array to restart binary file 
!------------------------------Arguments--------------------------------
    integer , intent(in) :: iu                     !input unit
    integer , intent(in) :: ndim1                  !dimension
    integer , intent(in) :: ndim2                  !dimension
    real, dimension(ndim1,ndim2,*), intent(in) :: arr      !output data
!---------------------------Local variables-----------------------------
    integer :: ierr              
#if ( defined OPENDAP )
    integer, dimension(ndim1,ndim2,lis%d%glbnch) :: buf3d
    integer, allocatable, dimension(:) :: buf
    integer :: i,j,k
!-----------------------------------------------------------------------
#if ( defined SPMD )
    if ( npes > 1 ) then
    
       if ( masterproc ) then
          allocate(buf(lis%d%glbnch))
       else
          allocate(buf(di_array(iam)))
       endif
       do i = 1, ndim1
          do j = 1, ndim2
             buf(1:di_array(iam)) = arr(i,j,1:di_array(iam))
             call MPI_GATHERV(buf, di_array(iam),                    &
                              MPI_REAL,buf,di_array,displs,MPI_REAL, &
                              0, MPI_COMM_WORLD, ierr)
   
             if ( masterproc ) then
                buf3d(i,j,:) = buf
             endif
          enddo
       enddo
    endif
#endif
    if (masterproc) then
       if ( npes > 1 ) then
          write (iu,iostat=ierr)(((buf3d(k,j,i),k=1,ndim1),j=1,ndim2),i=1,lis%d%glbnch)
       else
          write (iu,iostat=ierr) arr(1:ndim1,1:ndim2,1:lis%d%glbnch)
       endif
       if (ierr /= 0 ) then
          write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
          call endrun
       end if
    endif
    deallocate(buf)
#else 
! Write out array directly
    write (iu,iostat=ierr) arr(1:ndim1,1:ndim2,1:lis%d%glbnch)
    if (ierr /= 0 ) then
       write (6,*) 'WRTOUT ierror ',ierr,' on i/o unit = ',iu
       call endrun
    end if
#endif
    return
  end subroutine wrtout3d_real

!=======================================================================
! END GENERIC PROCEDURE DEFINITIONS
!=======================================================================

end module restFileMod


