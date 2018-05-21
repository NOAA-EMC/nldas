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

module initializeMod

    use spmdMod, only : masterproc,iam

!=======================================================================
CONTAINS
!=======================================================================

  subroutine initialize()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Land model initialization. 
!
! Method: 
! Initialization routine for land surface model. Initialize land surface
! variables for atmosphere model if not a continuation run.
! Model constants are set in block data subprograms. This subroutine: 
!   o Initializes run control variables via the [clmexp] namelist 
!   o Initializes active Fortran unit numbers 1 to 99, except 5 and 6, to false
!   o Reads surface data on [lsmlon] x [lsmlat] grid 
!   o Defines the multiple plant types and fraction areas for each surface type 
!   o Builds the appropriate subgrid <-> grid mapping indices and weights 
!   o Sets up parallel processing by dividing the "big" vector of [numpatch] 
!     subgrid patches into "little" vectors of [sizpatchlit] patches
!   o Assigns subgrid patches the appropriate time-constant data, i.e., 
!     initializes the time constant variables
!   o Initializes history file variables
!   o Initializes river basin data
!   o Initializes accumulation variables
! For a restart or branch run only:
!   o Reads restart data 
! For an initial run only:
!   o Reads initial data and initializes the time variant variables
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: initializeMod.F90,v 1.17 2004/11/24 23:32:24 jim Exp $
!-----------------------------------------------------------------------

!<debug>
!    use precision
!</debug>
    use lisdrv_module, only :lis, grid, gindex,tile
    use clm_varcon
!<debug this module is not showing in stack frame pane>
!    use clm_varpar     , only : lsmlon, lsmlat, maxpatch
!</debug>
    use clm_varmap     , only : numland, numpatch, begpatch, endpatch
!    use clm_varsur     , only : lsmedge, numlon
    use clm_varctl     , only : nsrest
!    use clm_varctl     , only : fsurdat, finidat, nsrest, irad, &
!                                mksrf_offline_fgrid, mksrf_offline_fnavyoro
    use controlMod     , only : control_init, control_print
    use fileutils      , only : lsmiou
!    use mksrfdatMod    , only : mksrfdat
!    use surfFileMod    , only : surfrd
    use pftcFileMod    , only : pftconrd
!    use DustEmissionMod, only : Dustini
!    use mvegFileMod    , only : interpMonthlyVeg
!    use histFileMod    , only : histini, histvar_ini
!    use restFileMod    , only : restrd
!     use atmdrvMod    
!    use lis_module     ! LDAS variables
!    use tile_module     ! Tile variables
!    use grid_module     ! Grid variables
   use lis_indices_module    
    implicit none

!     type (tiledec) tile(lis%d%nch)
! ------------------------ arguments ----------------------------------
    integer  :: i,j,k                 !indices
    integer  :: ncsec                 !current time of day [seconds]
!<debug>
!    integer  :: vegxy(lis%d%nc,lis%d%nr,maxpatch) !vegetation type
!    real(r8) :: wtxy(lis%d%nc,lis%d%nr,maxpatch)  !subgrid weights
!using save to force wtxy and vegxy off of the stack>
!    integer,save  :: vegxy(7200,3000,17) !vegetation type
!    real(r8),save :: wtxy(7200,3000,17)  !subgrid weights
    integer, allocatable :: vegxy(:,:,:)
    real(r8), allocatable :: wtxy(:,:,:)
    integer  :: ierr
!</debug>
    integer  :: nstep                        !current time step
    real*8   :: temptime
    real     :: tempgmt, tempmnf
    integer  :: tempdoy,nyr,nmth,nday,nhr,nmin,nsec,tempts
    integer  :: tempyr,tempmo,tempda,temphr,tempmn,tempss

! ----------------------------------------------------------------------
! echo initialization to standard output
       call header()

       if (masterproc) then
       write (6,*) 'Attempting to initialize the land model .....'
       write (6,*)
       endif
       
! ----------------------------------------------------------------------
! Initialize run control variables, time manager, timestep
! ----------------------------------------------------------------------

    call control_init (lis)
  
    if (masterproc) call control_print()
    
! ----------------------------------------------------------------------
! Initialize Fortran unit numbers 1 to 99 to inactive: except standard 
! input (5) and standard output (6) 
! ----------------------------------------------------------------------

    lsmiou(:) = .false.
    lsmiou(5) = .true.
    lsmiou(6) = .true.

    if (masterproc) then
       write (6,*) 'Preset Fortran unit numbers:'
       write (6,*) '   unit  5 = standard input'
       write (6,*) '   unit  6 = standard output'
    endif


! ----------------------------------------------------------------------
! Read list of PFTs and their corresponding parameter values
! ----------------------------------------------------------------------
    call pftconrd ()
!    allocate(vegxy(nc_dom,nr_dom,17))
!    allocate(wtxy(nc_dom,nr_dom,17))
    !if ( masterproc ) then
!       call tile2clm2(lis,grid,gindex,wtxy,vegxy)
    !endif
! ----------------------------------------------------------------------
! Use [wtxy] to build mapping indices and weights: 
! [lsmlon] x [lsmlat] grid <-> [numland] vector of land points <-> 
! [numpatch] vector of subgrid points
! ----------------------------------------------------------------------
    !if ( masterproc ) then
       call clm_map (gindex,lis_nc_working,lis_nr_working)
    !endif
!    deallocate(wtxy)
! ----------------------------------------------------------------------
! Initialize time invariant variables as subgrid vectors of length [numpatch] 
! ----------------------------------------------------------------------
    if (masterproc) then
       write (6,*) ('Attempting to initialize time invariant variables')
    endif

!=== LDAS modification: Passed in lis and grid modules
    !call iniTimeConst (vegxy, lis,grid,gindex)
    call iniTimeConst ()
!    deallocate(vegxy)
    if (masterproc) then
       write (6,*) ('Successfully initialized time invariant variables')
       write (6,*)
    endif

  end subroutine initialize

!=======================================================================

  subroutine header

!----------------------------------------------------------------------- 
! 
! Purpose: 
! echo and save model version number
!
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

    use precision
    use clm_varctl, only : version
    implicit none

    version = 'CLM MODEL version 2.0'

    if ( masterproc )then
      write (6,*) trim(version)
      write (6,*)
    end if

    return
  end subroutine header

!=======================================================================

end module initializeMod
