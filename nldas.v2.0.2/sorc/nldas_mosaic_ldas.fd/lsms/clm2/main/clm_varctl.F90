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

module clm_varctl

!<debug>
!  use precision
  use precision, only : r8
!  use clm_varpar
  use clm_varpar, only : maxhist, maxalflds
  use clm2drv_module
!</debug>
  implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Module of run control variables 
!
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: clm_varctl.F90,v 1.8 2004/11/24 22:57:11 jim Exp $
!-----------------------------------------------------------------------

! Run control variables

  character(len=16) :: caseid         !case id
  character(len=80) :: ctitle         !case title
  character(len=40) :: caseoutdir     !case output dir
  integer :: nsrest                   !0: initial run. 1: restart: 3: branch

! History file variables
  type (clmdrvdec) :: clmdrv  

  logical           :: hist_dov2xy(maxhist)     !true: grid-average history field. false: vector
  integer           :: hist_mfilt(maxhist)      !max number of time samples per history file
  character(len= 8) :: hist_fldadd(maxalflds)   !names of fields to change to active
  character(len= 8) :: hist_chntyp(2,maxalflds) !paired field names and type for namelist
  integer           :: hist_nhtfrq(maxhist)     !history interval (iterations) (0=> monthly average)
  character(len= 8) :: hist_crtinic             !if set to 'MONTHLY' or 'YEARLY', write initial cond. file

! Long term archive variables

  character(len=256):: archive_dir              !long term archive directory (can be mass store)
  character(len= 8) :: mss_wpass                !mass store write password for output files
  integer           :: mss_irt                  !mass store retention period

! Run input files

  character(len=256) :: finidat                 !initial conditions file name
  character(len=256) :: fsurdat                 !surface data file name
  character(len=256) :: fpftcon                 !ASCII data file with PFT physiological constants
  character(len=256) :: nrevsn                  !restart data file name for branch run
  character(len=256) :: frivinp_rtm             !RTM input data file name
  character(len=256) :: offline_atmdir          !directory for input offline model atm data forcing files (Mass Store ok)

! Files and logical variables for generating surface dataset

  character(len=256) :: mksrf_offline_fgrid    !land grid file name to use instead of generating grid
  character(len=256) :: mksrf_offline_fnavyoro !directory for 20 min navy orography dataset
  real(r8)           :: mksrf_offline_edgen    !northern edge of grid (degrees): >  -90 and < 90
  real(r8)           :: mksrf_offline_edgee    !eastern edge of grid (degrees) : see following notes
  real(r8)           :: mksrf_offline_edges    !southern edge of grid (degrees): >= -90 and <  90
  real(r8)           :: mksrf_offline_edgew    !western edge of grid (degrees) : see following notes
  character(len=256) :: mksrf_fvegtyp          !when making [fsurdat]: vegetation data file name
  character(len=256) :: mksrf_fsoitex          !when making [fsurdat]: soil texture data file name
  character(len=256) :: mksrf_fsoicol          !when making [fsurdat]: soil color data file name
  character(len=256) :: mksrf_flanwat          !when making [fsurdat]: inland water data file name
  character(len=256) :: mksrf_furban           !when making [fsurdat]: urban data file name
  character(len=256) :: mksrf_fglacier         !when making [fsurdat]: glacier data file name
  character(len=256) :: mksrf_flai             !when making [fsurdat]: lai data filename
  character(len=256) :: mksrf_firr             !obsolete

! Physics

  integer :: irad         !solar radiation frequency (iterations)
  logical :: conchk       !true => turn on error energy and water conservation check
  logical :: wrtdia       !true => write global average diagnostics to std out
  logical :: csm_doflxave !true => only communicate with flux coupler on albedo calc time steps

! Rtm control variables
 
  integer :: rtm_nsteps   !if > 1, average rtm over rtm_nsteps time steps  

! Derived variables (run, history and restart file)

  integer :: ncprec                      !output precision (set to nc_float or nc_double)
  character(len=256) :: locfnh(maxhist)  !local history file names                    
  character(len=256) :: rpntdir          !directory name for local restart pointer file
  character(len=256) :: rpntfil          !file name for local restart pointer file
  character(len=80) ::  version          !model version number                       

end module clm_varctl
