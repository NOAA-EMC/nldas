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
!BOP
! 
! !ROUTINE: vic_readrestarst.F90
!
! !DESCRIPTION:
!  This program reads restart files for VIC.  
!
! !REVISION HISTORY:
! 
! 17 Nov 2003; Justin Sheffield : Initial Version
! 21 Nov 2003; Sujay Kumar : Added the time manager restart. 
! 
! !INTERFACE:
subroutine vic_readrestart
! !USES:
  use lisdrv_module, only: lis,tile
  use vic_varder, only : vicdrv
  !use time_manager
  implicit none
  !integer :: curSec
!EOP
!BOC
  if(lis%o%startcode .eq. 1) then
!<kluge>
! We no longer write time info into the restart file
!     OPEN(40,FILE='time.rst',FORM='unformatted')
!     write(*,*) "read time.rst"
!     call timemgr_read_restart(40)
!     write(*,*) "read time.rst"
!     call timemgr_restart()
!     write(*,*) "read time.rst"
!</kluge>
!     call get_curr_date(lis%t%yr,lis%t%mo,lis%t%da,curSec)
!     call sec2time(curSec,lis%t%hr,lis%t%mn,lis%t%ss)
!     call updatetime(lis%t) 
     call read_initial_model_state(trim(vicdrv%VIC_RFILE), & 
          len(trim(vicdrv%VIC_RFILE)), &
          vicdrv%vic_snowband, vicdrv%vic_nlayer, vicdrv%vic_nnode)

    call initialize_model_state2(vicdrv%vic_nlayer,     &
                                 vicdrv%vic_snowband,   &
                                 vicdrv%vic_nnode,      &
                                 vicdrv%vic_quick_flux, &
                                 vicdrv%vic_grnd_flux,  &
                                 vicdrv%vic_frozen_soil)
!<kluge>
! We no longer write time info into the restart file
!     close(40)
!</kluge>
  endif
!EOC
end subroutine vic_readrestart
