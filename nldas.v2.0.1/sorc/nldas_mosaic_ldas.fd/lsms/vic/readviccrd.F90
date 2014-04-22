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
! !ROUTINE: readviccrd.F90
!
! !DESCRIPTION:
!  Routine to read Vic specific parameters from the card file. 
!
! !REVISION HISTORY:
! 14 Oct 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readviccrd(vicdrv)
! !USES:
  use vicdrv_module
!EOP
  implicit none
  integer :: lsm
  type(vicdrvdec) :: vicdrv
  namelist /vic/vicdrv
!BOC
  open(11,file='lis.crd',form='formatted',status='old')
  read(unit=11,NML=vic)
  PRINT*,'Running VIC LSM:'
  PRINT*,'VIC Soil File: ', vicdrv%VIC_SFILE
  vicdrv%vic_quick_flux = 1 !default
  vicdrv%vic_full_energy  = 1
!  vicdrv%vic_frozen_soil  = 1 moved into card file
  if(vicdrv%vic_full_energy==1) then
     vicdrv%vic_quick_flux = 1
     vicdrv%vic_grnd_flux = 1
  endif
  if(vicdrv%vic_frozen_soil==1) then
     vicdrv%vic_quick_flux = 0
     vicdrv%vic_grnd_flux = 1
  endif
! Justin: number of nodes needs to be reduced if the following are true
  if(vicdrv%vic_frozen_soil==0 .and. vicdrv%vic_quick_flux == 1) then
     vicdrv%vic_nnode = 3
  endif
! Justin: end of changes
  vicdrv%vicopen = 0


  close(11)
!EOC
end subroutine readviccrd
