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
! !MODULE: clm2drv_module.F90
!
! !DESCRIPTION:
!  Module for runtime specific CLM2 variables
!
! !REVISION HISTORY:
!
! 14 Oct 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
module clm2drv_module
! !ARGUMENTS:
  type clmdrvdec
     integer :: numout            !Counts number of output times for CLM2
     integer :: clm2open          ! Keeps track of opening files
     integer :: varid(29)
     character*40 :: clm2_rfile   !CLM2 Active restart file
     character*40 :: clm2_vfile   !CLM2 Vegetation Tile Specification File
     character*40 :: clm2_chtfile !CLM2 Canopy Heights File
     real :: clm2_ism             !CLM2 intial soil moisture
     real :: clm2_it              !CLM2 initial soil temperature
     real :: clm2_iscv            !CLM2 Initial snow mass
     real :: writeintc2           !CLM2 Output Interval (hours)	
  end type clmdrvdec
!EOP  
end MODULE clm2drv_module
