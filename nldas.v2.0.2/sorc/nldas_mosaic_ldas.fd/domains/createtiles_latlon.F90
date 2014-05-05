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
! !ROUTINE: createtiles_latlon
!
! !DESCRIPTION:
!  This primary goal of this routine is to determine tile space for 
!  a lat/lon domain
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Major F90 and major structure revision
!  3  Jan 2000: Minor T=0 bug fix, should have no effect on output
!  8  Mar 2000: Brian Cosgrove; Initialized FGRD to 0 For Dec Alpha Runs
! 22  Aug 2000: Brian Cosgrove; Altered code for US/Mexico/Canada Mask
! 04  Feb 2001: Jon Gottschalck; Added option to read and use Koster tile space
! 17  Oct 2003: Sujay Kumar ; Initial version of subsetting code
!
! !INTERFACE:
subroutine createtiles_latlon()
! !USES:
  use lisdrv_module, only: lis
  use grid_module
  use spmdMod
!EOP
  IMPLICIT NONE
  real, allocatable :: elevdiff(:, :)
!=== Local Variables =====================================================
  integer ::c,r,t,i,j,count  ! Loop counters
  real, allocatable :: VEG(:,:,:) !Temporary vegetation processing variable
  real :: isum
  real, allocatable :: tmpelev(:)
  integer :: landnveg
  real,allocatable :: tsum(:,:)  !Temporary processing variable
  real,allocatable :: fgrd(:,:,:)

  integer :: ios1,mlat,mlon,line,glnc,glnr
  integer :: line1,line2
  integer :: nc_dom

  integer :: ierr
  integer :: gnc, gnr
  integer :: cindex, rindex
  
  real, allocatable :: localmask(:,:)
  real :: locallat
  real :: locallon
  
!=== End Variable Definition =============================================
!BOC

	print *,'here brian 11'
  if ( masterproc ) then

     allocate(localmask(lis%d%lnc,lis%d%lnr))
     call readlandmask(lis%d%landcover, localmask)

     allocate(elevdiff(lis%d%lnc, lis%d%lnr), stat=ierr)
     call check_error(ierr,'Error allocating elev diff.',iam)
     
     call readelevdiff(lis%d%elev, elevdiff)

     allocate(fgrd(lis%d%lnc,lis%d%lnr,lis%p%nt), stat=ierr)
     call check_error(ierr,'Error allocating fgrd.',iam)

     call readlandcover(lis%d%landcover, fgrd)

     allocate(tsum(lis%d%lnc, lis%d%lnr), stat=ierr)
     call check_error(ierr,'Error allocating tsum.',iam)
     tsum = 0.0
     
     call calculate_domveg(fgrd, tsum)
     	print *,'calling vegtilespace' 
     call create_vegtilespace(fgrd, tsum, localmask, elevdiff)
	print *,'called vegtilsapce brian'

     deallocate(elevdiff)

     deallocate(localmask, stat=ierr)
     call check_error(ierr,'Error allocating glbmask',iam)      
     deallocate(fgrd, stat=ierr)
     call check_error(ierr,'Error allocating glbfgrd',iam)
     deallocate(tsum, stat=ierr)
     call check_error(ierr,'Error allocating glbtsum.',iam)
     
     write(*,*) 'msg: createtiles_latlon -- actual number of tiles:', & 
          lis%d%glbnch,' (',iam,')'
     write(*,*)
     
     write(*,*) 'msg: createtiles_latlon -- size of grid dimension:', & 
          lis%d%glbngrid,' (',iam,')'
      
  endif
  print*,'MSG: createtiles_latlon -- done',' (',iam,')'   
  return
!EOC
end subroutine createtiles_latlon
