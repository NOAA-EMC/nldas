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
! !ROUTINE: recmwfmask.F90: 
!
! !DESCRIPTION:
! Reads in land-sea mask for 0.5 degree Reanalysis ECMWF data set,
! changes grid from ECMWF convention to GLIS convention by calling
! Recmwfgrid_2_gldasgrid, and transforms grid array into 1D vector 
! for later use in ret_reanlecmwf.F90.  
!
! !REVISION HISTORY:
!  16 Apr 2002: Urszula Jambor; Initial Code
!  24 Nov 2003: Sujay Kumar; Included in LIS
! 
! !INTERFACE:
subroutine readbergmask()
! !USES: 
  use lisdrv_module, only :lis    ! LIS non-model-specific 1-D variables
  use bergdomain_module, only : bergdrv
!EOP
  implicit none
  
  integer :: i, j, c
  integer :: istat
  integer :: lmask(bergdrv%ncold,bergdrv%nrold)

  !----------------------------------------------------------------
  !Open land-sea mask file
  !----------------------------------------------------------------
  print*, 'opening mask file ',bergdrv%emaskfile
  open(88,file=bergdrv%emaskfile, form='formatted', iostat=istat)
  if (istat/=0) then
     print *, 'Problem opening file: ', trim(bergdrv%emaskfile)
  else
     do i = 1,bergdrv%nrold
        read(88,'(720(2x,i1))') (lmask(j,i),j=1,bergdrv%ncold)
     end do
  end if
  close(88)

  !----------------------------------------------------------------
  !Change grid from ECMWF convention to GLDAS convention
  !----------------------------------------------------------------

  call berggrid_2_gldasgrid(bergdrv%ncold, bergdrv%nrold, lmask)

  !----------------------------------------------------------------
  !Transfer 2D grid to 1D vector
  !----------------------------------------------------------------

  c = 0
  do i=1, bergdrv%nrold
     do j=1, bergdrv%ncold
	c = c + 1
        bergdrv%remask1d(c) = lmask(j,i)
     enddo
  enddo

end subroutine readbergmask

