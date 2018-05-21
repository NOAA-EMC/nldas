!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! recmwfmask.F90: 
!
! DESCRIPTION:
! Reads in land-sea mask for 0.5 degree Reanalysis ECMWF data set,
! changes grid from ECMWF convention to GLDAS convention by calling
! Recmwfgrid_2_gldasgrid, and transforms grid array into 1D vector 
! for later use in ret_reanlecmwf.F90.  
!
! REVISION HISTORY:
!  16 Apr 2002: Urszula Jambor; Initial Code
!=========================================================================

subroutine recmwfmask(ldas)

  use ldas_module      ! LDAS non-model-specific 1-D variables

  implicit none
  
  type (ldasdec) LDAS 

!=== Local Variables ===================================================

  integer :: i, j, c
  integer :: istat
  integer :: lmask(ldas%ncold,ldas%nrold)

!=== End Variable Definition ===========================================

  !----------------------------------------------------------------
  !Open land-sea mask file
  !----------------------------------------------------------------

  open(88,file=ldas%emaskfile, form='formatted', iostat=istat)
  if (istat/=0) then
     print *, 'Problem opening file: ', trim(ldas%emaskfile)
  else
     do i = 1,ldas%nrold
        read(88,'(720(2x,i1))') (lmask(j,i),j=1,ldas%ncold)
     end do
  end if
  close(88)

  !----------------------------------------------------------------
  !Change grid from ECMWF convention to GLDAS convention
  !----------------------------------------------------------------

  call Recmwfgrid_2_gldasgrid(ldas%ncold, ldas%nrold, lmask)

  open(55,file='../testmask.bin',form='unformatted')
  write(55) lmask*500.0
  close(55)

  !----------------------------------------------------------------
  !Transfer 2D grid to 1D vector
  !----------------------------------------------------------------

  c = 0
  do i=1, ldas%nrold
     do j=1, ldas%ncold
	c = c + 1
        ldas%remask1d(c) = lmask(j,i)
     enddo
  enddo

end subroutine recmwfmask

