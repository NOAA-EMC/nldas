!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine fill_land(ldas, grid, c, r, array)
!
! DESCRIPTION:
! Averages above-zero grid values surrounding grid point (c,r).
! Called by retglbSW to modify SW fluxes of zero over the equator  
! incorrectly computed by ipolates for the 2.0 by 2.5 degree domain.
! This is a coarse fix to avoid making changes to weights used in
! the interpolation source code.
!
! INPUT:
! c,r         Grid column and row location
! array       Full data array
!
! OUTPUT:
! array(c,r)  Original zero-value overwritten with average of 
!             surrounding values.
!
! REVISION HISTORY:
! 2002.02.27  U. Jambor; original code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fill_land(ldas, grid, c, r, array)

  use ldas_module      	! LDAS non-model-specIFic 1-D variables
  use grid_module      	! LDAS non-model-specIFic grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc, ldas%nr)

  integer :: c, r, npts

  real :: array(ldas%nc,ldas%nr)
  real :: fillvalue

!!!end of declarations

  npts = 0
  fillvalue = 0.0
  if (grid(c,r)%fimask==1) then !land point
    if (array(c-1,r-1)>0.0) then
      npts=npts+1
      fillvalue=fillvalue+array(c-1,r-1)
    endif
    if (array(c,r-1)>0.0) then
      npts=npts+1
      fillvalue=fillvalue+array(c,r-1)
    endif
    if (array(c+1,r-1)>0.0) then
      npts=npts+1
      fillvalue=fillvalue+array(c+1,r-1)
    endif
    if (array(c-1,r)>0.0) then
      npts=npts+1
      fillvalue=fillvalue+array(c-1,r)
    endif
    if (array(c+1,r)>0.0) then
      npts=npts+1
      fillvalue=fillvalue+array(c+1,r)
    endif
    if (array(c-1,r+1)>0.0) then
      npts=npts+1
      fillvalue=fillvalue+array(c-1,r+1)
    endif
    if (array(c,r+1)>0.0) then
      npts=npts+1
      fillvalue=fillvalue+array(c,r+1)
    endif
    if (array(c+1,r+1)>0.0) then
      npts=npts+1
      fillvalue=fillvalue+array(c+1,r+1)
    endif

!!!Compute new grid value based on surrounding points
    if (npts > 0) then
      fillvalue = fillvalue / npts
      array(c,r) = fillvalue
!      print*, 'npts, fillvalue', npts, fillvalue 
    end if

  end if !land point
 
end subroutine fill_land
