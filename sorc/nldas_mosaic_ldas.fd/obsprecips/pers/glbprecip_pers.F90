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
! !ROUTINE: glbprecip_pers.F90
!
! !DESCRIPTION:
!  Includes reading routines for global PERSIANN precipitation product
!  Used instead of GDAS/GEOS precipitation forcing
!
! !REVISION HISTORY:
!  17 Jul 2001: Jon Gottschalck; Initial code
!  04 Feb 2002: Jon Gottschalck; Added necessary code to use global precip
!               observations with domain 3 (2x2.5)
!  30 Jul 2002: Jon Gottschalck; Added code to use Huffman and Persiann precip data
!
! !INTERFACE:
subroutine glbprecip_pers ( ld, gindex, name_pers, ferror_pers )
! !USES:
  use lis_module 
  use obsprecipforcing_module, only: obsprecip
!EOP
  implicit none

  type (lisdec) ld
  integer :: gindex(ld%d%lnc, ld%d%lnr)
  integer :: index,row(243003),col(243003)

!==== Local Variables=======================

  integer :: c,r,ferror_pers,ios,lrec,offsetx,offsety      ! Loop indicies and error flags
  integer :: i,j,xd,yd,gyd,ibad,mm,nn,xdc,ydc,offset       ! Domain specific parameters
  parameter(xd=1440,yd=480)                                ! Dimension of original PERSIANN data
  parameter(ibad=-9999.0)                                  ! Bad (missing data) value
  real :: rr(xd,yd)                                        ! Original precip array
  real :: precip(xd,yd),testout(xd,yd)                     ! Reconfigured original precip array
  real, pointer :: precip_regrid(:,:)                      ! Interpolated precip array
  character(len=80) :: fname, name_pers                    ! Filename variables

!=== End Variable Definition =======================
!BOC
  fname = name_pers
  obsprecip = -1.0
!----------------------------------------------------------------------
! Determine offset in number of rows from 60 S 
! since PERSIANN starts at 50 S
!----------------------------------------------------------------------
  open(unit=10,file=fname, status='old',access='direct', &
      form='unformatted',recl=xd*yd*4,iostat=ios)

  if (ios .eq. 0) then
    read (10,rec=1) precip

    do i = 1,yd
       do j = 1,xd
          index = gindex(j,i)
          col(index) = j
          row(index) = i
          if (index .ne. -1) then
             obsprecip(index) = precip(j,i)
          endif
       enddo
    enddo
    ferror_pers = 1
    close(10)
    print*, "Obtained PERSIANN precipitation data ", fname
 else
    print*, "Missing PERSIANN precipitation data ", fname
    ferror_pers = 0
 endif
!EOC
end subroutine glbprecip_pers
