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
! !ROUTINE: glbprecip_huff.F90
!
! !DESCRIPTION:
!  Includes reading routines for global HUFFMAN precipitation product
!  Used instead of GDAS/GEOS precipitation forcing
!
! !REVISION HISTORY: 
!  17 Jul 2001: Jon Gottschalck; Initial code
!  04 Feb 2002: Jon Gottschalck; Added necessary code to use global precip
!               observations with domain 3 (2x2.5)
!  30 Jul 2002: Jon Gottschalck; Added code to use Huffman and Persiann precip data
!
! !INTERFACE:
subroutine glbprecip_huff (name_huff, ferror_huff )
! !USES:
 use lisdrv_module, only : gindex 
 use obsprecipforcing_module, only: obsprecip
!EOP
  implicit none

  integer :: index

!==== Local Variables=======================
       
 integer :: c,r,ferror_huff,ios,lrec,offsetx,offsety  ! Loop indicies and error flags
 integer :: i,j,xd,yd,gyd,ibad,mm,nn,xdc,ydc          ! Domain specific parameters
 parameter(xd=1440,yd=480)                            ! Dimension of original HUFFMAN data
 parameter(ibad=-31999)                               ! Bad (missing data) value
 integer*2 :: rr(xd,yd)                               ! Original 2 byte integer HUFFMAN data
 real :: precip(xd,yd)                                ! Original real precipitation array
 real, pointer :: precip_regrid(:,:)                  ! Interpolated precipitation array
 character(len=80) :: fname, name_huff                ! Filename variables
 character(len=2880) ::  head                         ! Header variable
       
!=== End Variable Definition =======================
!BOC       
 fname = name_huff
!------------------------------------------------------------------------
! Fill necessary arrays to assure not using old HUFFMAN data
!------------------------------------------------------------------------
 precip = -1.0
 obsprecip = -1.0
!------------------------------------------------------------------------
! Find HUFFMAN precip data, read it in and assign to forcing precip array.
! Must reverse grid in latitude dimension to be consistent with LDAS grid
!------------------------------------------------------------------------
  open(unit=10,file=fname, status='old', &
   &          access='direct',recl=xd*yd*4, &
   &          form='unformatted',iostat=ios)

  if (ios .eq. 0) then
     read (10,rec=1) head (1:2880),&
          ( ( rr (i, j), i = 1, xd ), j = 1, yd)
   do i = 1,yd
    do j = 1,xd
     if (j .lt. 721) then
      precip(j,i) = float(rr(j+720,yd+1-i)) / 100.0
      if ( rr(j+720,yd+1-i) .eq. ibad ) precip(j,i) = -1.0
      if ( rr(j+720,yd+1-i) .lt. 0.0 )  precip(j,i) = -1.0
     else
      precip(j,i) = float(rr(j-720,yd+1-i)) / 100.0
      if ( rr(j-720,yd+1-i) .eq. ibad ) precip(j,i) = -1.0
      if ( rr(j-720,yd+1-i) .lt. 0.0 )  precip(j,i) = -1.0
     endif
    enddo
  enddo
!------------------------------------------------------------------------
! Interpolating to desired domain and resolution
! Global precip datasets not used currently to force NLDAS
!------------------------------------------------------------------------
  do i = 1,yd
   do j = 1,xd
       index = gindex(j,i)
       if (index .ne. -1) then
         obsprecip(index) = precip(j,i)
       endif
   enddo
  enddo

  ferror_huff = 1
 close(10)
  print*, "Obtained HUFFMAN precipitation data ", fname
 else
   print*, "Missing HUFFMAN precipitation data ", fname
   ferror_huff = 0
 endif
!EOC
end subroutine glbprecip_huff
