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
! moslairead.f90:
!
! DESCRIPTION:
!  This program reads in AVHRR LAI and DSAI data for MOSAIC
!
! REVISION HISTORY:
!  10 Dec 2001: Jon Gottschalck; Initial code
!  20 Feb 2002: Jon Gottschalck; Modified to use for 1/4 and 2x2.5 using 1/8 degree monthly data
!  01 Oct 2002: Jon Gottschalck; Modified to allow MODIS LAI
!=========================================================================
 
subroutine moslairead (yr1,mo1,yr2,mo2,time1,time2, &
     wt1,wt2)

  use lisdrv_module, only : lis, tile
  use mos_varder
  use precision
  implicit none

  integer            :: yr1,mo1,yr2,mo2                          ! temporary time variables
  real*8             :: time1,time2                              ! temporary time variables
  real(r8) :: lai1(lis%d%nch),lai2(lis%d%nch)
  real(r8) :: sai1(lis%d%nch),sai2(lis%d%nch)
  real :: green1(lis%d%nch),green2(lis%d%nch)
  real :: wt1, wt2
  integer :: t

  call readlai(lis%p%lai, lai1, lai2, wt1, wt2)
  call readsai(lis%p%lai, sai1, sai2, wt1, wt2)
!	wt1=0.5
!	wt2=0.5
!	lai1=1.0
!	lai2=1.0
!	sai1=0.5
!	sai2=0.5
  if(lis%p%laiflag .eq. 1) then 
     mos%lai1 = lai1
     mos%lai2 = lai2
  endif

  if(lis%p%saiflag .eq. 1) then 
     mos%sai1 = sai1
     mos%sai2 = sai2
  endif

  do t=1,lis%d%nch
     if (tile(t)%vegt .eq. 12) then
        mos(t)%lai1 = 0.001
        mos(t)%sai1  = 0.001
        mos(t)%lai2 = 0.001
        mos(t)%sai2  = 0.001
        green1(t)   = 0.001
        green2(t)   = 0.001
     else
        if ((mos(t)%lai1 + mos(t)%sai1) .ne. 0.0) then 
           green1(t) = mos(t)%lai1 / (mos(t)%lai1 + mos(t)%sai1)
        else
           green1(t) = 0.001
        endif
        if ((mos(t)%lai2 + mos(t)%sai2) .ne. 0.0) then
           green2(t) = mos(t)%lai2 / (mos(t)%lai2 + mos(t)%sai2)
        else
           green2(t) = 0.001
        endif
     endif
  end do
 
  do t=1,lis%d%nch
     mos(t)%lai   = wt1 * (mos(t)%lai1 + mos(t)%sai1) &
          + wt2 * (mos(t)%lai2 + mos(t)%sai2)
     mos(t)%dsai  = wt1 * mos(t)%sai1 + wt2 * mos(t)%sai2
     mos(t)%green = wt1 * green1(t)  + wt2 * green2(t)
     
  enddo
  
!      allocate(domlai(lis%d%lnc,lis%d%lnr))
!      domlai = -9999.0
!      if(domain.eq.8) then 
!         do t=1,lis%d%nch
!            print*, 'here ',t,tile(t)%row, tile(t)%col
!            domlai(tile(t)%col,tile(t)%row) = mos(t)%lai
!         enddo
!      endif
!      open(32,file="domlai.bin",form='unformatted')
!      write(32) domlai
!      close(32)
!      deallocate(domlai)

end subroutine moslairead
