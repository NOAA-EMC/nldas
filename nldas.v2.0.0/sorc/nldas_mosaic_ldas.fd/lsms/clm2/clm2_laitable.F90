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
module clm2_laitable
  real :: lai_table(13,12)

contains
  subroutine readlaitable
    use lisdrv_module, only : lis, tile
    use clm_varder, only :clm
    integer :: i,j,t

    open(92,file='./BCS/clm_parms/veg_lib.txt',status='old')
    do t=1,13
       read(92,*) (lai_table(t,j), j=1,12)
    enddo
    close(92)

   do t=1,lis%d%nch
      clm(t)%tlai=lai_table(tile(t)%vegt,lis%t%mo)
   enddo
  end subroutine readlaitable
end module clm2_laitable
