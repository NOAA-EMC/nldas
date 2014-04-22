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
! !ROUTINE: read_gswpsai
!
! !DESCRIPTION:
!  This program reads in GSWP-2 SAI data
!
! !REVISION HISTORY:
!  07 Jul 2005: James Geiger, Initial Specification
! 
! !INTERFACE: 
subroutine read_gswpsai(sai1, sai2, wt1, wt2)
! !USES:
   use lisdrv_module, only : lis
   use time_manager
   use clm_varder
!EOP
   implicit none

!=== Arguments ===========================================================
   real, dimension(lis%d%nch), intent(out) :: sai1, sai2
   real, intent(out)                        :: wt1, wt2

   integer :: yr1, mo1, yr2, mo2 ! Temporary Time variables
   real*8  :: time1, time2       ! Temporary Time variables
   integer :: doy1, doy2         ! Temporary Time variables
   real    :: gmt1, gmt2         ! Interpolation weights
   integer :: zeroi              ! Integer Number Holders

   real, dimension(16) :: sai_table
   integer :: t
!=== End Local variable list
!BOC
!------------------------------------------------------------------------
! Determine current time to find correct LAI files
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! Initialize LAI flag varaible
!------------------------------------------------------------------------
    lis%p%saiflag = 0
   
   mo1 = lis%t%mo
   yr1 = lis%t%yr 
   mo2 = lis%t%mo+1
   if ( mo2 == 13 ) then 
      mo2 = 1
      yr2 = lis%t%yr + 1
   endif

   call date2time(time1,doy1,gmt1,yr1,mo1,lis%t%da,zeroi,zeroi,zeroi)
   call date2time(time2,doy2,gmt2,yr2,mo2,1,zeroi,zeroi,zeroi)

!------------------------------------------------------------------------   
! Check to see if need new LAI data
!------------------------------------------------------------------------
   if (time2 > lis%p%saitime) then 
      lis%p%saiflag = 1
   else
      lis%p%saiflag = 0
   endif
   
!------------------------------------------------------------------------
! Determine weights between months
!------------------------------------------------------------------------
!   wt1 = (time2-lis%t%time)/(time2-time1)
!   wt2 = (lis%t%time-time1)/(time2-time1)
   wt1 = 1.0
   wt2 = 0.0

!------------------------------------------------------------------------   
! Get new SAI data if required
!------------------------------------------------------------------------
   if ( lis%p%saiflag == 1 ) then
      lis%p%saitime = time2

      !Overwrite with the table-based sai for IGBP
      sai_table(1) = 2.0
      sai_table(2) = 2.0
      sai_table(3) = 2.0
      sai_table(4) = 2.0
      sai_table(5) = 2.0
      sai_table(6) = 2.0
      sai_table(7) = 2.0
      sai_table(8) = 2.0
      sai_table(9) = 2.4352
      sai_table(10) = 4.0
      sai_table(11) = 0.0
      sai_table(12) = 0.95765
      sai_table(13) = 0.0 !2.0
      sai_table(14) = 0.95765
      sai_table(15) = 0.0
      sai_table(16) = 0.0

      do t=1,lis%d%nch
         sai1(t) = sai_table(clm(t)%itypveg)
      enddo

      sai2 = 0.0

   endif

end subroutine read_gswpsai

