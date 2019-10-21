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
#include "misc.h"
!BOP
!
! !ROUTINE: read_gswplai
!
! !DESCRIPTION:
!  This program reads in GSWP-2 LAI data
!
! !REVISION HISTORY:
!  07 Jul 2005: James Geiger, Initial Specification
! 
! !INTERFACE: 
subroutine read_gswplai(lai1, lai2, wt1, wt2)
! !USES:
   use lisdrv_module, only : lis
   use time_manager
   use gswp_module,   only : getgswp_monindex
#if ( defined USE_NETCDF )
   use netcdf
#endif
!EOP
   implicit none

!=== Arguments ===========================================================
   real, dimension(lis%d%nch), intent(out) :: lai1, lai2
   real, intent(out)                        :: wt1, wt2

   integer :: yr1, mo1, yr2, mo2 ! Temporary Time variables
   real*8  :: time1, time2       ! Temporary Time variables
   integer :: doy1, doy2         ! Temporary Time variables
   real    :: gmt1, gmt2         ! Interpolation weights
   integer :: zeroi, numi        ! Integer Number Holders

   integer           :: status, laiid, ncid, index
!=== End Local variable list
!BOC
!------------------------------------------------------------------------
! Determine current time to find correct LAI files
!------------------------------------------------------------------------
   if (lis%t%tscount .eq. 0) then
      lis%t%yr = lis%t%syr
      lis%t%mo = lis%t%smo
      lis%t%da = lis%t%sda
      lis%t%mn = lis%t%smn
      lis%t%ss = lis%t%sss
   endif
 
   call date2time(lis%t%time,lis%t%doy,lis%t%gmt,lis%t%yr, &
        lis%t%mo,lis%t%da,lis%t%hr,lis%t%mn,lis%t%ss)
!------------------------------------------------------------------------
! Initialize LAI flag varaible
!------------------------------------------------------------------------
   lis%p%laiflag = 0

!------------------------------------------------------------------------   
! Determine Monthly data Times (Assume Monthly 
! value valid at DA=16 HR=00Z)
!------------------------------------------------------------------------
   zeroi=0
   numi=16
   if (lis%t%da .lt. 16) then
      mo1 = lis%t%mo-1
      yr1 = lis%t%yr
      if (mo1 .eq. 0) then
         mo1 = 12
         yr1 = lis%t%yr - 1
      endif
      mo2 = lis%t%mo
      yr2 = lis%t%yr
   else
      mo1 = lis%t%mo
      yr1 = lis%t%yr
      mo2 = lis%t%mo+1
      yr2 = lis%t%yr
      if (mo2 .eq. 13) then
         mo2 = 1
         yr2 = lis%t%yr + 1
      endif
   endif
   
   call date2time(time1,doy1,gmt1,yr1,mo1,numi,zeroi,zeroi,zeroi)
   call date2time(time2,doy2,gmt2,yr2,mo2,numi,zeroi,zeroi,zeroi) 
   
!------------------------------------------------------------------------   
! Check to see if need new LAI data
!------------------------------------------------------------------------
   if (time2 > lis%p%laitime) then 
      lis%p%laiflag = 1
   else
      lis%p%laiflag = 0
   endif
   
!------------------------------------------------------------------------
! Determine weights between months
!------------------------------------------------------------------------
   wt1 = (time2-lis%t%time)/(time2-time1)
   wt2 = (lis%t%time-time1)/(time2-time1)
!   wt1 = 1.0
!   wt2 = 0.0

!------------------------------------------------------------------------   
! Get new LAI data if required
!------------------------------------------------------------------------
   if ( lis%p%laiflag == 1 ) then
      lis%p%laitime = time2
#if ( defined USE_NETCDF )
      call getgswp_monindex(lis%t%yr,lis%t%mo, index)

      print*, 'MSG: read_gswplai -- Reading LAI file: ',trim(lis%p%gswplai), &
              ' for ',index

      status = nf90_open(path=trim(lis%p%gswplai), mode=nf90_nowrite, ncid=ncid)
      status = nf90_inq_varid(ncid, "LAI", laiid)

      status = nf90_get_var(ncid, laiid, lai1,       &
                            start=(/1,index/),       &
                            count=(/lis%d%glbnch,1/))

      status = nf90_get_var(ncid, laiid, lai2,       &
                            start=(/1,index+1/),     &
                            count=(/lis%d%glbnch,1/))

      status = nf90_close(ncid)
      print*, 'MSG: read_gswplai -- Read LAI data ',status
#else
      call lis_log_msg("ERR: read_gswplai -- Don't know how to read LAI")
      call endrun
#endif
   endif

end subroutine read_gswplai

