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
! !ROUTINE: noah_gswp_gfrac
!
! !DESCRIPTION:
! This routine reads in GSWP-2 greenness fraction data.
!
! !REVISION HISTORY:
!
!  28 Apr 2002: K. Arsenault; Added NOAH LSM to LDAS, initial code
!
! !INTERFACE:
subroutine noah_gswp_gfrac
! !USES:
#if ( defined USE_NETCDF )
  use netcdf
  use noah_varder      
  use time_manager
  use lisdrv_module, only : grid,tile,lis,gindex
  use gswp_module,   only : getgswp_monindex
  use gswpdomain_module, only : gswpdrv
#if ( defined OPENDAP )
  use opendap_module
#endif
!EOP
  implicit none
!=== Local Variables =====================================================
  integer :: ncid, status
  integer :: grnid
  real,allocatable :: gfracout(:,:)
  real, allocatable, dimension(:) :: gfrac1, gfrac2
  integer :: index
  integer :: cindex, rindex
  integer :: i,t,c,r              ! Loop counters
  real*8  :: time1,time2          ! Temporary Time variables
  integer :: yr1,mo1,yr2,mo2      ! Temporary Time variables
  integer :: doy1,doy2            ! Temporary Time variables
  integer :: zeroi,numi           ! Integer Number Holders
  real :: wt,gmt1,gmt2            ! Interpolation weights
#if ( defined OPENDAP )
  real :: value1(parm_nc,1+nroffset:parm_nr+nroffset) ! Temporary value holder for MO1
  real :: value2(parm_nc,1+nroffset:parm_nr+nroffset) ! Temporary value holder for MO2
#else
  real,allocatable :: value1(:,:)   ! Temporary value holder for MO1
  real,allocatable :: value2(:,:)   ! Temporary value holder for MO2
#endif
  integer :: gfrac_flag 
#if ( ! defined OPENDAP )
  integer :: tnroffset = 0
#endif

!=== End Variable Definition =============================================
!BOC

  noahdrv%noah_gflag = 0

!------------------------------------------------------------------------
! Determine Monthly data Times (Assume Monthly value valid at DA=16)
!------------------------------------------------------------------------
  zeroi=0
  numi=16

  if(lis%t%da.lt.16)then
     mo1=lis%t%mo-1
     yr1=lis%t%yr 
     if(mo1.eq.0)then
        mo1=12
        yr1=lis%t%yr-1
     endif
     mo2=lis%t%mo
     yr2=lis%t%yr
  else
     mo1=lis%t%mo
     yr1=lis%t%yr
     mo2=lis%t%mo+1
     yr2=lis%t%yr
     if(mo2.eq.13)then
        mo2=1
        yr2=lis%t%yr+1
     endif
  endif

  call date2time(time1,doy1,gmt1,yr1,mo1,numi,zeroi,zeroi,zeroi)
  call date2time(time2,doy2,gmt2,yr2,mo2,numi,zeroi,zeroi,zeroi)

!------------------------------------------------------------------------  
!  Weight to be used to interpolate greenness fraction values.
!------------------------------------------------------------------------
  wt = (lis%t%time-time1)/(time2-time1)

!------------------------------------------------------------------------
!  Determine if GFRAC files need to be updated
!------------------------------------------------------------------------
  if(time2 .gt. noahdrv%noah_gfractime) then 
     gfrac_flag = 1
  else 
     gfrac_flag = 0
  endif

  if ( gfrac_flag == 1 ) then 
!     noahdrv%noah_gfractime = time2 !will be set in noah_alb
     noahdrv%noah_gflag = 1
!------------------------------------------------------------------------
! Open greenness fraction dataset of months corresponding to   
! time1 and time2 for selected LDAS domain and read data.
!------------------------------------------------------------------------

     allocate(value1(lis%d%lnc,lis%d%lnr))
     allocate(value2(lis%d%lnc,lis%d%lnr))
     allocate(gfrac1(lis%d%glbnch))
     allocate(gfrac2(lis%d%glbnch))

     call getgswp_monindex(yr1, mo1, index)

     call lis_log_msg('MSG: noah_gswp_gfrac -- Reading GFRAC file: ' &
                      //trim(gswpdrv%gfrac))

     status = nf90_open(path=trim(gswpdrv%gfrac), mode=nf90_nowrite, &
                        ncid=ncid)
     status = nf90_inq_varid(ncid, "grnFrac", grnid)
     status = nf90_get_var(ncid, grnid, gfrac1,     &
                           start=(/1,index/),       &
                           count=(/lis%d%glbnch,1/))
     status = nf90_get_var(ncid, grnid, gfrac2,     &
                           start=(/1,index+1/),     &
                           count=(/lis%d%glbnch,1/))
     status = nf90_close(ncid)

     do c=1,lis%d%lnc
        do r=1,lis%d%lnr
           rindex = 150-r+1
           cindex = c
           if(gindex(cindex,rindex).ne.-1) then 
              value1(cindex,rindex) = gfrac1(gindex(cindex,rindex))
              value2(cindex,rindex) = gfrac2(gindex(cindex,rindex))
           endif
        enddo
     enddo

!------------------------------------------------------------------------     
! Assign MONTHLY vegetation greenness fractions to each tile. 
!------------------------------------------------------------------------
     do i=1,lis%d%nch
        if((value1(tile(i)%col, tile(i)%row-tnroffset) .ne. -9999.000) & 
         .and.(value2(tile(i)%col, tile(i)%row-tnroffset).ne.-9999.000)) &
         then
           noah(i)%vegmp1=value1(tile(i)%col, tile(i)%row-tnroffset) 
           noah(i)%vegmp2=value2(tile(i)%col, tile(i)%row-tnroffset) 
        endif
     end do

     deallocate(value1)
     deallocate(value2)
     deallocate(gfrac1)
     deallocate(gfrac2)

   endif

!------------------------------------------------------------------------
!  Interpolate greenness fraction values once daily
!------------------------------------------------------------------------

   if (noahdrv%noah_gfracdchk .ne. lis%t%da) then 

      noahdrv%noah_gflag = 1

      do i=1,lis%d%nch
         noah(i)%vegip = wt*noah(i)%vegmp2 + (1.0-wt)*noah(i)%vegmp1
      end do

      noahdrv%noah_gfracdchk = lis%t%da

      if(lis%o%wparam.eq.1) then
         allocate(gfracout(lis%d%lnc,lis%d%lnr))
         gfracout = -9999.0
#if ( defined OPENDAP )
         do i=1,lis%d%nch      
            gfracout(tile(i)%col,tile(i)%row) = noah(i)%vegip*1.0
         enddo
#else
         do i=1,lis%d%nch      
            if(grid(i)%lat.ge.lis%d%gridDesc(4).and. & 
               grid(i)%lat.le.lis%d%gridDesc(7).and. & 
               grid(i)%lon.ge.lis%d%gridDesc(5).and. & 
               grid(i)%lon.le.lis%d%gridDesc(8)) then
               rindex=tile(i)%row-nint((lis%d%gridDesc(4)-lis%d%gridDesc(44)) &
                      /lis%d%gridDesc(9))
               cindex=tile(i)%col-nint((lis%d%gridDesc(5)-lis%d%gridDesc(45)) &
                      /lis%d%gridDesc(10))
               gfracout(cindex,rindex) = noah(i)%vegip*1.0
            endif
         enddo
#endif
     
         open(32,file="gfracout.bin",form='unformatted')
         write(32) gfracout
         close(32)   
         deallocate(gfracout)
      endif
   endif
   return
!EOC
#endif
end subroutine noah_gswp_gfrac
