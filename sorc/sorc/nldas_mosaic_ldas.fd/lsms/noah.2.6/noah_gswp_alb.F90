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
! !ROUTINE: noah_gswp_alb
!
! !DESCRIPTION:
! This routine reads in GSWP-2 albedo data.
!
! !REVISION HISTORY:
! 14 Jun 2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah_gswp_alb
! !USES:     
#if ( defined USE_NETCDF )
  use netcdf
  use noah_varder      ! NOAH tile variables
  use time_manager
  use lisdrv_module, only : grid,tile,lis,gindex
  use gswp_module,   only : getgswp_monindex
  use gswpdomain_module, only : gswpdrv
  implicit none
!EOP
!=== Local Variables =====================================================
  integer :: index
  real, allocatable :: albout(:,:)
  integer :: cindex, rindex
  real*8  :: time1,time2           ! Temporary Time variables
  integer :: yr1,mo1,yr2,mo2       ! Temporary Time variables
  INTEGER :: I,J,c,r               ! Loop counters
  integer :: doy1,doy2             ! Temporary Time variables
  integer :: zeroi,numi            ! Integer Number Holders
  integer :: alb_flag              ! Flag to update albedo values
  real :: gmt1,gmt2                ! GMT Values 
  real :: wt
  real,allocatable :: value1(:,:)  ! Temporary value holder for QQ1
  real,allocatable :: value2(:,:)  ! Temporary value holder for QQ2
  integer :: tnroffset = 0
  real, allocatable, dimension(:) :: alb1, alb2
  integer :: ncid, status
  integer :: albid
!=== End Variable Definition =============================================
!BOC

  noahdrv%noah_aflag = 0

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

  wt = ( lis%t%time - time1 ) / ( time2 - time1 )

  if(time2 .gt. noahdrv%noah_gfractime) then 
     alb_flag = 1
  else 
     alb_flag = 0
  endif

  if ( alb_flag == 1 ) then 
     noahdrv%noah_gfractime = time2
!-------------------------------------------------------------------------  
!  Open the needed two quarterly snow-free albedo files   
!-------------------------------------------------------------------------  
     allocate(value1(lis%d%lnc,lis%d%lnr))
     allocate(value2(lis%d%lnc,lis%d%lnr))

     allocate(alb1(lis%d%glbnch))
     allocate(alb2(lis%d%glbnch))

     call getgswp_monindex(yr1, mo1, index)

     call lis_log_msg('MSG: noah_gswp_alb -- Reading ALBEDO file: ' &
                      //trim(gswpdrv%albedo))

     status = nf90_open(path=trim(gswpdrv%albedo), mode=nf90_nowrite, &
                        ncid=ncid)
     status = nf90_inq_varid(ncid, "Albedo", albid)
     status = nf90_get_var(ncid, albid, alb1,       &
                           start=(/1,index/),       &
                           count=(/lis%d%glbnch,1/))
     status = nf90_get_var(ncid, albid, alb2,       &
                           start=(/1,index+1/),     &
                           count=(/lis%d%glbnch,1/))
     status = nf90_close(ncid)

     do c=1,lis%d%lnc
        do r=1,lis%d%lnr
           rindex = 150-r+1
           cindex = c
           if(gindex(cindex,rindex).ne.-1) then 
              value1(cindex,rindex) = alb1(gindex(cindex,rindex))
              value2(cindex,rindex) = alb2(gindex(cindex,rindex))
           endif
        enddo
     enddo

     do i=1,lis%d%nch
        if((value1(tile(i)%col, tile(i)%row-tnroffset).ne.-9999.000)) then 
           noah(i)%albsf1 = value1(tile(i)%col, tile(i)%row-tnroffset)
           noah(i)%albsf2 = value2(tile(i)%col, tile(i)%row-tnroffset)
        endif
     enddo

     deallocate(value1)
     deallocate(value2)
     deallocate(alb1)
     deallocate(alb2)

   endif

!-------------------------------------------------------------------------  
! Assign albedo fractions to each tile and interpolate daily.
!-------------------------------------------------------------------------  
  if (noahdrv%noah_albdchk .ne. lis%t%da) then 
     noahdrv%noah_aflag = 1     
     do i=1,lis%d%nch
        if (noah(i)%albsf1 .ne. -9999.000) then 
           noah(i)%albsf = wt*noah(i)%albsf2 + (1.0-wt)*noah(i)%albsf1
        endif
     end do

     noahdrv%noah_albdchk=lis%t%da

     if(lis%o%wparam.eq.1) then 
        allocate(albout(lis%d%lnc,lis%d%lnr))
        albout = -9999.0
#if ( defined OPENDAP )
        do i=1,lis%d%nch      
           albout(tile(i)%col,tile(i)%row) = noah(i)%albsf
        enddo
#else
        do i=1,lis%d%nch      
           if(grid(i)%lat.ge.lis%d%gridDesc(4).and. & 
                grid(i)%lat.le.lis%d%gridDesc(7).and. & 
                grid(i)%lon.ge.lis%d%gridDesc(5).and. & 
                grid(i)%lon.le.lis%d%gridDesc(8)) then
              rindex = tile(i)%row - (lis%d%gridDesc(4)-lis%d%gridDesc(44)) &
                       /lis%d%gridDesc(9)
              cindex = tile(i)%col - (lis%d%gridDesc(5)-lis%d%gridDesc(45)) &
                       /lis%d%gridDesc(10)
              albout(cindex,rindex) = noah(i)%albsf
           endif
        enddo
#endif

        open(32,file="albout.bin",form='unformatted')
        write(32) albout
        close(32)  
        deallocate(albout)
     end if
  endif ! End daily interpolation
  return
!EOC
#endif
end subroutine noah_gswp_alb
