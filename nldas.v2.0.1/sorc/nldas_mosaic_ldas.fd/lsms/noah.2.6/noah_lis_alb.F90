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
! !ROUTINE: noah_lis_alb
!
! !DESCRIPTION:
!  This subroutine takes quarterly surface albedo (snow-free) data and  
!  day to interpolate and determine the actual value of the albedo 
!  for that date.  This actual value is then returned to the main
!  program.  The assumption is that the data point is valid for the 
!  dates of January 31, April 30, July 31, and October 31.
!
! !REVISION HISTORY:
!
!  28 Apr 2002: K. Arsenault; Added NOAH LSM to LDAS, initial code
!
! !INTERFACE:
subroutine noah_lis_alb
! !USES:     
  use time_manager
  use noah_varder      ! NOAH tile variables
  use time_manager
  use lisdrv_module, only : grid,tile,lis
  use lis_openfileMod
  use lis_indices_module
#if ( defined OPENDAP )
  use opendap_module
#endif
  implicit none
!EOP
!=== Local Variables =====================================================
  real, allocatable :: albout(:,:)
  integer :: cindex, rindex
  integer :: line,line1,line2,glnc,glnr
  integer :: ios1
  integer :: i,j,c,r              ! loop counters
  integer :: janda,janmo          ! january 31 
  integer :: aprda,aprmo          ! april 30
  integer :: julda,julmo          ! july 31
  integer :: octda,octmo          ! october 31 
  integer :: yr                   ! year of run  
  integer :: doy1                 ! temporary time variables
  integer :: zeroi                ! integer number holders
  integer :: albflag              ! flag to update albedo values
  real*8 :: time                  ! current model time variable
  real*8 :: jan31,apr30,jul31,oct31 ! dates of quarterly albedo files
  real*8 :: qdif                  ! difference between q1 and q2 times
  real*8 :: timdif                ! difference between time and q1 time 
  real :: gmt1,gmt2               ! gmt values 
  real :: value1(lis_nc_data,lis_nr_data) ! Temporary value holder for QQ1
  real :: value2(lis_nc_data,lis_nr_data) ! Temporary value holder for QQ2

  real :: valdif(lis%d%nch)       ! Difference of QQ2 and QQ1 albedo
  character*2 :: qq1,qq2       ! Filename places for quarter values 
  character*50 :: file1, file2
!=== End Variable Definition =============================================
!BOC
  zeroi=0
  noahdrv%noah_aflag = 0
!-------------------------------------------------------------------------  
! Determine Dates of the quarters in terms of Year (e.g., 1999.3) 
!-------------------------------------------------------------------------  
  time=lis%t%time
  yr=lis%t%yr
!-------------------------------------------------------------------------  
!  January 31
!-------------------------------------------------------------------------  
  janda=31
  janmo=01
  call date2time(jan31,doy1,gmt1,yr,janmo,&
       janda,zeroi,zeroi,zeroi)
!-------------------------------------------------------------------------  
!  April 30
!-------------------------------------------------------------------------  
  aprda=30
  aprmo=04
  call date2time(apr30,doy1,gmt1,yr,aprmo,&
       aprda,zeroi,zeroi,zeroi)
!-------------------------------------------------------------------------  
!  July 31
!-------------------------------------------------------------------------  
  julda=31
  julmo=07
  call date2time(jul31,doy1,gmt1,yr,julmo,&
       julda,zeroi,zeroi,zeroi)
!-------------------------------------------------------------------------  
!  October 31 
!-------------------------------------------------------------------------  
  octda=31
  octmo=10
  call date2time(oct31,doy1,gmt1,yr,octmo,&
       octda,zeroi,zeroi,zeroi)
!-------------------------------------------------------------------------  
! Determine which two quarterly albedo files book-end model time.
!-------------------------------------------------------------------------  

  if ( time.ge.jan31 .and. time.le.apr30 ) then
     qq1="01"
     qq2="02"
     qdif = apr30-jan31
     timdif = time-jan31
     albflag = 1
  elseif ( time.ge.apr30 .and. time.le.jul31 ) then
     qq1="02"
     qq2="03"
     qdif = jul31-apr30
     timdif = time-apr30
     albflag = 2
  elseif ( time.ge.jul31 .and. time.le.oct31 ) then
     qq1="03"
     qq2="04"
     qdif = oct31-jul31
     timdif = time-jul31
     albflag = 3
  elseif ( time.ge.oct31 ) then
     qq1="04"
     qq2="01"
     qdif = (jan31+1.0)-oct31
     timdif = time-oct31
     albflag = 4
  elseif ( time.lt.jan31) then
     qq1="04"
     qq2="01"
     oct31=oct31-1.0
     qdif = jan31-oct31
     timdif = time-oct31
     albflag = 5
  endif
  
  if(noahdrv%noah_albtime .ne. albflag) then 
     noahdrv%noah_albtime = albflag
     noahdrv%noah_aflag = 1
!-------------------------------------------------------------------------  
!  Open the needed two quarterly snow-free albedo files   
!-------------------------------------------------------------------------  
     file1 = trim(noahdrv%noah_albfile)//'alb_'//QQ1//'.1gd4r'
     call lis_open_file(10, file=file1, status='old', form='unformatted', &
                        access='direct', recl=4, script='getalbedo.pl',   &
                        time_offset=qq1)
     file2 = trim(noahdrv%noah_albfile)//'alb_'//QQ2//'.1gd4r'
     call lis_open_file(11, file=file2, status='old', form='unformatted', &
                        access='direct', recl=4, script='getalbedo.pl',   &
                        time_offset=qq2)
     call lis_read_file(10,value1)
     call lis_read_file(11,value2)
     close(10)
     close(11)

!-------------------------------------------------------------------------  
! Assign quarterly albedo fractions to each tile. 
!-------------------------------------------------------------------------  
     do i=1,lis%d%nch
        if((value1(tile(i)%col, tile(i)%row-lis_tnroffset).ne.-9999.000) &
             .and. (value2(tile(i)%col,tile(i)%row-lis_tnroffset)&
             .ne.-9999.000)) then 
           noah(i)%albsf1= value1(tile(i)%col, tile(i)%row-lis_tnroffset)
           noah(i)%albsf2= value2(tile(i)%col, tile(i)%row-lis_tnroffset)
        endif
     enddo
  endif    ! End albflag selection

!-------------------------------------------------------------------------  
! Assign albedo fractions to each tile and interpolate daily.
!-------------------------------------------------------------------------  
  if (noahdrv%noah_albdchk .ne. lis%t%da) then 
     noahdrv%noah_aflag = 1     
     do i=1,lis%d%nch
        if (noah(i)%albsf1 .ne. -9999.000) then 
           valdif(i) = noah(i)%albsf2 - noah(i)%albsf1
           noah(i)%albsf = (timdif*valdif(i)/qdif)+noah(i)%albsf1
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
              rindex = tile(i)%row-nint((lis%d%gridDesc(4)-lis%d%gridDesc(44))&
                       /lis%d%gridDesc(9))
              cindex = tile(i)%col-nint((lis%d%gridDesc(5)-lis%d%gridDesc(45))&
                       /lis%d%gridDesc(10))
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
end subroutine noah_lis_alb
