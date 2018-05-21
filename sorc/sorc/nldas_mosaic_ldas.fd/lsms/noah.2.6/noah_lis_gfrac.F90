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
! !ROUTINE: noah_lis_gfrac
!
! !DESCRIPTION:
!  This subroutine takes vegetation greenness fraction data and the date to 
!  interpolate and determine the actual value of the greenness fraction 
!  for that date.  This actual value is then returned to the main
!  program.  The assumption is that the data point is valid for the 16th
!  of the given month, at 00Z.
!
! !REVISION HISTORY:
!
!  28 Apr 2002: K. Arsenault; Added NOAH LSM to LDAS, initial code
!
! !INTERFACE:
subroutine noah_lis_gfrac
! !USES:
  use noah_varder      
  use time_manager
  use lisdrv_module, only : grid,tile,lis
  use lis_openfileMod
  use lis_indices_module
!EOP
  implicit none
!=== Local Variables =====================================================
  real,allocatable :: gfracout(:,:)
  integer :: cindex, rindex
  integer :: line,line1,line2,glnc,glnr
  integer :: ios1
  integer :: i,t,c,r              ! loop counters
  real*8  :: time1,time2          ! temporary time variables
  integer :: yr1,mo1,yr2,mo2      ! temporary time variables
  integer :: doy1,doy2            ! temporary time variables
  integer :: zeroi,numi           ! integer number holders
  real :: wt1,wt2,gmt1,gmt2       ! interpolation weights
  real :: value1(lis_nc_data,lis_nr_data) ! temporary value holder for mo1
  real :: value2(lis_nc_data,lis_nr_data) ! temporary value holder for mo2

  character*2 :: mm1,mm2        ! Filename places for integer. MO1, MO2
  integer :: gfrac_flag 
  character*50 :: file1, file2
!=== End Variable Definition =============================================
!BOC
  noahdrv%noah_gflag = 0
  zeroi=0
  numi=16
!------------------------------------------------------------------------
! Determine Monthly data Times (Assume Monthly value valid at DA=16)
!------------------------------------------------------------------------
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
  
  call date2time(time1,doy1,gmt1,yr1,mo1,&
       numi,zeroi,zeroi,zeroi)
  call date2time(time2,doy2,gmt2,yr2,mo2,&
       numi,zeroi,zeroi,zeroi)
!------------------------------------------------------------------------  
!  Weights to be used to interpolate greenness fraction values.
!------------------------------------------------------------------------
  wt1= (time2-lis%t%time)/(time2-time1)
  wt2= (lis%t%time-time1)/(time2-time1)
!------------------------------------------------------------------------
!  Determine if GFRAC files need to be updated
!------------------------------------------------------------------------
  if(time2 .gt. noahdrv%noah_gfractime) then 
     gfrac_flag = 1
  else 
     gfrac_flag = 0
  endif
  
  if(gfrac_flag .eq. 1) then 
     noahdrv%noah_gfractime = time2
     noahdrv%noah_gflag = 1
!------------------------------------------------------------------------
! Open greenness fraction dataset of months corresponding to   
! time1 and time2 for selected LDAS domain and read data.
!------------------------------------------------------------------------
     write(mm1,3) mo1
     write(mm2,3) mo2
3    format(i2.2)

     file1 = trim(noahdrv%noah_mgfile)//'gfrac_'//mm1//'.1gd4r'
     call lis_open_file(10, file=file1, status='old', form='unformatted', &
                        access='direct', recl=4, script='getgfrac.pl',    &
                        time_offset=mm1)
     file2 = trim(noahdrv%noah_mgfile)//'gfrac_'//mm2//'.1gd4r'
     call lis_open_file(11, file=file2, status='old', form='unformatted', &
                        access='direct', recl=4, script='getgfrac.pl',    &
                        time_offset=mm2)
     
     call lis_read_file(10,value1)
     call lis_read_file(11,value2)
     close(10)
     close(11)

!------------------------------------------------------------------------     
! Assign MONTHLY vegetation greenness fractions to each tile. 
!------------------------------------------------------------------------
     do i=1,lis%d%nch
        if((value1(tile(i)%col, tile(i)%row-lis_tnroffset) .ne. -9999.000) & 
         .and.(value2(tile(i)%col, tile(i)%row-lis_tnroffset).ne.-9999.000)) &
         then
           noah(i)%vegmp1=value1(tile(i)%col, tile(i)%row-lis_tnroffset) 
           noah(i)%vegmp2=value2(tile(i)%col, tile(i)%row-lis_tnroffset) 
        endif
     end do
  endif

!------------------------------------------------------------------------
!  Interpolate greenness fraction values once daily
!------------------------------------------------------------------------

  if (noahdrv%noah_gfracdchk .ne. lis%t%da) then 
     noahdrv%noah_gflag = 1
     do i=1,lis%d%nch
        noah(i)%vegip = (wt1*noah(i)%vegmp1)+(wt2*noah(i)%vegmp2)
     end do
      
     noahdrv%noah_gfracdchk = lis%t%da
     print*, 'Done noah_gfrac',' (',iam,')'

     if(lis%o%wparam.eq.1) then
        allocate(gfracout(lis%d%lnc,lis%d%lnr))
        gfracout = -9999.0
#if ( defined OPENAP )
        do i=1,lis%d%nch      
           gfracout(tile(i)%col,tile(i)%row) = noah(i)%vegip*1.0
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
              gfracout(cindex,rindex) = noah(i)%vegip*1.0
           endif
        enddo
#endif
     
        open(32,file="gfracout.bin",form='unformatted')
        write(32) gfracout
        close(32)   
        deallocate(gfracout)
     endif
  end if
  return
!EOC
end subroutine noah_lis_gfrac
