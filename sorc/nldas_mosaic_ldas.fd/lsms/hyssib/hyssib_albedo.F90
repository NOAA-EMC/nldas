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
! !ROUTINE: hyssib_albedo.F90: 
!
! !DESCRIPTION:
!  This subroutine takes quarterly surface albedo (snow-free) data and  
!  day to interpolate and determine the actual value of the albedo 
!  for that date.  This actual value is then returned to the main
!  program.  The assumption is that the data point is valid for the 
!  dates of January 31, April 30, July 31, and October 31.
!
! !REVISION HISTORY:
! 28 Apr 2002: K. Arsenault, Added NOAH LSM to LDAS, initial code
! 21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
!
! !INTERFACE:
      subroutine hyssib_albedo
! !USES:
      use time_manager
      use hyssib_varder           ! HY-SSiB tile variables
      use time_manager
      use lisdrv_module, only : grid,tile,lis
#if ( defined OPENDAP )
      use opendap_module
#endif
!EOP
      implicit none
!=== Local Variables =====================================================
      real, allocatable :: albout(:,:)
      integer :: cindex, rindex

      INTEGER :: I,J,c,r        ! Loop counters
      INTEGER :: JANDA,JANMO    ! January 31
      INTEGER :: APRDA,APRMO    ! April 30
      INTEGER :: JULDA,JULMO    ! July 31
      INTEGER :: OCTDA,OCTMO    ! October 31
      INTEGER :: YR             ! Year of run
      INTEGER :: DOY1           ! Temporary Time variables
      INTEGER :: ZEROI          ! Integer Number Holders
      INTEGER :: ALBFLAG        ! Flag to update albedo values

      REAL*8 :: TIME            ! Current Model Time variable
      REAL*8 :: JAN31,APR30,JUL31,OCT31 ! Dates of quarterly albedo files
      REAL*8 :: QDIF            ! Difference between Q1 and Q2 times
      REAL*8 :: TIMDIF          ! Difference between TIME and Q1 time
      REAL :: GMT1,GMT2         ! GMT Values
#if ( defined OPENDAP )
      REAL :: VALUE1(parm_nc,parm_nr) ! Temporary value holder for QQ1
      REAL :: VALUE2(parm_nc,parm_nr) ! Temporary value holder for QQ2
#else
      REAL :: VALUE1(LIS%D%GNC,LIS%D%GNR) ! Temporary value holder for QQ1
      REAL :: VALUE2(LIS%D%GNC,LIS%D%GNR) ! Temporary value holder for QQ2
#endif
      REAL :: VALDIF(LIS%D%NCH) ! Difference of QQ2 and QQ1 albedo

      CHARACTER*2 :: QQ1,QQ2    ! Filename places for quarter values
#if ( ! defined OPENDAP )
      integer :: tnroffset = 0
#endif

!=== End Variable Definition =============================================
!BOC
      zeroi=0
      hyssibdrv%hyssib_aflag = 0
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
      call date2time(jan31,doy1,gmt1,yr,janmo,janda,zeroi,zeroi,zeroi)
!-------------------------------------------------------------------------
!  April 30
!-------------------------------------------------------------------------
      aprda=30
      aprmo=04
      call date2time(apr30,doy1,gmt1,yr,aprmo,aprda,zeroi,zeroi,zeroi)
!-------------------------------------------------------------------------
!  July 31
!-------------------------------------------------------------------------
      julda=31
      julmo=07
      call date2time(jul31,doy1,gmt1,yr,julmo,julda,zeroi,zeroi,zeroi)
!-------------------------------------------------------------------------
!  October 31
!-------------------------------------------------------------------------
      octda=31
      octmo=10
      call date2time(oct31,doy1,gmt1,yr,octmo,octda,zeroi,zeroi,zeroi)
!-------------------------------------------------------------------------
! Determine which two quarterly albedo files book-end model time.
!-------------------------------------------------------------------------
      if ((time.ge.jan31).and.(time.le.apr30)) then
         qq1="01"
         qq2="02"
         qdif = apr30-jan31
         timdif = time-jan31
         albflag = 1
      elseif ((time.ge.apr30).and.(time.le.jul31)) then
         qq1="02"
         qq2="03"
         qdif = jul31-apr30
         timdif = time-apr30
         albflag = 2
      elseif ((time.ge.jul31).and.(time.le.oct31)) then
         qq1="03"
         qq2="04"
         qdif = oct31-jul31
         timdif = time-jul31
         albflag = 3
      elseif (time.ge.oct31) then
         qq1="04"
         qq2="01"
         qdif = (jan31+1.0)-oct31
         timdif = time-oct31
         albflag = 4
      elseif (time.lt.jan31) then
         qq1="04"
         qq2="01"
         oct31=oct31-1.0
         qdif = jan31-oct31
         timdif = time-oct31
         albflag = 5
      endif

      if (hyssibdrv%hyssib_albtime.ne.albflag) then 
         hyssibdrv%hyssib_albtime = albflag
         hyssibdrv%hyssib_aflag = 1
!-------------------------------------------------------------------------
!  Open the needed two quarterly snow-free albedo files
!-------------------------------------------------------------------------

#if ( defined OPENDAP )
         print*,'MSG: hyssib_alb -- Retrieving ALBEDO file ',&
                trim(hyssibdrv%hyssib_albfile)//'albedo_'//QQ1// &
                '_'//trim(hyssibdrv%hyssib_galbres)//'.bfsa',' (',iam,')'
         call system("opendap_scripts/getalbedo.pl "//ciam//" "// &
                     trim(hyssibdrv%hyssib_albfile)//'albedo_'//QQ1// &
                     '_'//trim(hyssibdrv%hyssib_galbres)//'.bfsa'     &
                     //" "//cparm_slat//" "//cparm_nlat           &
                     //" "//cparm_wlon//" "//cparm_elon//" "//QQ1)
         print*,'MSG: hyssib_alb -- Retrieving ALBEDO file ', &
                trim(hyssibdrv%hyssib_albfile)//'albedo_'//QQ2// &
                '_'//trim(hyssibdrv%hyssib_galbres)//'.bfsa',' (',iam,')'
         call system("opendap_scripts/getalbedo.pl "//ciam//" "// &
                     trim(hyssibdrv%hyssib_albfile)//'albedo_'//QQ2// &
                     '_'//trim(hyssibdrv%hyssib_galbres)//'.bfsa'     &
                     //" "//cparm_slat//" "//cparm_nlat           &
                     //" "//cparm_wlon//" "//cparm_elon//" "//QQ2)
#endif
         print*,'MSG: hyssib_alb -- Retrieving ALBEDO file ', &
                trim(hyssibdrv%hyssib_albfile)//'albedo_'//QQ1// &
                '_'//trim(hyssibdrv%hyssib_galbres)//'.bfsa',' (',iam,')'
         OPEN(10,FILE=trim(HYSSIBDRV%HYSSIB_ALBFILE)//'albedo_'//QQ1// &
              '_'//trim(hyssibdrv%hyssib_galbres)//'.bfsa',&
              STATUS='OLD',FORM='UNFORMATTED')
         print*,'MSG: hyssib_alb -- Retrieving ALBEDO file ', &
                trim(hyssibdrv%hyssib_albfile)//'albedo_'//QQ2// &
                '_'//trim(hyssibdrv%hyssib_galbres)//'.bfsa',' (',iam,')'
         OPEN(11,FILE=trim(HYSSIBDRV%HYSSIB_ALBFILE)//'albedo_'//QQ2// &
              '_'//trim(hyssibdrv%hyssib_galbres)//'.bfsa',&
              STATUS='OLD',FORM='UNFORMATTED')

         read(10) value1
         read(11) value2
         close(10)
         close(11)
!-------------------------------------------------------------------------
! Assign quarterly albedo fractions to each tile.
!-------------------------------------------------------------------------
         do i=1,lis%d%nch
            if ((value1(tile(i)%col,tile(i)%row-tnroffset).ne.-9999.000) &
               .and.(value2(tile(i)%col,tile(i)%row-tnroffset)           &
               .ne.-9999.000)) then
               hyssib(i)%albsf1= value1(tile(i)%col,tile(i)%row-tnroffset)
               hyssib(i)%albsf2= value2(tile(i)%col,tile(i)%row-tnroffset)
            endif
         enddo
      endif                     ! End albflag selection
!-------------------------------------------------------------------------
! Assign albedo fractions to each tile and interpolate daily.
!-------------------------------------------------------------------------
      if (hyssibdrv%hyssib_albdchk.ne.lis%t%da) then
         hyssibdrv%hyssib_aflag = 1
         do i=1,lis%d%nch
            if (hyssib(i)%albsf1.ne.-9999.000) then
               valdif(i) = hyssib(i)%albsf2 - hyssib(i)%albsf1
               hyssib(i)%albsf = (timdif*valdif(i)/qdif)+hyssib(i)%albsf1
            endif
         enddo
         hyssibdrv%hyssib_albdchk=lis%t%da

         if (lis%o%wparam.eq.1) then
            allocate(albout(lis%d%lnc,lis%d%lnr))
            do i=1,lis%d%nch
               if (grid(i)%lat*1000.ge.lis%d%gridDesc(4).and. &
                   grid(i)%lat*1000.le.lis%d%gridDesc(7).and. &
                   grid(i)%lon*1000.ge.lis%d%gridDesc(5).and. &
                   grid(i)%lon*1000.le.lis%d%gridDesc(8)) then
                  rindex = tile(i)%row - (lis%d%gridDesc(4)-lis%d%gridDesc(44)) &
                           /lis%d%gridDesc(9)
                  cindex = tile(i)%col - (lis%d%gridDesc(5)-lis%d%gridDesc(45)) &
                           /lis%d%gridDesc(10)
                  albout(cindex,rindex) = hyssib(i)%albsf
               endif
            enddo
            open(32,file="albout.bin",form='unformatted')
            write(32) albout
            close(32)
            deallocate(albout)
         endif
      endif                     ! End daily interpolation
      return
!EOC
      end subroutine hyssib_albedo

