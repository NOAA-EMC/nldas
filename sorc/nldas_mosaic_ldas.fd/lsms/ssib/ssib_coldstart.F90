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
! !ROUTINE: ssib_coldstart.F90
!
! !DESCRIPTION:
!  Routine for SSiB initialization from cold start
!
! !REVISION HISTORY:
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
! 18 May 2004: David Mocko, made compatible with SiB-lings
!
! !INTERFACE:
      subroutine ssib_coldstart
! !USES:
      use lisdrv_module, only : lis, grid, tile,gindex
      use ssib_varder      ! SSIB tile variables
      use time_manager
      use tile_spmdMod
!EOP
      implicit none

      INTEGER :: RW              ! 1=read restart, 2=write restart
      INTEGER :: C,R,T,I,J,L,N,F ! Loop counters
      INTEGER :: FOUND           ! Counting variable
      INTEGER :: YR,MO,DA,HR,MN,SS  !Time variables
      INTEGER :: VCLASS,lnc,lnr,NCH,io

!   NOTE: Next four variables originally not included LIS
      INTEGER, ALLOCATABLE :: TMP_COL(:),TMP_ROW(:),TMP_VEGT(:)
      INTEGER, ALLOCATABLE :: SIB_VEGT(:)
      REAL, ALLOCATABLE :: TMP_FGRD(:)
      REAL :: RHOICE=917.0           ! Density of ice
      REAL :: WT1,WT2                ! Weights for soil wetness initialization

      CHARACTER*80 FILEN,MKFYRMO
      CHARACTER*1  FNAME(80),FBASE(80),FSUBS(80),FMKDIR(80)
      CHARACTER*1  FTIME(10),FYRMODIR(80)
      INTEGER K
      INTEGER NOAAIC   ! 0=Use IC's from card file, 1=Use NOAA IC's from

      integer :: curSec
      real, allocatable :: tmptile(:)
      real varfield(lis%d%lnc,lis%d%lnr)
      real vecvarfield(5,lis%d%nch)
      PARAMETER (NOAAIC=0)
!BOC
      if (ssibdrv%SSIB_FLGRES.eq.0) then
      if (lis%o%startcode.eq.2) then
         print*,'MSG: ssib_coldstart -- cold-starting ssib', &
                '...using ics from card file',' (', iam, ')'
         print*,'DBG: ssib_coldstart -- nch',lis%d%nch,' (', iam, ')'
!      write(*,*)lis%d%nch,lis%d%lnc*lis%d%lnr
         varfield=0.0
!      open(1,file='/u2/gustavo/LIS2.3/LIS/tmp/INITIAL/init.txt',&
!      status='unknown')
         do io=1,5
!         read(1,*) varfield
            do c = 1,lis%d%lnc
               do r = 1,lis%d%lnr
                  if (gindex(c,r).ne.-1) then 
                     vecvarfield(io,gindex(c,r)) = varfield(c,r)
                     if (vecvarfield(io,gindex(c,r)).le.9.9989996E+20) then 
                        if(io.eq.1)then
                           vecvarfield(io,gindex(c,r)) = ssibdrv%SSIB_ISM
                        endif
                        if(io.eq.2)then
                           vecvarfield(io,gindex(c,r)) = ssibdrv%SSIB_ISM
                        endif
                        if(io.eq.3)then
                           vecvarfield(io,gindex(c,r)) = ssibdrv%SSIB_IT
                        endif
                        if(io.eq.4)then
                           vecvarfield(io,gindex(c,r)) = ssibdrv%SSIB_IT
                        endif
                        if(io.eq.5)then
                           vecvarfield(io,gindex(c,r)) = ssibdrv%SSIB_IT
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
!      CLOSE(1)

         do t=1,lis%d%nch
            SSIB(T)%TMINI       = vecvarfield(5,t) ! TC=TA=TM
            SSIB(T)%TCINI       = vecvarfield(5,t) ! TC=TA=TM
            SSIB(T)%TAINI       = vecvarfield(5,t) ! TC=TA=TM
            SSIB(T)%HTINI       = 0.0 ! no need to specify
            SSIB(T)%QAINI       = 0.0 ! no need to specify
            SSIB(T)%TGSINI      = vecvarfield(3,t) ! top soil temperature
            SSIB(T)%TDINI       = vecvarfield(4,t) ! deep soil temperature
            SSIB(T)%WWWINI(1)   = vecvarfield(1,t) ! top soil moisture
            SSIB(T)%WWWINI(2)   = vecvarfield(2,t) ! 2 layer soil moisture
            SSIB(T)%WWWINI(3)   = vecvarfield(2,t) ! 3 layer soil moisture
            SSIB(T)%CAPACINI(1) = 0.0
            SSIB(T)%CAPACINI(2) = 0.0
         enddo                  ! end of tile loop for card options

         lis%t%yr=lis%t%syr
         lis%t%mo=lis%t%smo
         lis%t%da=lis%t%sda
         lis%t%hr=lis%t%shr
         lis%t%mn=lis%t%smn
         lis%t%ss=lis%t%sss

         call date2time(lis%t%time,lis%t%doy,lis%t%gmt,lis%t%yr, &
                        lis%t%mo,lis%t%da,lis%t%hr,lis%t%mn,lis%t%ss)
         write(*,*) 'MSG: ssib_coldstart -- Using lis.crd start time ',&
                    lis%t%time, ' (', iam, ')'
      endif
      endif

      if (ssibdrv%SSIB_FLGRES.eq.1) then
         if (LIS%O%STARTCODE.EQ.2) then
            allocate(tmptile(lis%d%nch))
            OPEN(40,FILE=ssibdrv%SSIB_RFILE,FORM='unformatted')

            call timemgr_read_restart(40)
!        call timemgr_restart()
!        call get_curr_date(lis%t%yr,lis%t%mo,lis%t%da,curSec)
!        call sec2time(curSec,lis%t%hr,lis%t%mn,lis%t%ss)
!        call updatetime(lis%t) !Updates LIS variables.

            WRITE(*,*)'SSIB Restart File Used: ',ssibdrv%SSIB_RFILE

            READ(40) VCLASS,lnc,lnr,NCH !Time, veg class, no. tiles

!   Check for Vegetation Class Conflict
            IF (VCLASS.NE.LIS%P%VCLASS) THEN
               WRITE(*,*) ssibdrv%SSIB_RFILE, &
                          'Vegetation class conflict'
               STOP
            ENDIF

!   Check for Grid Space Conflict
            IF (lnc.NE.LIS%D%lnc.OR.lnr.NE.LIS%D%lnr) THEN
               WRITE(*,*) ssibdrv%SSIB_RFILE, &
                          'Grid space mismatch - SSiB HALTED'
               STOP
            ENDIF

!-- Transfer Restart tile space to LIS tile space
            IF (NCH.NE.LIS%D%NCH) THEN
               WRITE(*,*)'Restart Tile Space Mismatch, Halting..'
               stop
            endif

!-- Allocate temporary arrays.
            allocate(TMP_COL(NCH))
            allocate(TMP_ROW(NCH))
            allocate(TMP_FGRD(NCH))
            allocate(TMP_VEGT(NCH))
            allocate(SIB_VEGT(NCH))

!-- NOTE: Next four read statements originally removed from LIS
            READ(40) TMP_COL    !Grid Col of Tile
            READ(40) TMP_ROW    !Grid Row of Tile
            READ(40) TMP_FGRD   !Fraction of Grid covered by tile
            READ(40) TMP_VEGT   !Vegetation Type of Tile

!*** TESTING
!        open(unit=31, file='test.arrays.ssibrst.txt', &
!	     form='formatted')
!	write(31,*) '  COL TCOL  ROW TROW  VEGT TVEGT'
!	do t=1,nch
!	  write(31,310) TILE(t)%COL,TMP_COL(t),TILE(t)%ROW,TMP_ROW(t),&
!	        TILE(t)%VEGT,TMP_VEGT(t)
!	enddo
!  310   format(6i5)
!  	close(31)

!-- Convert maketiles UMD vegetation classes to Noah SiB classes before
!--  testing for tile definition conflict.
            SIB_VEGT = TILE%VEGT
            do t=1,nch
               call SSIB_MAPVEGC(SIB_VEGT(t))
            enddo

!-- Check for tile definition conflict
            if ((sum(TMP_COL - TILE%COL) .ne. 0) .or. &
                (sum(TMP_ROW - TILE%ROW) .ne. 0) .or. &
                (sum(TMP_VEGT - SIB_VEGT) .ne. 0)) then
               write(*,*) sum(TMP_vegt),sum(sib_vegt)
               write(*,*) 'Restart tile definition mismatch, halting.'
!               stop
            endif
            READ(40) SSIB%TCINI
            READ(40) SSIB%TGSINI
            READ(40) SSIB%TDINI
            READ(40) SSIB%TAINI
            READ(40) SSIB%TMINI
            READ(40) SSIB%HTINI
            READ(40) SSIB%QAINI
            DO L=1,3
               READ(40) TMPTILE
               SSIB%WWWINI(L) = TMPTILE
            ENDDO
            DO L=1,2
               READ(40) TMPTILE
               SSIB%CAPACINI(L) = TMPTILE
            ENDDO

            CLOSE(40)

            deallocate(tmptile)
            deallocate(TMP_COL)
            deallocate(TMP_ROW)
            deallocate(TMP_FGRD)
            deallocate(TMP_VEGT)
            deallocate(SIB_VEGT)
         endif
      endif
      return
!EOC
      end subroutine ssib_coldstart

