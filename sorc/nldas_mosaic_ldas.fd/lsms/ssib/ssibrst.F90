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
! !ROUTINE: ssibrst.F90
!
! !DESCRIPTION:
!  This program reads restart files for SSiB.  This includes all
!  relevant water/energy storages, tile information, and time information.
!  It also rectifies changes in the tile space.  
!
! !REVISION HISTORY:
!  1 Oct 1999: Jared Entin, Initial Code
! 15 Oct 1999: Paul Houser, Significant F90 Revision
! 05 Sep 2001: Brian Cosgrove, Modified code to use Dag Lohmann's NOAA
!              initial conditions if necessary.  This is controlled with
!              local variable NOAAIC.  Normally set to 0 in this subroutine
!              but set to 1 if want to use Dag's NOAA IC's.  Changed output
!              directory structure, and commented out if-then check so that
!              directory is always made.
! 28 Apr 2002: Kristi Arsenault, Added SSIB LSM into LDAS
! 28 May 2002: Kristi Arsenault, For STARTCODE=4, corrected SNEQV values
!              and put SMC, SH2O, STC limit for GDAS and GEOS forcing.
! 30 Oct 2003: Matt Rodell, Added back COL,ROW,FGRD,VEGT to restart files
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
!  5 May 2004: David Mocko, made compatible with SiB-lings
!
! RESTART FILE FORMAT(fortran sequential binary):
!  YR,MO,DA,HR,MN,SS,VCLASS,NCH !Restart time,Veg class,no.tiles, no.soil lay
!  TILE(NCH)%COL        !Grid Col of Tile
!  TILE(NCH)%ROW        !Grid Row of Tile
!  TILE(NCH)%FGRD       !Fraction of Grid covered by Tile
!  TILE(NCH)%VEGT       !Vegetation Type of Tile
!  SSIB(NCH)%STATES     !Model States in Tile Space
!
! !INTERFACE:
      SUBROUTINE SSIBRST
! !USES:
      use lisdrv_module, only : lis, grid, tile
      use ssib_varder, only : ssibdrv
      use ssib_varder      ! SSiB tile variables
      use time_manager
      use tile_spmdMod
!EOP
      implicit none

      INTEGER :: RW                 ! 1=read restart, 2=write restart
      INTEGER :: C,R,T,I,J,L,N,F    ! Loop counters
      INTEGER :: FOUND              ! Counting variable
      INTEGER :: YR,MO,DA,HR,MN,SS  ! Time variables
      INTEGER :: VCLASS,lnc,lnr,NCH

!   NOTE: Next four variables originally not included LIS
      INTEGER, ALLOCATABLE :: TMP_COL(:),TMP_ROW(:),TMP_VEGT(:)
      INTEGER, ALLOCATABLE :: SIB_VEGT(:)
      REAL, ALLOCATABLE :: TMP_FGRD(:)

      REAL :: RHOICE=917.0          ! Density of ice
      REAL :: WT1,WT2               ! Weights for soil wetness initialization

      CHARACTER*80 FILEN,MKFYRMO
      CHARACTER*1  FNAME(80),FBASE(80),FSUBS(80),FMKDIR(80)
      CHARACTER*1  FTIME(10),FYRMODIR(80)
      INTEGER K
      INTEGER NOAAIC   ! 0=Use IC's from card file, 1=Use NOAA IC's from

      integer :: curSec
      real, allocatable :: tmptile(:)
      PARAMETER (NOAAIC=0)
!BOC
!=== End Variable Definition =============================================

!=== Read Active Archive File ============================================
      print*,'DBG: ssibrst -- in ssibrst',' (',iam,')'
!<kluge don't do ssibrst -- see ssib_coldstart>
!</kluge don't do ssibrst -- see ssib_coldstart>
!      write(*,*) masterproc

      if (masterproc) then
         if (LIS%O%STARTCODE.EQ.1) then
            allocate(tmptile(lis%d%nch))
            OPEN(40,FILE=ssibdrv%SSIB_RFILE,FORM='unformatted')

            call timemgr_read_restart(40)
!            call timemgr_restart()
!            call get_curr_date(lis%t%yr,lis%t%mo,lis%t%da,curSec)
!            call sec2time(curSec,lis%t%hr,lis%t%mn,lis%t%ss)
!            call updatetime(lis%t) !Updates LIS variables.

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
               write(*,*) 'Restart tile definition mismatch, halting.'
               stop
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
      end subroutine ssibrst

