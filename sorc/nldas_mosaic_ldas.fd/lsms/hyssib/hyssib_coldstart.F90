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
! !ROUTINE: hyssib_coldstart.F90
!
! !DESCRIPTION:
! Routine for HY-SSiB initialization from cold start
! 
! !REVISION HISTORY:
!    Dec 2003: Luis-Gustavo Goncalves, Initial version
! 29 Apr 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      subroutine hyssib_coldstart()
! !USES:
      use lisdrv_module, only : lis, grid, tile
      use hyssib_varder      ! HY-SSIB tile variables
      use time_manager
      use tile_spmdMod
!EOP
      implicit none

      INTEGER :: RW              ! 1=read restart, 2=write restart
      INTEGER :: C,R,T,I,J,L,N,F ! Loop counters
      INTEGER :: FOUND           ! Counting variable
      INTEGER :: YR,MO,DA,HR,MN,SS  !Time variables
      INTEGER :: VCLASS,lnc,lnr,NCH

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
      PARAMETER (NOAAIC=0)
!BOC
      if (lis%o%startcode.eq.2) then
         print*,'MSG: hyssib_coldstart -- cold-starting hyssib', &
                '...using ics from card file',' (', iam, ')'
         print*,'DBG: hyssib_coldstart -- nch',lis%d%nch,' (', iam, ')'
         do t=1,lis%d%nch
            HYSSIB(T)%TC          =  hyssibdrv%HYSSIB_IT
            HYSSIB(T)%TG          =  hyssibdrv%HYSSIB_IT
            HYSSIB(T)%TSN         =  hyssibdrv%HYSSIB_IT
            HYSSIB(T)%TD          =  hyssibdrv%HYSSIB_IT
            HYSSIB(T)%WWW(1)      =  hyssibdrv%HYSSIB_ISM
            HYSSIB(T)%WWW(2)      =  hyssibdrv%HYSSIB_ISM
            HYSSIB(T)%WWW(3)      =  hyssibdrv%HYSSIB_ISM
            HYSSIB(T)%CAPAC(1)    =  0.0
            HYSSIB(T)%CAPAC(2)    =  0.0
            HYSSIB(T)%SNOW(1)     =  0.0
            HYSSIB(T)%SNOW(2)     =  0.0
            HYSSIB(T)%SGFG        =  0.0
            HYSSIB(T)%SDENS       =  0.0
         enddo                  ! end of tile loop for card options

         lis%t%yr=lis%t%syr
         lis%t%mo=lis%t%smo
         lis%t%da=lis%t%sda
         lis%t%hr=lis%t%shr
         lis%t%mn=lis%t%smn
         lis%t%ss=lis%t%sss

         call date2time(lis%t%time,lis%t%doy,lis%t%gmt,lis%t%yr, &
                        lis%t%mo,lis%t%da,lis%t%hr,lis%t%mn,lis%t%ss)
         write(*,*) 'MSG: hyssib_coldstart -- Using lis.crd start time ', &
                    lis%t%time, ' (', iam, ')'
      endif

      return
!EOC
      end subroutine hyssib_coldstart

