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
! !ROUTINE: mosrst.F90
!
! !DESCRIPTION:
!  This program reads restart files for Mosaic.  This
!   includes all relevant water/energy storages, tile information,
!   and time information.  It also rectifies changes in the tile space.  
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  19 Jan 2001: Brian Cosgrove; Added CLOSE statement
!  15 Mar 2001: Jon Gottschalck; Added option 4 (initialize Mosaic with GDAS forcing)
!  17 Apr 2001: Jon Gottschalck; Modified option 4 to work with GEOS forcing
!  05 Sep 2001: Brian Cosgrove; Modified code to use Dag Lohmann's NOAA
!               initial conditions if necessary.  This is controlled with
!               local variable NOAAIC.  Normally set to 0 in this subroutine
!               but set to 1 if want to use Dag's NOAA IC's.  Changed output
!               directory structure, and commented out if-then check so that
!               directory is always made.
!  18 Sep 2001: Brian Cosgrove; Changed use of Dag's ICs so that third mosaic
!               soil layer is set to layer 2's soil moisture since layer
!               3 doesn't have a wilting point
!  04 Feb 2002: Jon Gottschalck; Added section to set all Koster tiles of veg
!               type 9 (land ice) to have tons of snow so it never melts since the soil/veg
!               parameters are set for bare soil (Koster -- personal communication).
!  05 Feb 2002: Brian Cosgrove; Added a few diagnostic vars to the NOAAIC=1 option
!  20 Nov 2002: Jon Radakovich; Updated for PSAS temperature assimilation and BC so
!               the forecast bias is included in the restart file.  Added
!               conditionals based on startcode=5 (Using spun-up restart so model time comes
!               from card file) and startcode=6 (Restarting a bias correction run).
! 12 Dec. 2002: Brian Cosgrove; Fixed usage of Wiltpoint variable.  Before,
!               Wiltpoint1 and Wiltpoint2 were used in calculation of
!               root zone soil moisture availability...now, only Wiltpoint2
!               is used since wiltpoint1 is not the correct wilting point
!               needed for the calculation.
!  14 Jan 2003: Urszula Jambor; Added deallocation statements near end of routine
!               and changed pointer variables to allocatable.
!  23 Jan 2003: Urszula Jambor; Switch index order of GEOS forcing
!               array.  Snow is 12, soil wetness is 13.
!=========================================================================
! RESTART FILE FORMAT(fortran sequential binary):
!  VCLASS,NC,NR,NCH     !Veg class,no. columns, no. rows, no.tiles
!  MOS(NCH)%STATES      !Model States in Tile Space
!=========================================================================
! 
! !INTERFACE:
subroutine MOSrst
! !USES:
  use lisdrv_module, only : lis, grid, tile
  use mos_varder, only : mosdrv
  USE mos_varder      ! MOSAIC tile variables
  use time_manager
  use tile_spmdMod
!EOP
  IMPLICIT NONE      
  
  INTEGER :: RW              ! 1=read restart, 2=write restart
  INTEGER :: C,R,T,I,J,L,N,F ! Loop counters
  INTEGER :: FOUND           ! Counting variable
  
  INTEGER :: YR,MO,DA,HR,MN,SS  !Time variables
  INTEGER :: VCLASS,NC,NR,NCH
  
  REAL :: RHOICE=917.0              ! Density of ice
  REAL :: WT1,WT2                   ! Weights for soil wetness initialization

  CHARACTER*80 FILEN,MKFYRMO
  CHARACTER*1  FNAME(80),FBASE(80),FSUBS(80),FMKDIR(80)
  CHARACTER*1  FTIME(10),FYRMODIR(80)
  INTEGER K
  INTEGER NOAAIC   ! 0=Use IC's from card file, 1=Use NOAA IC's from 
  
  integer :: curSec
  real, allocatable :: tmptile(:)
      REAL  VGZDEX1(lis%d%nch),VGZDEX2(lis%d%nch),VGZDEX3(lis%d%nch)
      REAL  POR1(lis%d%nch),POR2(lis%d%nch),POR3(lis%d%nch)
      REAL  VGWMAX1(lis%d%nch),VGWMAX2(lis%d%nch),VGWMAX3(lis%d%nch)
      REAL  VGPH1X(lis%d%nch),VGPH2X(lis%d%nch),VGPH3X(lis%d%nch)
      REAL  VGBEEX(lis%d%nch),VGPSAX(lis%d%nch)
      REAL  Wiltpoint1(lis%d%nch),Wiltpoint2(lis%d%nch)
      REAL  Wiltpoint3(lis%d%nch)
      REAL SOILMOIST(lis%d%lnc,lis%d%lnr),SOILTEMP(lis%d%lnc,lis%d%lnr)

      REAL SOILWW,SOILWM

  PARAMETER (NOAAIC=0)

!=== End Variable Definition =============================================
!BOC
!-------------------------------------------------------------------------
! Read Active Archive File
!-------------------------------------------------------------------------
	print *,'in mosrst, lis%d%nch=',lis%d%nch
  print*,'DBG: mosrst -- in mosrst',' (',iam,')',lis%o%startcode
  if(masterproc) then 
     IF(LIS%O%STARTCODE.EQ.1)THEN
        allocate(tmptile(lis%d%nch))
        OPEN(40,FILE=mosdrv%MOS_RFILE,FORM='unformatted')
        
        !call timemgr_read_restart(40)
!        call timemgr_restart()
!        call get_curr_date(lis%t%yr,lis%t%mo,lis%t%da,curSec)
!        call sec2time(curSec,lis%t%hr,lis%t%mn,lis%t%ss)
!        call updatetime(lis%t) !Updates LIS variables.
        WRITE(*,*)'MOSAIC Restart File Used: ',mosdrv%MOS_RFILE
        !READ(40) VCLASS,NC,NR,NCH  !Veg class, no. columns, no. rows, no. tiles    
        READ(40) NC,NR,NCH  !Veg class, no. columns, no. rows, no. tiles
	print *,'vclass,nc,nr,nch=',vclass,nc,nr,nch
!------------------------------------------------------------------------
!   Check for Vegetation Class Conflict 
!------------------------------------------------------------------------
        !IF(VCLASS.NE.LIS%P%VCLASS)THEN
        !   WRITE(*,*)mosdrv%MOS_RFILE,' Vegetation class conflict'
        !   call endrun
        !ENDIF
!------------------------------------------------------------------------
!   Check for Grid Space Conflict 
!------------------------------------------------------------------------
        IF(NC.NE.LIS%D%LNC.OR.NR.NE.LIS%D%LNR)THEN
           WRITE(*,*) mosdrv%MOS_RFILE,'Grid space mismatch - MOSAIC HALTED'
           WRITE(*,*) nc,nr,lis%d%lnc,lis%d%lnr
           call endrun
        ENDIF
!------------------------------------------------------------------------
! Transfer Restart tile space to LIS tile space
!------------------------------------------------------------------------
	print *,'lis%d%nch=',lis%d%nch
        IF(NCH.NE.LIS%D%NCH)THEN           
           WRITE(*,*)'Restart Tile Space Mismatch, Halting..'
	   print *,'nch,lis%d%nch=',nch,lis%d%nch
           call endrun
        endif
        read(40) mos%ct         !MOSAIC Canopy/Soil Temperature
	read(40) mos%qa         !MOSAIC Canopy Humidity
        read(40) mos%ics        !MOSAIC Interception Canopy Storage
        read(40) mos%snow       !MOSAIC Snow Depth
	read(40) mos%SoT        !MOSAIC Deep Soil Temperaure
        do l=1,3
	  read(40) tmptile  !MOSAIC Soil Wetness (3 layers)
	  mos%SoWET(l)=tmptile
        enddo 
        close(40) 
       
	if (NOAAIC.eq.1) then
!       READ IN INITIAL SOIL TEMP AND SOIL WETNESS FROM NOAA's initial
!       conditions from Dag
	print *,'reading DAGs conditions'
        print *,'reading DAGs conditions'
        print *,'reading DAGs conditions'
        print *,'reading DAGs conditions'
        print *,'reading DAGs conditions'
        print *,'reading DAGs conditions'
        print *,'reading DAGs conditions'
        print *,'reading DAGs conditions'

        open (unit=22,file= &
         './input/BCS/1_8deg/1996093012.moist_scale.bin', &
       form='unformatted')
        read(22) ((soilmoist(i,j),i=1,464),j=1,224)
        close (22)
        open (unit=22,file= &
         './input/BCS/1_8deg/1996093012.temp2.bin', &
       form='unformatted')
        read(22)((soiltemp(i,j),i=1,464),j=1,224)
        close(22)

        DO T=1,lis%d%nch
        VGZDEX1(T)=  MOS(T)%SOILP(4)
        VGZDEX2(T)=  MOS(T)%SOILP(5)
        VGZDEX3(T)=  MOS(T)%SOILP(6)
        VGWMAX1(T)=  MOS(T)%SOILP(8)
        VGWMAX2(T)=  MOS(T)%SOILP(9)
        VGWMAX3(T)=  MOS(T)%SOILP(10)

        POR1(T)=VGWMAX1(T)/ &
                  (VGZDEX1(T)*1000.)
        POR2(T)=VGWMAX2(T)/ &
                  (VGZDEX2(T)*1000.)
        POR3(T)=VGWMAX3(T)/ &
                  (VGZDEX3(T)*1000.)


        VGPH1X(T)=MOS(T)%VEGP(10)
        VGPH2X(T)=MOS(T)%VEGP(11)
        VGPSAX(T)=MOS(T)%SOILP(2)
        VGBEEX(T)=MOS(T)%SOILP(1)
        Wiltpoint1(T)=((VGPH1X(T)/VGPSAX(T)) ** (1.0 / (-VGBEEX(T))))
        Wiltpoint2(T)=((VGPH2X(T)/VGPSAX(T)) ** (1.0 / (-VGBEEX(T))))
        Wiltpoint3(T)=0.0
        ENDDO

       DO T=1,lis%d%nch
        MOS(T)%CT=MOSDRV%MOS_IT
        MOS(T)%QA=0.0
        MOS(T)%ICS=0.0
        MOS(T)%SNOW=0.0

         MOS(T)%SoT=SOILTEMP(TILE(T)%COL,TILE(T)%ROW)
         MOS(T)%SoWET(1)= ((soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR1(T)) &
       -(soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR1(T)*Wiltpoint2(T))+ &
       (POR1(T)*Wiltpoint2(T)))/POR1(T)
         MOS(T)%SoWET(2)= ((soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR2(T)) &
       -(soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR2(T)*Wiltpoint2(T))+ &
       (POR2(T)*Wiltpoint2(T)))/POR2(T)
!         MOS(T)%SoWET(3)= ((soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR3(T))
!     &  -(soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR3(T)*Wiltpoint2(T))+
!     &  (POR3(T)*Wiltpoint2(T)))/POR3(T)
!       Wiltpoint of layer 3 is 0 so need to set layer 3 sowet to layer 2 sowet
          MOS(T)%SoWET(3)=MOS(T)%SoWET(2)

        if (t.eq.61444) then
        print *,'soilwet 61444',MOS(T)%SoWET(1),MOS(T)%SoWET(2), &
       MOS(T)%SoWET(3)
        endif
        if (t.eq.61447) then
        print *,'soilwet 61447',MOS(T)%SoWET(1),MOS(T)%SoWET(2), &
        MOS(T)%SoWET(3)
        endif

       SOILWM=0.0
       SOILWW=0.0
       SOILWM=(POR1(T)-(POR1(T)*WILTPOINT2(T)))*VGZDEX1(T)
       SOILWM=SOILWM+((POR2(T)-(POR2(T)*WILTPOINT2(T)))*VGZDEX2(T))

       SOILWW=((MOS(T)%SoWET(1)*POR1(T))-(POR1(T)*WILTPOINT2(T)))* &
       VGZDEX1(T)
       SOILWW=SOILWW+(((MOS(T)%SoWET(2)*POR2(T))- &
             (POR2(T)*WILTPOINT2(T)))*VGZDEX2(T))

        if (t.eq.61444) print *,'61444 vgzdex1,2,3', &
       VGZDEX1(T),VGZDEX2(T),VGZDEX3(T)
        if (t.eq.61444) print *,'61444 vgwmax1,2,3', &
       VGWMAX1(T),VGWMAX2(T),VGWMAX3(T)
        if (t.eq.61444) print *,'61444 POR1,2,3', &
       POR1(T),POR2(T),POR3(T)
        if (t.eq.61444) print *,'61444 VGPH1X,VGPH2X,VGPSAX', &
       VGPH1X(T),VGPH2X(T),VGPSAX(T)
        if (t.eq.61444) print *,'61444 VGBEEX,Wilt1,Wilt2,Wilt3',&
       VGBEEX(T),Wiltpoint1(T),Wiltpoint2(T),Wiltpoint3(T)

        if (t.eq.61447) print *,'61447 vgzdex1,2,3',&
       VGZDEX1(T),VGZDEX2(T),VGZDEX3(T)
        if (t.eq.61447) print *,'61447 vgwmax1,2,3',&
       VGWMAX1(T),VGWMAX2(T),VGWMAX3(T)
        if (t.eq.61447) print *,'61447 POR1,2,3',&
       POR1(T),POR2(T),POR3(T)
        if (t.eq.61447) print *,'61447 VGPH1X,VGPH2X,VGPSAX',&
       VGPH1X(T),VGPH2X(T),VGPSAX(T)
        if (t.eq.61447) print *,'61447 VGBEEX,Wilt1,Wilt2,Wilt3',&
       VGBEEX(T),Wiltpoint1(T),Wiltpoint2(T),Wiltpoint3(T)
	ENDDO
	ENDIF !endif reading in Dags initial conditions

        deallocate(tmptile)
     endif
  endif

!EOC
end subroutine mosrst
