!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! mosrst.f: 
!
! DESCRIPTION:
!  This program reads and writes restart files for MOSAIC.  This
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
!               the forecast bias is included in the restart file.  Added conditionals
!               based on startcode=5 (Using spun-up restart so model time comes
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
!  03 Jun 2003: Jon Gottschalck; Modified model initialization section      
!=========================================================================
! RESTART FILE FORMAT(fortran sequential binary):
!  YR,MO,DA,HR,MN,SS,VCLASS,NCH !Restart time,Veg class,no.tiles, no.soil lay 
!  TILE(NCH)%COL        !Grid Col of Tile   
!  TILE(NCH)%ROW        !Grid Row of Tile
!  TILE(NCH)%FGRD       !Fraction of Grid covered by tile
!  TILE(NCH)%VEGT       !Vegetation Type of Tile
!  MOS(NCH)%STATES      !Model States in Tile Space
!=========================================================================

      SUBROUTINE MOSRST(RW,LDAS,TILE,GRID,MOS)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      USE grid_module      ! LDAS non-model-specific grid variables
      USE mos_module       ! Mosaic tile variables
      IMPLICIT NONE      
      type (ldasdec) LDAS              
      type (tiledec) TILE(LDAS%NCH)
      type (griddec) GRID(LDAS%NC,LDAS%NR)   
      type (mosdec)  MOS(LDAS%NCH)

!=== Local Variables =====================================================
      INTEGER :: RW              ! 1=read restart, 2=write restart
      INTEGER :: C,R,T,I,J,L,N,F ! Loop counters
      INTEGER :: FOUND           ! Counting variable

!=== Temporary tile space transfer files (different than in ldas_module)
      INTEGER :: YR,MO,DA,HR,MN,SS  !Time variables
      INTEGER :: VCLASS,NC,NR,NCH
      INTEGER, allocatable :: COL(:)    ! Column
      INTEGER, allocatable :: ROW(:)    ! Row
      INTEGER, allocatable :: VEGT(:)   ! Tile veg type
      
      
      REAL, allocatable :: FGRD(:)      ! Grid Fraction of Tile      
      REAL, allocatable :: CT(:)        ! Mosaic Canopy/Soil Temperature
      REAL, allocatable :: QA(:)        ! Mosaic Canopy Humidity
      REAL, allocatable :: ICS(:)       ! Mosaic Interception
      REAL, allocatable :: SNOW(:)      ! Mosaic Snow Depth
      REAL, allocatable :: SoT(:)       ! Mosaic Deep Soil Temperature
      REAL, allocatable :: SoWET(:,:)   ! Mosaic Soil Wetness
      REAL, allocatable :: TMPTILEO(:)  ! Temporary Transfer Array   
      REAL :: TMPTILEN(LDAS%NCH)    ! Temporary Transfer Array   
      REAL :: TMPG(LDAS%NC,LDAS%NR) ! Temporary Transfer Array   

      REAL :: G_CT(LDAS%NC,LDAS%NR)      ! Mosaic Canopy/Soil Temperature
      REAL :: G_QA(LDAS%NC,LDAS%NR)      ! Mosaic Canopy Humidity
      REAL :: G_ICS(LDAS%NC,LDAS%NR)     ! Mosaic Interception
      REAL :: G_SNOW(LDAS%NC,LDAS%NR)    ! Mosaic Snow Depth
      REAL :: G_SoT(LDAS%NC,LDAS%NR)     ! Mosaic Deep Soil Temperature
      REAL :: G_FBIAS(LDAS%NC,LDAS%NR)   ! Mosaic Bias Term
      REAL :: G_SoWET(LDAS%NC,LDAS%NR,3) ! Mosaic Soil Wetness
      REAL :: WT1,WT2                    ! Weights for soil wetness initialization
      REAL :: RHOICE=917.0               ! Density of ice

      CHARACTER*80 FILEN,MKFYRMO
      CHARACTER*1  FNAME(80),FBASE(80),FSUBS(80),FMKDIR(80)
      CHARACTER*1  FTIME(10),FYRMODIR(80)
      INTEGER K
      INTEGER NOAAIC   ! 0=Use IC's from card file, 1=Use NOAA IC's from 
                       !                            Dag's data files 
      REAL  VGZDEX1(LDAS%NCH),VGZDEX2(LDAS%NCH),VGZDEX3(LDAS%NCH)
      REAL  POR1(LDAS%NCH),POR2(LDAS%NCH),POR3(LDAS%NCH)
      REAL  VGWMAX1(LDAS%NCH),VGWMAX2(LDAS%NCH),VGWMAX3(LDAS%NCH)
      REAL  VGPH1X(LDAS%NCH),VGPH2X(LDAS%NCH),VGPH3X(LDAS%NCH)
      REAL  VGBEEX(LDAS%NCH),VGPSAX(LDAS%NCH)
      REAL  Wiltpoint1(LDAS%NCH),Wiltpoint2(LDAS%NCH)
      REAL  Wiltpoint3(LDAS%NCH)
      REAL SOILMOIST(LDAS%NC,LDAS%NR),SOILTEMP(LDAS%NC,LDAS%NR)
      REAL VGZDEX(LDAS%NCH,3),VGWMAX(LDAS%NCH,3),POR(LDAS%NCH,3)     
 
      REAL SOILWW,SOILWM
      PARAMETER (NOAAIC=0)

!=== End Variable Definition =============================================

!=== Read Active Archive File ============================================
      IF((RW.EQ.1.AND.LDAS%MOS_IC.EQ.1).OR.
     &   (RW.EQ.1.AND.LDAS%STARTCODE.EQ.1).OR.
     &   (RW.EQ.1.AND.LDAS%STARTCODE.EQ.5).OR.
     &   (RW.EQ.1.AND.LDAS%STARTCODE.EQ.6))THEN

       OPEN(40,FILE=LDAS%MOS_RFILE,FORM='unformatted')

       READ(40) YR,MO,DA,HR,MN,SS,VCLASS,NC,NR,NCH  !Time, veg class, no. tiles
       
       ALLOCATE (COL(NCH),ROW(NCH),FGRD(NCH),VEGT(NCH))
       ALLOCATE (CT(NCH),QA(NCH),ICS(NCH),SNOW(NCH),SoT(NCH))
       ALLOCATE (SoWET(NCH,3),TMPTILEO(NCH))

       READ(40) COL        !Grid Col of Tile   
       READ(40) ROW        !Grid Row of Tile
       READ(40) FGRD       !Fraction of Grid covered by tile
       READ(40) VEGT       !Vegetation Type of Tile
       READ(40) CT         !MOSAIC Canopy/Soil Temperature 
       READ(40) QA         !MOSAIC Canopy Humidity
       READ(40) ICS        !MOSAIC Interception Canopy Storage
       READ(40) SNOW       !MOSAIC Snow Depth
       READ(40) SoT        !MOSAIC Deep Soil Temperaure
       IF(LDAS%RBIAS.EQ.1.AND.LDAS%STARTCODE.EQ.6)THEN
        READ(40) G_FBIAS      !MOSAIC Forecast Bias
        READ(40) GRID%MOSDELT  !MOSAIC Analysis Increment
        DO R=1,LDAS%NR
         DO C=1,LDAS%NC
         GRID(C,R)%MOSBETAK(1)=G_FBIAS(C,R)
         ENDDO
        ENDDO
       ENDIF
       DO L=1,3
        READ(40) TMPTILEO  !MOSAIC Soil Wetness (3 layers)
        DO T=1,NCH
         SoWET(T,L)=TMPTILEO(T)
        ENDDO
       ENDDO

       CLOSE(40)

       WRITE(*,*)'MOSAIC Restart File Read: ',LDAS%MOS_RFILE
       WRITE(*,*)'MOSAIC Restart File Time: ',YR,MO,DA,HR,MN,SS

       WRITE(79,*)'MOSAIC Restart File Read: ',LDAS%MOS_RFILE
       WRITE(79,*)'MOSAIC Restart File Time: ',YR,MO,DA,HR,MN,SS

       IF(LDAS%STARTCODE.EQ.1.OR.LDAS%STARTCODE.EQ.6)THEN
! Establish Model Restart Time
        LDAS%YR=YR
        LDAS%MO=MO 
        LDAS%DA=DA
        LDAS%HR=HR
        LDAS%MN=MN
        LDAS%SS=SS
        CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,YR,MO,DA,
     &                 HR,MN,SS) 
        LDAS%MOS_STIME = LDAS%TIME
        WRITE(*,*)'MOS Restart File Time Used: ',LDAS%MOS_RFILE
       ENDIF

!=== Using spun-up IC, must used card file time
      IF(RW.EQ.1.AND.LDAS%STARTCODE.EQ.5)THEN
       LDAS%YR=LDAS%SYR
       LDAS%MO=LDAS%SMO
       LDAS%DA=LDAS%SDA
       LDAS%HR=LDAS%SHR
       LDAS%MN=LDAS%SMN
       LDAS%SS=LDAS%SSS
       CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,LDAS%YR,
     &                LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS)
       WRITE(*,*)'Using ldas.crd start time for spun-up IC ',LDAS%TIME
       WRITE(79,*)'Using ldas.crd start time for spun-up IC ',LDAS%TIME
      ENDIF


!=== Rectify Restart Tile Space to LDAS Tile Space =====================
       IF(LDAS%MOS_IC.EQ.1)THEN

! Check for Vegetation Class Conflict 
        IF(VCLASS.NE.LDAS%VCLASS)THEN
         WRITE(*,*)LDAS%MOS_RFILE,' Vegetation class conflict'
         WRITE(79,*)LDAS%MOS_RFILE,' Vegetation class conflict'
         STOP
        ENDIF

! Check for Grid Space Conflict 
        IF(NC.NE.LDAS%NC.OR.NR.NE.LDAS%NR)THEN
         WRITE(*,*)LDAS%MOS_RFILE,'Grid space mismatch - MOS HALTED'
         STOP
        ENDIF

! Transfer Restart tile space to LDAS tile space
        IF(NCH.NE.LDAS%NCH)THEN
         WRITE(*,*)'Restart Tile Space Mismatch-Transfer in Progress'
         WRITE(79,*)'Restart Tile Space Mismatch-Transfer in Progress'

!  Start by finding grid averages
         CALL T2GR(CT,G_CT,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
         CALL T2GR(QA,G_QA,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
         CALL T2GR(ICS,G_ICS,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
         CALL T2GR(SNOW,G_SNOW,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
         CALL T2GR(SoT,G_SoT,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
         DO L=1,3
          CALL T2GR(SoWet(:,L),G_SoWet(:,:,L),
     1              LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
         ENDDO

! Perform state transfer
         C=0
         DO 555 T=1,LDAS%NCH 
          IF(AMOD(FLOAT(T),10000.0).EQ.0.0)WRITE(*,23)'  Transferred ',
     1     100.0*FLOAT(T)/FLOAT(LDAS%NCH),' Percent of Tiles'

          IF(AMOD(FLOAT(T),10000.0).EQ.0.0)WRITE(79,23)'  Transferred ',
     1     100.0*FLOAT(T)/FLOAT(LDAS%NCH),' Percent of Tiles'

 23       FORMAT(A14,F5.2,A17)
          FOUND=0
          DO N=1,NCH
           IF(TILE(T)%VEGT.EQ.VEGT(N).AND.
     1       TILE(T)%COL .EQ. COL(N).AND.
     2       TILE(T)%ROW .EQ. ROW(N))THEN
            MOS(T)%CT=CT(N)
            MOS(T)%QA=Qa(N)
            MOS(T)%ICS=ICS(N)
            MOS(T)%SNOW=SNOW(N)
            MOS(T)%SoT=SoT(N)
            DO L=1,3
             MOS(T)%SoWET(L)=SoWET(N,L)
            ENDDO
            FOUND=1
            GOTO 555 
           ENDIF
          ENDDO
          IF(FOUND.EQ.0)THEN        
           MOS(T)%CT=G_CT(TILE(T)%COL,TILE(T)%ROW)
           MOS(T)%QA=G_QA(TILE(T)%COL,TILE(T)%ROW)
           MOS(T)%ICS=G_ICS(TILE(T)%COL,TILE(T)%ROW)
           MOS(T)%SNOW=G_SNOW(TILE(T)%COL,TILE(T)%ROW)
           MOS(T)%SoT=G_SoT(TILE(T)%COL,TILE(T)%ROW)
           DO L=1,3
            MOS(T)%SoWET(L)=G_SoWET(TILE(T)%COL,TILE(T)%ROW,L)
           ENDDO
           C=0
          ENDIF
 555     CONTINUE
         WRITE(*,*)'Tile Space Transfer Complete'
         WRITE(*,*)'Mosaic Restart NCH:',NCH,'Current NCH:',LDAS%NCH
         WRITE(*,*) C, ' Tiles not found in old Mosaic restart'        
         WRITE(*,*)  

         WRITE(79,*)'Tile Space Transfer Complete'
         WRITE(79,*)'Mosaic Restart NCH:',NCH,'Current NCH:',LDAS%NCH
         WRITE(79,*) C, ' Tiles not found in old Mosaic restart'
         WRITE(79,*)

       ELSE !The number of tiles is a match
         DO T=1,LDAS%NCH
          MOS(T)%CT=CT(T)
          MOS(T)%QA=QA(T)
          MOS(T)%ICS=ICS(T)
          MOS(T)%SNOW=SNOW(T)
          MOS(T)%SoT=SoT(T)
          DO L=1,3
           MOS(T)%SoWET(L)=SoWET(T,L)
          ENDDO
         ENDDO	
        ENDIF   
       ENDIF
      ENDIF !RW option 1

!=== In the case of no previously existing restart, the user can 
!===  specify the restart condifions in the card file
      IF(RW.EQ.1.AND.LDAS%MOS_IC.EQ.3)THEN !User restart

	IF (NOAAIC.EQ.1) THEN
c       READ IN INITIAL SOIL TEMP AND SOIL WETNESS FROM NOAA's initial
C       conditions from Dag
        open (unit=22,file=
     &    './BCS/N0.125/1996093012.moist_scale.bin',
     &  form='unformatted')
        read(22) ((soilmoist(i,j),i=1,464),j=1,224)
        close (22)
        open (unit=22,file=
     &    './BCS/N0.125/1996093012.temp2.bin',
     &  form='unformatted')
        read(22)((soiltemp(i,j),i=1,464),j=1,224)
        close(22)

        DO T=1,LDAS%NCH
        VGZDEX1(T)=  MOS(T)%SOILP(4)
        VGZDEX2(T)=  MOS(T)%SOILP(5)
        VGZDEX3(T)=  MOS(T)%SOILP(6)
        VGWMAX1(T)=  MOS(T)%SOILP(8)
        VGWMAX2(T)=  MOS(T)%SOILP(9)
        VGWMAX3(T)=  MOS(T)%SOILP(10)

        POR1(T)=VGWMAX1(T)/
     1             (VGZDEX1(T)*1000.)
        POR2(T)=VGWMAX2(T)/
     1             (VGZDEX2(T)*1000.)
        POR3(T)=VGWMAX3(T)/
     1             (VGZDEX3(T)*1000.)


        VGPH1X(T)=MOS(T)%VEGP(10)
        VGPH2X(T)=MOS(T)%VEGP(11)
        VGPSAX(T)=MOS(T)%SOILP(2)
        VGBEEX(T)=MOS(T)%SOILP(1)
        Wiltpoint1(T)=((VGPH1X(T)/VGPSAX(T)) ** (1.0 / (-VGBEEX(T))))
        Wiltpoint2(T)=((VGPH2X(T)/VGPSAX(T)) ** (1.0 / (-VGBEEX(T))))
        Wiltpoint3(T)=0.0
        ENDDO
        ENDIF

       DO T=1,LDAS%NCH
        MOS(T)%CT=LDAS%MOS_IT
        MOS(T)%QA=0.0
        MOS(T)%ICS=0.0
        MOS(T)%SNOW=0.0

        IF (NOAAIC.EQ.1) THEN
         MOS(T)%SoT=SOILTEMP(TILE(T)%COL,TILE(T)%ROW)
         MOS(T)%SoWET(1)= ((soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR1(T))
     &  -(soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR1(T)*Wiltpoint2(T))+
     &  (POR1(T)*Wiltpoint2(T)))/POR1(T)
         MOS(T)%SoWET(2)= ((soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR2(T))
     &  -(soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR2(T)*Wiltpoint2(T))+
     &  (POR2(T)*Wiltpoint2(T)))/POR2(T)
c         MOS(T)%SoWET(3)= ((soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR3(T))
c     &  -(soilmoist(TILE(T)%COL,TILE(T)%ROW)*POR3(T)*Wiltpoint2(T))+
c     &  (POR3(T)*Wiltpoint2(T)))/POR3(T)
C	Wiltpoint of layer 3 is 0 so need to set layer 3 sowet to layer 2 sowet
          MOS(T)%SoWET(3)=MOS(T)%SoWET(2)
       SOILWM=0.0
       SOILWW=0.0
       SOILWM=(POR1(T)-(POR1(T)*WILTPOINT2(T)))*VGZDEX1(T)
       SOILWM=SOILWM+((POR2(T)-(POR2(T)*WILTPOINT2(T)))*VGZDEX2(T))

       SOILWW=((MOS(T)%SoWET(1)*POR1(T))-(POR1(T)*WILTPOINT2(T)))*
     &  VGZDEX1(T)
       SOILWW=SOILWW+(((MOS(T)%SoWET(2)*POR2(T))-
     &        (POR2(T)*WILTPOINT2(T)))*VGZDEX2(T))


        ENDIF

        IF (NOAAIC.EQ.0) THEN
         MOS(T)%SoT=LDAS%MOS_IT
         DO L=1,3
           MOS(T)%SoWET(L)=LDAS%MOS_ISM
         ENDDO

!=== For Koster tile space runs set snow depths (liquid equivalent) to 100000 kg/m2
!=== for vegetation type 9 (land ice) so that it never melts since vegetation type 9 set
!=== to bare soil parameters (communication with Randy Koster -- Nov 2001)

         IF (LDAS%KOSTER.EQ.1.AND.TILE(T)%VEGT.EQ.9) MOS(T)%SNOW=100000
        
        ENDIF

       ENDDO

        IF (NOAAIC.EQ.1) THEN
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'USING DAGS IC DATA, NOT CARD!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'USING DAGS IC DATA, NOT CARD!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(79,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	ENDIF

       IF (NOAAIC.EQ.1) THEN
         WRITE(*,*)'MOSAIC Initialized SM/T=',LDAS%MOS_ISM,LDAS%MOS_IT
         WRITE(79,*)'MOSAIC Initialized SM/T=',LDAS%MOS_ISM,LDAS%MOS_IT
       ENDIF
      ENDIF

!=== Set starttime to ldas.crd Stime 
      IF(RW.EQ.1.AND.LDAS%STARTCODE.EQ.3)THEN 
       LDAS%YR=LDAS%SYR
       LDAS%MO=LDAS%SMO 
       LDAS%DA=LDAS%SDA
       LDAS%HR=LDAS%SHR
       LDAS%MN=LDAS%SMN
       LDAS%SS=LDAS%SSS
       CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,LDAS%YR,
     &                LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS) 
       WRITE(*,*)'Using ldas.crd start time ',LDAS%TIME
       WRITE(79,*)'Using ldas.crd start time ',LDAS%TIME
      ENDIF

!=== In the case of no previously existing restart, the user can
!=== initialize land surface state using model forcing
      IF(RW.EQ.1.AND.LDAS%MOS_IC.EQ.4)THEN !User restart with model forcing

!=== Uses specified start time set in ldas.crd similar to user defined restart
        LDAS%YR=LDAS%SYR
        LDAS%MO=LDAS%SMO
        LDAS%DA=LDAS%SDA
        LDAS%HR=LDAS%SHR
        LDAS%MN=LDAS%SMN
        LDAS%SS=LDAS%SSS
        CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,LDAS%YR,
     &                     LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS)
        WRITE(*,*)'Using ldas.crd model set init start time ',LDAS%TIME
        WRITE(79,*)'Using ldas.crd model set init start time ',LDAS%TIME

!===Determine the model forcing
        IF (LDAS%FGDAS.EQ.1) CALL GETGDAS(LDAS,GRID)
        IF (LDAS%FGEOS.EQ.1) CALL GETGEOS(LDAS,GRID)

!===Convert model interpolated forcing data from grid to tile space
        DO F=1,LDAS%NMIF
          DO T=1,LDAS%NCH
            TILE(T)%FORCING(F)=GRID(TILE(T)%COL,TILE(T)%ROW)%FORCING(F)
          ENDDO
        ENDDO

!===Initialize Mosaic parameters with model interpolated data

        DO T=1,LDAS%NCH
          MOS(T)%CT=TILE(T)%FORCING(1)          ! Model 2 m air temperature (K)
          MOS(T)%QA=TILE(T)%FORCING(2)          ! Model 2 m air specific humidity (kg/kg)
          MOS(T)%ICS=0.0                        ! Set to zero
          IF(LDAS%FGDAS .EQ. 1) THEN            ! For GDAS 
            MOS(T)%SNOW=TILE(T)%FORCING(15)     ! Model liq. equiv. snow depth (kg/m2)
            MOS(T)%SoT=TILE(T)%FORCING(1)       ! Model 10-200 cm soil layer temperature (K)
          ELSE                                                  ! For GEOS
            MOS(T)%SNOW=(TILE(T)%FORCING(12)/1000.0)*RHOICE     ! Model liq. equiv. snow depth (mm)/1000*density of ice
            MOS(T)%SoT=TILE(T)%FORCING(1)       ! Model 2 m air temperature (K)
          ENDIF

          VGZDEX(T,1)=  MOS(T)%SOILP(4)
          VGZDEX(T,2)=  MOS(T)%SOILP(5)
          VGZDEX(T,3)=  MOS(T)%SOILP(6)
          VGWMAX(T,1)=  MOS(T)%SOILP(8)
          VGWMAX(T,2)=  MOS(T)%SOILP(9)
          VGWMAX(T,3)=  MOS(T)%SOILP(10)

          POR(T,1)=VGWMAX(T,1)/(VGZDEX(T,1)*1000.)
          POR(T,2)=VGWMAX(T,2)/(VGZDEX(T,2)*1000.)
          POR(T,3)=VGWMAX(T,3)/(VGZDEX(T,3)*1000.)
          
!===Initialize soil wetness MOSAIC layers based on a weighted average between 2 model layers
!===MOS(T)%SOILP(4:6) are soil layer depths (layer 1 through 3) (m)          
         IF(LDAS%FGDAS.EQ.1) THEN                   ! For GDAS
          DO L=1,3
           IF (L.EQ.1) THEN
            WT1=(MOS(T)%SOILP(L+3) - 0.0) / MOS(T)%SOILP(L+3)
            WT2=(MOS(T)%SOILP(L+3) - 0.10) / MOS(T)%SOILP(L+3)
           ELSE
            IF ((MOS(T)%SOILP(L+3)-MOS(T)%SOILP(L+3-1)).EQ.0.0) THEN
             WT1=(0.10-MOS(T)%SOILP(L+3-1))/((MOS(T)%SOILP(L+3)+0.00001)
     &            - MOS(T)%SOILP(L+3-1))
             WT2=(MOS(T)%SOILP(L+3)-0.10)/((MOS(T)%SOILP(L+3)+0.00001)
     &            - MOS(T)%SOILP(L+3-1))
            ELSE
             WT1=(0.10-MOS(T)%SOILP(L+3-1))/(MOS(T)%SOILP(L+3)
     &            - MOS(T)%SOILP(L+3-1))
             WT2=(MOS(T)%SOILP(L+3)-0.10)/(MOS(T)%SOILP(L+3)
     &            - MOS(T)%SOILP(L+3-1))
            ENDIF
           ENDIF
            
           IF (WT1.LT.0.0) WT1=0.0
           IF (WT2.LT.0.0) WT2=0.0
           IF (WT1.GT.1.0) WT1=1.0
           IF (WT2.GT.1.0) WT2=1.0
            
!!!       MOS(T)%SoWET(L)=WT1*TILE(T)%FORCING(11)+WT2*TILE(T)%FORCING(12)
       MOS(T)%SoWET(L)=WT1*(TILE(T)%FORCING(11)/POR(T,L))
     &       + WT2*(TILE(T)%FORCING(12)/POR(T,L))

          ENDDO
        ELSE                                         ! For GEOS
!=== Assign GEOS percent capacity soil moisture  
!=== MOS(T)%SOILP(8:10) are moisture holding capacity of soil layers 1 through 3 (kg/m2)
          DO L=1,3
            MOS(T)%SoWET(L) = TILE(T)%FORCING(13) 
!!!         MOS(T)%SoWET(L) = TILE(T)%FORCING(13) *
!!!     &   (MOS(T)%SOILP(L+7) / (MOS(T)%SOILP(L+3)*1000))
          ENDDO
        ENDIF
       ENDDO

        WRITE(*,*)'MOSAIC Initialized with model data'
        WRITE(79,*)'MOSAIC Initialized with model data'

!=== For Koster tile space runs set snow depths (liquid equivalent) to 100000 kg/m2
!=== for vegetation type 9 (land ice) so that it never melts since vegetation type 9 set
!=== to bare soil parameters (communication with Randy Koster -- Nov 2001)

        IF (LDAS%KOSTER .EQ. 1) THEN
          DO T=1,LDAS%NCH
              IF (TILE(T)%VEGT .EQ. 9) MOS(T)%SNOW = 100000
          ENDDO
        ENDIF

!=== Write out model land surface state as initialized with model data
!=== for a future restart if requested
        OPEN(40,FILE=LDAS%MOS_MFILE,FORM='unformatted')
        
        WRITE(40) LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS,
     1      LDAS%VCLASS,LDAS%NC,LDAS%NR,LDAS%NCH  !Veg class, no tiles
        WRITE(40) TILE%COL       !Grid Col of Tile
        WRITE(40) TILE%ROW       !Grid Row of Tile
        WRITE(40) TILE%FGRD      !Fraction of Grid covered by tile
        WRITE(40) TILE%VEGT      !Vegetation Type of Tile
        WRITE(40) MOS%CT         !MOSAIC Canopy/Soil Temperature
        WRITE(40) MOS%QA         !MOSAIC Canopy Humidity
        WRITE(40) MOS%ICS        !MOSAIC Interception Canopy Storage
        WRITE(40) MOS%SNOW       !MOSAIC Snow Depth
        WRITE(40) MOS%SoT        !MOSAIC Deep Soil Temperaure
        IF(LDAS%RBIAS.EQ.1)THEN
         WRITE(40) GRID%MOSBETAK(1)  !MOSAIC Forecast BIAS
         WRITE(40) GRID%MOSDELT  !MOSAIC Analysis Increment
        ENDIF
        DO L=1,3
          DO T=1,LDAS%NCH
            TMPTILEN(T)=MOS(T)%SoWET(L)
          ENDDO
          WRITE(40) TMPTILEN      !MOSAIC Soil Wetness (3 layers)
        ENDDO

        CLOSE(40)

        DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF (LDAS%DOMAIN .EQ. 2) THEN
              DEALLOCATE(GRID(C,R)%GLBDATA1)
              DEALLOCATE(GRID(C,R)%GLBDATA2)
              ALLOCATE(GRID(C,R)%GLBDATA1(LDAS%NF))
              ALLOCATE(GRID(C,R)%GLBDATA2(LDAS%NF))
            ENDIF
            DEALLOCATE(GRID(C,R)%FORCING)
            ALLOCATE(GRID(C,R)%FORCING(LDAS%NF))
          ENDDO
        ENDDO

        DO T=1,LDAS%NCH
           DEALLOCATE(TILE(T)%FORCING)
           ALLOCATE(TILE(T)%FORCING(LDAS%NF))
        ENDDO

      ENDIF 

!=== Restart Writing (2 file are written - active and archive)
      IF(RW.EQ.2)THEN
       OPEN(40,FILE=LDAS%MOS_RFILE,FORM='unformatted') !Active archive restart

       WRITE(40) LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS,
     1           LDAS%VCLASS,LDAS%NC,LDAS%NR,LDAS%NCH  !Veg class, no tiles       
       WRITE(40) TILE%COL       !Grid Col of Tile   
       WRITE(40) TILE%ROW       !Grid Row of Tile
       WRITE(40) TILE%FGRD      !Fraction of Grid covered by tile
       WRITE(40) TILE%VEGT      !Vegetation Type of Tile
       WRITE(40) MOS%CT         !MOSAIC Canopy/Soil Temperature 
       WRITE(40) MOS%QA         !MOSAIC Canopy Humidity
       WRITE(40) MOS%ICS        !MOSAIC Interception Canopy Storage
       WRITE(40) MOS%SNOW       !MOSAIC Snow Depth
       WRITE(40) MOS%SoT        !MOSAIC Deep Soil Temperaure
       IF(LDAS%RBIAS.EQ.1)THEN
        WRITE(40) GRID%MOSBETAK(1) !MOSAIC Forecast BIAS
        WRITE(40) GRID%MOSDELT  !MOSAIC Analysis Increment
       ENDIF
       DO L=1,3
        DO T=1,LDAS%NCH
         TMPTILEN(T)=MOS(T)%SoWET(L)
        ENDDO
        WRITE(40) TMPTILEN      !MOSAIC Soil Wetness (3 layers)
       ENDDO

       CLOSE(40)
   
       WRITE(*,*)'MOSAIC Active Restart Written: ',LDAS%MOS_RFILE
       WRITE(79,*)'MOSAIC Active Restart Written: ',LDAS%MOS_RFILE

!=== Now write the archived restart file
 91    FORMAT(A4,I3,A5,I4,A1,I4,I2,I2,A7,I3,A1)
 92    FORMAT(80A1)
 93    FORMAT(A80)
 94    FORMAT(I4,I2,I2,I2)
 95    FORMAT(10A1)
 96    FORMAT(A40)
 97    FORMAT(A4,I3,A4) 
 98    FORMAT(A4,I3,A5,I4,A1,I4,I2,I2)
 99    FORMAT(A4,I3)
100    FORMAT(A9)
101    FORMAT(A7)
       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)
       WRITE(90,94,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
       READ(90,95,REC=1)FTIME
       DO I=1,10
        IF(FTIME(I).EQ.(' '))FTIME(I)='0'
       ENDDO
       WRITE(90,91,REC=1) '/EXP',LDAS%EXPCODE,'/MOS/',LDAS%YR,
     &  '/',LDAS%YR,LDAS%MO,
     &  LDAS%DA,'/LDAS.E',LDAS%EXPCODE,'.'
       READ(90,92,REC=1) (FNAME(I),I=1,36)
       DO I=1,36
        IF(FNAME(I).EQ.(' '))FNAME(I)='0'
       ENDDO

       WRITE(90,100,REC=1)'mkdir -p '
       READ(90,92,REC=1)(FMKDIR(I),I=1,9)
       WRITE(90,98,REC=1)'/EXP',LDAS%EXPCODE,'/MOS/',
     &  LDAS%YR,'/',LDAS%YR,LDAS%MO,LDAS%DA
       READ(90,92,REC=1) (FYRMODIR(I),I=1,25)
       DO I=1,25
        IF(FYRMODIR(I).EQ.(' '))FYRMODIR(I)='0'
       ENDDO

       WRITE(90,101,REC=1)'.MOSrst'
       READ(90,92,REC=1) (FSUBS(I),I=1,7)

       WRITE(90,96,REC=1) LDAS%ODIR                       
       READ(90,92,REC=1) (FBASE(I),I=1,80)
       C=0
       DO I=1,80
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
       ENDDO
       WRITE(90,92,REC=1)(FBASE(I),I=1,C),(FNAME(I),I=1,36),
     &                   (FTIME(I),I=1,10),(FSUBS(I),I=1,7 ) 
       READ(90,93,REC=1)FILEN
 
       WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &   (FYRMODIR(I),I=1,25)
       READ(90,93,REC=1)MKFYRMO
       CLOSE(90)
!== Archive File Name Generation Complete
!== Make the directories for the MOS restart file
c       IF(LDAS%GMT.LE.LDAS%WRITEINTM)THEN             
        CALL SYSTEM(MKFYRMO)
c       ENDIF

!== Archive File Name Generation Complete
       OPEN(40,FILE=FILEN,status='unknown',FORM='unformatted')
       WRITE(40) LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS,
     1           LDAS%VCLASS,LDAS%NC,LDAS%NR,LDAS%NCH  !Veg class, no tiles       
       WRITE(40) TILE%COL       !Grid Col of Tile   
       WRITE(40) TILE%ROW       !Grid Row of Tile
       WRITE(40) TILE%FGRD      !Fraction of Grid covered by tile
       WRITE(40) TILE%VEGT      !Vegetation Type of Tile
       WRITE(40) MOS%CT         !MOSAIC Canopy/Soil Temperature 
       WRITE(40) MOS%QA         !MOSAIC Canopy Humidity
       WRITE(40) MOS%ICS        !MOSAIC Interception Canopy Storage
       WRITE(40) MOS%SNOW       !MOSAIC Snow Depth
       WRITE(40) MOS%SoT        !MOSAIC Deep Soil Temperaure
       IF(LDAS%RBIAS.EQ.1)THEN
        WRITE(40) GRID%MOSBETAK(1) !MOSAIC Forecast BIAS
        WRITE(40) GRID%MOSDELT  !MOSAIC Analysis Increment
       ENDIF
       DO L=1,3
        DO T=1,LDAS%NCH
         TMPTILEN(T)=MOS(T)%SoWET(L)
        ENDDO
        WRITE(40) TMPTILEN      !MOSAIC Soil Wetness (3 layers)
       ENDDO

       CLOSE(40)
   
       WRITE(*,*)'MOSAIC Archive Restart Written: ',FILEN
       WRITE(79,*)'MOSAIC Archive Restart Written: ',FILEN

      ENDIF

!== Deallocate arrays if necessary
      IF (ALLOCATED(COL))   DEALLOCATE(COL)
      IF (ALLOCATED(ROW))   DEALLOCATE(ROW)
      IF (ALLOCATED(FGRD))  DEALLOCATE(FGRD)
      IF (ALLOCATED(VEGT))  DEALLOCATE(VEGT)
      IF (ALLOCATED(CT))    DEALLOCATE(CT)
      IF (ALLOCATED(QA))    DEALLOCATE(QA)
      IF (ALLOCATED(ICS))   DEALLOCATE(ICS)
      IF (ALLOCATED(SNOW))  DEALLOCATE(SNOW)
      IF (ALLOCATED(SoT))   DEALLOCATE(SoT)
      IF (ALLOCATED(SoWET)) DEALLOCATE(SoWET)
      IF (ALLOCATED(TMPTILEO)) DEALLOCATE(TMPTILEO)

      Return
      End






