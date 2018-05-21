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
! noahrst.f: 
!
! DESCRIPTION:
!  This program reads and writes restart files for NOAH.  This
!   includes all relevant water/energy storages, tile information,
!   and time information.  It also rectifies changes in the tile space.  
!
! *NOTE:
!   LDAS%STARTCODE = NOAH_IC:  
!     1=restart file, 2=realtime, 3=defined (cold start), 4=model init
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  05 Sep 2001: Brian Cosgrove; Modified code to use Dag Lohmann's NOAA
!               initial conditions if necessary.  This is controlled with
!               local variable NOAAIC.  Normally set to 0 in this subroutine
!               but set to 1 if want to use Dag's NOAA IC's.  Changed output
!               directory structure, and commented out if-then check so that
!               directory is always made.
!  28 Apr 2002: Kristi Arsenault; Added NOAH LSM into LDAS
!  28 May 2002: Kristi Arsenault; For STARTCODE=4, corrected SNEQV values  
!                and put SMC, SH2O, STC limit for GDAS and GEOS forcing.
!  04 Nov 2002: Kristi Arsenault; Added correction for liquid soil moisture
!                initialization and calls new subroutine, SH2O_INIT 
!  07 Jan 2003: Kristi Arsenault; Set POR4 for soil layers 
!  14 Jan 2003: Urszula Jambor; Added deallocation statements near end of routine
!               and changed pointer variables to allocatable.
!=========================================================================
! RESTART FILE FORMAT(fortran sequential binary):
!  YR,MO,DA,HR,MN,SS,VCLASS,NCH !Restart time,Veg class,no.tiles, no.soil lay 
!  TILE(NCH)%COL        !Grid Col of Tile   
!  TILE(NCH)%ROW        !Grid Row of Tile
!  TILE(NCH)%FGRD       !Fraction of Grid covered by tile
!  TILE(NCH)%VEGT       !Vegetation Type of Tile
!  NOAH(NCH)%STATES     !Model States in Tile Space
!=========================================================================

      SUBROUTINE NOAHRST(RW,LDAS,TILE,GRID,NOAH)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      USE grid_module      ! LDAS non-model-specific grid variables
      USE noah_module      ! NOAH tile variables
      IMPLICIT NONE      
      type (ldasdec) LDAS              
      type (tiledec) TILE(LDAS%NCH)
      type (griddec) GRID(LDAS%NC,LDAS%NR)   
      type (noahdec) NOAH(LDAS%NCH)

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

      REAL, allocatable :: T1(:)        ! NOAH Skin Temperature (K)
      REAL, allocatable :: CMC(:)       ! NOAH Canopy Water Content
      REAL, allocatable :: SNOWH(:)     ! NOAH Actual Snow Depth (m)
      REAL, allocatable :: SNEQV(:)     ! NOAH Water Equivalent Snow Depth (m)
      REAL, allocatable :: STC(:,:)     ! NOAH Soil Layer Temperature (K)
      REAL, allocatable :: SMC(:,:)     ! NOAH Soil Layer Total Moisture (liq+frzn) 
      REAL, allocatable :: SH2O(:,:)    ! NOAH Soil Layer Liquid Moisture
      REAL, allocatable :: CH(:)        ! NOAH Sfc Exchange Coef. for Heat/Moisture
      REAL, allocatable :: CM(:)        ! NOAH Sfc Exchange Coef. for Momentum

      REAL  SOILMOIST(LDAS%NC,LDAS%NR),SOILTEMP(LDAS%NC,LDAS%NR)
      REAL  POR1(LDAS%NCH),POR2(LDAS%NCH),POR3(LDAS%NCH),POR4(LDAS%NCH)
      REAL  SMCWLT(LDAS%NCH)        ! NOAH Wilting point for soil moisture (vol.)
      REAL  SLDPTH1(LDAS%NCH)       ! Soil layer thicknesses (m)
      REAL  SLDPTH2(LDAS%NCH), SLDPTH3(LDAS%NCH), SLDPTH4(LDAS%NCH)
      REAL  ZDPTH1(LDAS%NCH)       ! Soil layer depths (m)
      REAL  ZDPTH2(LDAS%NCH), ZDPTH3(LDAS%NCH), ZDPTH4(LDAS%NCH)

      REAL, allocatable :: TMPTILE1(:)  ! Temporary Transfer Array
      REAL, allocatable :: TMPTILE2(:)  ! Temporary Transfer Array
      REAL, allocatable :: TMPTILE3(:)  ! Temporary Transfer Array
      REAL :: TMPTILEN1(LDAS%NCH)   ! Temporary Transfer Array
      REAL :: TMPTILEN2(LDAS%NCH)   ! Temporary Transfer Array
      REAL :: TMPTILEN3(LDAS%NCH)   ! Temporary Transfer Array
      REAL :: TMPG(LDAS%NC,LDAS%NR) ! Temporary Transfer Array

      REAL :: G_T1(LDAS%NC,LDAS%NR)     ! NOAH Skin Temperature (K) 
      REAL :: G_CMC(LDAS%NC,LDAS%NR)    ! NOAH Canopy Water Content 
      REAL :: G_SNOWH(LDAS%NC,LDAS%NR)  ! NOAH Actual Snow Depth (m) 
      REAL :: G_SNEQV(LDAS%NC,LDAS%NR)  ! NOAH Water Equivalent Snow Depth (m)
      REAL :: G_STC(LDAS%NC,LDAS%NR,4)  ! NOAH Soil Layer Temperature (K) 
      REAL :: G_SMC(LDAS%NC,LDAS%NR,4)  ! NOAH Soil Total Moisture (liq+frzn) 
      REAL :: G_SH2O(LDAS%NC,LDAS%NR,4) ! NOAH Soil Layer Liquid Moisture 
      REAL :: G_CH(LDAS%NC,LDAS%NR)     ! NOAH Sfc Exchange Coef. for Heat/Moisture
      REAL :: G_CM(LDAS%NC,LDAS%NR)     ! NOAH Sfc Exchange Coef. for Momentum

      REAL :: RHOICE=917.0              ! Density of ice
      REAL :: WT1,WT2                   ! Weights for soil wetness initialization

      CHARACTER*80 FILEN,MKFYRMO
      CHARACTER*1  FNAME(80),FBASE(80),FSUBS(80),FMKDIR(80)
      CHARACTER*1  FTIME(10),FYRMODIR(80)
      INTEGER K
      INTEGER NOAAIC   ! 0=Use IC's from card file, 1=Use NOAA IC's from 
                       !                              Dag's data files 

      PARAMETER (NOAAIC=0)

!=== End Variable Definition =============================================

!=== Read Active Archive File ============================================

      IF((RW.EQ.1.AND.LDAS%NOAH_IC.EQ.1).OR.
     &   (RW.EQ.1.AND.LDAS%STARTCODE.EQ.1))THEN

        OPEN(40,FILE=LDAS%NOAH_RFILE,FORM='unformatted')

        READ(40) YR,MO,DA,HR,MN,SS,VCLASS,NC,NR,NCH  !Time, veg class, no. tiles
       
        ALLOCATE (COL(NCH),ROW(NCH),FGRD(NCH),VEGT(NCH))
        ALLOCATE (T1(NCH),CMC(NCH),SNOWH(NCH),SNEQV(NCH),STC(NCH,4))
        ALLOCATE (SMC(NCH,4),SH2O(NCH,4),CH(NCH),CM(NCH)) 
        ALLOCATE (TMPTILE1(NCH),TMPTILE2(NCH),TMPTILE3(NCH))

        READ(40) COL        !Grid Col of Tile   
        READ(40) ROW        !Grid Row of Tile
        READ(40) FGRD       !Fraction of Grid covered by tile
        READ(40) VEGT       !Vegetation Type of Tile
        READ(40) T1         !NOAH Skin Temperature (K) 
        READ(40) CMC        !NOAH Canopy Water Content 
        READ(40) SNOWH      !NOAH Actual Snow Depth (m) 
        READ(40) SNEQV      !NOAH Water Equivalent Snow Depth (m) 
        DO L=1,4
          READ(40) TMPTILE1  !NOAH Soil Layer Temp (4 layers)
          DO T=1,NCH
            STC(T,L)=TMPTILE1(T)
          ENDDO
        ENDDO
        DO L=1,4
          READ(40) TMPTILE2  !NOAH Total soil moist. (4 layers)
          DO T=1,NCH
            SMC(T,L)=TMPTILE2(T)
          ENDDO
        ENDDO
        DO L=1,4
          READ(40) TMPTILE3  !NOAH Liquid-only soil moist. (4 layers)
          DO T=1,NCH
            SH2O(T,L)=TMPTILE3(T)
          ENDDO
        ENDDO
        READ(40) CH          !NOAH Sfc Exchange Coef. for Heat/Moisture
        READ(40) CM          !NOAH Sfc Exchange Coef. for Momentum

        CLOSE(40)

        WRITE(*,*)'NOAH Restart File Read: ',LDAS%NOAH_RFILE
        WRITE(*,*)'NOAH Restart File Time: ',YR,MO,DA,HR,MN,SS

        WRITE(79,*)'NOAH Restart File Read: ',LDAS%NOAH_RFILE
        WRITE(79,*)'NOAH Restart File Time: ',YR,MO,DA,HR,MN,SS

!=== Establish Model Restart Time
        IF(LDAS%STARTCODE.EQ.1)THEN
          LDAS%YR=YR
          LDAS%MO=MO 
          LDAS%DA=DA
          LDAS%HR=HR
          LDAS%MN=MN
          LDAS%SS=SS
          CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,YR,MO,DA,
     &                  HR,MN,SS) 
          LDAS%NOAH_STIME = LDAS%TIME
          WRITE(*,*)'NOAH Restart File Time Used: ',LDAS%NOAH_RFILE
        ENDIF

!=== Rectify Restart Tile Space to LDAS Tile Space =====================
       IF(LDAS%NOAH_IC.EQ.1)THEN

!   Check for Vegetation Class Conflict 
         IF(VCLASS.NE.LDAS%VCLASS)THEN
           WRITE(*,*)LDAS%NOAH_RFILE,' Vegetation class conflict'
           WRITE(79,*)LDAS%NOAH_RFILE,' Vegetation class conflict'
           STOP
         ENDIF

!   Check for Grid Space Conflict 
         IF(NC.NE.LDAS%NC.OR.NR.NE.LDAS%NR)THEN
           WRITE(*,*)LDAS%NOAH_RFILE,'Grid space mismatch - NOAH HALTED'
           STOP
         ENDIF

!-- Transfer Restart tile space to LDAS tile space

         IF(NCH.NE.LDAS%NCH)THEN

           WRITE(*,*)'Restart Tile Space Mismatch-Transfer in Progress'
           WRITE(79,*)'Restart Tile Space Mismatch-Transfer in Progress'

!   Start by finding grid averages
           CALL T2GR(T1,G_T1,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
           CALL T2GR(CMC,G_CMC,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
           CALL T2GR(SNOWH,G_SNOWH,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
           CALL T2GR(SNEQV,G_SNEQV,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
           DO L=1,4
             CALL T2GR(STC(:,L),G_STC(:,:,L),
     1               LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
           END DO
           DO L=1,4
             CALL T2GR(SMC(:,L),G_SMC(:,:,L),
     1               LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
           ENDDO
           DO L=1,4
             CALL T2GR(SH2O(:,L),G_SH2O(:,:,L),
     1               LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
           ENDDO
           CALL T2GR(CH,G_CH,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)
           CALL T2GR(CM,G_CM,LDAS%NC,LDAS%NR,NCH,FGRD,COL,ROW)

!   Perform state transfer
           C=0
           DO 555 T=1,LDAS%NCH 
             IF(AMOD(FLOAT(T),10000.0).EQ.0.0)
     1         WRITE(*,23)'  Transferred ',
     2          100.0*FLOAT(T)/FLOAT(LDAS%NCH),' Percent of Tiles'

             IF(AMOD(FLOAT(T),10000.0).EQ.0.0)
     1         WRITE(79,23)'  Transferred ',
     2          100.0*FLOAT(T)/FLOAT(LDAS%NCH),' Percent of Tiles'

 23           FORMAT(A14,F5.2,A17)
              FOUND=0
              DO N=1,NCH
               IF(TILE(T)%VEGT.EQ.VEGT(N).AND.
     1           TILE(T)%COL .EQ. COL(N).AND.
     2           TILE(T)%ROW .EQ. ROW(N))THEN
                  NOAH(T)%T1=T1(N)
                  NOAH(T)%CMC=CMC(N)
                  NOAH(T)%SNOWH=SNOWH(N)
                  NOAH(T)%SNEQV=SNEQV(N)
                  DO L=1,4
                    NOAH(T)%STC(L)=STC(N,L)
                  ENDDO
                  DO L=1,4
                    NOAH(T)%SMC(L)=SMC(N,L)
                  ENDDO
                  DO L=1,4
                    NOAH(T)%SH2O(L)=SH2O(N,L)
                  ENDDO
                  NOAH(T)%CH=CH(N)
                  NOAH(T)%CM=CM(N)

                FOUND=1
                GOTO 555 
               ENDIF
              ENDDO

             IF(FOUND.EQ.0)THEN        
               NOAH(T)%T1=G_T1(TILE(T)%COL,TILE(T)%ROW)
               NOAH(T)%CMC=G_CMC(TILE(T)%COL,TILE(T)%ROW)
               NOAH(T)%SNOWH=G_SNOWH(TILE(T)%COL,TILE(T)%ROW)
               NOAH(T)%SNEQV=G_SNEQV(TILE(T)%COL,TILE(T)%ROW)
              DO L=1,4
                NOAH(T)%STC(L)=G_STC(TILE(T)%COL,TILE(T)%ROW,L)
              ENDDO
              DO L=1,4
                NOAH(T)%SMC(L)=G_SMC(TILE(T)%COL,TILE(T)%ROW,L)
              ENDDO
              DO L=1,4
                NOAH(T)%SH2O(L)=G_SH2O(TILE(T)%COL,TILE(T)%ROW,L)
              ENDDO
               NOAH(T)%CH=G_CH(TILE(T)%COL,TILE(T)%ROW)
               NOAH(T)%CM=G_CM(TILE(T)%COL,TILE(T)%ROW)
              C=0
             ENDIF

 555       CONTINUE

           WRITE(*,*)'Tile Space Transfer Complete'
           WRITE(*,*)'NOAH Restart NCH:',NCH,'Current NCH:',LDAS%NCH
           WRITE(*,*) C, ' Tiles not found in old NOAH restart'        
           WRITE(*,*)  

           WRITE(79,*)'Tile Space Transfer Complete'
           WRITE(79,*)'NOAH Restart NCH:',NCH,'Current NCH:',LDAS%NCH
           WRITE(79,*) C, ' Tiles not found in old NOAH restart'
           WRITE(79,*)

         ELSE  !The number of tiles is a match

           DO T=1,LDAS%NCH
             NOAH(T)%T1=T1(T)
             NOAH(T)%CMC=CMC(T)
             NOAH(T)%SNOWH=SNOWH(T)
             NOAH(T)%SNEQV=SNEQV(T)
             DO L=1,4
              NOAH(T)%STC(L)=STC(T,L)
             ENDDO
             DO L=1,4
              NOAH(T)%SMC(L)=SMC(T,L)
             ENDDO
             DO L=1,4
              NOAH(T)%SH2O(L)=SH2O(T,L)
             ENDDO
              NOAH(T)%CH=CH(T)
              NOAH(T)%CM=CM(T)
           ENDDO	

         ENDIF !End for restart tile space transfer  
       ENDIF   !End of matching restart tile space to LDAS  

      ENDIF    !RW option 1

!=== In the case of NO previously existing restart, the user can 
!===  specify the restart condifions in the card file.

      IF(RW.EQ.1.AND.LDAS%NOAH_IC.EQ.3) THEN   !COLD-START NOAH 

	IF (NOAAIC.EQ.1) THEN      !Using Dag's data files

!---  READ IN INITIAL SOIL TEMP AND SOIL WETNESS FROM NOAA's 
!---   INITIAL CONDITIONS (FROM DAG).
!      Soil Moisture file
          open (unit=22,file=
     &      './BCS/N0.125/1996093012.moist_scale.bin',
     &     form='unformatted')
           read(22) ((soilmoist(i,j),i=1,464),j=1,224)
          close (22)
!      Soil Temperaure file
          open (unit=22,file=
     &      './BCS/N0.125/1996093012.temp2.bin',
     &     form='unformatted')
           read(22)((soiltemp(i,j),i=1,464),j=1,224)
          close(22)

          DO T=1,LDAS%NCH    ! Tile loop for Dag's data files

             NOAH(T)%T1   =LDAS%NOAH_IT
             NOAH(T)%CMC  =0.0004
             NOAH(T)%SNOWH=0.0
             NOAH(T)%SNEQV=0.0
c             NOAH(T)%CH   =0.0001
c             NOAH(T)%CM   =0.0001
             NOAH(T)%CH   =0.0150022404
             NOAH(T)%CM   =0.0205970779

             DO L=1,4
               NOAH(T)%STC(L)=SOILTEMP(TILE(T)%COL,TILE(T)%ROW)
             ENDDO

! ***  The soil moisture values from NOAA dataset range from 
!       wilting point (0) to fully saturated or porosity (1)  

             SLDPTH1(T) = 0.1     
             SLDPTH2(T) = 0.3
             SLDPTH3(T) = 0.6
             SLDPTH4(T) = 1.0

             SMCWLT(T) = NOAH(T)%SOILP(7)        

             POR1(T) = NOAH(T)%SOILP(1)     ! For now SMCMAX parameter is used 
             POR2(T) = NOAH(T)%SOILP(1)
             POR3(T) = NOAH(T)%SOILP(1)
             POR4(T) = NOAH(T)%SOILP(1)

             NOAH(T)%SMC(1)= (soilmoist(TILE(T)%COL,TILE(T)%ROW)
     &        *POR1(T))-(soilmoist(TILE(T)%COL,TILE(T)%ROW)*SMCWLT(T))
     &        +(SMCWLT(T))
             NOAH(T)%SMC(2)= (soilmoist(TILE(T)%COL,TILE(T)%ROW)
     &        *POR2(T))-(soilmoist(TILE(T)%COL,TILE(T)%ROW)*SMCWLT(T))
     &        +(SMCWLT(T))
             NOAH(T)%SMC(3)= (soilmoist(TILE(T)%COL,TILE(T)%ROW)
     &        *POR3(T))-(soilmoist(TILE(T)%COL,TILE(T)%ROW)*SMCWLT(T))
     &        +(SMCWLT(T))
             NOAH(T)%SMC(4)= (soilmoist(TILE(T)%COL,TILE(T)%ROW)
     &        *POR4(T))-(soilmoist(TILE(T)%COL,TILE(T)%ROW)*SMCWLT(T))
     &        +(SMCWLT(T))

             NOAH(T)%SH2O(1)= (soilmoist(TILE(T)%COL,TILE(T)%ROW)
     &        *POR1(T))-(soilmoist(TILE(T)%COL,TILE(T)%ROW)*SMCWLT(T))
     &        +(SMCWLT(T))
             NOAH(T)%SH2O(2)= (soilmoist(TILE(T)%COL,TILE(T)%ROW)
     &        *POR2(T))-(soilmoist(TILE(T)%COL,TILE(T)%ROW)*SMCWLT(T))
     &        +(SMCWLT(T))
             NOAH(T)%SH2O(3)= (soilmoist(TILE(T)%COL,TILE(T)%ROW)
     &        *POR3(T))-(soilmoist(TILE(T)%COL,TILE(T)%ROW)*SMCWLT(T))
     &        +(SMCWLT(T))
             NOAH(T)%SH2O(4)= (soilmoist(TILE(T)%COL,TILE(T)%ROW)
     &        *POR4(T))-(soilmoist(TILE(T)%COL,TILE(T)%ROW)*SMCWLT(T))
     &        +(SMCWLT(T))

          ENDDO

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

           WRITE(*,*)'NOAH Init. SM/T=',LDAS%NOAH_ISM,LDAS%NOAH_IT
           WRITE(79,*)'NOAH Init. SM/T=',LDAS%NOAH_ISM,LDAS%NOAH_IT

        ENDIF   !End use of NOAA IC's from Dag's data files

        IF (NOAAIC.EQ.0) THEN      ! Use IC's from card file

          print *, 'USING ICs FROM CARD FILE TO INIT NOAH LSM!!' 

          DO T=1,LDAS%NCH

            NOAH(T)%T1=LDAS%NOAH_IT
            NOAH(T)%T1=280.0

            NOAH(T)%CMC=0.0004
            NOAH(T)%SNOWH=0.0
            NOAH(T)%SNEQV=0.0
c            NOAH(T)%CH=0.0001
c            NOAH(T)%CM=0.0001
            NOAH(T)%CH=0.0150022404
            NOAH(T)%CM=0.0205970779

            DO L=1,4
              NOAH(T)%STC(L)=LDAS%NOAH_IT
            ENDDO
!             NOAH(T)%STC(1)=266.0
!             NOAH(T)%STC(2)=273.0
!             NOAH(T)%STC(3)=276.0
!             NOAH(T)%STC(4)=280.0

             NOAH(T)%SMC(1)=0.3252287
             NOAH(T)%SMC(2)=0.3194746
             NOAH(T)%SMC(3)=0.3172167
             NOAH(T)%SMC(4)=0.3078052

             NOAH(T)%SH2O(1)=0.1660042
             NOAH(T)%SH2O(2)=0.2828006
             NOAH(T)%SH2O(3)=0.3172163
             NOAH(T)%SH2O(4)=0.3078025

c            DO L=1,4
c              NOAH(T)%SMC(L)=LDAS%NOAH_ISM
c            ENDDO
c            DO L=1,4
c              NOAH(T)%SH2O(L)=LDAS%NOAH_ISM
c            ENDDO

          ENDDO  ! End of Tile loop for card options
        ENDIF    ! End use of IC's from card file

      ENDIF    ! COLD-START OPTION FOR NOAH

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

!=== In the case of NO previously existing restart, the user can
!=== initialize land surface state using model forcing.

      IF(RW.EQ.1.AND.LDAS%NOAH_IC.EQ.4) THEN  !Use restart with model forcing

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

!===Initialize NOAH parameters with global model-interpolated data

         DO T=1,LDAS%NCH

           NOAH(T)%T1 =TILE(T)%FORCING(1)       ! Model 2 m air temperature (K)
           NOAH(T)%CMC=0.0001       
           NOAH(T)%CH =0.0150022404
           NOAH(T)%CM =0.0205970779

           IF(LDAS%FGDAS .EQ. 1) THEN   ! For GDAS 
             NOAH(T)%SNOWH=(TILE(T)%FORCING(15)/1000.0)*RHOICE  ! Model liq. equiv. snow depth (kg/m2)
             NOAH(T)%SNEQV=TILE(T)%FORCING(15)/1000.0           ! Model liq. equiv. snow depth (kg/m2)

           ELSE     ! For GEOS
             NOAH(T)%SNOWH=(TILE(T)%FORCING(13)/1000.0)*RHOICE  ! Model liq. equiv. snow depth (mm)
                                                                !  /1000*density of ice
             NOAH(T)%SNEQV=TILE(T)%FORCING(13)/1000.0           ! Model liq. equiv. snow depth (mm)
             DO L=1,4
               NOAH(T)%STC(L)=TILE(T)%FORCING(1)     ! Model 2 m air temperature (K)
               IF(NOAH(T)%STC(L).EQ.0.0)  NOAH(T)%STC(L)=285.0
             ENDDO
           ENDIF
          
!=== GDAS: Initialize soil moisture and temperature in 4 NOAH layers based   
!===       values from GDAS soil layers at depths: 0-0.1 m and 0.1-2.0 m.

           IF(LDAS%FGDAS.EQ.1) THEN    ! For GDAS

            DO L=1,4

             IF (L.EQ.1) THEN          ! LAYER 1

               NOAH(T)%SMC(L)=TILE(T)%FORCING(11)  ! 0-10 cm soil moisture layer
                IF(NOAH(T)%SMC(L).EQ.0.0)  NOAH(T)%SMC(L)=0.15
               NOAH(T)%STC(L)=TILE(T)%FORCING(13)  ! 0-10 cm soil layer temperature (K)
                IF(NOAH(T)%STC(L).EQ.0.0)  NOAH(T)%STC(L)=276.0

               IF (NOAH(T)%STC(L) .GE. 273.149) THEN   ! Soil Temp above or at freezing point
                  NOAH(T)%SH2O(L)=NOAH(T)%SMC(L)
               ELSE                                    ! Soil Temp below freezing point

                  CALL SH2OINIT(NOAH(T)%SMC(L),NOAH(T)%STC(L),
     &              NOAH(T)%SOILP(1),NOAH(T)%SOILP(2),NOAH(T)%SOILP(4),
     &              NOAH(T)%SH2O(L))   
               ENDIF

             ELSE                      ! LAYERS 2-4

               NOAH(T)%SMC(L)=TILE(T)%FORCING(12)      ! 10-200 cm soil moisture  
                IF(NOAH(T)%SMC(L).EQ.0.0)  NOAH(T)%SMC(L)=0.15
               NOAH(T)%STC(L)=TILE(T)%FORCING(14)      ! 10-200 cm soil layer temp (K)  
                IF(NOAH(T)%STC(L).EQ.0.0)  NOAH(T)%STC(L)=276.0

               IF (NOAH(T)%STC(L) .GE. 273.149) THEN   ! Soil Temp above or at freezing point
                  NOAH(T)%SH2O(L)=NOAH(T)%SMC(L)
               ELSE                                    ! Soil Temp below freezing point

                  CALL SH2OINIT(NOAH(T)%SMC(L),NOAH(T)%STC(L),
     &              NOAH(T)%SOILP(1),NOAH(T)%SOILP(2),NOAH(T)%SOILP(4),
     &              NOAH(T)%SH2O(L)) 
               ENDIF

             ENDIF
            ENDDO  ! End layer loop for GDAS forcing

!=== Need to convert GEOS percent field capacity soil moisture to NOAH volumetric  
!===  soil moisture for layers 1 through 4. 

          ELSE     ! For GEOS
            DO L=1,4
              NOAH(T)%SMC(L) = TILE(T)%FORCING(12) * NOAH(T)%SOILP(1)
              NOAH(T)%SH2O(L) = TILE(T)%FORCING(12) * NOAH(T)%SOILP(1)
               IF(NOAH(T)%SMC(L).EQ.0.0)  NOAH(T)%SMC(L)=0.15
               IF(NOAH(T)%SH2O(L).EQ.0.0)  NOAH(T)%SH2O(L)=0.15
            ENDDO

          ENDIF  ! End GDAS/GEOS option
         ENDDO   ! End initializing NOAH parameters with model interpolated data

         WRITE(*,*)'NOAH Initialized with model data'
         WRITE(79,*)'NOAH Initialized with model data'

!==============================================================================
!=== Write out model land surface state as initialized with model data
!=== for a future restart if requested

        OPEN(40,FILE=LDAS%NOAH_MFILE,FORM='unformatted')
        
        WRITE(40) LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS,
     1      LDAS%VCLASS,LDAS%NC,LDAS%NR,LDAS%NCH  !Veg class, no tiles
        WRITE(40) TILE%COL     !Grid Col of Tile
        WRITE(40) TILE%ROW     !Grid Row of Tile
        WRITE(40) TILE%FGRD    !Fraction of Grid covered by tile
        WRITE(40) TILE%VEGT    !Vegetation Type of Tile
        WRITE(40) NOAH%T1      !NOAH Skin Temperature (K) 
        WRITE(40) NOAH%CMC     !NOAH Canopy Water Content
        WRITE(40) NOAH%SNOWH   !NOAH Actual Snow Depth
        WRITE(40) NOAH%SNEQV   !NOAH Water Equivalent Snow Depth
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN1(T)=NOAH(T)%STC(L)
          ENDDO
          WRITE(40) TMPTILEN1  !NOAH Soil Temperature (4 layers)
        ENDDO
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN2(T)=NOAH(T)%SMC(L)
          ENDDO
          WRITE(40) TMPTILEN2  !NOAH Total Soil Moist. (4 layers)
        ENDDO
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN3(T)=NOAH(T)%SH2O(L)
          ENDDO
          WRITE(40) TMPTILEN3  !NOAH Liquid Soil Moist. (4 layers)
        ENDDO
        WRITE(40) NOAH%CH      !NOAH Heat/Moisture Sfc Exchange Coef. 
        WRITE(40) NOAH%CM      !NOAH Momentum Sfc Exchange Coef.  

        CLOSE(40)

        DO C=1,LDAS%NC
          DO R=1,LDAS%NR
            IF (LDAS%DOMAIN .EQ. 2) THEN         ! GLDAS
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

!=== Restart Writing (2 files are written = active and archive)

      IF(RW.EQ.2)THEN
        OPEN(40,FILE=LDAS%NOAH_RFILE,FORM='unformatted') !Active archive restart

        WRITE(40) LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS,
     1            LDAS%VCLASS,LDAS%NC,LDAS%NR,LDAS%NCH  !Veg class, no tiles       
        WRITE(40) TILE%COL       !Grid Col of Tile   
        WRITE(40) TILE%ROW       !Grid Row of Tile
        WRITE(40) TILE%FGRD      !Fraction of Grid covered by tile
        WRITE(40) TILE%VEGT      !Vegetation Type of Tile
        WRITE(40) NOAH%T1        !NOAH Skin Temperature (K)
        WRITE(40) NOAH%CMC       !NOAH Canopy Water Content
        WRITE(40) NOAH%SNOWH     !NOAH Actual Snow Depth
        WRITE(40) NOAH%SNEQV     !NOAH Water Equivalent Snow Depth
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN1(T)=NOAH(T)%STC(L)
          ENDDO
          WRITE(40) TMPTILEN1  !NOAH Soil Temperature (4 layers)
        ENDDO
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN2(T)=NOAH(T)%SMC(L)
          ENDDO
          WRITE(40) TMPTILEN2  !NOAH Total Soil Moist. (4 layers)
        ENDDO
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN3(T)=NOAH(T)%SH2O(L)
          ENDDO
          WRITE(40) TMPTILEN3  !NOAH Liquid Soil Moist. (4 layers)
        ENDDO
        WRITE(40) NOAH%CH      !NOAH Heat/Moisture Sfc Exchange Coef.
        WRITE(40) NOAH%CM      !NOAH Momentum Sfc Exchange Coef.

        CLOSE(40)
   
        WRITE(*,*)'NOAH Active Restart Written: ',LDAS%NOAH_RFILE
        WRITE(79,*)'NOAH Active Restart Written: ',LDAS%NOAH_RFILE

!=== Now write the archived restart file
 91     FORMAT(A4,I3,A6,I4,A1,I4,I2,I2,A7,I3,A1)
 92     FORMAT(80A1)
 93     FORMAT(A80)
 94     FORMAT(I4,I2,I2,I2)
 95     FORMAT(10A1)
 96     FORMAT(A40)
 97     FORMAT(A4,I3,A4) 
 98     FORMAT(A4,I3,A6,I4,A1,I4,I2,I2)
 99     FORMAT(A4,I3)
100     FORMAT(A9)
101     FORMAT(A8)

       OPEN(90,FILE='temp',FORM='FORMATTED',ACCESS='DIRECT',RECL=80)

        WRITE(90,94,REC=1)LDAS%YR,LDAS%MO,LDAS%DA,LDAS%HR
        READ(90,95,REC=1)FTIME
        DO I=1,10
          IF(FTIME(I).EQ.(' '))FTIME(I)='0'
        ENDDO
        WRITE(90,91,REC=1) '/EXP',LDAS%EXPCODE,'/NOAH/',LDAS%YR,
     &    '/',LDAS%YR,LDAS%MO,
     &    LDAS%DA,'/LDAS.E',LDAS%EXPCODE,'.'
        READ(90,92,REC=1) (FNAME(I),I=1,37)
        DO I=1,73
          IF(FNAME(I).EQ.(' '))FNAME(I)='0'
        ENDDO

        WRITE(90,100,REC=1)'mkdir -p '
        READ(90,92,REC=1)(FMKDIR(I),I=1,9)
        WRITE(90,98,REC=1)'/EXP',LDAS%EXPCODE,'/NOAH/',
     &    LDAS%YR,'/',LDAS%YR,LDAS%MO,LDAS%DA
        READ(90,92,REC=1) (FYRMODIR(I),I=1,26)
        DO I=1,26
          IF(FYRMODIR(I).EQ.(' '))FYRMODIR(I)='0'
        ENDDO

        WRITE(90,101,REC=1)'.NOAHrst'
        READ(90,92,REC=1) (FSUBS(I),I=1,8)

        WRITE(90,96,REC=1) LDAS%ODIR                       
        READ(90,92,REC=1) (FBASE(I),I=1,80)
        C=0
        DO I=1,80
          IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
        ENDDO
        WRITE(90,92,REC=1)(FBASE(I),I=1,C),(FNAME(I),I=1,37),
     &                   (FTIME(I),I=1,10),(FSUBS(I),I=1,8) 
        READ(90,93,REC=1)FILEN
 
        WRITE(90,92,REC=1)(FMKDIR(I),I=1,9),(FBASE(I),I=1,C),
     &    (FYRMODIR(I),I=1,26)
        READ(90,93,REC=1)MKFYRMO

       CLOSE(90)

!== Archive File Name Generation Complete
!== Make the directories for the NOAH restart file
c       IF(LDAS%GMT.LE.LDAS%WRITEINTN)THEN             
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
        WRITE(40) NOAH%T1        !NOAH Skin Temperature (K)
        WRITE(40) NOAH%CMC       !NOAH Canopy Water Content
        WRITE(40) NOAH%SNOWH     !NOAH Actual Snow Depth
        WRITE(40) NOAH%SNEQV     !NOAH Water Equivalent Snow Depth
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN1(T)=NOAH(T)%STC(L)
          ENDDO
          WRITE(40) TMPTILEN1  !NOAH Soil Temperature (4 layers)
        ENDDO
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN2(T)=NOAH(T)%SMC(L)
          ENDDO
          WRITE(40) TMPTILEN2  !NOAH Total Soil Moist. (4 layers)
        ENDDO
        DO L=1,4
          DO T=1,LDAS%NCH
            TMPTILEN3(T)=NOAH(T)%SH2O(L)
          ENDDO
          WRITE(40) TMPTILEN3  !NOAH Liquid Soil Moist. (4 layers)
        ENDDO
        WRITE(40) NOAH%CH        !NOAH Heat/Moisture Sfc Exchange Coef.
        WRITE(40) NOAH%CM        !NOAH Momentum Sfc Exchange Coef.

       CLOSE(40)
   
        WRITE(*,*)'NOAH Archive Restart Written: ',FILEN
        WRITE(79,*)'NOAH Archive Restart Written: ',FILEN

      ENDIF  !End restart writing (active and archive files)

!== Deallocate arrays if necessary
      IF (ALLOCATED(COL))   DEALLOCATE(COL)
      IF (ALLOCATED(ROW))   DEALLOCATE(ROW)
      IF (ALLOCATED(FGRD))  DEALLOCATE(FGRD)
      IF (ALLOCATED(VEGT))  DEALLOCATE(VEGT)
      IF (ALLOCATED(T1))    DEALLOCATE(T1)
      IF (ALLOCATED(CMC))   DEALLOCATE(CMC)
      IF (ALLOCATED(SNOWH)) DEALLOCATE(SNOWH)
      IF (ALLOCATED(SNEQV)) DEALLOCATE(SNEQV)
      IF (ALLOCATED(STC))   DEALLOCATE(STC)
      IF (ALLOCATED(SMC))   DEALLOCATE(SMC)
      IF (ALLOCATED(SH2O))  DEALLOCATE(SH2O)
      IF (ALLOCATED(CH))    DEALLOCATE(CH)
      IF (ALLOCATED(CM))    DEALLOCATE(CM)
      IF (ALLOCATED(TMPTILE1)) DEALLOCATE(TMPTILE1)
      IF (ALLOCATED(TMPTILE2)) DEALLOCATE(TMPTILE2)
      IF (ALLOCATED(TMPTILE3)) DEALLOCATE(TMPTILE3)

      Return
      End

