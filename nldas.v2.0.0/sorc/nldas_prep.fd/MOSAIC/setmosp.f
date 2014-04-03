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
! SETMOSP.f: 
!
! DESCRIPTION:
!  This subroutine retrieves MOSAIC parameters - Significant F90 revisions
!   below this subroutine will be requiret in the future.  
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Code
!  31 Jul 2001: Matt Rodell; Updated for new soil parameter definition.
!  14 Feb 2002: Jon Gottschalck; Added allocated space for AVHRR LAI/DSAI     
!  07 Mar 2002: Brian Cosgrove; Corrected declaration of TEX1 var from real to int
!  14 Jan 2003: Urszula Jambor; Added conditional to check if need exists 
!               to allocate for AVHRR LAI/DSAI variables.
!  05 Aug 2003: Brian Cosgrove; Added new soil scheme choice...option 5 now
!               uses the same parameters as option 4 (Cosby/Rawls)
!               but uses vegetation-dependent soil layer thicknesses 
!=========================================================================

      SUBROUTINE SETMOSP(LDAS,TILE,MOS)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      USE mos_module       ! Mosaic tile variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (tiledec) TILE(LDAS%NCH)
      type (mosdec)  MOS(LDAS%NCH)

!+++ Local Parameters for new soil definition ++++++++++++++++++++++++++++
!+++  Layer thicknesss correspond to soil property map files - do not change!
C       Thicknesses of soil layers for any simulations
C       These are compatible with the FIRST set of Matt Rodell's 
C       soil files, i.e. sand_nldas3.5.bfsa

      REAL, PARAMETER :: D1 = 0.02	! thickness of soil layer 1, m
      REAL, PARAMETER :: D2 = 1.48	! thickness of soil layer 2, m
      REAL, PARAMETER :: D3 = 2.00	! thickness of soil layer 3, m

C	thicknesss of soil layers (semi-official) for NLDAS simulations
C	These are compatible with either Yun Duans's soil maps
C       or the SECOND set of Matt Rodell's soil files, i.e., sand_nldas2m.bfsa
      REAL, PARAMETER :: DNLDAS1 = 0.1    ! thickness of soil layer 1, m
      REAL, PARAMETER :: DNLDAS2 = 0.3    ! thickness of soil layer 2, m
      REAL, PARAMETER :: DNLDAS3 = 1.6    ! thickness of soil layer 3, m


!+++  Minimum values of sin(theta) based on vegetation type.
      REAL, PARAMETER :: S8 = 0.57787	! Closed shrubland
      REAL, PARAMETER :: S9 = 0.95504	! Open shrubland
      REAL, PARAMETER :: S12 = 0.1736	! Bare soil
      REAL, PARAMETER :: S0 = 0.05	! All others 
!+++  Maximum allowable porosity.
      REAL, PARAMETER :: PORMAX=0.70

!=== Local Variables =====================================================
      INTEGER :: C,R,T,I,J,K,JJ,II     ! Loop counters
      INTEGER :: SOILTEXT(11,LDAS%NC,LDAS%NR)
      INTEGER :: TEX1(LDAS%NC,LDAS%NR)
      REAL :: VEGP(LDAS%NCH,LDAS%MOS_NVEGP)   
      REAL :: VEGMP(LDAS%NCH,LDAS%MOS_NMVEGP,12)   
      REAL :: SOILP(LDAS%NCH,LDAS%MOS_NSOILP)   
      REAL :: VALUE(LDAS%NT,LDAS%MOS_NVEGP)
      REAL :: VALUEMON(LDAS%NT,LDAS%MOS_NMVEGP,12)
      REAL :: BASICSET(LDAS%NT,LDAS%MOS_NSOILP)
      REAL :: SAND1(LDAS%NC,LDAS%NR)
      REAL :: SILT1(LDAS%NC,LDAS%NR)
      REAL :: CLAY1(LDAS%NC,LDAS%NR)
      REAL :: POR1(LDAS%NC,LDAS%NR)
      REAL :: POR2(LDAS%NC,LDAS%NR)
      REAL :: POR3(LDAS%NC,LDAS%NR)
      REAL :: POR1A,POR2A,POR3A
      REAL :: SLOPE(LDAS%NC,LDAS%NR)
      REAL :: KSAT1(LDAS%NC,LDAS%NR)
      REAL :: PSI1(LDAS%NC,LDAS%NR)
      REAL :: B1(LDAS%NC,LDAS%NR)
      REAL :: HYCON(12)
      REAL :: SOILPOT(12)
      REAL :: POROS(12)
      REAL :: BPARAM(12)     

      REAL :: POR1TILE(LDAS%NCH)
      REAL :: POR2TILE(LDAS%NCH)
      REAL :: POR3TILE(LDAS%NCH)
!=== End Variable Definition =============================================

!=== Allocate Memory for Mosaic Variables ================================    
      DO T=1,LDAS%NCH
       ALLOCATE (MOS(T)%VEGP(LDAS%MOS_NVEGP))
       ALLOCATE (MOS(T)%VEGMP(LDAS%MOS_NMVEGP,12))
       ALLOCATE (MOS(T)%VEGIP(LDAS%MOS_NMVEGP))
       ALLOCATE (MOS(T)%SOILP(LDAS%MOS_NSOILP))
       IF (LDAS%LAI.GE.2) THEN
        ALLOCATE (MOS(T)%LAI_T1_F)
        ALLOCATE (MOS(T)%LAI_T2_F)
        ALLOCATE (MOS(T)%DSAI_T1)
        ALLOCATE (MOS(T)%DSAI_T2)
        ALLOCATE (MOS(T)%GREEN1)
        ALLOCATE (MOS(T)%GREEN2)
        ALLOCATE (MOS(T)%LAI)
        ALLOCATE (MOS(T)%DSAI)
        ALLOCATE (MOS(T)%GREEN)
       ENDIF
      ENDDO

!=== Get Vegetation Parameters for Mosaic Model in Tile Space

!=== Read in the Mosaic Static and Monthly Vegetation Parameter Files
      OPEN(UNIT=15,FILE=LDAS%MOS_VFILE,STATUS='OLD')
      OPEN(UNIT=16,FILE=LDAS%MOS_MVFILE,STATUS='OLD')

      DO J=1,LDAS%MOS_NMVEGP
       DO K=1,12
        READ(16,*)(VALUEMON(I,J,K),I=1,LDAS%NT)
       ENDDO !K
      ENDDO !J
      CLOSE(16)

      DO J=1,LDAS%MOS_NVEGP
       READ(15,*)(VALUE(I,J),I=1,LDAS%NT)
      ENDDO 
      CLOSE(15)

!=== Assign STATIC vegetation parameters to each tile based on the
!=== type of vegetation present in that tile.
!=== These parameters will be stored in one long array--structured
!=== as follows: Tile 1, all the parameters (1 through numparam)
!=== then Tile 2, all the parameters. 
!=== Then Tile 3, all the parameters etc.
      DO I=1,LDAS%NCH
       DO J=1,LDAS%MOS_NVEGP
        MOS(I)%VEGP(J)=VALUE(TILE(I)%VEGT,J)				
       ENDDO !J
      ENDDO !I

!=== Assign MONTHLY vegetation parameters to each tile based on the
!=== type of vegetation present in that tile.
!=== These parameters will be stored in one long array--structured
!=== as follows: Tile 1, parameter 1 Jan. through paramter 1 Dec.,
!=== Tile 1, parameter 2 Jan. through Dec., and so on until all
!=== months of all the paramters are cycled through..then the
!=== process repeats for Tile 2 and so forth until all the
!=== tiles are cycled through
      DO I=1,LDAS%NCH
       DO J=1,LDAS%MOS_NMVEGP
        DO K=1,12
         MOS(I)%VEGMP(J,K)=VALUEMON(TILE(I)%VEGT,J,K)
        ENDDO !K
       ENDDO !J
      ENDDO !I

!=== Assign and Allocate Number of soil layers (Always 3 for Mosaic)
      DO T=1,LDAS%NCH
       MOS(T)%NSLAY=3
       ALLOCATE (MOS(T)%SoWet(MOS(T)%NSLAY))
      ENDDO
	
!+++ Vegetation-based soil parameterizations.
!    Option 1 (original)
!    Option 5 (soil layer depths are veg based, while other pararams are not)
      if ((ldas%soil .eq. 1).or.(ldas%soil.eq.5)) then

!===   Get Soil Parameters (Based on Vegetation) for Mosaic Model in Tile Space
!===   Read in Soil Parameter Data
	    
        OPEN(10,FILE=LDAS%MOS_SFILE,STATUS='OLD',
     &   ACCESS='SEQUENTIAL')
     
        DO I=1,LDAS%MOS_NSOILP
         READ(10,*)(BASICSET(JJ,I),JJ=1,LDAS%NT)
        ENDDO
        CLOSE(10)
    
        DO I=1,LDAS%NCH
         K=TILE(I)%VEGT
          DO J=1,LDAS%MOS_NSOILP
           MOS(I)%SOILP(J)=BASICSET(K,J)                  
          ENDDO !J 
        ENDDO !I	 

        end if	!soil=1/5

!    New soil parameterization.
      if ((ldas%soil.eq.2).or.(ldas%soil.eq.3)) then

!+++   Open soil files (sand, silt, clay, porosity, and slope).
        OPEN(11,FILE=LDAS%SAFILE,FORM='UNFORMATTED',STATUS='OLD')
        OPEN(12,FILE=LDAS%SIFILE,FORM='UNFORMATTED',STATUS='OLD')
        OPEN(13,FILE=LDAS%CLFILE,FORM='UNFORMATTED',STATUS='OLD')
        OPEN(14,FILE=LDAS%POFILE,FORM='UNFORMATTED',STATUS='OLD')
        OPEN(15,FILE=LDAS%SLFILE,FORM='UNFORMATTED',STATUS='OLD')

!+++   Read soil properties and compute soil parameters in grid space.
!+++   Note that the first 4 files each contain 3 records, one
!+++    for each layer thickness (0-2, 2-150, 150-350 cm or 0-10, 10-40, 40-200 cm
!+++    depending on what soil scheme was chosen).  At this time only
!+++    the top layer data are used to calculate soil parameters - except
!+++    for the water holding capacity, which is layer-specific, based on
!+++    porosity and layer thickness.
        READ(11) SAND1
        READ(12) SILT1
        READ(13) CLAY1
        READ(14) POR1
        READ(14) POR2
        READ(14) POR3
        READ(15) SLOPE

        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
        CLOSE(14)
        CLOSE(15)

!+++    Determine USDA texture class (local variable).
        CALL TEXTURE (LDAS%NC,LDAS%NR,SAND1,SILT1,CLAY1,TEX1)

!+++    Determine saturated hydraulic conductivity at the surface.
        CALL HYDCON_TEX (LDAS%NC,LDAS%NR,TEX1,KSAT1)

!+++    Determine saturated soil potential based on texture class.
        CALL PSIS (LDAS%NC,LDAS%NR,TEX1,PSI1)

!+++    Determine parameter b based on texture class.
        CALL PARAMB (LDAS%NC,LDAS%NR,TEX1,B1)

!+++    Assign soil parameters in tile space.
        write(*,*) 'Assigning soil parameters in tile space.'
        DO I=1,LDAS%NCH
         MOS(I)%SOILP(1) = B1(TILE(I)%COL,TILE(I)%ROW)
	  if (MOS(I)%SOILP(1) .le. 0.0) then 
	    write(*,*) 'COL,ROW,TILE,B1',TILE(I)%COL,TILE(I)%ROW,
     &	     I,MOS(I)%SOILP(1)
	    stop
	  end if
         MOS(I)%SOILP(2) = PSI1(TILE(I)%COL,TILE(I)%ROW)
	  if (MOS(I)%SOILP(2) .ge. 0.0) then 
	    write(*,*) 'COL,ROW,TILE,PSI1',TILE(I)%COL,TILE(I)%ROW,
     &	     I,MOS(I)%SOILP(2)
	    stop
	  end if
         MOS(I)%SOILP(3) = KSAT1(TILE(I)%COL,TILE(I)%ROW)
	  if (MOS(I)%SOILP(3) .le. 0.0) then 
	    write(*,*) 'COL,ROW,TILE,KSAT1',TILE(I)%COL,TILE(I)%ROW,
     &	     I,MOS(I)%SOILP(3)
	    stop
	  end if
C	Set layer thicknesses based on soil scheme chosen in card file
	IF (LDAS%SOIL.EQ.2) THEN
         MOS(I)%SOILP(4) = D1
         MOS(I)%SOILP(5) = D2
         MOS(I)%SOILP(6) = D3
        ENDIF
        IF (LDAS%SOIL.EQ.3) THEN
         MOS(I)%SOILP(4) = DNLDAS1
         MOS(I)%SOILP(5) = DNLDAS2
         MOS(I)%SOILP(6) = DNLDAS3
        ENDIF

!+++     Special for slope: use vegetation-based values as minima for
!+++      desert-like vegetation classes.
         SELECT CASE (TILE(I)%VEGT)
	  CASE (8)	! Closed shrubland
	   MOS(I)%SOILP(7) = MAX(S8, SIN(SLOPE(TILE(I)%COL,
     &	    TILE(I)%ROW)))
	  CASE (9)	! Open shrubland
	   MOS(I)%SOILP(7) = MAX(S9, SIN(SLOPE(TILE(I)%COL,
     &	    TILE(I)%ROW)))
	  CASE (12)	! Bare ground
	   MOS(I)%SOILP(7) = MAX(S12, SIN(SLOPE(TILE(I)%COL,
     &	    TILE(I)%ROW)))
	  CASE DEFAULT	! All other vegetation classes
	   MOS(I)%SOILP(7) = MAX(S0, SIN(SLOPE(TILE(I)%COL,
     &	    TILE(I)%ROW)))
         END SELECT
	  if (MOS(I)%SOILP(7) .le. 0.0) then 
	    write(*,*) 'COL,ROW,TILE,SIN(THETA)',TILE(I)%COL,TILE(I)%ROW,
     &	     I,MOS(I)%SOILP(7)
	  end if

         POR1A=AMIN1(POR1(TILE(I)%COL,TILE(I)%ROW),PORMAX)
         POR2A=AMIN1(POR2(TILE(I)%COL,TILE(I)%ROW),PORMAX)
         POR3A=AMIN1(POR3(TILE(I)%COL,TILE(I)%ROW),PORMAX)


        IF (LDAS%SOIL.EQ.2) THEN
         MOS(I)%SOILP(8) = POR1A * D1 * 1000.0
	ENDIF
        IF (LDAS%SOIL.EQ.3) THEN
         MOS(I)%SOILP(8) = POR1A * DNLDAS1 * 1000.0
	ENDIF
	  if (MOS(I)%SOILP(8) .le. 0.0) then 
	    write(*,*) 'COL,ROW,TILE,WSAT1',TILE(I)%COL,TILE(I)%ROW,
     &	     I,MOS(I)%SOILP(8)
	    stop
	  end if

        IF (LDAS%SOIL.EQ.2) THEN
         MOS(I)%SOILP(9) = POR2A * D2 * 1000.0
	ENDIF
        IF (LDAS%SOIL.EQ.3) THEN
         MOS(I)%SOILP(9) = POR2A * DNLDAS2 * 1000.0
	ENDIF

	  if (MOS(I)%SOILP(9) .le. 0.0) then 
	    write(*,*) 'COL,ROW,TILE,WSAT2',TILE(I)%COL,TILE(I)%ROW,
     &	     I,MOS(I)%SOILP(9)
	    stop
	  end if
        IF (LDAS%SOIL.EQ.2) THEN
         MOS(I)%SOILP(10) = POR3A * D3 * 1000.0
	ENDIF
        IF (LDAS%SOIL.EQ.3) THEN
         MOS(I)%SOILP(10) = POR3A * DNLDAS3 * 1000.0
	ENDIF
	  if (MOS(I)%SOILP(10) .le. 0.0) then 
	    write(*,*) 'COL,ROW,TILE,WSAT3',TILE(I)%COL,TILE(I)%ROW,
     &	     I,MOS(I)%SOILP(10)
	    stop
	  end if

        ENDDO !I	 

      end if	!soil=2 or 3

!    N-LDAS COSBY/RAWLS soil parameterization.

	if ((ldas%soil.eq.4).or.(ldas%soil.eq.5)) then
        DATA POROS/0.339,0.421,0.434,0.476,0.476,0.439,0.404,0.464,
     &             0.465,0.406,0.468,0.457/
        DATA HYCON/3.75E-05,1.81E-05,6.36E-06,2.14E-06,2.14E-06,
     &             1.44E-06, 
     &             9.72E-07,
     &             1.19E-06,3.33E-07,2.78E-07,4.44E-07,5.28E-07/
        DATA SOILPOT/-0.0692,-0.0363,-0.1413,-0.7586,-0.7586,-0.3548,
     &               -0.1349,-0.6166,-0.263,-0.0977,-0.3236,-0.4677/
        DATA BPARAM/2.79,4.26,4.74,5.33,5.33,5.25,6.77,
     &               8.72,8.17,10.73,10.39,11.55/
	endif

      if (ldas%soil .eq. 4) then

!+++   Open soil file (slope).
        OPEN(14,FILE=LDAS%TXFILE,STATUS='OLD')
        OPEN(15,FILE=LDAS%SLFILE,FORM='UNFORMATTED',STATUS='OLD')

!+++    Read soil texture and compute soil parameters in grid space.
!+++    At this time only
!+++    the top layer data are used to calculate soil parameters - except
!+++    for the water holding capacity, which is layer-specific, based on
!+++    porosity and layer depth.  Porosity values for 3 mosaic layers
!++     were derived from wieghted average of 11 layers of LDAS-Yun data
!++     Soil parameters come mostly from Cosby, but Sand class parameters and ksat
!++     values for all classes come from Rawls

        do k=1,11
        do j=1,224
        do i=1,464
        read (14,'(I3,1X,I3,1X,3X,I2)')
     &  ii,jj,soiltext(k,i,j)
        if (soiltext(k,i,j).eq.13) soiltext(k,i,j)=12
        enddo
        enddo
        enddo

        print *,'   '
        print *,'SETTING SOIL CLASS 13 to 12...since have no params for 
     &  organic class'
        print *,'WARNING!!!!!!!!!!!!!!!!!!'
        print *,'   '

        READ(15) SLOPE
        CLOSE(14)
        CLOSE(15)

!+++    Determine Porosity at 3 Levels. HARDWIRED FOR 0-10,10-40,40-200cm
!++     Porosity data for each of 12 soil textures
c        DATA POROS/0.339,0.421,0.434,0.476,0.476,0.439,0.404,0.464,
c     &             0.465,0.406,0.468,0.457/

        do j=1,224
        do i=1,464
        POR1(i,j)=(0.5*poros(soiltext(1,i,j)))+(0.5*
     &             poros(soiltext(2,i,j)))
        POR2(i,j)=((1.0/3.0)*poros(soiltext(3,i,j)))+
     &            ((1.0/3.0)*poros(soiltext(4,i,j)))+
     &            ((1.0/3.0)*poros(soiltext(5,i,j)))
        POR3(i,j)=(0.125*poros(soiltext(6,i,j)))+(0.125*
     &             poros(soiltext(7,i,j)))+
     &            (0.125*poros(soiltext(8,i,j)))+(0.3125*
     &            poros(soiltext(9,i,j)))+
     &            (0.3125*poros(soiltext(10,i,j)))
        enddo
        enddo

!+++    Determine saturated hydraulic conductivity at the surface.
c        DATA HYCON/3.75E-05,1.81E-05,6.36E-06,2.14E-06,2.14E-06,
c     &             1.44E-06,
c     &             9.72E-07,
c     &             1.19E-06,3.33E-07,2.78E-07,4.44E-07,5.28E-07/
        do j=1,224
        do i=1,464
        KSAT1(i,j)=HYCON(soiltext(1,i,j))
        enddo
        enddo
!+++    Determine saturated soil potential based on texture class.
c        DATA SOILPOT/-0.0692,-0.0363,-0.1413,-0.7586,-0.7586,-0.3548,
c     &               -0.1349,-0.6166,-0.263,-0.0977,-0.3236,-0.4677/
        do j=1,224
        do i=1,464
        PSI1(i,j)=SOILPOT(soiltext(1,i,j))
        enddo
        enddo

!+++    Determine parameter b based on texture class.
c        DATA BPARAM/2.79,4.26,4.74,5.33,5.33,5.25,6.77,
c     &               8.72,8.17,10.73,10.39,11.55/
        do j=1,224
        do i=1,464
        B1(i,j)=BPARAM(soiltext(1,i,j))
        enddo
        enddo


!+++    Assign soil parameters in tile space.
        write(*,*) 'Assigning soil parameters in tile space.'
        DO I=1,LDAS%NCH
         MOS(I)%SOILP(1) = B1(TILE(I)%COL,TILE(I)%ROW)
          if (MOS(I)%SOILP(1) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,B1',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(1)
            stop
          end if
         MOS(I)%SOILP(2) = PSI1(TILE(I)%COL,TILE(I)%ROW)
          if (MOS(I)%SOILP(2) .ge. 0.0) then
            write(*,*) 'COL,ROW,TILE,PSI1',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(2)
            stop
          end if
         MOS(I)%SOILP(3) = KSAT1(TILE(I)%COL,TILE(I)%ROW)
          if (MOS(I)%SOILP(3) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,KSAT1',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(3)
            stop
          end if
         MOS(I)%SOILP(4) = DNLDAS1
         MOS(I)%SOILP(5) = DNLDAS2
         MOS(I)%SOILP(6) = DNLDAS3
!+++     Special for slope: use vegetation-based values as minima for
!+++      desert-like vegetation classes.
         SELECT CASE (TILE(I)%VEGT)
          CASE (8)      ! Closed shrubland
           MOS(I)%SOILP(7) = MAX(S8, SIN(SLOPE(TILE(I)%COL,
     &      TILE(I)%ROW)))
          CASE (9)      ! Open shrubland
           MOS(I)%SOILP(7) = MAX(S9, SIN(SLOPE(TILE(I)%COL,
     &      TILE(I)%ROW)))
          CASE (12)     ! Bare ground
           MOS(I)%SOILP(7) = MAX(S12, SIN(SLOPE(TILE(I)%COL,
     &      TILE(I)%ROW)))
          CASE DEFAULT  ! All other vegetation classes
           MOS(I)%SOILP(7) = MAX(S0, SIN(SLOPE(TILE(I)%COL,
     &      TILE(I)%ROW)))
         END SELECT
          if (MOS(I)%SOILP(7) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,SIN(THETA)',TILE(I)%COL,
     &      TILE(I)%ROW,
     &       I,MOS(I)%SOILP(7)
          end if

         POR1A=AMIN1(POR1(TILE(I)%COL,TILE(I)%ROW),PORMAX)
         POR2A=AMIN1(POR2(TILE(I)%COL,TILE(I)%ROW),PORMAX)
         POR3A=AMIN1(POR3(TILE(I)%COL,TILE(I)%ROW),PORMAX)
         MOS(I)%SOILP(8) = POR1A * DNLDAS1 * 1000.0
          if (MOS(I)%SOILP(8) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,WSAT1',TILE(I)%COL,
     &  TILE(I)%ROW,
     &       I,MOS(I)%SOILP(8)
            stop
          end if
         MOS(I)%SOILP(9) = POR2A * DNLDAS2 * 1000.0
          if (MOS(I)%SOILP(9) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,WSAT2',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(9)
            stop
          end if
         MOS(I)%SOILP(10) = POR3A * DNLDAS3 * 1000.0
          if (MOS(I)%SOILP(10) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,WSAT3',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(10)
            stop
          end if

        ENDDO !I

      end if    !soil=4


!    N-LDAS COSBY/RAWLS soil parameterization WITH RANDY's VEG-DEPENDENT SOIL LAYER THICKNESSES.
      if (ldas%soil .eq. 5) then
 
!+++   Open soil file (slope).
        OPEN(14,FILE=LDAS%TXFILE,STATUS='OLD')
        OPEN(15,FILE=LDAS%SLFILE,FORM='UNFORMATTED',STATUS='OLD')
 
!+++    Read soil texture and compute soil parameters in grid space.
!+++    At this time only
!+++    the top layer data are used to calculate soil parameters - except
!+++    for the water holding capacity, which is layer-specific, based on
!+++    porosity and layer depth.  Porosity values for 3 mosaic layers
!++     were derived from wieghted average of 11 layers of LDAS-Yun data
!++     Soil parameters come mostly from Cosby, but Sand class parameters and ksat
!++     values for all classes come from Rawls
 
        do k=1,11
        do j=1,224
        do i=1,464
        read (14,'(I3,1X,I3,1X,3X,I2)')
     &  ii,jj,soiltext(k,i,j)
        if (soiltext(k,i,j).eq.13) soiltext(k,i,j)=12
        enddo
        enddo
        enddo
 
        print *,'   '
        print *,'SETTING SOIL CLASS 13 to 12...since have no params for
     &  organic class'
        print *,'WARNING!!!!!!!!!!!!!!!!!!'
        print *,'   '
 
        READ(15) SLOPE
        CLOSE(14)
        CLOSE(15)
 
!+++    Determine Porosity at 3 Levels.  Hardwired to Randy's veg-dependent depths
!++     Porosity data for each of 12 soil textures, varies by veg tile
c        DATA POROS/0.339,0.421,0.434,0.476,0.476,0.439,0.404,0.464,
c     &             0.465,0.406,0.468,0.457/


	DO I=1,LDAS%NCH

        if ((tile(i)%vegt.ge.1).and.(tile(i)%vegt.le.5)) then
	por1tile(i)=poros(soiltext(1,tile(i)%col,tile(i)%row))
	por2tile(i)=(
     &  (0.02027*poros(soiltext(1,tile(i)%col,tile(i)%row)))+
     &  (0.033784*poros(soiltext(2,tile(i)%col,tile(i)%row)))+
     &  (0.067568*poros(soiltext(3,tile(i)%col,tile(i)%row)))+
     &  (0.067568*poros(soiltext(4,tile(i)%col,tile(i)%row)))+
     &  (0.067568*poros(soiltext(5,tile(i)%col,tile(i)%row)))+
     &  (0.135135*poros(soiltext(6,tile(i)%col,tile(i)%row)))+
     &  (0.135135*poros(soiltext(7,tile(i)%col,tile(i)%row)))+
     &  (0.135135*poros(soiltext(8,tile(i)%col,tile(i)%row)))+
     &  (0.337838*poros(soiltext(9,tile(i)%col,tile(i)%row))))
	 por3tile(i)=(
     &  (0.25*poros(soiltext(10,tile(i)%col,tile(i)%row)))+
     &  (0.75*poros(soiltext(11,tile(i)%col,tile(i)%row))))
	endif

        if (tile(i)%vegt.eq.6) then
        por1tile(i)=poros(soiltext(1,tile(i)%col,tile(i)%row))
        por2tile(i)=(      
     &  (0.02027*poros(soiltext(1,tile(i)%col,tile(i)%row)))+
     &  (0.023851*poros(soiltext(2,tile(i)%col,tile(i)%row)))+
     &  (0.039752*poros(soiltext(3,tile(i)%col,tile(i)%row)))+
     &  (0.079504*poros(soiltext(4,tile(i)%col,tile(i)%row)))+
     &  (0.079504*poros(soiltext(5,tile(i)%col,tile(i)%row)))+
     &  (0.079504*poros(soiltext(6,tile(i)%col,tile(i)%row)))+
     &  (0.159008*poros(soiltext(7,tile(i)%col,tile(i)%row)))+
     &  (0.159008*poros(soiltext(8,tile(i)%col,tile(i)%row)))+
     &  (0.220862*poros(soiltext(9,tile(i)%col,tile(i)%row))))
         por3tile(i)=(      
     &  (0.124831*poros(soiltext(9,tile(i)%col,tile(i)%row)))+
     &  (0.280899*poros(soiltext(10,tile(i)%col,tile(i)%row)))+
     &  (0.59427*poros(soiltext(11,tile(i)%col,tile(i)%row))))
        endif

        if (tile(i)%vegt.eq.7) then
        por1tile(i)=poros(soiltext(1,tile(i)%col,tile(i)%row))
        por2tile(i)=(
     &  (0.036345*poros(soiltext(1,tile(i)%col,tile(i)%row)))+
     &  (0.060575*poros(soiltext(2,tile(i)%col,tile(i)%row)))+
     &  (0.121151*poros(soiltext(3,tile(i)%col,tile(i)%row)))+
     &  (0.121151*poros(soiltext(4,tile(i)%col,tile(i)%row)))+
     &  (0.121151*poros(soiltext(5,tile(i)%col,tile(i)%row)))+
     &  (0.242301*poros(soiltext(6,tile(i)%col,tile(i)%row)))+
     &  (0.242301*poros(soiltext(7,tile(i)%col,tile(i)%row)))+
     &  (0.055025*poros(soiltext(8,tile(i)%col,tile(i)%row))))
         por3tile(i)=(       
     &  (0.114344*poros(soiltext(8,tile(i)%col,tile(i)%row)))+
     &  (0.36985*poros(soiltext(9,tile(i)%col,tile(i)%row)))+
     &  (0.36985*poros(soiltext(10,tile(i)%col,tile(i)%row)))+
     &  (0.145957*poros(soiltext(11,tile(i)%col,tile(i)%row))))
        endif

        if (tile(i)%vegt.eq.8) then
        por1tile(i)=poros(soiltext(1,tile(i)%col,tile(i)%row))
        por2tile(i)=(
     &  (0.141963*poros(soiltext(1,tile(i)%col,tile(i)%row)))+
     &  (0.204232*poros(soiltext(2,tile(i)%col,tile(i)%row)))+
     &  (0.408464*poros(soiltext(3,tile(i)%col,tile(i)%row)))+
     &  (0.245342*poros(soiltext(4,tile(i)%col,tile(i)%row))))
         por3tile(i)=(
     &  (0.057728*poros(soiltext(4,tile(i)%col,tile(i)%row)))+
     &  (0.144553*poros(soiltext(5,tile(i)%col,tile(i)%row)))+
     &  (0.289105*poros(soiltext(6,tile(i)%col,tile(i)%row)))+
     &  (0.289105*poros(soiltext(7,tile(i)%col,tile(i)%row)))+
     &  (0.21951*poros(soiltext(8,tile(i)%col,tile(i)%row))))
        endif

        if (tile(i)%vegt.eq.9) then
        por1tile(i)=poros(soiltext(1,tile(i)%col,tile(i)%row))
        por2tile(i)=(
     &  (0.531749*poros(soiltext(1,tile(i)%col,tile(i)%row)))+
     &  (0.468251*poros(soiltext(2,tile(i)%col,tile(i)%row))))
         por3tile(i)=( 
     &  (0.038294*poros(soiltext(2,tile(i)%col,tile(i)%row)))+
     &  (0.24717*poros(soiltext(3,tile(i)%col,tile(i)%row)))+
     &  (0.24717*poros(soiltext(4,tile(i)%col,tile(i)%row)))+
     &  (0.24717*poros(soiltext(5,tile(i)%col,tile(i)%row)))+
     &  (0.220196*poros(soiltext(6,tile(i)%col,tile(i)%row))))
        endif
	
        if ((tile(i)%vegt.eq.10).or.(tile(i)%vegt.eq.11)) then
        por1tile(i)=poros(soiltext(1,tile(i)%col,tile(i)%row))
        por2tile(i)=(
     &  (0.06383*poros(soiltext(1,tile(i)%col,tile(i)%row)))+
     &  (0.106383*poros(soiltext(2,tile(i)%col,tile(i)%row)))+
     &  (0.212766*poros(soiltext(3,tile(i)%col,tile(i)%row)))+
     &  (0.212766*poros(soiltext(4,tile(i)%col,tile(i)%row)))+
     &  (0.212766*poros(soiltext(5,tile(i)%col,tile(i)%row)))+
     &  (0.191489*poros(soiltext(6,tile(i)%col,tile(i)%row))))
         por3tile(i)=(
     &  (0.11*poros(soiltext(6,tile(i)%col,tile(i)%row)))+
     &  (0.2*poros(soiltext(7,tile(i)%col,tile(i)%row)))+
     &  (0.2*poros(soiltext(8,tile(i)%col,tile(i)%row)))+
     &  (0.49*poros(soiltext(9,tile(i)%col,tile(i)%row))))
        endif

        if (tile(i)%vegt.eq.12) then
        por1tile(i)=poros(soiltext(1,tile(i)%col,tile(i)%row))
        por2tile(i)=(
     &  (1.0*poros(soiltext(1,tile(i)%col,tile(i)%row))))
        por3tile(i)=(
     &  (0.105333*poros(soiltext(1,tile(i)%col,tile(i)%row)))+
     &  (0.166667*poros(soiltext(2,tile(i)%col,tile(i)%row)))+
     &  (0.333333*poros(soiltext(3,tile(i)%col,tile(i)%row)))+
     &  (0.333333*poros(soiltext(4,tile(i)%col,tile(i)%row)))+
     &  (0.061333*poros(soiltext(5,tile(i)%col,tile(i)%row))))
        endif

        if (tile(i)%vegt.eq.13) then
        por1tile(i)=poros(soiltext(1,tile(i)%col,tile(i)%row))
        por2tile(i)=(
     &  (0.051115*poros(soiltext(1,tile(i)%col,tile(i)%row)))+
     &  (0.082654*poros(soiltext(2,tile(i)%col,tile(i)%row)))+
     &  (0.165308*poros(soiltext(3,tile(i)%col,tile(i)%row)))+
     &  (0.165308*poros(soiltext(4,tile(i)%col,tile(i)%row)))+
     &  (0.165308*poros(soiltext(5,tile(i)%col,tile(i)%row)))+
     &  (0.330616*poros(soiltext(6,tile(i)%col,tile(i)%row)))+
     &  (0.039691*poros(soiltext(7,tile(i)%col,tile(i)%row))))
        por3tile(i)=(
     &  (0.157402*poros(soiltext(7,tile(i)%col,tile(i)%row)))+
     &  (0.178876*poros(soiltext(8,tile(i)%col,tile(i)%row)))+
     &  (0.447191*poros(soiltext(9,tile(i)%col,tile(i)%row)))+
     &  (0.21653*poros(soiltext(10,tile(i)%col,tile(i)%row))))
        endif
	ENDDO

!+++    Determine saturated hydraulic conductivity at the surface.
c        DATA HYCON/3.75E-05,1.81E-05,6.36E-06,2.14E-06,2.14E-06,
c     &             1.44E-06,
c     &             9.72E-07,
c     &             1.19E-06,3.33E-07,2.78E-07,4.44E-07,5.28E-07/
        do j=1,224
        do i=1,464
        KSAT1(i,j)=HYCON(soiltext(1,i,j))
        enddo
        enddo
!+++    Determine saturated soil potential based on texture class.
c        DATA SOILPOT/-0.0692,-0.0363,-0.1413,-0.7586,-0.7586,-0.3548,
c     &               -0.1349,-0.6166,-0.263,-0.0977,-0.3236,-0.4677/
        do j=1,224
        do i=1,464
        PSI1(i,j)=SOILPOT(soiltext(1,i,j))
        enddo
        enddo
 
!+++    Determine parameter b based on texture class.
c        DATA BPARAM/2.79,4.26,4.74,5.33,5.33,5.25,6.77,
c     &               8.72,8.17,10.73,10.39,11.55/
        do j=1,224
        do i=1,464
        B1(i,j)=BPARAM(soiltext(1,i,j))
        enddo
        enddo
 
 
!+++    Assign soil parameters in tile space.
        write(*,*) 'Assigning soil parameters in tile space.'
        DO I=1,LDAS%NCH
         MOS(I)%SOILP(1) = B1(TILE(I)%COL,TILE(I)%ROW)
          if (MOS(I)%SOILP(1) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,B1',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(1)
            stop
          end if
         MOS(I)%SOILP(2) = PSI1(TILE(I)%COL,TILE(I)%ROW)
          if (MOS(I)%SOILP(2) .ge. 0.0) then
            write(*,*) 'COL,ROW,TILE,PSI1',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(2)
            stop
          end if
         MOS(I)%SOILP(3) = KSAT1(TILE(I)%COL,TILE(I)%ROW)
          if (MOS(I)%SOILP(3) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,KSAT1',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(3)
            stop
          end if
C	Following lines are commented out so that soil depths for each
C       tile are kept as assigned at start of program, and are 
C       veg-related
c         MOS(I)%SOILP(4) = DNLDAS1
c         MOS(I)%SOILP(5) = DNLDAS2
c         MOS(I)%SOILP(6) = DNLDAS3
!+++     Special for slope: use vegetation-based values as minima for
!+++      desert-like vegetation classes.
         SELECT CASE (TILE(I)%VEGT)
          CASE (8)      ! Closed shrubland
           MOS(I)%SOILP(7) = MAX(S8, SIN(SLOPE(TILE(I)%COL,
     &      TILE(I)%ROW)))
          CASE (9)      ! Open shrubland
           MOS(I)%SOILP(7) = MAX(S9, SIN(SLOPE(TILE(I)%COL,
     &      TILE(I)%ROW)))
          CASE (12)     ! Bare ground
           MOS(I)%SOILP(7) = MAX(S12, SIN(SLOPE(TILE(I)%COL,
     &      TILE(I)%ROW)))
          CASE DEFAULT  ! All other vegetation classes
           MOS(I)%SOILP(7) = MAX(S0, SIN(SLOPE(TILE(I)%COL,
     &      TILE(I)%ROW)))
         END SELECT
          if (MOS(I)%SOILP(7) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,SIN(THETA)',TILE(I)%COL,
     &      TILE(I)%ROW,
     &       I,MOS(I)%SOILP(7)
          end if
 
         POR1A=AMIN1(POR1TILE(I),PORMAX)
         POR2A=AMIN1(POR2TILE(I),PORMAX)
         POR3A=AMIN1(POR3TILE(I),PORMAX)
         MOS(I)%SOILP(8) = POR1A * mos(i)%soilp(4) * 1000.0
          if (MOS(I)%SOILP(8) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,WSAT1',TILE(I)%COL,
     &  TILE(I)%ROW,
     &       I,MOS(I)%SOILP(8)
            stop
          end if
         MOS(I)%SOILP(9) = POR2A * mos(i)%soilp(5) * 1000.0
          if (MOS(I)%SOILP(9) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,WSAT2',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(9)
            stop
          end if
         MOS(I)%SOILP(10) = POR3A * mos(i)%soilp(6) * 1000.0
          if (MOS(I)%SOILP(10) .le. 0.0) then
            write(*,*) 'COL,ROW,TILE,WSAT3',TILE(I)%COL,TILE(I)%ROW,
     &       I,MOS(I)%SOILP(10)
            stop
          end if
 
        ENDDO !I
 
      end if    !soil=5

      RETURN
      END

