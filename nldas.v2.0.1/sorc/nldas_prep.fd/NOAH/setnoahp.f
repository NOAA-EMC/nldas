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
! SETNOAHP.f: 
!
! DESCRIPTION:
!  This subroutine retrieves NOAH parameters - Significant F90 revisions
!   below this subroutine will be required in the future.  
!
! REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added NOAH LSM, Initial Code
!  04 Nov 2002: Kristi Arsenault;  Added new resolutions for MXSNALB
!                                  and incorporated bottom temp fields
!                                  with elevation adjustments
!  07 Jan 2003: Kristi Arsenault;  Corrected MXSNALB and TBOT to work
!                                  properly with multiple tiles
!  21 May 2003: Kristi Arsenault;  Enabled Soil Options 2-4 (corresponding
!                                    to card file)
!  24 Jun 2003: Kristi Arsenault;  Added and initialize time flags for 
!                                  noah_gfrac and noah_alb
!=========================================================================

      SUBROUTINE SETNOAHP(LDAS,TILE,NOAH)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      USE grid_module      ! LDAS non-model-specific grid variables
      USE noah_module      ! NOAH tile variables

      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (tiledec) TILE(LDAS%NCH)
      type (griddec) GRID(LDAS%NC,LDAS%NR)
      type (noahdec) NOAH(LDAS%NCH)

!+++  Layer depths correspond to soil property map files 
!      Thicknesses of soil layers (set originally for Mosaic) 
!      These are compatible with the FIRST set of Matt Rodell's
!      soil files, i.e. sand_nldas3.5.bfsa

      REAL, PARAMETER :: D1 = 0.02      ! thickness of soil layer 1, m
      REAL, PARAMETER :: D2 = 1.48      ! thickness of soil layer 2, m
      REAL, PARAMETER :: D3 = 2.00      ! thickness of soil layer 3, m

!     Thicknesss of soil layers (semi-official) for NLDAS simulations
!      These are compatible with either Yun Duans's soil maps
!      or the SECOND set of Matt Rodell's soil files, i.e., sand_nldas2m.bfsa
      REAL, PARAMETER :: DNLDAS1 = 0.1    ! thickness of soil layer 1, m
      REAL, PARAMETER :: DNLDAS2 = 0.3    ! thickness of soil layer 2, m
      REAL, PARAMETER :: DNLDAS3 = 1.6    ! thickness of soil layer 3, m

!** NOAH LSM Layer Thicknesses:
!      REAL, PARAMETER :: D1 = 0.10      ! Thickness of soil layer 1, m
!      REAL, PARAMETER :: D2 = 0.30      ! Thickness of soil layer 2, m
!      REAL, PARAMETER :: D3 = 0.60      ! Thickness of soil layer 3, m
!      REAL, PARAMETER :: D4 = 1.00      ! Thickness of soil layer 4, m 
!
! - NOTE: Lowest boundary layer associated with temperature, TBOT, 
!         is specified in NOAH code as ZBOT = 8.0 m

!+++  Maximum allowable porosity.
      REAL, PARAMETER :: PORMAX=0.70

!=== Local Variables =====================================================

      INTEGER :: N,C,R,T,I,J,K,II,JJ               !Loop counters
      INTEGER :: LINE
        
      INTEGER :: SIBVEG(LDAS%NCH)                  !SIB2 Vegetation Classes  
      REAL :: VEGP(LDAS%NCH,LDAS%NOAH_NVEGP)       !Static Veg. Parameters
      REAL :: VEGMP(LDAS%NC,LDAS%NR)               !Monthly GFRAC data

      INTEGER :: SOILTEXT(11,LDAS%NC,LDAS%NR)      !Yun's 11-layer soil classes
      INTEGER :: SOILTYP(LDAS%NC,LDAS%NR)
      INTEGER :: PLACESLTYP(LDAS%NCH)
      INTEGER :: ZOBSOIL(LDAS%NCH)                 !Zobler Soil Classes 
      REAL :: SOILP(LDAS%NCH,LDAS%NOAH_NSOILP)   
      REAL :: VALUE(LDAS%NT,LDAS%NOAH_NVEGP)
      REAL :: VALUEMON(LDAS%NT,LDAS%NOAH_NMVEGP,12)
      REAL :: BASICSET(LDAS%NOAH_ZST,LDAS%NOAH_NSOILP)
      REAL :: SAND1(LDAS%NC,LDAS%NR)
      REAL :: CLAY1(LDAS%NC,LDAS%NR)

c      REAL :: SILT1(LDAS%NC,LDAS%NR)
      REAL :: POR1(LDAS%NC,LDAS%NR)
      REAL :: POR2(LDAS%NC,LDAS%NR)
      REAL :: POR3(LDAS%NC,LDAS%NR)
      REAL :: POR1A,POR2A,POR3A
c      REAL :: SLOPE(LDAS%NC,LDAS%NR)
c      REAL :: KSAT1(LDAS%NC,LDAS%NR)
c      REAL :: TEX1(LDAS%NC,LDAS%NR)
c      REAL :: PSI1(LDAS%NC,LDAS%NR)
c      REAL :: B1(LDAS%NC,LDAS%NR)

      REAL :: TMPALB(LDAS%NC,LDAS%NR)
      REAL :: MXSNALB(LDAS%NCH)            !Max. Snow Albedo
      REAL :: TEMPBOT(LDAS%NCH)            !Bottom Temperature fields
      REAL :: PLACETBOT(LDAS%NC,LDAS%NR) 
      REAL, ALLOCATABLE, DIMENSION (:) :: ELEVDIFF

!=== End Variable Definition =============================================

!=== Allocate Memory for NOAH Variables ================================  
  
      DO T=1,LDAS%NCH
       ALLOCATE (NOAH(T)%VEGP(LDAS%NOAH_NVEGP))
       ALLOCATE (NOAH(T)%VEGIP(1))
       ALLOCATE (NOAH(T)%VEGMP1(1))
       ALLOCATE (NOAH(T)%VEGMP2(1))
       ALLOCATE (NOAH(T)%SOILP(LDAS%NOAH_NSOILP))
       ALLOCATE (NOAH(T)%ZOBSOIL(1))
       ALLOCATE (NOAH(T)%MXSNALB(1))
       ALLOCATE (NOAH(T)%ALBSF(1))
       ALLOCATE (NOAH(T)%ALBSF1(1))
       ALLOCATE (NOAH(T)%ALBSF2(1))
       ALLOCATE (NOAH(T)%TEMPBOT(1))
      ENDDO

!=== Initialize time flags for noah_gfrac and noah_alb
       ldas%noah_gfractime=0.0
       ldas%noah_albtime=0
       ldas%noah_albdchk=0
       ldas%noah_gfracdchk=0

!=== Convert UMD Classes to SIB Classes for Each Tile
     
       print *, 'Calling MAPVEGC to convert UMD to SIB' 

       DO N=1,LDAS%NCH
         CALL MAPVEGC(N,TILE(N)%VEGT)
!         write(56,*) N, TILE(N)%VEGT
       ENDDO   !N

!=== Get Vegetation Parameters for NOAH Model in Tile Space

!=== Read in the NOAH Static Vegetation Parameter Files

      OPEN(UNIT=11,FILE=LDAS%NOAH_VFILE,STATUS='OLD')

      DO J=1,LDAS%NOAH_NVEGP
        READ(11,*)(VALUE(I,J),I=1,LDAS%NT)
      ENDDO 
      CLOSE(11)

!=== Assign STATIC vegetation parameters to each tile based on the
!=== type of vegetation present in that tile.
!=== These parameters will be stored in one long array--structured
!=== as follows: Tile 1, all the parameters (1 through numparam)
!=== then Tile 2, all the parameters. 
!=== Then Tile 3, all the parameters etc.

      DO I=1,LDAS%NCH
       DO J=1,LDAS%NOAH_NVEGP
        NOAH(I)%VEGP(J)=VALUE(TILE(I)%VEGT,J)
       ENDDO !J
      ENDDO !I

!=== Read in bottom temperature fields and adjust for elevation differnce
!=== with either Eta (NLDAS) or GDAS (GLDAS) correction datasets. 

      SELECT CASE (LDAS%DOMAIN)

!  NLDAS (0.125 degree domain)
      CASE(1)
        print *, 'Opening NLDAS TBOT File'
        OPEN(UNIT=12, FILE=LDAS%NOAH_NTBOT,
     &    ACCESS ='DIRECT', RECL=4)
!  GLDAS domains
      CASE(2)    ! 0.25 degree domain
        print *, 'Opening GLDAS (0.25) TBOT File'
        OPEN(UNIT=12, FILE=LDAS%NOAH_GTBOT,
     &    ACCESS ='DIRECT', RECL=4)
      CASE(3)    ! 2x2.5 degree domain
        print *, 'Opening GLDAS (2x2.5) TBOT File'
        OPEN(UNIT=12, FILE=LDAS%NOAH_STBOT,
     &    ACCESS ='DIRECT', RECL=4)
      CASE(4)    ! 1.0 degree domain
        print *, 'Opening GLDAS (1.0) TBOT File'
        OPEN(UNIT=12, FILE=LDAS%NOAH_OTBOT,
     &    ACCESS ='DIRECT', RECL=4)
      CASE(5)    ! 0.5 degree domain
        print *, 'Opening GLDAS (0.5) TBOT File'
        OPEN(UNIT=12, FILE=LDAS%NOAH_ETBOT,
     &    ACCESS ='DIRECT', RECL=4)

      END SELECT

      LINE=0
      DO J=1,LDAS%NR
        DO I=1,LDAS%NC
         LINE=LINE+1
         READ(12,REC=LINE) PLACETBOT(I,J)
        ENDDO !I
      ENDDO !J
      CLOSE(12)
!  Convert TBOT from grid to tile space.
      DO I=1,LDAS%NCH
       IF(PLACETBOT(TILE(I)%COL,TILE(I)%ROW).NE.-9999.00) THEN
         NOAH(I)%TEMPBOT=PLACETBOT(TILE(I)%COL,TILE(I)%ROW)
       ENDIF 
      ENDDO

!=== The MAX SNOW ALBEDO field is opened and read in here:

      SELECT CASE (LDAS%DOMAIN)

!  NLDAS (0.125 degree domain)
      CASE(1)  
        print *, 'Opening NLDAS MAXSNALB File'
        OPEN(UNIT=14,FILE=LDAS%NOAH_NMXSNAL,
     &    ACCESS='DIRECT', RECL=4) 
!  GLDAS domains
      CASE(2)    ! 0.25 degree domain       
        print *, 'Opening GLDAS (0.25) MAXSNALB File'
        OPEN(UNIT=14,FILE=LDAS%NOAH_GMXSNAL,
     &    ACCESS='DIRECT', RECL=4)
      CASE(3)    ! 2x2.5 degree domain
        print *, 'Opening GLDAS (2x2.5) MAXSNALB File'
        OPEN(UNIT=14,FILE=LDAS%NOAH_SMXSNAL,
     &    ACCESS='DIRECT', RECL=4)
      CASE(4)    ! 1.0 degree domain
        print *, 'Opening GLDAS (1.0) MAXSNALB File'
        OPEN(UNIT=14,FILE=LDAS%NOAH_OMXSNAL,
     &    ACCESS='DIRECT', RECL=4)
      CASE(5)    ! 0.5 degree domain
        print *, 'Opening GLDAS (0.5) MAXSNALB File'
        OPEN(UNIT=14,FILE=LDAS%NOAH_EMXSNAL,
     &    ACCESS='DIRECT', RECL=4)

      END SELECT 

      LINE=0
      DO J=1,LDAS%NR
        DO I=1,LDAS%NC
         LINE=LINE+1
         READ(14,REC=LINE) TMPALB(I,J)
        ENDDO !I
      ENDDO !J
      CLOSE(14)
!  Convert max. snow albedo from grid to tile space.
      DO I=1,LDAS%NCH
        IF(TMPALB(TILE(I)%COL,TILE(I)%ROW).NE.-9999.00) THEN
          NOAH(I)%MXSNALB=TMPALB(TILE(I)%COL,TILE(I)%ROW)
        ENDIF
      END DO

!=== Assign and Allocate Number of Soil Layers (Currently, 4 for NOAH)

      DO T=1,LDAS%NCH
        NOAH(T)%NSLAY=4
        ALLOCATE (NOAH(T)%SMC(NOAH(T)%NSLAY))
      ENDDO
      DO T=1,LDAS%NCH
        NOAH(T)%NSLAY=4
        ALLOCATE (NOAH(T)%STC(NOAH(T)%NSLAY))
      ENDDO
      DO T=1,LDAS%NCH
        NOAH(T)%NSLAY=4
        ALLOCATE (NOAH(T)%SH2O(NOAH(T)%NSLAY))
      ENDDO

!=== Original vegetation-based soil parameterization - not set for Noah. 
!      if (ldas%soil .eq. 1) then

!=== NOAH-Derived Zobler Soil Classes and Associated Parameters =======

!=== Reynolds Soils options
      IF ((LDAS%SOIL.EQ.2).OR.(LDAS%SOIL.EQ.3)) THEN

!=== Open soil files (sand, clay, and porosity).
       OPEN(15,FILE=LDAS%SAFILE,FORM='UNFORMATTED',STATUS='OLD')
       OPEN(16,FILE=LDAS%CLFILE,FORM='UNFORMATTED',STATUS='OLD')
!       OPEN(17,FILE=LDAS%POFILE,FORM='UNFORMATTED',STATUS='OLD')

!=== Read soil properties and convert to Zobler soil classes.
!     Note that the 3 files each contain 3 global records, one
!     for each layer depth (0-2, 2-150, 150-350 cm).  At this
!     time only the top layer data are used to map soil 
!     parameters.  Since NOAH has 4 layers, the third layer of 
!     porosity is temporarily used for layer 4 as well.

        READ(15) SAND1
        READ(16) CLAY1
!        READ(17) POR1
!        READ(17) POR2
!        READ(17) POR3

        CLOSE(15)
        CLOSE(16)
!        CLOSE(17)

!+++ Determine Zobler-equivalent soil classes derived from
!    percentages of sand and clay, used in NOAH.

        CALL SOILTYPE(LDAS%NC,LDAS%NR,SAND1,CLAY1,SOILTYP)

!=== Assign soil parameters in tile space.
!       write(*,*) 'Assigning soil parameters in tile space.'

!=== Read in the NOAH Soil Parameter File
        OPEN(UNIT=18,FILE=LDAS%NOAH_SFILE,STATUS='OLD',
     &   ACCESS='SEQUENTIAL')
         DO I=1,LDAS%NOAH_NSOILP
           READ(18,*)(BASICSET(JJ,I),JJ=1,LDAS%NOAH_ZST)
         ENDDO
        CLOSE(18)

!=== Convert grid space to tile space for soil type values.
        DO I=1,LDAS%NCH         
          PLACESLTYP(I)=SOILTYP(TILE(I)%COL,TILE(I)%ROW)
          NOAH(I)%ZOBSOIL=PLACESLTYP(I)
!           write(44,*) I,NOAH(I)%ZOBSOIL
        END DO

!=== Assign SOIL Parameters to each tile based on the
!    type of Zobler soil class present in that tile.

        DO I=1,LDAS%NCH             !Tile loop
          K=PLACESLTYP(I)           !Soil type
          DO J=1,LDAS%NOAH_NSOILP   !Soil parameter loop
           NOAH(I)%SOILP(J)=BASICSET(K,J)
          ENDDO !J
        ENDDO !I

!=== FOLLOWING CODE IS ONLY TO REPLACE ABOVE POROSITY PARAMETER
!     WITH MULTIPLE LAYER POROSITY VALUES.
!         POR1A=AMIN1(POR1(TILE(I)%COL,TILE(I)%ROW),PORMAX)
!         POR2A=AMIN1(POR2(TILE(I)%COL,TILE(I)%ROW),PORMAX)
!         POR3A=AMIN1(POR3(TILE(I)%COL,TILE(I)%ROW),PORMAX)
!
!         NOAH(I)%SOILP(1) = POR1A * D1 * 1000.0
!         if (NOAH(I)%SOILP(1) .le. 0.0) then
!           write(*,*) 'COL,ROW,TILE,WSAT1',TILE(I)%COL,TILE(I)%ROW,
!     &      I,NOAH(I)%SOILP(1)
!           stop
!         end if
!
!        ENDDO !I
!---------------------------------------------------------------------
      ENDIF   !  End option for Reynolds Soil Types 

!=== NLDAS COSBY/RAWLS (Yun - 2m) Soil Parameterization ============== 
      IF (LDAS%SOIL.EQ.4) THEN

!+++   Open soil file.
        OPEN(14,FILE=LDAS%TXFILE,STATUS='OLD')

!+++    Read soil texture and map to Zobler soil classes.
!+++    At this time only the top layer data are used to assign 
!+++    soil parameters in grid space.    

        do k=1,11
         do j=1,224
          do i=1,464
           read (14,'(I3,1X,I3,1X,3X,I2)')
     &     ii,jj,soiltext(k,i,j)                                 
           if (soiltext(k,i,j).eq.13) soiltext(k,i,j)=12         
          enddo                                                 
         enddo                                                 
        enddo                                                 

        print *,'   '                                         
        print *,'SETTING SOIL CLASS 13 to 12...since have no params for
     & organic class'                                        
        print *,'WARNING!!!!!!!!!!!!!!!!!!'                   
        print *,'   '                      

!=== Map the soil classes to Zobler classes 

        DO J=1,LDAS%NR
         DO I=1,LDAS%NC
            IF (soiltext(1,I,J) .EQ. 1) SOILTYP(I,J) = 1
            IF (soiltext(1,I,J) .EQ. 2) SOILTYP(I,J) = 1
            IF (soiltext(1,I,J) .EQ. 3) SOILTYP(I,J) = 4
            IF (soiltext(1,I,J) .EQ. 4) SOILTYP(I,J) = 2
            IF (soiltext(1,I,J) .EQ. 5) SOILTYP(I,J) = 2
            IF (soiltext(1,I,J) .EQ. 6) SOILTYP(I,J) = 2
            IF (soiltext(1,I,J) .EQ. 7) SOILTYP(I,J) = 7
            IF (soiltext(1,I,J) .EQ. 8) SOILTYP(I,J) = 2
            IF (soiltext(1,I,J) .EQ. 9) SOILTYP(I,J) = 6
            IF (soiltext(1,I,J) .EQ. 10) SOILTYP(I,J) = 5       
            IF (soiltext(1,I,J) .EQ. 11) SOILTYP(I,J) = 5        
            IF (soiltext(1,I,J) .EQ. 12) SOILTYP(I,J) = 3        
            IF (soiltext(1,I,J) .EQ. 13) SOILTYP(I,J) = 8        
            IF (soiltext(1,I,J) .GT. 13) THEN             
              SOILTYP(I,J) = 2
!              WRITE(*,*) I, J                           
            END IF                                         
         END DO                                               
        END DO                                                

!=== Assign soil parameters in tile space.                    
                                                              
!        write(*,*) 'Assigning soil parameters in tile space.'
                                                              
!=== Read in the NOAH Soil Parameter File                     
                                                              
        OPEN(UNIT=18,FILE=LDAS%NOAH_SFILE,STATUS='OLD',       
     &   ACCESS='SEQUENTIAL')                                 
                                                              
        DO I=1,LDAS%NOAH_NSOILP                               
          READ(18,*)(BASICSET(JJ,I),JJ=1,LDAS%NOAH_ZST)       
        ENDDO                                                 
        CLOSE(18)                                             
                                                              
!=== Convert grid space to tile space for soil type values.   
        DO I=1,LDAS%NCH                                       
          PLACESLTYP(I)=SOILTYP(TILE(I)%COL,TILE(I)%ROW)      
          NOAH(I)%ZOBSOIL=PLACESLTYP(I)                       
        END DO                                                
                                                              
!=== Assign SOIL Parameters to each tile based on the         
!    type of Zobler soil class present in that tile.          
                                                              
        DO I=1,LDAS%NCH             !Tile loop                
          K=PLACESLTYP(I)           !Soil type                
          DO J=1,LDAS%NOAH_NSOILP   !Soil parameter loop      
                                                              
           NOAH(I)%SOILP(J)=BASICSET(K,J)                     
                                                              
          ENDDO !J                                            
        ENDDO !I       

      ENDIF  ! End of Yun's soil option  

!=== END OF SOIL OPTION 4 ================================= 

      RETURN
      END

