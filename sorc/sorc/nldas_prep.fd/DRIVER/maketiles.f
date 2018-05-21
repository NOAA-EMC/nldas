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
! maketiles.f: 
!
! DESCRIPTION:
!  This primary goal of this routine is to determine tile space.
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Major F90 and major structure revision
!  3  Jan 2000: Minor T=0 bug fix, should have no effect on output
!  8  Mar 2000: Brian Cosgrove; Initialized FGRD to 0 For Dec Alpha Runs
! 22  Aug 2000: Brian Cosgrove; Altered code for US/Mexico/Canada Mask
! 04  Feb 2001: Jon Gottschalck; Added option to read and use Koster tile space
!=========================================================================
      Subroutine maketiles(LDAS,GRID,TILE)
     
      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE tile_module      ! LDAS non-model-specific tile variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (tiledec) TILE(LDAS%NCH)
      type (griddec) GRID(LDAS%NC,LDAS%NR)

!=== Local Variables =====================================================
      INTEGER :: C,R,T,I,J     ! Loop counters
      INTEGER :: VEG(LDAS%NT)  !Temporary vegetation processing variable
      INTEGER :: WATER         !Temporary vegetation processing variable
      INTEGER :: ISUM          !Temporary vegetation processing variable
      INTEGER :: NKT = 14856
      INTEGER :: KVEG, J2, LANDNVEG
!      REAL :: VEG(LDAS%NT)
!      REAL :: WATER
!      REAL :: ISUM
!      REAL :: F1TEST
      REAL :: TPGRID
      REAL    :: RSUM          !Temporary vegetation processing variable
      REAL    :: FVT(LDAS%NT)  !Temporary vegetation processing variable
      REAL    :: MAX           !Temporary vegetation processing variable
      REAL    :: TSUM(LDAS%NC,LDAS%NR)  !Temporary processing variable
      INTEGER :: NCHP          !Number of tiles use for array size
      REAL    :: XLAT, XLON    !Latitude and Longitude for soil read in

!=== End Variable Definition =============================================
      NCHP = LDAS%NCH
      TSUM = 0.0

!=== Get the Land/water mask
      OPEN(30,file=LDAS%MFILE,FORM='formatted',STATUS='old')
      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
        READ(30,*) I,J,GRID(I,J)%LAT,GRID(I,J)%LON,GRID(I,J)%IMASK
        GRID(I,J)%MASK=FLOAT(GRID(I,J)%IMASK)
       ENDDO
      ENDDO
      CLOSE(30)

      OPEN(30,file=LDAS%FMFILE,FORM='formatted',STATUS='old')
      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
        READ(30,*) I,J,GRID(I,J)%LAT,GRID(I,J)%LON,GRID(I,J)%FIMASK
        GRID(I,J)%FMASK=FLOAT(GRID(I,J)%FIMASK)
       ENDDO
      ENDDO
      CLOSE(30)

!=== Read in Old-type Soils Data (values are non-varying).
      OPEN(40,FILE=LDAS%SFILE,form='formatted',STATUS='old')
      DO C=1,LDAS%NC
       DO R=1,LDAS%NR
        READ(40,*) I,J,XLAT,XLON,GRID(I,J)%SOILT
       ENDDO
      ENDDO
      CLOSE(40)

!=== Read in Vegetation Data
      DO R=1,LDAS%NR  !rows
       DO C=1,LDAS%NC  !columns
         ALLOCATE (GRID(C,R)%FGRD(LDAS%NT))
       ENDDO !R
      ENDDO !C 

!=== Select which tile-space veg. info to use (UMD or Koster)

      IF (LDAS%KOSTER .LT. 1) THEN             ! Use original UMD indexing
      OPEN(98,FILE=LDAS%VFILE,FORM='FORMATTED')
      DO R=1,LDAS%NR  !rows
       DO C=1,LDAS%NC  !columns
	READ(98,*)I,J,GRID(I,J)%LAT,GRID(I,J)%LON,ISUM,WATER,
     1   (VEG(T),T=1,LDAS%NT)
!         READ(98,*)I,J,GRID(I,J)%LAT,GRID(I,J)%LON,ISUM,F1TEST,WATER,
!     1   (VEG(T),T=1,LDAS%NT)
!NOTE: Current logic assumes all water or all land in a grid - this can
!  be easily changed in this part of the code if a lake model is included.
        ISUM=0  
        DO T=1,LDAS%NT 
         ISUM=ISUM+VEG(T)  !recompute ISUM without water points
        ENDDO
        DO T=1,LDAS%NT 
c brian added next line
         GRID(I,J)%FGRD(T)=0.0
         IF(ISUM.GT.0)GRID(I,J)%FGRD(T)=VEG(T)/FLOAT(ISUM)
!         IF(ISUM.GT.0)GRID(I,J)%FGRD(T)=VEG(T)/ISUM
        ENDDO !T
       ENDDO !R
      ENDDO !C 
      CLOSE(98)
      
      ELSE                               ! Use Koster tile space

       DO R=1,LDAS%NR  ! rows       
        DO C=1,LDAS%NC ! cols
         DO T=1,LDAS%NT ! veg types
           GRID(C,R)%FGRD(T)=0.0
           TSUM(C,R)=0.0
         ENDDO        
        ENDDO
       ENDDO
      
      OPEN(98,FILE=LDAS%KVFILE,FORM='FORMATTED')
      DO T=1,NKT  !no. of koster tiles
        READ(98,15)I,J,KVEG,TPGRID
        J2=J-15 ! to properly index latitude dimension (1-76)
        IF (KVEG .NE. 0) 
     &   GRID(I,J2)%FGRD(KVEG) = TPGRID
!        READ(98,15)I,J,KVEG,GRID(I,J)%FGRD(KVEG)
 15     FORMAT(6x,i3,10x,i2,18x,i1,19x,f8.5)
      ENDDO !T
      CLOSE(98)

      ENDIF

*** Exclude tiles with MINA (minimum tile grid area),  
***   normalize remaining tiles to 100%
      DO R=1,LDAS%NR  !rows
       DO C=1,LDAS%NC  !columns         

        RSUM=0.0
        DO T=1,LDAS%NT
         IF(GRID(C,R)%FGRD(T).LT.LDAS%MINA)THEN
          GRID(C,R)%FGRD(T)=0.0              ! impose area percent cutoff
         ENDIF
         RSUM=RSUM+GRID(C,R)%FGRD(T)
        ENDDO

* Renormalize veg fractions within a grid to 1
        IF(RSUM.GT.0.0) THEN  
         DO T=1,LDAS%NT  !Renormalize SUMT back to 1.0
          IF(RSUM.GT.0.0)GRID(C,R)%FGRD(T)=GRID(C,R)%FGRD(T)/RSUM
         ENDDO
 
         RSUM=0.0
         DO T=1,LDAS%NT
          RSUM=RSUM+GRID(C,R)%FGRD(T)  !Recalculate RSUM to check 
         ENDDO
 
         IF(RSUM.LT.0.9999.or.RSUM.GT.1.0001)THEN  !Check Renormalization
          WRITE(*,*) 'ERROR1 IN VEGETATION TILES',RSUM,C,R
          WRITE(79,*) 'ERROR1 IN VEGETATION TILES',RSUM,C,R
         ENDIF
        ENDIF

       ENDDO 
      ENDDO

*** Exclude tiles with MAXT (Maximum Tiles per grid), 
***   normalize remaining tiles to 100%
*** Determine the grid predominance order of the tiles
***  PVEG(NT) will contain the predominance order of tiles
      DO R=1,LDAS%NR  !rows
       DO C=1,LDAS%NC  !columns
        ALLOCATE (GRID(C,R)%PVEG(LDAS%NT))
       ENDDO
      ENDDO

      DO R=1,LDAS%NR  !rows
       DO C=1,LDAS%NC  !columns
        DO T=1,LDAS%NT
         FVT(T)=GRID(C,R)%FGRD(T)  !FVT= temp fgrd working array
         GRID(C,R)%PVEG(T)=0
        ENDDO
        DO I=1,LDAS%NT  !Loop through predominance level
         MAX=0.0
         T=0
         DO J=1,LDAS%NT
          IF(FVT(J).GT.MAX)THEN
           IF(GRID(C,R)%FGRD(J).GT.0) THEN
            MAX=FVT(J)
            T=J
           ENDIF
          ENDIF
         ENDDO
         IF(T.GT.0) THEN
          GRID(C,R)%PVEG(T)=I
          FVT(T)=-999.0       !eliminate chosen from next search 
         ENDIF
        ENDDO
       ENDDO !IR
      ENDDO !IC

*** Impose MAXT Cutoff
      DO R=1,LDAS%NR  !rows
       DO C=1,LDAS%NC  !columns         
        RSUM=0.0
        DO T=1,LDAS%NT
         IF(GRID(C,R)%PVEG(T).LT.1) THEN
          GRID(C,R)%FGRD(T)=0.0    
          GRID(C,R)%PVEG(T)=0  
         ENDIF
         IF(GRID(C,R)%PVEG(T).GT.LDAS%MAXT) THEN
          GRID(C,R)%FGRD(T)=0.0              ! impose MAXT cutoff
          GRID(C,R)%PVEG(T)=0  
         ENDIF
         RSUM=RSUM+GRID(C,R)%FGRD(T)
        ENDDO

* Renormalize veg fractions within a grid to 1
        IF(RSUM.GT.0.0) THEN  
         DO T=1,LDAS%NT  !Renormalize SUMT back to 1.0
          IF(RSUM.GT.0.0)GRID(C,R)%FGRD(T)=GRID(C,R)%FGRD(T)/RSUM
         ENDDO
 
         RSUM=0.0
         DO T=1,LDAS%NT
          RSUM=RSUM+GRID(C,R)%FGRD(T)  !Recalculate RSUM to check 
         ENDDO
         TSUM(C,R)=RSUM
 
         IF(RSUM.LT.0.9999.or.RSUM.GT.1.0001)THEN  !Check Renormalization
          WRITE(*,*) 'ERROR2 IN VEGETATION TILES',RSUM,C,R
          WRITE(79,*) 'ERROR2 IN VEGETATION TILES',RSUM,C,R
         ENDIF
        ENDIF

       ENDDO 
      ENDDO

*** Account for if land but no vegetation assigned, use mixed type for UMD, grassland for Koster
      LANDNVEG = 5 ! UMD specific, changed below if using Koster tilespace
      IF (LDAS%KOSTER .EQ. 1) LANDNVEG = 4

*** Make Tile Space
      LDAS%NCH=0
      DO T=1,LDAS%NT  !loop through each tile type
       DO R=1,LDAS%NR  !loop through rows
         DO C=1,LDAS%NC  !loop through columns
         IF(GRID(C,R)%MASK.gt.0.99.and.
     1       GRID(C,R)%MASK.lt.3.01)THEN  !we have land
             IF(GRID(C,R)%FGRD(T).GT.0.0)THEN
            LDAS%NCH=LDAS%NCH+1  !count the number of tiles (or chips)
            TILE(LDAS%NCH)%ROW=R    !keep track of tile row
            TILE(LDAS%NCH)%COL=C    !keep track of tile column
            TILE(LDAS%NCH)%VEGT=T    !Keep track of tile surface type
            TILE(LDAS%NCH)%FGRD=GRID(C,R)%FGRD(T)   !keep track of tile fraction
            TILE(LDAS%NCH)%LON=GRID(C,R)%LON       !Longitude of tile
            TILE(LDAS%NCH)%LAT=GRID(C,R)%LAT       !latitude of tile
            TILE(LDAS%NCH)%SOILT=GRID(C,R)%SOILT   !Soil Type of Tile
            TILE(LDAS%NCH)%PVEG=GRID(C,R)%PVEG(T)  !Predominance of vegetation class in grid
          ENDIF

** What if we we have land without vegetation assigned! Use mixed type
          IF(TSUM(C,R).eq.0.0.and.T.eq.LANDNVEG)THEN  !VEGETATION TYPE SPECIFIC
           LDAS%NCH=LDAS%NCH+1  !count the number of tiles (or chips)
           TILE(LDAS%NCH)%ROW=R    !keep track of tile row
           TILE(LDAS%NCH)%COL=C    !keep track of tile column
           TILE(LDAS%NCH)%VEGT=T    !Keep track of tile surface type
           TILE(LDAS%NCH)%FGRD=1.0   !keep track of tile fraction
           TILE(LDAS%NCH)%LON=GRID(C,R)%LON       !Longitude of tile
           TILE(LDAS%NCH)%LAT=GRID(C,R)%LAT       !latitude of tile
           TILE(LDAS%NCH)%SOILT=GRID(C,R)%SOILT   !Soil Type of Tile
           TILE(LDAS%NCH)%PVEG=GRID(C,R)%PVEG(T)  !Predominance of vegetation class in grid
         ENDIF
        ENDIF
        ENDDO
       ENDDO
      ENDDO

      WRITE(*,*) 'Size of Mosaic Tile Dimension:',NCHP
      WRITE(*,*) 'Actual Number of Mosaic Tiles:',LDAS%NCH
      WRITE(*,*)

      WRITE(79,*) 'Size of Mosaic Tile Dimension:',NCHP
      WRITE(79,*) 'Actual Number of Mosaic Tiles:',LDAS%NCH
      WRITE(79,*)

      RETURN
      END




























