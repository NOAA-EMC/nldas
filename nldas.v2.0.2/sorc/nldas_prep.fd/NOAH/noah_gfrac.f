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
! noah_gfrac.f: 
!
! DESCRIPTION:
!  This subroutine takes vegetation greenness fraction data and the date to 
!  interpolate and determine the actual value of the greenness fraction 
!  for that date.  This actual value is then returned to the main
!  program.  The assumption is that the data point is valid for the 16th
!  of the given month, at 00Z.
!
! REVISION HISTORY:
!  28 Apr 2002: K. Arsenault; Added NOAH LSM to LDAS, initial code
!  04 Nov 2002: K. Arsenault; Added 0.5 and 1.0 deg resolution fields
!  24 Jun 2003: K. Arsenault; GFRAC flags added to update once daily
!=========================================================================

      SUBROUTINE NOAH_GFRAC(LDAS,NOAH,TILE)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE noah_module      ! NOAH tile variables
      USE tile_module	   ! LDAS non-model-specific tile variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (noahdec) NOAH(LDAS%NCH)
      type (tiledec) TILE(LDAS%NCH)

!=== Local Variables =====================================================
      INTEGER :: I,J,P,T              ! Loop counters
      INTEGER :: LINE1,LINE2          ! Read file counters  
      REAL*8 :: TIME1,TIME2           ! Temporary Time variables
      INTEGER :: YR1,MO1,YR2,MO2      ! Temporary Time variables
      INTEGER :: DOY1,DOY2            ! Temporary Time variables
      INTEGER :: ZEROI,NUMI           ! Integer Number Holders
      INTEGER :: GFRAC_FLAG           ! Flag to update gfrac 
      REAL :: WT1,WT2,GMT1,GMT2       ! Interpolation weights
      REAL :: VALUE1(LDAS%NC,LDAS%NR) ! Temporary value holder for MO1
      REAL :: VALUE2(LDAS%NC,LDAS%NR) ! Temporary value holder for MO2
      REAL :: VEGMP1(LDAS%NCH)        ! Month 1 GFRAC in tile space
      REAL :: VEGMP2(LDAS%NCH)        ! Month 2 GFRAC in tile space

      CHARACTER*2 :: MM1,MM2        ! Filename places for integer. MO1, MO2

!=== End Variable Definition =============================================
	  
C   Initialize Numbers
      ZEROI=0
      NUMI=16

!=== Determine Monthly data Times (Assume Monthly value valid at DA=16)
      IF(LDAS%DA.LT.16)THEN
        MO1=LDAS%MO-1
        YR1=LDAS%YR 
        IF(MO1.EQ.0)THEN
          MO1=12
          YR1=LDAS%YR-1
        ENDIF
        MO2=LDAS%MO
        YR2=LDAS%YR
      ELSE
        MO1=LDAS%MO
        YR1=LDAS%YR
        MO2=LDAS%MO+1
        YR2=LDAS%YR
        IF(MO2.EQ.13)THEN
          MO2=1
          YR2=LDAS%YR+1
        ENDIF
      ENDIF

      CALL DATE2TIME(TIME1,DOY1,GMT1,YR1,MO1,
     &  NUMI,ZEROI,ZEROI,ZEROI)
      CALL DATE2TIME(TIME2,DOY2,GMT2,YR2,MO2,
     &  NUMI,ZEROI,ZEROI,ZEROI)

!--  Weights to be used to interpolate greenness fraction values.  
       WT1= (TIME2-LDAS%TIME)/(TIME2-TIME1)
       WT2= (LDAS%TIME-TIME1)/(TIME2-TIME1)

!--  Determine if GFRAC files need to be updated
      IF (TIME2 .GT. LDAS%NOAH_GFRACTIME) THEN
        GFRAC_FLAG = 1
      ELSE
        GFRAC_FLAG = 0
      ENDIF

      IF (GFRAC_FLAG .EQ. 1) THEN   !Open needed gfrac files
       LDAS%NOAH_GFRACTIME = TIME2

!=== Open greenness fraction dataset of months corresponding to   
!    time1 and time2 for selected LDAS domain and read data.

       WRITE(MM1,3) MO1
       WRITE(MM2,3) MO2
 3     FORMAT(I2.2)

       SELECT CASE (LDAS%DOMAIN)

       CASE (1)    ! NLDAS (0.125)
        OPEN (10, FILE='BCS/NOAH/gfrac_'//MM1//'_0.125.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4) 
        OPEN (11, FILE='BCS/NOAH/gfrac_'//MM2//'_0.125.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4) 

       CASE (2)    ! GLDAS (0.25)
        OPEN (10, FILE='BCS/NOAH/gfrac_'//MM1//'_0.25.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4)
        OPEN (11, FILE='BCS/NOAH/gfrac_'//MM2//'_0.25.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4) 

       CASE (3)    ! GLDAS (2x2.5)
        OPEN (10, FILE='BCS/NOAH/gfrac_'//MM1//'_2-2.5.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4)
        OPEN (11, FILE='BCS/NOAH/gfrac_'//MM2//'_2-2.5.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4)

       CASE (4)    ! GLDAS (1.0)
        OPEN (10, FILE='BCS/NOAH/gfrac_'//MM1//'_1.0.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4)
        OPEN (11, FILE='BCS/NOAH/gfrac_'//MM2//'_1.0.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4)

       CASE (5)    ! GLDAS (0.5)
        OPEN (10, FILE='BCS/NOAH/gfrac_'//MM1//'_0.5.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4)
        OPEN (11, FILE='BCS/NOAH/gfrac_'//MM2//'_0.5.bin',
     &    STATUS='OLD', ACCESS='DIRECT', RECL=4)

       END SELECT

         LINE1=0
         LINE2=0
         DO J=1,LDAS%NR
          DO I=1,LDAS%NC
            LINE1=LINE1+1
            LINE2=LINE2+1
             READ(10,REC=LINE1) VALUE1(I,J)
             READ(11,REC=LINE2) VALUE2(I,J)
          ENDDO !I
         ENDDO !J
        CLOSE(10)
        CLOSE(11)
     
!=== Assign MONTHLY vegetation greenness fractions to each tile
       DO I=1,LDAS%NCH          !Tile loop
        IF((VALUE1(TILE(I)%COL,TILE(I)%ROW).NE.-9999.000)
     &  .AND.(VALUE2(TILE(I)%COL,TILE(I)%ROW).NE.-9999.000))THEN
          NOAH(I)%VEGMP1=VALUE1(TILE(I)%COL,TILE(I)%ROW)
          NOAH(I)%VEGMP2=VALUE2(TILE(I)%COL,TILE(I)%ROW)
        ENDIF
       END DO

      ENDIF   ! End of gfrac file flag 

!===  Interpolate greenness fraction values once daily 

      IF (LDAS%NOAH_GFRACDCHK .NE. LDAS%DA) THEN

       DO I=1,LDAS%NCH          !Tile loop
         NOAH(I)%VEGIP = (WT1*NOAH(I)%VEGMP1)+(WT2*NOAH(I)%VEGMP2)
       END DO     ! End tile DO Loop	 

       LDAS%NOAH_GFRACDCHK = LDAS%DA

      ENDIF  ! End daily interpolation
    
      RETURN
      END

