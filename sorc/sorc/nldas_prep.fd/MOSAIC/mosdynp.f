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
! mosdynp.f: 
!
! DESCRIPTION:
!  This subroutine take all the monthly varying parameters
!  and the date and determine the actual value of the parameter for that date
!  this actual value is returned to the main program
!  The assumption is that the data point is valid for the 16th
!  of the given month at 00hr
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  8  Mar 2000: Brian Cosgrove; Added Integer Number Holders For Dec Alpha Runs
!  6  Apr 2001: Matt Rodell; Assign veg parameters based on N/S hemisphere
!  11 Feb 2002: Jon Gottschalck; Added use of AVHRR derived LAI/Greenness
!  01 Oct 2002: Jon Gottschalck; Modified to allow for MODIS LAI      
!=========================================================================

      SUBROUTINE MOSDYNP(LDAS,MOS,TILE)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE mos_module       ! Mosaic tile variables
      USE tile_module	   ! LDAS non-model-specific tile variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (mosdec)  MOS(LDAS%NCH)
      type (tiledec) TILE(LDAS%NCH)

!=== Local Variables =====================================================
      INTEGER :: P,T                 ! Loop counters
      REAL*8  :: TIME1,TIME2         ! Temporary Time variables
      INTEGER :: YR1,MO1,YR2,MO2     ! Temporary Time variables
      INTEGER :: NHMO1,NHMO2	     ! Temp var - N. Hemisphere equiv month
      INTEGER :: DOY1,DOY2           ! Temporary Time variables
      REAL    :: WT1,WT2,GMT1,GMT2   ! Interpolation weights
      INTEGER :: ZEROI,NUMI          ! Integer Number Holders
!=== End Variable Definition =============================================
	  
C	Initialize Numbers
	ZEROI=0
	NUMI=16
!=== Determine Monthly data Times (Assume Monthly value valid at DA=16 HR=00Z)
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
      WT1= (TIME2-LDAS%TIME)/(TIME2-TIME1)
      WT2= (LDAS%TIME-TIME1)/(TIME2-TIME1)

      IF (LDAS%LAI .GE. 2) CALL MOSLAIREAD(LDAS,MOS,TILE,YR1,MO1,
     &                           YR2,MO2,TIME1,TIME2,WT1,WT2)

      DO T=1,LDAS%NCH

!      If tile is in the southern hemisphere, use veg parameters from the 
!	opposite (6 months away) time of year.
       IF (TILE(T)%LAT .lt. 0.0) then
	NHMO1 = MOD((MO1+6),12)
	IF (NHMO1 .EQ. 0) NHMO1=12
	NHMO2 = MOD((MO2+6),12)
	IF (NHMO2 .EQ. 0) NHMO2=12
       ELSE
	NHMO1 = MO1
	NHMO2 = MO2
       END IF

       IF (LDAS%LAI .GE. 2) THEN

        MOS(T)%VEGIP(1) = MOS(T)%GREEN
        MOS(T)%VEGIP(2) = MOS(T)%LAI
        IF (TILE(T)%VEGT .EQ. 12) THEN
         MOS(T)%VEGIP(1) = 0.001
         MOS(T)%VEGIP(2) = 0.001
        ENDIF

        DO P=3,LDAS%MOS_NMVEGP
          MOS(T)%VEGIP(P) = WT1*MOS(T)%VEGMP(P,NHMO1) +
     &                      WT2*MOS(T)%VEGMP(P,NHMO2)
        ENDDO

      ELSE

       DO P=1,LDAS%MOS_NMVEGP
        MOS(T)%VEGIP(P) = WT1*MOS(T)%VEGMP(P,NHMO1) + 
     &	 WT2*MOS(T)%VEGMP(P,NHMO2)
        ENDDO
	 
      ENDIF

      ENDDO
    
      RETURN
      END
 
