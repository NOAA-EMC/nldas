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
! GEOGFILL.f: 
!
! DESCRIPTION:
!  Fill in grid points that are invalid in LDAS due to differences in
!  geography between forcing and land surface (sub surface parameters).
!
!
! REVISION HISTORY:
!  12 Apr 2001: Jon Gottschalck; Initial Code
!  17 Apr 2002: Urszula Jambor; Modified slightly for use with the 
!               0.5deg reanal.ECMWF data set.
!  20 Sep 2002: Urszula Jambor; Reverted to March 1, 2002 version of 
!               geogfill.F90 and created geogfill2.F90, which has
!               Reanalysis forcing data specifically in mind.
!  13 Nov 2002: Urszula Jambor; Switch in parameter numbers
!               for snow and soil moisture.  Snow is now 12, and 
!               soil moisture is now 13.
!  29 Jan 2003: Urszula Jambor; Removed LDAS & GRID modules from list of
!               arguments, instead pass only needed variables directly
!               (nc,nr,fimask).
!=========================================================================

      SUBROUTINE GEOGFILL(NC,NR,FIMASK,GEOGMASK,DATA,V)

      IMPLICIT NONE

!=== Local Variables =====================================================

      INTEGER :: I,J,KK,KT,II,JJ,V           ! Loop counters
      INTEGER :: NC, NR, FIMASK(NC,NR)       ! Dimensions and Forcing Mask
      INTEGER :: SMASK(NC,NR)                ! Bitmap of points that are LDAS land 
                                             ! but have a forcing bitmap of
                                             ! .FALSE. (sub-surface parameter bitmap)
      REAL :: DATA(NC,NR)                    ! GLDAS field
      REAL :: SSUM,SAVG                      ! Sum and average, used to fill in points
      INTEGER :: SCNT                        ! Count, used to fill in points
      LOGICAL*1 :: GEOGMASK(NC,NR) ! Forcing bitmap

!=== End Variable Definition =============================================

!=== Initializing sub-surface parameter bitmap ===========================
     DO J=1,NR
        DO I=1,NC
	  SMASK(I,J) = 0
	ENDDO
     ENDDO

!=== Defining sub-surface parameter bitmap ===============================
      DO J=1,NR
        DO I=1,NC
          IF (FIMASK(I,J) .GT. 0 .AND. &
               .NOT.(GEOGMASK(I,J))          ) THEN
            SMASK(I,J) = 1
          ENDIF
        ENDDO
      ENDDO

!=== Filling in points where sub-surface parameter bitmap is equal to 1
      DO J = 1, NR
        DO I = 1, NC
          SSUM = 0
          SCNT = 0
          IF (SMASK(I,J) .GT. 0) THEN             ! Is it one of these problem land points
            IF (J .LE. 2 .OR. J .GE. NR-1) THEN
              DATA(I,J) = DATA(I,J)
            ELSEIF (I .LE. 2) THEN                 ! Points for first two columns
              IF(I .EQ. 1) THEN
                KT = NC-2
              ELSE
                KT = NC-1
              ENDIF
              DO KK = 1, 5
                KT = KT + 1
                DO JJ = J-2, J+2
                  IF(V .EQ. 13) THEN    ! Soil moisture/temp, do not allow 0's
                    IF (DATA(KT,JJ).NE.0 .AND. DATA(KT,JJ).NE.-1) THEN
                      SSUM = SSUM + DATA(KT,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ELSE                                                     ! Snow (0 may be valid data)
                    IF (DATA(KT,JJ).NE.-1 .AND. SMASK(KT,JJ).NE.1) THEN
                      SSUM = SSUM + DATA(KT,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ENDIF
                ENDDO
                IF (KT .EQ. NC) KT = 0
              ENDDO
              IF (SCNT .NE. 0) THEN
                SAVG = SSUM / SCNT     ! Set data point to average of points (1/2 degree box)
              ELSE
                SAVG = -1               ! Set to a flag to be acted on later
              ENDIF
              DATA(I,J) = SAVG                 ! Assign returned array to average
            ELSEIF (I .GE. NC-1) THEN           ! Points for last two columns
              KT = I-3
              DO KK = 1, 5
                KT = KT + 1
                DO JJ = J-2, J+2
                  IF(V .EQ. 13) THEN  ! Soil moisture/temp, do not allow 0's
                    IF (DATA(KT,JJ).NE.0 .AND. DATA(KT,JJ).NE.-1) THEN
                      SSUM = SSUM + DATA(KT,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ELSE                                                  ! Snow (0 may be valid data)
                    IF (DATA(KT,JJ).NE.-1 .AND. SMASK(KT,JJ).NE.1) THEN
                      SSUM = SSUM + DATA(KT,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ENDIF
                ENDDO
                IF (KT .EQ. NC) KT = 0
              ENDDO
              IF (SCNT .NE. 0) THEN
                SAVG = SSUM / SCNT
              ELSE
                SAVG = -1
              ENDIF
              DATA(I,J) = SAVG
            ELSE                           ! All other points
              DO II = I-2, I+2
                DO JJ = J-2, J+2
                  IF(V .EQ. 13) THEN  ! Soil moisture, do not allow 0's
                    IF (DATA(II,JJ).NE.0 .AND. DATA(II,JJ).NE.-1) THEN
                      SSUM = SSUM + DATA(II,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ELSE                                                  ! Snow (0 may be valid data)
                    IF (DATA(II,JJ).NE.-1 .AND. SMASK(II,JJ).NE.1) THEN
                      SSUM = SSUM + DATA(II,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
              IF (SCNT .NE. 0) THEN
                SAVG = SSUM / SCNT
              ELSE
                SAVG = -1
              ENDIF
              DATA(I,J) = SAVG
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
      
