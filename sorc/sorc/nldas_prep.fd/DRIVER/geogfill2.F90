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
! GEOGFILL2.F90: 
!
! DESCRIPTION:
!  Fill in grid points that are invalid in LDAS due to differences in
!  geography between forcing and land surface.
!
!  Based on original geogfill.F90 for use with Reanalysis forcing 
!  data, to account for differences in land masks.  
!  NOTE** for Reanalysis data, v=6 is the V-component of wind,
!  which is always set to ZERO since the U-component (v=5) is assigned 
!  the absolute value of wind speed.
!
!  For v=1:temperature, v=2:specific humidity, v=4:LW radiation,
!  v=5:wind, and v=7:pressure, data values of zero are not allowed to 
!  contribute to average fill-value.
!  For v=3:SW radiation, v=8,9:precipitation, use the mask defined by 
!  valid temperature data points to establish whether to seek fill-value
!  for given column and row, rather than include vs exclude zeroes.
!
! REVISION HISTORY:
!  20 Sep 2002: Urszula Jambor; Modified original geogfill routine
!               for use with Reanalysis forcing sets prepared by
!               Aaron Berg via NSIPP
!  29 Jan 2003: Urszula Jambor; Removed LDAS & GRID modules from list of
!               arguments, instead pass only needed variables directly
!               (nc,nr,fimask).
!=========================================================================

      SUBROUTINE GEOGFILL2(NC,NR,FIMASK,GEOGMASK,DATA,V,TMASK)

      IMPLICIT NONE

!=== Local Variables =====================================================

      INTEGER :: I,J,KK,KT,II,JJ,V           ! Loop counters
      INTEGER :: NC,NR,FIMASK(NC,NR)         ! Dimensions and Forcing Mask
      INTEGER :: SMASK(NC,NR)                ! Grid points where LDAS
                                             ! defines land but
                                             ! forcing bitmap is.FALSE.
      REAL :: DATA(NC,NR)                    ! Interpolated Forcing field
      REAL :: SSUM,SAVG                      ! Sum, average--fill in points
      INTEGER :: SCNT                        ! Counter--fill in points
      LOGICAL*1 :: TMASK(NC,NR)              ! Valid temperature points
      LOGICAL*1 :: GEOGMASK(NC,NR)           ! Interpolated Forcing bitmap

!=== End Variable Definition =============================================

!=== Initializing sub-surface parameter bitmap ===========================
     SMASK = 0

!=== Defining sub-surface parameter bitmap ===============================
      DO J=1,NR
        DO I=1,NC
          IF (FIMASK(I,J) .GT. 0 .AND. &
               .NOT.(GEOGMASK(I,J))          ) THEN
            SMASK(I,J) = 1
          ENDIF
        ENDDO
      ENDDO
!      print*, 'number of smask points is ', v, sum(smask)

!=== Filling in points where sub-surface parameter bitmap is equal to 1
      DO J = 1, NR
        DO I = 1, NC
          SSUM = 0
          SCNT = 0
          IF (SMASK(I,J) .GT. 0) THEN              ! Problem land point
            IF (J .LE. 2 .OR. J .GE. NR-1) THEN
              DATA(I,J) = DATA(I,J)
            ELSEIF (I .LE. 2) THEN                 ! First two columns
              IF(I .EQ. 1) THEN
                KT = NC-2
              ELSE
                KT = NC-1
              ENDIF
              DO KK = 1, 5
                KT = KT + 1
                DO JJ = J-2, J+2
                  IF ( (v.ne.3) .and. (v.lt.8) ) THEN
                    IF ((DATA(KT,JJ).GT.0.0) .AND. DATA(KT,JJ).NE.-1) THEN
                      SSUM = SSUM + DATA(KT,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ELSE IF ( TMASK(I,J) .AND. DATA(KT,JJ).NE.-1) THEN
                    IF (SMASK(KT,JJ).EQ.0) THEN
                      SSUM = SSUM + DATA(KT,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ENDIF
                ENDDO
                IF (KT .EQ. NC) KT = 0
              ENDDO
              IF (SCNT .NE. 0) THEN
                SAVG = SSUM / SCNT      ! Average of surrounding points
              ELSE
                SAVG = -1               ! Set flag for later
              ENDIF
              DATA(I,J) = SAVG                 ! Assign fill-in data point
            ELSEIF (I .GE. NC-1) THEN     ! Last two columns
              KT = I-3
              DO KK = 1, 5
                KT = KT + 1
                DO JJ = J-2, J+2
                  IF ( (v.ne.3) .and. (v.lt.8) ) THEN
                    IF ((DATA(KT,JJ).GT.0.0) .AND. DATA(KT,JJ).NE.-1) THEN
                      SSUM = SSUM + DATA(KT,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ELSE IF ( TMASK(I,J) .AND. DATA(KT,JJ).NE.-1 ) THEN
                    IF (SMASK(KT,JJ).EQ.0) THEN
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
            ELSE                                   ! All other points
              DO II = I-2, I+2
                DO JJ = J-2, J+2
                  IF ( (v.ne.3) .and. (v.lt.8) ) THEN
                    IF ((DATA(II,JJ).GT.0.0) .AND. DATA(II,JJ).NE.-1) THEN
                      SSUM = SSUM + DATA(II,JJ)
                      SCNT = SCNT + 1
                    ENDIF
                  ELSE IF ( TMASK(I,J) .AND. DATA(II,JJ).NE.-1 ) THEN
                    IF (SMASK(II,JJ).EQ.0) THEN
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

      IF (v==1) THEN ! Assign TMASK for future use
        DO J=1,NR
          DO I=1,NC
            IF (DATA(I,J) > 0.0) THEN
              tmask(i,j) = .true.
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
      
