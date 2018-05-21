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
! noah_alb.f: 
!
! DESCRIPTION:
!  This subroutine takes quarterly surface albedo (snow-free) data and  
!  day to interpolate and determine the actual value of the albedo 
!  for that date.  This actual value is then returned to the main
!  program.  The assumption is that the data point is valid for the 
!  dates of January 31, April 30, July 31, and October 31.
!
! REVISION HISTORY:
!  28 Apr 2002: K. Arsenault; Added NOAH LSM to LDAS, initial code
!  04 Nov 2002: K. Arsenault; Added 1.0 and 0.5 deg GLDAS fields
!  07 Jan 2003: K. Arsenault; Make corrections to time dimensions 
!                             and variables
!  20 Jan 2003: K. Arsenault; Removed a call to date2time subroutine 
!  24 Jun 2003: K. Arsenault; ALB flags added to interpolate once daily
!  16 Oct 2003: K. Arsenault; Removed explicit leap year features
!=========================================================================

      SUBROUTINE NOAH_ALB(LDAS,NOAH,TILE)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE noah_module      ! NOAH tile variables
      USE tile_module      ! LDAS non-model-specific tile variables
      IMPLICIT NONE
      type (ldasdec) LDAS              
      type (noahdec) NOAH(LDAS%NCH)
      type (tiledec) TILE(LDAS%NCH)

!=== Local Variables =====================================================
      INTEGER :: I,J,T                ! Loop counters
      INTEGER :: LINE1,LINE2          ! Read file counters
      INTEGER :: JANDA,JANMO          ! January 31 
      INTEGER :: APRDA,APRMO          ! April 30
      INTEGER :: JULDA,JULMO          ! July 31
      INTEGER :: OCTDA,OCTMO          ! October 31 
      INTEGER :: LPDA,LPMO,DEC31      ! Used to determine last day of year 
      INTEGER :: YR                   ! Year of run  
      INTEGER :: DOY1                 ! Temporary Time variables
      INTEGER :: ZEROI                ! Integer Number Holders
      INTEGER :: ALBFLAG              ! Flag to update albedo files

      REAL*8 :: TIME                  ! Current Model Time variable
      REAL*8 :: JAN31,APR30,JUL31,OCT31 ! Dates of quarterly albedo files
      REAL*8 :: LP                    ! Used to account for end of year
      REAL*8 :: QDIF                  ! Difference between Q1 and Q2 times
      REAL*8 :: TIMDIF                ! Difference between TIME and Q1 time
      REAL :: GMT1,GMT2               ! Interpolation weights
      REAL :: VALUE1(LDAS%NC,LDAS%NR) ! Temporary value holder for QQ1
      REAL :: VALUE2(LDAS%NC,LDAS%NR) ! Temporary value holder for QQ2
      REAL :: VALDIF(LDAS%NCH)        ! Difference of QQ2 and QQ1 albedo

      CHARACTER*2 :: QQ1,QQ2          ! Filename places for quarter values 

!=== End Variable Definition =============================================
  
!   Initialize Numbers
      ZEROI=0

!=== Determine Dates of the quarters in terms of Year (e.g., 1999.3) 
      TIME=LDAS%TIME
      YR=LDAS%YR

!  January 31
      JANDA=31
      JANMO=01
      CALL DATE2TIME(JAN31,DOY1,GMT1,YR,JANMO,
     &  JANDA,ZEROI,ZEROI,ZEROI)

!  April 30
      APRDA=30
      APRMO=04
      CALL DATE2TIME(APR30,DOY1,GMT1,YR,APRMO,
     &  APRDA,ZEROI,ZEROI,ZEROI)

!  July 31
      JULDA=31
      JULMO=07
      CALL DATE2TIME(JUL31,DOY1,GMT1,YR,JULMO,
     &  JULDA,ZEROI,ZEROI,ZEROI)

!  October 31 
      OCTDA=31
      OCTMO=10
      CALL DATE2TIME(OCT31,DOY1,GMT1,YR,OCTMO,
     &  OCTDA,ZEROI,ZEROI,ZEROI)

!--- Determine which two quarterly albedo files book-end model time.

      IF ( TIME.GE.JAN31 .and. TIME.LE.APR30 ) THEN
         QQ1="01"
         QQ2="02"
         QDIF = APR30-Jan31
         TIMDIF = TIME-JAN31
          print *, 'JAN31=',JAN31,' APR30=',APR30
         ALBFLAG=1 

      ELSEIF ( TIME.GE.APR30 .and. TIME.LE.JUL31 ) THEN
         QQ1="02"
         QQ2="03"
         QDIF = JUL31-APR30
         TIMDIF = TIME-APR30
!          print *, 'APR30=',APR30,' JUL31=',JUL31
         ALBFLAG=2
 
      ELSEIF ( TIME.GE.JUL31 .and. TIME.LE.OCT31 ) THEN
         QQ1="03"
         QQ2="04"
         QDIF = OCT31-JUL31
         TIMDIF = TIME-JUL31
!          print *, 'JUL31=',JUL31,' OCT31=',OCT31
         ALBFLAG=3

      ELSEIF ( TIME.GE.OCT31 ) THEN
         QQ1="04"
         QQ2="01"
         QDIF = (JAN31+1.0)-OCT31
         TIMDIF = TIME-OCT31
!          print *, 'OCT31=',OCT31,' JAN31=',JAN31
         ALBFLAG=4

      ELSEIF ( TIME.LT.JAN31 ) THEN
         QQ1="04"
         QQ2="01"
         OCT31=OCT31-1.0
         QDIF = JAN31-OCT31
         TIMDIF = TIME-OCT31
!          print *, 'OCT31=',OCT31,' JAN31=',JAN31
         ALBFLAG=4

      ENDIF

!===  Open the needed two quarterly snow-free albedo files   
      IF (LDAS%NOAH_ALBTIME .NE. ALBFLAG) THEN
        LDAS%NOAH_ALBTIME = ALBFLAG

       SELECT CASE (LDAS%DOMAIN)

       CASE (1)    ! NLDAS (0.125 deg)  
        OPEN (10, FILE='BCS/NOAH/albedo_'//QQ1//'_0.125.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4) 
        OPEN (11, FILE='BCS/NOAH/albedo_'//QQ2//'_0.125.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4) 

       CASE (2)    ! GLDAS (0.25 deg) 
         OPEN (10, FILE='BCS/NOAH/albedo_'//QQ1//'_0.25.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4) 
         OPEN (11, FILE='BCS/NOAH/albedo_'//QQ2//'_0.25.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4) 

       CASE (3)    ! GLDAS (2x2.5 deg) 
         OPEN (10, FILE='BCS/NOAH/albedo_'//QQ1//'_2-2.5.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4) 
         OPEN (11, FILE='BCS/NOAH/albedo_'//QQ2//'_2-2.5.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4)

       CASE (4)    ! GLDAS (1.0 deg)
         OPEN (10, FILE='BCS/NOAH/albedo_'//QQ1//'_1.0.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4)
         OPEN (11, FILE='BCS/NOAH/albedo_'//QQ2//'_1.0.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4)

       CASE (5)    ! GLDAS (0.5 deg)
         OPEN (10, FILE='BCS/NOAH/albedo_'//QQ1//'_0.5.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4)
         OPEN (11, FILE='BCS/NOAH/albedo_'//QQ2//'_0.5.bin',
     &    STATUS='OLD',ACCESS='DIRECT', RECL=4)

       END SELECT   ! End domain selection 

!=== Read albedo dataset of quarters corresponding to   
!    bookends of TIME for selected LDAS domain.

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

!=== Assign albedo fractions to each tile and interpolate daily.
       DO I=1,LDAS%NCH          !Tile loop
        IF((VALUE1(TILE(I)%COL,TILE(I)%ROW).NE.-9999.000)
     &  .AND.(VALUE2(TILE(I)%COL,TILE(I)%ROW).NE.-9999.000))THEN
          NOAH(I)%ALBSF1(1) = VALUE1(TILE(I)%COL,TILE(I)%ROW)  
          NOAH(I)%ALBSF2(1) = VALUE2(TILE(I)%COL,TILE(I)%ROW)
        ENDIF 
       END DO
       print *, ' Finished opening new albedo files '

      ENDIF   ! End albflag selection

!=== Assign albedo fractions to each tile and interpolate daily. 

      IF (LDAS%NOAH_ALBDCHK .NE. LDAS%DA) THEN  
          print *, 'ALBDCHCK, DA: ',LDAS%NOAH_ALBDCHK,LDAS%DA

       DO I=1,LDAS%NCH          
         VALDIF(I)=NOAH(I)%ALBSF2(1)-NOAH(I)%ALBSF1(1)
         NOAH(I)%ALBSF = (TIMDIF*VALDIF(I)/QDIF)+NOAH(I)%ALBSF1(1)

        if (i.eq.2000) then
         print*,' (2) ALBSF1 & ALBSF2(i=2000):',NOAH(i)%ALBSF1,
     &           NOAH(i)%ALBSF2
         print *, ' (2) VALDIF, ALBSF = ', VALDIF(I),noah(i)%albsf
         print *, '     qdif,timdif     = ',qdif,timdif
        endif

       END DO

       LDAS%NOAH_ALBDCHK=LDAS%DA
 
      ENDIF   ! End daily interpolation

      RETURN
      END
