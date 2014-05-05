!-------------------------------------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
! Released October 2005
!
! See SOFTWARE DISTRIBUTION POLICY for software distribution policies
!
! The LIS source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See COPYRIGHT.TXT for copyright details.
!
!-------------------------------------------------------------------------
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

subroutine mosdynp(ld,lt,lp)
      
  use lisdrv_module, only : tile, grid, gindex
  use lis_module      ! ldas non-model-specific 1-d variables
  use mos_varder      ! noah tile variables
  use time_manager
  use tile_module
  implicit none
  type (lisdomain) ld
  type (listime) LT
  type (lisparameters) LP
      
!=== Local Variables =====================================================
      integer :: index,i,j,k
      INTEGER :: P,T                 ! Loop counters
      REAL*8  :: TIME1,TIME2         ! Temporary Time variables
      INTEGER :: YR1,MO1,YR2,MO2     ! Temporary Time variables
      INTEGER :: NHMO1,NHMO2	     ! Temp var - N. Hemisphere equiv month
      INTEGER :: DOY1,DOY2           ! Temporary Time variables
      REAL    :: WT1,WT2,GMT1,GMT2   ! Interpolation weights
      INTEGER :: ZEROI,NUMI          ! Integer Number Holders
      REAL :: VALUEMON(LP%NT,MOSDRV%MOS_NMVEGP,12)
      REAL :: VEGMP(LD%NCH,MOSDRV%MOS_NMVEGP,12)   
!=== End Variable Definition =============================================
!	Initialize Numbers
	ZEROI=0
	NUMI=16
!=== Determine Monthly data Times (Assume Monthly value valid at DA=16 HR=00Z)
      IF(LT%DA.LT.16)THEN
       MO1=LT%MO-1
       YR1=LT%YR 
       IF(MO1.EQ.0)THEN
        MO1=12
        YR1=LT%YR-1
       ENDIF
       MO2=LT%MO
       YR2=LT%YR
      ELSE
       MO1=LT%MO
       YR1=LT%YR
       MO2=LT%MO+1
       YR2=LT%YR
       IF(MO2.EQ.13)THEN
        MO2=1
        YR2=LT%YR+1
       ENDIF
      ENDIF
      CALL DATE2TIME(TIME1,DOY1,GMT1,YR1,MO1, &
       NUMI,ZEROI,ZEROI,ZEROI)
      CALL DATE2TIME(TIME2,DOY2,GMT2,YR2,MO2, &
       NUMI,ZEROI,ZEROI,ZEROI)
      WT1= (TIME2-LT%TIME)/(TIME2-TIME1)
      WT2= (LT%TIME-TIME1)/(TIME2-TIME1)
      IF (LP%LAI .GE. 2) CALL MOSLAIREAD(YR1,MO1, &
                                YR2,MO2,TIME1,TIME2,WT1,WT2)
      
      OPEN(UNIT=16,FILE=MOSDRV%MOS_MVFILE,STATUS='OLD')
      DO J=1,MOSDRV%MOS_NMVEGP
       DO K=1,12
        READ(16,*)(VALUEMON(I,J,K),I=1,LP%NT)
       ENDDO !K
      ENDDO !J
      CLOSE(16)

      DO I=1,LD%NCH
         DO J=1,MOSDRV%MOS_NMVEGP
            DO K=1,12
               VEGMP(I,J,K)=VALUEMON(TILE(I)%VEGT,J,K)
            ENDDO !K
         ENDDO !J
      ENDDO !I

!!$      open(unit=41,file='lai1old.bin',form="unformatted")
!!$      write(41) mos%lai1
!!$      write(41) mos%lai2
!!$      write(41) mos%sai1
!!$      write(41) mos%sai2
!!$      close(41)
!!$
!!$      print*,"WTF: LAI ",wt1,wt2
!!$      print*,mos(200000)%lai1
!!$      print*,mos(200000)%lai2
!!$      print*,mos(200000)%sai1
!!$      print*,mos(200000)%sai2
!!$      print*,mos(200000)%lai
!!$      print*,mos(200000)%dsai
!!$      print*,mos(200000)%green

      DO T=1,LD%NCH

!      If tile is in the southern hemisphere, use veg parameters from the 
!	opposite (6 months away) time of year.
!         index = gindex(tile(t)%col, tile(t)%row)
         index = tile(t)%index
       IF (grid(index)%LAT .lt. 0.0) then
	NHMO1 = MOD((MO1+6),12)
	IF (NHMO1 .EQ. 0) NHMO1=12
	NHMO2 = MOD((MO2+6),12)
	IF (NHMO2 .EQ. 0) NHMO2=12
       ELSE
	NHMO1 = MO1
	NHMO2 = MO2
       END IF

       IF (LP%LAI .GE. 2) THEN

        MOS(T)%VEGIP(1) = MOS(T)%GREEN
        MOS(T)%VEGIP(2) = MOS(T)%LAI
        IF (TILE(T)%VEGT .EQ. 12) THEN
         MOS(T)%VEGIP(1) = 0.001
         MOS(T)%VEGIP(2) = 0.001
        ENDIF

        DO P=3,MOSDRV%MOS_NMVEGP
          MOS(T)%VEGIP(P) = WT1*VEGMP(T,P,NHMO1) + &
                           WT2*VEGMP(T,P,NHMO2)
        ENDDO

      ELSE

       DO P=1,MOSDRV%MOS_NMVEGP
        MOS(T)%VEGIP(P) = WT1*VEGMP(T,P,NHMO1) +  &
     	 WT2*VEGMP(T,P,NHMO2)
        ENDDO
	 
      ENDIF

      ENDDO
    
      RETURN
    END subroutine mosdynp
 
