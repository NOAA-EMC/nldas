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
! clm1lairead.f90:
!
! DESCRIPTION:
!  This program reads in AVHRR LAI data for clm1
!
! REVISION HISTORY:
!  27 Nov 2001: Jon Gottschalck; Initial code
!  20 Feb 2002: Jon Gottschalck; Modified to use for 1/4 and 2x2.5 using 1/8 degree monthly data
!  01 Oct 2002: Jon Gottschalck; Modified to add MODIS LAI data
!=========================================================================
 
  SUBROUTINE CLM1LAIREAD (LDAS,DRV,TILE,CLM1)

  use ldas_module         ! LDAS parameters
  use drv_module          ! 1-D Land Model Driver variables
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! 1-D CLM variables
  
  IMPLICIT NONE

!=== Arguments ===========================================================

  TYPE (LDASDEC)     :: LDAS
  TYPE (DRVDEC)      :: DRV
  TYPE (CLM_TILEDEC) :: TILE(DRV%NCH)
  TYPE (clm11D)       :: clm1 (DRV%NCH)

!=== Local variables

  INTEGER            :: LINE,T,V,MLAT,MLON,L
  CHARACTER*1        :: LAI_1(14),LAI_2(14),SAI_1(14),SAI_2(14)  ! LAI/SAI for all vegetation types
  REAL               :: LAT1,LON1,LAT2,LON2                      ! Lat/Lon to determine specific data for each tile in direct access files
  REAL*8             :: TIME1,TIME2                              ! Temporary Time variables
  INTEGER            :: YR1,MO1,YR2,MO2                          ! Temporary Time variables
  INTEGER            :: DOY1,DOY2                                ! Temporary Time variables
  REAL               :: WT1,WT2,GMT1,GMT2                        ! Interpolation weights
  INTEGER            :: ZEROI,NUMI                               ! Integer Number Holders
  INTEGER            :: LAIFLAG,IOS1,IOS2                        ! Flag to read in new LAI data, file error variables 
  CHARACTER (LEN=4)  :: CYR1,CYR2                                ! Filename variables
  CHARACTER (LEN=2)  :: CMO1,CMO2                                ! Filename variables
  CHARACTER (LEN=80) :: NAME
  INTEGER            :: FLAG1, FLAG2
  INTEGER            :: CNT1,CNT2,CNT3,CNT4,II8,J8
  REAL               :: SUM1,SUM2,SUM3,SUM4
  INTEGER            :: D_START_NR,D_START_NC,START_8TH_NR,START_8TH_NC
  INTEGER            :: END_8TH_NR,END_8TH_NC,K,MM
  character(len=80) :: name9,  name10, name11, name12
  character(len=80) :: name13, name14, name15, name16

!=== End Local variable list

!=== Determine current time to find correct LAI files
  IF (LDAS%TSCOUNT .EQ. 0) THEN
   DRV%YR = DRV%SYR
   DRV%MO = DRV%SMO
   DRV%DA = DRV%SDA
   DRV%MN = DRV%SMN
   DRV%SS = DRV%SSS
  ELSE
   DRV%YR = DRV%YR
   DRV%MO = DRV%MO
   DRV%DA = DRV%DA
   DRV%MN = DRV%MN
   DRV%SS = DRV%SS
  ENDIF
  
  CALL DRV_DATE2TIME(DRV%TIME,DRV%DOY,DRV%DAY,DRV%GMT,DRV%YR, &
&                DRV%MO,DRV%DA,DRV%HR,DRV%MN,DRV%SS)

!=== Initialize LAI flag varaible
   LAIFLAG = 0
 
!=== Initialize Numbers
     ZEROI=0
     NUMI=16

!=== Determine Monthly data Times (Assume Monthly value valid at DA=16 HR=00Z)
      IF (DRV%DA .LT. 16) THEN
       MO1 = DRV%MO-1
       YR1 = DRV%YR
       IF (MO1 .EQ. 0) THEN
        MO1 = 12
        YR1 = DRV%YR - 1
       ENDIF
       MO2 = DRV%MO
       YR2 = DRV%YR
      ELSE
       MO1 = DRV%MO
       YR1 = DRV%YR
       MO2 = DRV%MO+1
       YR2 = DRV%YR
       IF (MO2 .EQ. 13) THEN
        MO2 = 1
        YR2 = DRV%YR + 1
       ENDIF
      ENDIF

      CALL DATE2TIME(TIME1,DOY1,GMT1,YR1,MO1,NUMI,ZEROI,ZEROI,ZEROI)
      CALL DATE2TIME(TIME2,DOY2,GMT2,YR2,MO2,NUMI,ZEROI,ZEROI,ZEROI) 

!=== Check to see if need new LAI data
  IF (TIME2 .GT. LDAS%LAITIME) THEN 
    LAIFLAG = 1
  ELSE
    LAIFLAG = 0
  ENDIF

!=== Get new LAI data if required
  IF (LAIFLAG .EQ. 1) THEN
  
     open(unit=90, file='temp', form='formatted', access='direct', &
     recl=80)
     write(90, 96, rec=1) yr1, mo1
     read (90, 92, rec=1) cyr1, cmo1
     close(90)
     open(unit=90, file='temp', form='formatted', access='direct', recl=80)
     write(90, 96, rec=1) yr2,  mo2  
     read (90, 92, rec=1) cyr2, cmo2
     close(90)
 96  format(i4,i2.2)
 92  format(a4,a2)
 
 LDAS%LAITIME = TIME2
 
  SELECT CASE (LDAS%LAI)
   CASE(2)
    CALL AVHRR_FILE(NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
                   LDAS%AVHRRDIR,CYR1,CYR2,CMO1,CMO2)		   
   CASE(3)
    CALL MODIS_FILE(NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
                  LDAS%MODISDIR,CYR1,CYR2,CMO1,CMO2)
   CASE DEFAULT
    PRINT*, "NOT A VALID LAI OPTION"
    STOP
  END SELECT
		   
!=== Open Satellite LAI files (assumes realtime monthly files are present first
!=== then uses climatology files)

!=== Assume realtime monthly files are present as default
  
  FLAG1 = 0
  FLAG2 = 0
   OPEN(10,FILE=NAME9,STATUS='OLD',FORM='UNFORMATTED',&
&  ACCESS='DIRECT',RECL=24,IOSTAT=IOS1)
   OPEN(11,FILE=NAME10,STATUS='OLD',FORM='UNFORMATTED',&
&  ACCESS='DIRECT',RECL=24,IOSTAT=IOS2)
   OPEN(12,FILE=NAME13,STATUS='OLD',FORM='UNFORMATTED',&
&  ACCESS='DIRECT',RECL=24,IOSTAT=IOS1)
   OPEN(13,FILE=NAME14,STATUS='OLD',FORM='UNFORMATTED',&
&  ACCESS='DIRECT',RECL=24,IOSTAT=IOS2)

  print*, "Using 1/8 Satellite LAI/DSAI data for month 1 ", &
  NAME9
  print*, "Using 1/8 Satellite LAI/DSAI data for month 2 ", &
  NAME10

  IF (IOS1 .NE. 0) THEN
   CLOSE(10)
   OPEN(10,FILE=NAME11,STATUS='OLD',FORM='UNFORMATTED',&
&   ACCESS='DIRECT',RECL=24)
   CLOSE(12)
   OPEN(12,FILE=NAME15,STATUS='OLD',FORM='UNFORMATTED',&
&   ACCESS='DIRECT',RECL=24)

   print*, "No realtime monthly data for month 1"
   print*, "Using 1/8 Satellite LAI/DSAI data for month 1 ", &
   NAME11

  FLAG1 = 1

  ENDIF
  IF (IOS2 .NE. 0) THEN
   CLOSE(11)
   OPEN(11,FILE=NAME12,STATUS='OLD',FORM='UNFORMATTED',&
&   ACCESS='DIRECT',RECL=24)
    CLOSE(13)
   OPEN(13,FILE=NAME16,STATUS='OLD',FORM='UNFORMATTED',&
&   ACCESS='DIRECT',RECL=24)
   FLAG2 = 1
   print*, "No realtime monthly data for month 2"
   print*, "Using 1/8 Satellite LAI/DSAI data for month 2 ", &
   NAME12
  ENDIF
  
!=== Loop through tiles to assign Satellite LAI values for each tile

  DO T=1,DRV%NCH

  IF (LDAS%DOMAIN .NE. 1) THEN
  
!=== Select either 1/4, 1/2, 1, or 2x2.5 aggregation
!=== Locates latitude and longitude of tile and determines the rows and columns in
!=== 1/8 grid space in which to read records and compute sums

   SELECT CASE (LDAS%DOMAIN)
   CASE (2)
    D_START_NR  = ((clm1(T)%LATDEG - (-59.875)) / 0.25) + 1
    START_8TH_NR = ((D_START_NR - 1) * 2) + 1
    END_8TH_NR   =   START_8TH_NR + 1
    D_START_NC  = ((clm1(T)%LONDEG - (-179.875)) / 0.25) + 1
    START_8TH_NC = ((D_START_NC - 1) * 2) + 1
    END_8TH_NC   =   START_8TH_NC + 1
   CASE (3)
    D_START_NR  = ((clm1(T)%LATDEG - (-60)) / 2.0) + 1
    START_8TH_NR = ((D_START_NR - 1) * 16) + 1
    END_8TH_NR   =   START_8TH_NR + 15
    D_START_NC  = ((clm1(T)%LONDEG - (-180)) / 2.5) + 1
    START_8TH_NC = ((D_START_NC - 1) * 20) + 1
    END_8TH_NC   =   START_8TH_NC + 19
   CASE (4)
    D_START_NR  = ((clm1(T)%LATDEG - (-59.500)) / 1.00) + 1
    START_8TH_NR = ((D_START_NR - 1) * 8) + 1
    END_8TH_NR   =   START_8TH_NR + 7
    D_START_NC  = ((clm1(T)%LONDEG - (-179.500)) / 1.00) + 1
    START_8TH_NC = ((D_START_NC - 1) * 8) + 1
    END_8TH_NC   =   START_8TH_NC + 7
   CASE (5)
    D_START_NR  = ((clm1(T)%LATDEG - (-59.750)) / 0.50) + 1
    START_8TH_NR = ((D_START_NR - 1) * 4) + 1
    END_8TH_NR   =   START_8TH_NR + 3
    D_START_NC  = ((clm1(T)%LONDEG - (-179.750)) / 0.50) + 1
    START_8TH_NC = ((D_START_NC - 1) * 4) + 1
    END_8TH_NC   =   START_8TH_NC + 3
   CASE DEFAULT
    PRINT*, "IMPROPER DOMAIN SELECTION"
    STOP
  END SELECT
  
!=== Initilaize sums for LAI month 1, LAI month 2, DSAI month 1, DSAI month 2

      SUM1 = 0.0
      SUM2 = 0.0
      SUM3 = 0.0
      SUM4 = 0.0
      CNT1 = 0
      CNT2 = 0
      CNT3 = 0
      CNT4 = 0
					
!=== Looping over 1/8 grid space that relates to 1/4 or 2x2.5 domains

     DO II8 = START_8TH_NR,END_8TH_NR
      DO J8 = START_8TH_NC,END_8TH_NC
        LINE = (II8 - 1)*2880 + J8

!=== Reading in record depending on type of data month 1 

	 READ(10,REC=LINE) LAT1, LON1, LAI_1
	 READ(12,REC=LINE) LAT1, LON1, SAI_1
	
!=== Reading in record depending on type of data month2 

	 READ(11,REC=LINE) LAT2, LON2, LAI_2
	 READ(13,REC=LINE) LAT2, LON2, SAI_2
	
!=== Convert 4 byte integers or 1 byte characters to real values for use in LDAS
!=== Summing over the 1/8 domain points, month 1 and month 2

     SELECT CASE (LDAS%LAI)
     
     CASE(2)     ! AVHRR LAI
      IF (ICHAR(LAI_1(TILE(T)%VEGT+1)) .NE. 251 .AND. ICHAR(LAI_1(TILE(T)%VEGT+1)) .NE. 0) THEN
       SUM1 = SUM1 + (ICHAR(LAI_1(TILE(T)%VEGT+1))) * 0.04
       CNT1 = CNT1 + 1
      ENDIF
      IF (ICHAR(SAI_1(TILE(T)%VEGT+1)) .NE. 251 .AND. ICHAR(SAI_1(TILE(T)%VEGT+1)) .NE. 0) THEN
       SUM3 = SUM3 + (ICHAR(SAI_1(TILE(T)%VEGT+1))) * 0.04
       CNT3 = CNT3 + 1
      ENDIF   
      IF (ICHAR(LAI_2(TILE(T)%VEGT+1)) .NE. 251 .AND. ICHAR(LAI_2(TILE(T)%VEGT+1)) .NE. 0) THEN
       SUM2 = SUM2 + (ICHAR(LAI_2(TILE(T)%VEGT+1))) * 0.04
       CNT2 = CNT2 + 1
      ENDIF
      IF (ICHAR(SAI_2(TILE(T)%VEGT+1)) .NE. 251 .AND. ICHAR(SAI_2(TILE(T)%VEGT+1)) .NE. 0) THEN
       SUM4 = SUM4 + (ICHAR(SAI_2(TILE(T)%VEGT+1))) * 0.04
       CNT4 = CNT4 + 1
      ENDIF
      
      CASE(3)     ! MODIS LAI
      IF (ICHAR(LAI_1(TILE(T)%VEGT+1)) .LT. 200) THEN
       SUM1 = SUM1 + (ICHAR(LAI_1(TILE(T)%VEGT+1))) * 0.10
       CNT1 = CNT1 + 1
      ENDIF
      IF (ICHAR(SAI_1(TILE(T)%VEGT+1)) .LT. 200) THEN
       SUM3 = SUM3 + (ICHAR(SAI_1(TILE(T)%VEGT+1))) * 0.10
       CNT3 = CNT3 + 1
      ENDIF
      IF (ICHAR(LAI_2(TILE(T)%VEGT+1)) .LT. 200) THEN
       SUM2 = SUM2 + (ICHAR(LAI_2(TILE(T)%VEGT+1))) * 0.10
       CNT2 = CNT2 + 1
      ENDIF
      IF (ICHAR(SAI_2(TILE(T)%VEGT+1)) .LT. 200) THEN
       SUM4 = SUM4 + (ICHAR(SAI_2(TILE(T)%VEGT+1))) * 0.10
       CNT4 = CNT4 + 1
      ENDIF				   
      CASE DEFAULT
       PRINT*, "Not a valid LAI Domain"
       STOP
     
     END SELECT
				     
      ENDDO
     ENDDO
     
!=== Compute averages for the vegetation type represented by tile

     IF (CNT1 .NE. 0) THEN
       clm1(T)%LAI_T1_F = SUM1 / CNT1
     ELSE
       clm1(T)%LAI_T1_F = 0
     ENDIF
     IF (CNT2 .NE. 0) THEN
       clm1(T)%LAI_T2_F = SUM2 / CNT2
     ELSE
       clm1(T)%LAI_T2_F = 0       
     ENDIF
     IF (CNT3 .NE. 0) THEN
       clm1(T)%SAI_T1_F = SUM3 / CNT3
     ELSE
       clm1(T)%SAI_T1_F = 0
     ENDIF
     IF (CNT4 .NE. 0) THEN
       clm1(T)%SAI_T2_F = SUM4 / CNT4
     ELSE
       clm1(T)%SAI_T2_F = 0
     ENDIF
							  
!=== Do NLDAS section (no aggregation)

   ELSE

     MLAT = (clm1(T)%LATDEG - (-59.9375)) / 0.125 + 1
     MLON = (clm1(T)%LONDEG - (-179.9375)) / 0.125 + 1
     LINE = (MLAT - 1)*2880 + MLON
     
!=== Read in data for both LAI and DSAI
!=== IF/ENDIF due to the fact climatology has different format then the actual monthly files
      READ (10,REC=LINE) LAT1, LON1, LAI_1
      READ (12,REC=LINE) LAT1, LON1, SAI_1
      READ (11,REC=LINE) LAT2, LON2, LAI_2
      READ (13,REC=LINE) LAT2, LON2, SAI_2

!=== Scale to real physical values
      SELECT CASE(LDAS%LAI)
       CASE(2)   ! AVHRR LAI
        clm1(T)%LAI_T1_F = ICHAR(LAI_1(TILE(T)%VEGT+1)) * 0.04
        clm1(T)%SAI_T1_F = ICHAR(SAI_1(TILE(T)%VEGT+1)) * 0.04
        clm1(T)%LAI_T2_F = ICHAR(LAI_2(TILE(T)%VEGT+1)) * 0.04
        clm1(T)%SAI_T2_F = ICHAR(SAI_2(TILE(T)%VEGT+1)) * 0.04
       CASE(3)
	clm1(T)%LAI_T1_F = ICHAR(LAI_1(TILE(T)%VEGT+1)) * 0.10
	clm1(T)%SAI_T1_F = ICHAR(SAI_1(TILE(T)%VEGT+1)) * 0.10
	clm1(T)%LAI_T2_F = ICHAR(LAI_2(TILE(T)%VEGT+1)) * 0.10
	clm1(T)%SAI_T2_F = ICHAR(SAI_2(TILE(T)%VEGT+1)) * 0.10
       CASE DEFAULT
	PRINT*, "Invalid domain for LAI data"
	STOP
      END SELECT
     
   ENDIF

  ENDDO
  
  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
  
  ENDIF
  
!=== Determine weights between months
      WT1= (TIME2-DRV%TIME)/(TIME2-TIME1)
      WT2= (DRV%TIME-TIME1)/(TIME2-TIME1)

!=== Assign interpolated LAI and DSAI values to the CLM variable names used in CLM main
     DO T=1,DRV%NCH

       clm1(T)%TLAI = WT1 * clm1(T)%LAI_T1_F + WT2 * clm1(T)%LAI_T2_F
       clm1(T)%TSAI = WT1 * clm1(T)%SAI_T1_F + WT2 * clm1(T)%SAI_T2_F
       IF (TILE(T)%VEGT .EQ. 12) THEN
         clm1(T)%TLAI=0.0
	 clm1(T)%TSAI=0.0
       ENDIF

     ENDDO
     
END SUBROUTINE CLM1LAIREAD

!==============================================================
!
!  DESCRIPTION: This subroutine puts together AVHRR file name
!==============================================================

subroutine avhrr_file (NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
                       avhrrdir, cyr1, cyr2, cmo1, cmo2 )

  implicit none

!==== Local Variables=======================

  character(len=80) :: name9,  name10, name11, name12, avhrrdir
  character(len=80) :: name13, name14, name15, name16
  character(len=4)  :: cyr1, cyr2
  character(len=2)  :: cmo1, cmo2
  integer :: i, c, flag
  character*1 :: fbase(80), fbase_2(80), fbase_3(80), fbase_4(80)
  character*1 :: fdir(7), fdir_2(7), fdir_3(7), fdir_4(7)
  character*1 :: fsubsn(15), fsubsn_2(15), fsubsn_3(15), fsubsn_4(15)
  character*1 :: fsubsn_5(15), fsubsn_6(15), fsubsn_7(15), fsubsn_8(15)

!=== End Variable Definition ===============
!=== formats for filename segments
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, i2)
95 format (10a1)
96 format (a40)
97 format (a9)
67 format (a15)
98 format (a1, a4, a2)
66 format (a5,a2)
99 format (7a1)		
			      
  open(unit=90, file='temp', form='formatted', access='direct', recl=80)
  open(unit=91, file='temp_2', form='formatted', access='direct', recl=80)
  open(unit=92, file='temp_3', form='formatted', access='direct', recl=80)
  open(unit=93, file='temp_4', form='formatted', access='direct', recl=80)
  write(90, 96, rec=1) avhrrdir
  write(91, 96, rec=1) avhrrdir
  write(92, 96, rec=1) avhrrdir
  write(93, 96, rec=1) avhrrdir
  read(90, 92, rec=1) (fbase(i), i=1,80)
  read(91, 92, rec=1) (fbase_2(i), i=1,80)
  read(92, 92, rec=1) (fbase_3(i), i=1,80)
  read(93, 92, rec=1) (fbase_4(i), i=1,80)

  write(90, 98, rec=1) '/', cyr1, cmo1
  read(90, 99, rec=1) fdir
  write(91, 98, rec=1) '/', cyr2, cmo2
  read(91, 99, rec=1) fdir_2
  write(92, 66, rec=1) '/CLIM', cmo1
  read(92, 99, rec=1) fdir_3
  write(93, 66, rec=1) '/CLIM', cmo2
  read(93, 99, rec=1) fdir_4
 
 do i = 1, 7
   if ( fdir(i) == ' ' ) fdir(i) = '0'
   if ( fdir_2(i) == ' ' ) fdir_2(i) = '0'
   if ( fdir_3(i) == ' ' ) fdir_3(i) = '0'
   if ( fdir_4(i) == ' ' ) fdir_4(i) = '0'
  enddo	

  write(90, 67, rec=1) '_AVHRRLAI_0.125'
  write(91, 67, rec=1) '_AVHRRLAI_0.125'
  write(92, 67, rec=1) '_AVHRRLAI_0.125'
  write(93, 67, rec=1) '_AVHRRLAI_0.125'
  read (90, 92, rec=1) (fsubsn(i), i=1,15)
  read (91, 92, rec=1) (fsubsn_2(i), i=1,15)
  read (92, 92, rec=1) (fsubsn_3(i), i=1,15)
  read (93, 92, rec=1) (fsubsn_4(i), i=1,15)
  write(90, 67, rec=1) '_AVHRRSAI_0.125'
  write(91, 67, rec=1) '_AVHRRSAI_0.125'
  write(92, 67, rec=1) '_AVHRRSAI_0.125'
  write(93, 67, rec=1) '_AVHRRSAI_0.125'
  read (90, 92, rec=1) (fsubsn_5(i), i=1,15)
  read (91, 92, rec=1) (fsubsn_6(i), i=1,15)
  read (92, 92, rec=1) (fsubsn_7(i), i=1,15)
  read (93, 92, rec=1) (fsubsn_8(i), i=1,15)

!sets c as the last character position of fbase
  c = 0
  do i = 1, 80
    if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
    if ( (fbase_2(i) == ' ') .and. (c == 0) ) c = i-1
    if ( (fbase_3(i) == ' ') .and. (c == 0) ) c = i-1
    if ( (fbase_4(i) == ' ') .and. (c == 0) ) c = i-1
  end do
  
  write(90, 92, rec=1) (fbase(i), i=1,c), (fdir(i), i=1,7),  &
                       (fsubsn(i), i=1,15)
  write(91, 92, rec=1) (fbase_2(i), i=1,c), (fdir_2(i), i=1,7),  &
                       (fsubsn_2(i), i=1,15)
  write(92, 92, rec=1) (fbase_3(i), i=1,c), (fdir_3(i), i=1,7),  &
                       (fsubsn_3(i), i=1,15)
  write(93, 92, rec=1) (fbase_4(i), i=1,c), (fdir_4(i), i=1,7),  &
                       (fsubsn_4(i), i=1,15)
  read(90, 93, rec=1) name9
  read(91, 93, rec=1) name10
  read(92, 93, rec=1) name11
  read(93, 93, rec=1) name12
  write(90, 92, rec=1) (fbase(i), i=1,c), (fdir(i), i=1,7),  &
                       (fsubsn_5(i), i=1,15)
  write(91, 92, rec=1) (fbase_2(i), i=1,c), (fdir_2(i), i=1,7),  &
                       (fsubsn_6(i), i=1,15)
  write(92, 92, rec=1) (fbase_3(i), i=1,c), (fdir_3(i), i=1,7),  &
                       (fsubsn_7(i), i=1,15)
  write(93, 92, rec=1) (fbase_4(i), i=1,c), (fdir_4(i), i=1,7),  &
                       (fsubsn_8(i), i=1,15)
  read(90, 93, rec=1) name13
  read(91, 93, rec=1) name14
  read(92, 93, rec=1) name15
  read(93, 93, rec=1) name16

  close(90)
  close(91)
  close(92)
  close(93)

  return
			       
 end subroutine avhrr_file

!==============================================================
!
!  DESCRIPTION: This subroutine puts together MODIS file name
!==============================================================

subroutine modis_file (NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
                       modisdir, cyr1, cyr2, cmo1, cmo2 )

  implicit none

!==== Local Variables=======================

  character(len=80) :: name9,  name10, name11, name12, modisdir
  character(len=80) :: name13, name14, name15, name16
  character(len=4)  :: cyr1, cyr2
  character(len=2)  :: cmo1, cmo2
  integer :: i, c, flag
  character*1 :: fbase(80), fbase_2(80), fbase_3(80), fbase_4(80)
  character*1 :: fdir(7), fdir_2(7), fdir_3(7), fdir_4(7)
  character*1 :: fsubsn(15), fsubsn_2(15), fsubsn_3(15), fsubsn_4(15)
  character*1 :: fsubsn_5(15), fsubsn_6(15), fsubsn_7(15), fsubsn_8(15)

!=== End Variable Definition ===============
!=== formats for filename segments
92 format (80a1)
93 format (a80)
94 format (i4, i2, i2, i2)
95 format (10a1)
96 format (a40)
97 format (a9)
67 format (a15)
98 format (a1, a4, a2)
66 format (a5,a2)
99 format (7a1)

  open(unit=90, file='temp', form='formatted', access='direct', recl=80)
  open(unit=91, file='temp_2', form='formatted', access='direct', recl=80)
  open(unit=92, file='temp_3', form='formatted', access='direct', recl=80)
  open(unit=93, file='temp_4', form='formatted', access='direct', recl=80)
  write(90, 96, rec=1) modisdir
  write(91, 96, rec=1) modisdir
  write(92, 96, rec=1) modisdir
  write(93, 96, rec=1) modisdir
  read(90, 92, rec=1) (fbase(i), i=1,80)
  read(91, 92, rec=1) (fbase_2(i), i=1,80)
  read(92, 92, rec=1) (fbase_3(i), i=1,80)
  read(93, 92, rec=1) (fbase_4(i), i=1,80)

  write(90, 98, rec=1) '/', cyr1, cmo1
  read(90, 99, rec=1) fdir
  write(91, 98, rec=1) '/', cyr2, cmo2
  read(91, 99, rec=1) fdir_2
  write(92, 66, rec=1) '/CLIM', cmo1
  read(92, 99, rec=1) fdir_3
  write(93, 66, rec=1) '/CLIM', cmo2
  read(93, 99, rec=1) fdir_4

  do i = 1, 7
   if ( fdir(i) == ' ' ) fdir(i) = '0'
   if ( fdir_2(i) == ' ' ) fdir_2(i) = '0'
   if ( fdir_3(i) == ' ' ) fdir_3(i) = '0'
   if ( fdir_4(i) == ' ' ) fdir_4(i) = '0'
  enddo

  write(90, 67, rec=1) '_MODISLAI_0.125'
  write(91, 67, rec=1) '_MODISLAI_0.125'
  write(92, 67, rec=1) '_MODISLAI_0.125'
  write(93, 67, rec=1) '_MODISLAI_0.125'
  read (90, 92, rec=1) (fsubsn(i), i=1,15)
  read (91, 92, rec=1) (fsubsn_2(i), i=1,15)
  read (92, 92, rec=1) (fsubsn_3(i), i=1,15)
  read (93, 92, rec=1) (fsubsn_4(i), i=1,15)
  write(90, 67, rec=1) '_MODISSAI_0.125'
  write(91, 67, rec=1) '_MODISSAI_0.125'
  write(92, 67, rec=1) '_MODISSAI_0.125'
  write(93, 67, rec=1) '_MODISSAI_0.125'
  read (90, 92, rec=1) (fsubsn_5(i), i=1,15)
  read (91, 92, rec=1) (fsubsn_6(i), i=1,15)
  read (92, 92, rec=1) (fsubsn_7(i), i=1,15)
  read (93, 92, rec=1) (fsubsn_8(i), i=1,15)

!sets c as the last character position of fbase
  c = 0
  do i = 1, 80
    if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
    if ( (fbase_2(i) == ' ') .and. (c == 0) ) c = i-1
    if ( (fbase_3(i) == ' ') .and. (c == 0) ) c = i-1
    if ( (fbase_4(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(90, 92, rec=1) (fbase(i), i=1,c), (fdir(i), i=1,7),  &
                       (fsubsn(i), i=1,15)
  write(91, 92, rec=1) (fbase_2(i), i=1,c), (fdir_2(i), i=1,7),  &
                       (fsubsn_2(i), i=1,15)
  write(92, 92, rec=1) (fbase_3(i), i=1,c), (fdir_3(i), i=1,7),  &
                       (fsubsn_3(i), i=1,15)
  write(93, 92, rec=1) (fbase_4(i), i=1,c), (fdir_4(i), i=1,7),  &
                       (fsubsn_4(i), i=1,15)
  read(90, 93, rec=1) name9
  read(91, 93, rec=1) name10
  read(92, 93, rec=1) name11
  read(93, 93, rec=1) name12
  write(90, 92, rec=1) (fbase(i), i=1,c), (fdir(i), i=1,7),  &
                       (fsubsn_5(i), i=1,15)
  write(91, 92, rec=1) (fbase_2(i), i=1,c), (fdir_2(i), i=1,7),  &
                       (fsubsn_6(i), i=1,15)
  write(92, 92, rec=1) (fbase_3(i), i=1,c), (fdir_3(i), i=1,7),  &
                       (fsubsn_7(i), i=1,15)
  write(93, 92, rec=1) (fbase_4(i), i=1,c), (fdir_4(i), i=1,7),  &
                       (fsubsn_8(i), i=1,15)
  read(90, 93, rec=1) name13
  read(91, 93, rec=1) name14
  read(92, 93, rec=1) name15
  read(93, 93, rec=1) name16

  close(90)
  close(91)
  close(92)
  close(93)

  return

 end subroutine modis_file
