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
! clmlairead.f90:
!
! DESCRIPTION:
!  This program reads in AVHRR LAI data for CLM
!
! REVISION HISTORY:
!  27 Nov 2001: Jon Gottschalck; Initial code
!  20 Feb 2002: Jon Gottschalck; Modified to use for 1/4 and 2x2.5 using 1/8 degree monthly data
!  01 Oct 2002: Jon Gottschalck; Modified to add MODIS LAI data
!=========================================================================
 
  SUBROUTINE CLM2LAIREAD (LDAS,TILE)

  use clm_varder
  use ldas_module         ! LDAS parameters
  use tile_module         ! LDAS Tile variables
  use clmtype             ! 1-D CLM variables
  
  IMPLICIT NONE

!=== Arguments ===========================================================

  TYPE (LDASDEC) :: LDAS
  TYPE (TILEDEC) :: TILE(LDAS%NCH)

!=== Local variables

  INTEGER            :: LINE,T,V,MLAT,MLON,L
  CHARACTER*1        :: LAI_1(14),LAI_2(14),SAI_1(14),SAI_2(14)  ! LAI/SAI for all vegetation types
  REAL               :: TOP_1(14),TOP_2(14),BOT_1(14),BOT_2(14)  ! TOP/BOT for all vegetation types
  REAL               :: LAT1,LON1,LAT2,LON2                      ! Lat/Lon to determine specific data for each tile in direct access files
  REAL*8             :: TIME1,TIME2                              ! Temporary Time variables
  INTEGER            :: YR1,MO1,YR2,MO2                          ! Temporary Time variables
  INTEGER            :: DOY1,DOY2                                ! Temporary Time variables
  REAL               :: WT1,WT2,GMT1,GMT2                        ! Interpolation weights
  INTEGER            :: ZEROI,NUMI                               ! Integer Number Holders
  INTEGER            :: LAIFLAG,IOS1,IOS2                        ! Flag to read in new LAI data, file error variables 
  CHARACTER (LEN=4)  :: CYR1,CYR2                                ! Filename variables
  CHARACTER (LEN=2)  :: CMO1,CMO2                                ! Filename variables
  CHARACTER (LEN=80) :: NAME, AVHRRDIR
  INTEGER            :: FLAG1, FLAG2
  INTEGER            :: CNT1,CNT2,CNT3,CNT4,CNT5,CNT6,CNT7,CNT8,II8,J8
  REAL               :: SUM1,SUM2,SUM3,SUM4,SUM5,SUM6,SUM7,SUM8
  INTEGER            :: D_START_NR,D_START_NC,START_8TH_NR,START_8TH_NC
  INTEGER            :: END_8TH_NR,END_8TH_NC,K,MM
  character(len=80) :: name2, name3, name4
  character(len=80) :: name5, name6, name7, name8
  character(len=80) :: name9,  name10, name11, name12
  character(len=80) :: name13, name14, name15, name16
  character(len=80) :: ntop1, ntop2, nbot1, nbot2
  real    :: pthtop(17),pthbot(17)  ! Canopy top and bottom height for the 13 UMD vegetation types
  integer :: p
 
!=== End Local variable list

!=== Determine current time to find correct LAI files
  IF (LDAS%TSCOUNT .EQ. 0) THEN
   LDAS%YR = LDAS%SYR
   LDAS%MO = LDAS%SMO
   LDAS%DA = LDAS%SDA
   LDAS%MN = LDAS%SMN
   LDAS%SS = LDAS%SSS
  ELSE
   LDAS%YR = LDAS%YR
   LDAS%MO = LDAS%MO
   LDAS%DA = LDAS%DA
   LDAS%MN = LDAS%MN
   LDAS%SS = LDAS%SS
  ENDIF
  
  CALL DATE2TIME(LDAS%TIME,LDAS%DOY,LDAS%GMT,LDAS%YR, &
&                LDAS%MO,LDAS%DA,LDAS%HR,LDAS%MN,LDAS%SS)

!=== Initialize LAI flag varaible
   LAIFLAG = 0
 
!=== Initialize Numbers
     ZEROI=0
     NUMI=16

!=== Determine Monthly data Times (Assume Monthly value valid at DA=16 HR=00Z)
      IF (LDAS%DA .LT. 16) THEN
       MO1 = LDAS%MO-1
       YR1 = LDAS%YR
       IF (MO1 .EQ. 0) THEN
        MO1 = 12
        YR1 = LDAS%YR - 1
       ENDIF
       MO2 = LDAS%MO
       YR2 = LDAS%YR
      ELSE
       MO1 = LDAS%MO
       YR1 = LDAS%YR
       MO2 = LDAS%MO+1
       YR2 = LDAS%YR
       IF (MO2 .EQ. 13) THEN
        MO2 = 1
        YR2 = LDAS%YR + 1
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

  AVHRRDIR = LDAS%AVHRRDIR

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
   CASE(1)
    CALL NCARLAI_CLM2(NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
              LDAS%AVHRRDIR,CYR1,CYR2,CMO1,CMO2)
   CASE(2)      
    CALL AVHRR_FILE_2(NAME,NAME2,NAME3,NAME4,NAME5,NAME6,NAME7,NAME8, &
                   NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
                   LDAS%AVHRRDIR,CYR1,CYR2,CMO1,CMO2)
   CASE(3)
    CALL MODIS_FILE_2(NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
                      LDAS%MODISDIR,CYR1,CYR2,CMO1,CMO2)
   CASE DEFAULT
    PRINT*, "NOT A VALID LAI OPTION"
    STOP
   END SELECT
	     
   CALL NCARCANHT_CLM2(NTOP1,NTOP2,NBOT1,NBOT2, &
                   AVHRRDIR,CYR1,CYR2,CMO1,CMO2)
		   
!=== Open AVHRR LAI files (assumes realtime monthly files are present first
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

  print*, "Using 1/8 AVHRR LAI/DSAI data for month 1 ", &
  NAME9
  print*, "Using 1/8 AVHRR LAI/DSAI data for month 2 ", &
  NAME10

  IF (IOS1 .NE. 0) THEN
   CLOSE(10)
   OPEN(10,FILE=NAME11,STATUS='OLD',FORM='UNFORMATTED',&
&   ACCESS='DIRECT',RECL=24)
   CLOSE(12)
   OPEN(12,FILE=NAME15,STATUS='OLD',FORM='UNFORMATTED',&
&   ACCESS='DIRECT',RECL=24)
   print*, "No realtime monthly data for month 1"
   print*, "Using 1/8 AVHRR LAI/DSAI data for month 1 ", &
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
   print*, "Using 1/8 AVHRR LAI/DSAI data for month 2 ", &
   NAME12
  ENDIF
  
!=== Opening up canopy top and bottom height data for both months  
  OPEN(50,FILE=NTOP1,STATUS='OLD',FORM='UNFORMATTED',&
       ACCESS='DIRECT',RECL=64)
  print*, "Using 1/8 NCAR CANTOP data for month 1 ", NTOP1
  OPEN(51,FILE=NTOP2,STATUS='OLD',FORM='UNFORMATTED',&
       ACCESS='DIRECT',RECL=64)
  print*, "Using 1/8 NCAR CANTOP data for month 2 ", NTOP2
  OPEN(52,FILE=NBOT1,STATUS='OLD',FORM='UNFORMATTED',&
       ACCESS='DIRECT',RECL=64)
  print*, "Using 1/8 NCAR CANBOT data for month 1 ", NBOT1
  OPEN(53,FILE=NBOT2,STATUS='OLD',FORM='UNFORMATTED',&
       ACCESS='DIRECT',RECL=64)
  print*, "Using 1/8 NCAR CANBOT data for month 2 ", NBOT2
  
!=== Loop through tiles to assign AVHRR LAI values for each tile

  DO T=1,LDAS%NCH

  IF (LDAS%DOMAIN .NE. 1) THEN
  
!=== Select either 1/4, 1/2, 1, or 2x2.5 aggregation
!=== Locates latitude and longitude of tile and determines the rows and columns in
!=== 1/8 grid space in which to read records and compute sums

   SELECT CASE (LDAS%DOMAIN)
   CASE (2)
    D_START_NR  = ((CLM(T)%LATDEG - (-59.875)) / 0.25) + 1
    START_8TH_NR = ((D_START_NR - 1) * 2) + 1
    END_8TH_NR   =   START_8TH_NR + 1
    D_START_NC  = ((CLM(T)%LONDEG - (-179.875)) / 0.25) + 1
    START_8TH_NC = ((D_START_NC - 1) * 2) + 1
    END_8TH_NC   =   START_8TH_NC + 1
   CASE (3)
    D_START_NR  = ((CLM(T)%LATDEG - (-60)) / 2.0) + 1
    START_8TH_NR = ((D_START_NR - 1) * 16) + 1
    END_8TH_NR   =   START_8TH_NR + 15
    D_START_NC  = ((CLM(T)%LONDEG - (-180)) / 2.5) + 1
    START_8TH_NC = ((D_START_NC - 1) * 20) + 1
    END_8TH_NC   =   START_8TH_NC + 19
   CASE (4)
    D_START_NR  = ((CLM(T)%LATDEG - (-59.500)) / 1.00) + 1
    START_8TH_NR = ((D_START_NR - 1) * 8) + 1
    END_8TH_NR   =   START_8TH_NR + 7
    D_START_NC  = ((CLM(T)%LONDEG - (-179.500)) / 1.00) + 1
    START_8TH_NC = ((D_START_NC - 1) * 8) + 1
    END_8TH_NC   =   START_8TH_NC + 7
   CASE (5)
    D_START_NR  = ((CLM(T)%LATDEG - (-59.750)) / 0.50) + 1
    START_8TH_NR = ((D_START_NR - 1) * 4) + 1
    END_8TH_NR   =   START_8TH_NR + 3
    D_START_NC  = ((CLM(T)%LONDEG - (-179.750)) / 0.50) + 1
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
      SUM5 = 0.0
      SUM6 = 0.0
      SUM7 = 0.0
      SUM8 = 0.0
      CNT5 = 0
      CNT6 = 0
      CNT7 = 0
      CNT8 = 0
											
!=== Looping over 1/8 grid space that relates to 1/4 or 2x2.5 domains

     DO II8 = START_8TH_NR,END_8TH_NR
      DO J8 = START_8TH_NC,END_8TH_NC
        LINE = (II8 - 1)*2880 + J8

!=== Reading in data for both months

	  READ(10,REC=LINE) LAT1, LON1, LAI_1
	  READ(12,REC=LINE) LAT1, LON1, SAI_1
	  READ(11,REC=LINE) LAT2, LON2, LAI_2
	  READ(13,REC=LINE) LAT2, LON2, SAI_2
	  READ(50,REC=LINE) LAT1, LON1, TOP_1
          READ(51,REC=LINE) LAT1, LON1, TOP_2
	  READ(52,REC=LINE) LAT1, LON1, BOT_1
	  READ(53,REC=LINE) LAT1, LON1, BOT_2
	 	 	 
!=== Convert 4 byte integers or 1 byte characters to real values for use in LDAS
!=== Summing over the 1/8 domain points for both months

       SELECT CASE (LDAS%LAI)
       
       CASE(2)     ! AVHRR LAI

       IF (ICHAR(LAI_1(CLM(T)%ITYPVEG+1)) .NE. 251 .AND. ICHAR(LAI_1(CLM(T)%ITYPVEG+1)) .NE. 0) THEN
         SUM1 = SUM1 + (ICHAR(LAI_1(CLM(T)%ITYPVEG+1))) * 0.04
         CNT1 = CNT1 + 1
       ENDIF
       IF (ICHAR(SAI_1(CLM(T)%ITYPVEG+1)) .NE. 251 .AND. ICHAR(SAI_1(CLM(T)%ITYPVEG+1)) .NE. 0) THEN
         SUM3 = SUM3 + (ICHAR(SAI_1(CLM(T)%ITYPVEG+1))) * 0.04
         CNT3 = CNT3 + 1
       ENDIF   
       IF (ICHAR(LAI_2(CLM(T)%ITYPVEG+1)) .NE. 251 .AND. ICHAR(LAI_2(CLM(T)%ITYPVEG+1)) .NE. 0) THEN
         SUM2 = SUM2 + (ICHAR(LAI_2(CLM(T)%ITYPVEG+1))) * 0.04
	 CNT2 = CNT2 + 1
       ENDIF
       IF (ICHAR(SAI_2(CLM(T)%ITYPVEG+1)) .NE. 251 .AND. ICHAR(SAI_2(CLM(T)%ITYPVEG+1)) .NE. 0) THEN
         SUM4 = SUM4 + (ICHAR(SAI_2(CLM(T)%ITYPVEG+1))) * 0.04
	 CNT4 = CNT4 + 1
       ENDIF
       
       CASE(3)     ! MODIS LAI
       IF (ICHAR(LAI_1(CLM(T)%ITYPVEG+1)) .LT. 200) THEN
         SUM1 = SUM1 + (ICHAR(LAI_1(CLM(T)%ITYPVEG+1))) * 0.10
	 CNT1 = CNT1 + 1
       ENDIF
       IF (ICHAR(SAI_1(CLM(T)%ITYPVEG+1)) .LT. 200) THEN
	 SUM3 = SUM3 + (ICHAR(SAI_1(CLM(T)%ITYPVEG+1))) * 0.10
         CNT3 = CNT3 + 1
       ENDIF
       IF (ICHAR(LAI_2(CLM(T)%ITYPVEG+1)) .LT. 200) THEN
         SUM2 = SUM2 + (ICHAR(LAI_2(CLM(T)%ITYPVEG+1))) * 0.10
	 CNT2 = CNT2 + 1
       ENDIF
       IF (ICHAR(SAI_2(CLM(T)%ITYPVEG+1)) .LT. 200) THEN
	 SUM4 = SUM4 + (ICHAR(SAI_2(CLM(T)%ITYPVEG+1))) * 0.10
         CNT4 = CNT4 + 1
       ENDIF
       CASE DEFAULT
	PRINT*, "Not a valid LAI Domain"
	STOP
				 
       END SELECT
				      
       SUM5 = SUM5 + TOP_1(CLM(T)%ITYPVEG+1)
       CNT5 = CNT5 + 1
       SUM6 = SUM6 + TOP_2(CLM(T)%ITYPVEG+1)
       CNT6 = CNT6 + 1
       SUM7 = SUM7 + BOT_1(CLM(T)%ITYPVEG+1)
       CNT7 = CNT7 + 1
       SUM8 = SUM8 + BOT_2(CLM(T)%ITYPVEG+1)
       CNT8 = CNT8 + 1
	      	      				     
      ENDDO
     ENDDO
     
!=== Compute averages for the vegetation type represented by tile

     IF (CNT1 .NE. 0) THEN
       CLM(T)%LAI_T1_F = SUM1 / CNT1
     ELSE
       CLM(T)%LAI_T1_F = 0
     ENDIF
     IF (CNT2 .NE. 0) THEN
       CLM(T)%LAI_T2_F = SUM2 / CNT2
     ELSE
       CLM(T)%LAI_T2_F = 0       
     ENDIF
     IF (CNT3 .NE. 0) THEN
       CLM(T)%SAI_T1_F = SUM3 / CNT3
     ELSE
       CLM(T)%SAI_T1_F = 0
     ENDIF
     IF (CNT4 .NE. 0) THEN
       CLM(T)%SAI_T2_F = SUM4 / CNT4
     ELSE
       CLM(T)%SAI_T2_F = 0
     ENDIF
     
     IF (CNT5 .NE. 0) THEN
       CLM(T)%TOP_T1_F = SUM5 / CNT5
     ELSE
       CLM(T)%TOP_T1_F = 0
     ENDIF
     IF (CNT6 .NE. 0) THEN
       CLM(T)%TOP_T2_F = SUM6 / CNT6
     ELSE
       CLM(T)%TOP_T2_F = 0
     ENDIF
     IF (CNT7 .NE. 0) THEN
       CLM(T)%BOT_T1_F = SUM7 / CNT7
     ELSE
       CLM(T)%BOT_T1_F = 0
     ENDIF
     IF (CNT8 .NE. 0) THEN
       CLM(T)%BOT_T2_F = SUM8 / CNT8
     ELSE
       CLM(T)%BOT_T2_F = 0
     ENDIF
			     						    							  
!=== Do NLDAS section (no aggregation)

   ELSE

     MLAT = (CLM(T)%LATDEG - (-59.9375)) / 0.125 + 1
     MLON = (CLM(T)%LONDEG - (-179.9375)) / 0.125 + 1
     LINE = (MLAT - 1)*2880 + MLON
     
!=== Read in data for both LAI and DSAI
!=== IF/ENDIF due to the fact climatology has different format then the actual monthly files
      READ(10,REC=LINE) LAT1, LON1, LAI_1
      READ(12,REC=LINE) LAT1, LON1, SAI_1
      READ(11,REC=LINE) LAT2, LON2, LAI_2
      READ(13,REC=LINE) LAT2, LON2, SAI_2
      READ(50,REC=LINE) LAT1, LON1, TOP_1
      READ(51,REC=LINE) LAT1, LON1, TOP_2
      READ(52,REC=LINE) LAT1, LON1, BOT_1
      READ(53,REC=LINE) LAT1, LON1, BOT_2

!=== Scale to real physical values 
      SELECT CASE(LDAS%LAI)
      CASE(2)   ! AVHRR LAI
       CLM(T)%LAI_T1_F = ICHAR(LAI_1(CLM(T)%ITYPVEG+1)) * 0.04
       CLM(T)%SAI_T1_F = ICHAR(SAI_1(CLM(T)%ITYPVEG+1)) * 0.04
       CLM(T)%LAI_T2_F = ICHAR(LAI_2(CLM(T)%ITYPVEG+1)) * 0.04
       CLM(T)%SAI_T2_F = ICHAR(SAI_2(CLM(T)%ITYPVEG+1)) * 0.04
      CASE(3)
       CLM(T)%LAI_T1_F = ICHAR(LAI_1(CLM(T)%ITYPVEG+1)) * 0.10
       CLM(T)%SAI_T1_F = ICHAR(SAI_1(CLM(T)%ITYPVEG+1)) * 0.10
       CLM(T)%LAI_T2_F = ICHAR(LAI_2(CLM(T)%ITYPVEG+1)) * 0.10
       CLM(T)%SAI_T2_F = ICHAR(SAI_2(CLM(T)%ITYPVEG+1)) * 0.10
      CASE DEFAULT
       PRINT*, "Invalid domain for LAI data"
       STOP
      END SELECT			    
      
      CLM(T)%TOP_T1_F = TOP_1(CLM(T)%ITYPVEG+1)
      CLM(T)%TOP_T2_F = TOP_2(CLM(T)%ITYPVEG+1)
      CLM(T)%BOT_T1_F = BOT_1(CLM(T)%ITYPVEG+1)
      CLM(T)%BOT_T2_F = BOT_2(CLM(T)%ITYPVEG+1)
           
   ENDIF

  ENDDO
  
  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
  CLOSE(50)
  CLOSE(51)
  CLOSE(52)
  CLOSE(53)
	
  ENDIF
  
      IF (LDAS%LAI .GE. 2) THEN
       CALL CANHTSET(pthtop,pthbot)
       DO T=1,LDAS%NCH
        IF (CLM(T)%TOP_T1_F .EQ. 0.0 .AND. (CLM(T)%LAI_T1_F .NE. 0.0 &
	.OR. CLM(T)%SAI_T1_F .NE. 0.0)) THEN
         CLM(T)%TOP_T1_F = pthtop(clm(t)%itypveg)
         CLM(T)%BOT_T1_F = pthbot(clm(t)%itypveg)
        ENDIF
        IF (CLM(T)%TOP_T2_F .EQ. 0.0 .AND. (CLM(T)%LAI_T2_F .NE. 0.0 &
	.OR. CLM(T)%SAI_T2_F .NE. 0.0)) THEN
         CLM(T)%TOP_T2_F = pthtop(clm(t)%itypveg)
         CLM(T)%BOT_T2_F = pthbot(clm(t)%itypveg)
        ENDIF
       ENDDO      
      ENDIF
  
!=== Determine weights between months
      WT1= (TIME2-LDAS%TIME)/(TIME2-TIME1)
      WT2= (LDAS%TIME-TIME1)/(TIME2-TIME1)

!=== Assign interpolated LAI and DSAI values to the CLM variable names used in CLM main
     DO T=1,LDAS%NCH

       CLM(T)%TLAI = WT1 * CLM(T)%LAI_T1_F + WT2 * CLM(T)%LAI_T2_F
       CLM(T)%TSAI = WT1 * CLM(T)%SAI_T1_F + WT2 * CLM(T)%SAI_T2_F
       CLM(T)%HTOP = WT1 * CLM(T)%TOP_T1_F + WT2 * CLM(T)%TOP_T2_F
       CLM(T)%HBOT = WT1 * CLM(T)%BOT_T1_F + WT2 * CLM(T)%BOT_T2_F
       IF (CLM(T)%ITYPVEG .EQ. 12) THEN
         CLM(T)%TLAI=0.0
	 CLM(T)%TSAI=0.0
	 CLM(T)%HTOP=0.0
	 CLM(T)%HBOT=0.0		  
       ENDIF

     ENDDO
     
END SUBROUTINE CLM2LAIREAD

!==============================================================
!
!  DESCRIPTION: This subroutine puts together AVHRR file name
!==============================================================

subroutine avhrr_file_2 (NAME,NAME2,NAME3,NAME4,NAME5,NAME6,NAME7,NAME8, &
                   NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
                   avhrrdir, cyr1, cyr2, cmo1, cmo2 )

  implicit none

!==== Local Variables=======================

  character(len=80) :: name, name2, name3, name4, avhrrdir 
  character(len=80) :: name5, name6, name7, name8
  character(len=80) :: name9,  name10, name11, name12
  character(len=80) :: name13, name14, name15, name16
  character(len=4)  :: cyr1, cyr2
  character(len=2)  :: cmo1, cmo2
  integer :: i, c, flag
  character*1 :: fbase(80), fbase_2(80), fbase_3(80), fbase_4(80)
  character*1 :: fdir(7), fdir_2(7), fdir_3(7), fdir_4(7)
  character*1 :: fsubs(9), fsubs_2(9), fsubs_3(9), fsubs_4(9)
  character*1 :: fsubs_5(9), fsubs_6(9), fsubs_7(9), fsubs_8(9)
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

!  print*, name9
!  print*, name10
!  print*, name11
!  print*, name12
!  print*, name13
!  print*, name14
!  print*, name15
!  print*, name16

  return
			       
 end subroutine avhrr_file_2

!==============================================================
!
!  DESCRIPTION: This subroutine puts together NCAR LAI/SAI/TOP/BOT file names
!==============================================================

subroutine ncarlai_clm2 (NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
                   avhrrdir, cyr1, cyr2, cmo1, cmo2 )

  implicit none

!==== Local Variables=======================

  character(len=80) :: name9,name10,name11,name12,name13
  character(len=80) :: avhrrdir,name14,name15,name16
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
                       
  write(90, 67, rec=1) '_NCLM2LAI_0.125'
  write(91, 67, rec=1) '_NCLM2LAI_0.125'
  write(92, 67, rec=1) '_NCLM2LAI_0.125'
  write(93, 67, rec=1) '_NCLM2LAI_0.125'
  read (90, 92, rec=1) (fsubsn(i), i=1,15)
  read (91, 92, rec=1) (fsubsn_2(i), i=1,15)
  read (92, 92, rec=1) (fsubsn_3(i), i=1,15)
  read (93, 92, rec=1) (fsubsn_4(i), i=1,15)
  write(90, 67, rec=1) '_NCLM2SAI_0.125'
  write(91, 67, rec=1) '_NCLM2SAI_0.125'
  write(92, 67, rec=1) '_NCLM2SAI_0.125'
  write(93, 67, rec=1) '_NCLM2SAI_0.125'
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
  
!  print*, name9
!  print*, name10
!  print*, name11
!  print*, name12
!  print*, name13
!  print*, name14
!  print*, name15
!  print*, name16

  return

 end subroutine ncarlai_clm2

!==============================================================
!
!  DESCRIPTION: This subroutine puts together NCAR TOP/BOT file names
!==============================================================

subroutine ncarcanht_clm2 (NTOP1,NTOP2,NBOT1,NBOT2, &
                   avhrrdir, cyr1, cyr2, cmo1, cmo2 )

  implicit none

!==== Local Variables=======================

  character(len=80) :: avhrrdir
  character(len=80) :: ntop1, nbot1, ntop2, nbot2
  character(len=4)  :: cyr1, cyr2
  character(len=2)  :: cmo1, cmo2
  integer :: i, c, flag
  character*1 :: fbase(80), fbase_2(80), fbase_3(80), fbase_4(80)
  character*1 :: fdir(7), fdir_2(7), fdir_3(7), fdir_4(7)
  character*1 :: fsubsn(15), fsubsn_2(15), fsubsn_3(15), fsubsn_4(15)

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

  write(90, 66, rec=1) '/CLIM', cmo1
  read(90, 99, rec=1) fdir
  write(91, 66, rec=1) '/CLIM', cmo2
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

  write(90, 67, rec=1) '_NCLM2TOP_0.125'
  write(91, 67, rec=1) '_NCLM2TOP_0.125'
  write(92, 67, rec=1) '_NCLM2BOT_0.125'
  write(93, 67, rec=1) '_NCLM2BOT_0.125'
  read (90, 92, rec=1) (fsubsn(i), i=1,15)
  read (91, 92, rec=1) (fsubsn_2(i), i=1,15)
  read (92, 92, rec=1) (fsubsn_3(i), i=1,15)
  read (93, 92, rec=1) (fsubsn_4(i), i=1,15)

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
  read(90, 93, rec=1) ntop1
  read(91, 93, rec=1) ntop2
  read(92, 93, rec=1) nbot1
  read(93, 93, rec=1) nbot2
  
!  print*, ntop1
!  print*, ntop2
!  print*, nbot1
!  print*, nbot2

  close(90)
  close(91)
  close(92)
  close(93)
 
  return

 end subroutine ncarcanht_clm2

!==============================================================
!
!  DESCRIPTION: This subroutine puts together MODIS file name
!==============================================================

subroutine modis_file_2 (NAME9,NAME10,NAME11,NAME12,NAME13,NAME14,NAME15,NAME16, &
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

 end subroutine modis_file_2

