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
! ldas_module.f: 
!
! DESCRIPTION:
!  Module for LDAS variable specification.  This file will contain no
!   tile space or grid space variables.  The LDAS, non-model-specific 
!   grid space variables will be in grid_module.f, and the tile space LDAS
!   non-model-specific variables will be in tile_module.f.  Model specific
!   variables will be in their own specific modules.
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial code
!  4  Apr 2000: Jeffrey Walker; Added some catcment model variables
!  11 Apr 2000: Brian Cosgrove; Elevation correction and Forcing Mask
!               variables added 
!  6  Jun 2000: Jon Radakovich; Updated for new version of CLM
!  23 Feb 2001: Urszula Jambor; Updated for GEOS & GDAS forcing in GLDAS
!  27 Feb 2001: Brian Cosgrove; Added Catchment forcing data variables
!  15 Mar 2001: Jon Gottschalck; Updated for GDAS initialization of Mosaic 
!  12 Apr 2001: Urszula Jambor; Added domain,lsm,& force namefile paramters     
!  30 Apr 2001: Jon Radakovich; Update for PSAS temperature assimilation
!  17 Jul 2001: Jon Gottschalck; Update for global precipitation variables
!  30 Jul 2001: Matt Rodell; Add new soil parameter variables
!  05 Sep 2001: Brian Cosgrove; Add variables for PAR and BRTTMP, remove 
!               1/4 to 1/8 interp variables
!  15 Oct 2001: Jesse Meng; Replace agrmet flag by agrmetsw and agrmetlw
!  27 Nov 2001: Jon Gottschalck; Added variables for AVHRR LAI data 
!  07 Dec 2001: Urszula Jambor; Added LDAS_KGDS array     
!  03 Feb 2002: Jon Gottschalck; Added Koster tilespace variables
!  05 Feb 2002: Brian Cosgrove; Added NLDAS 11 Layer soil class file from Yun Duan
!               and ftype variable to indicate NLDAS ETA forcing source used
!  08 Feb 2002: Urszula Jambor; Added latmax for AGRMET use.
!  15 Apr 2002: Urszula Jambor; Added ECMWF forcing options.
!  28 Apr 2002: Kristi Arsenault; Added NOAH LSM parameters and code
!  30 Jul 2002: Jon Gottschalck; Added variables to handle additional global observed precip data sources
!  25 Sep 2002: Jon Gottschalck; Added MODIS UMD vegetation classification
!  01 Oct 2002: Jon Gottschalck; Added comments for MODIS LAI data      
!  04 Nov 2002: Kristi Arsenault; Added new MAXSNALB and TBOT fields
!  10 Dec 2002: Urszula Jambor; Removed stat flags for AGRMET
!  05 Feb 2003: Jon Gottschalck; Added to include additions for CLM v. 2.0  
!  04 Mar 2003: Urszula Jambor; Added GRIDCHANGE flag to check for 
!               forcing grid changes.   
!  20 Jun 2003: Urszula Jambor; Added SW_TIMEs for 12-hrly ECMWF forcing   
!  24 Jun 2003: Kristi Arsenault; Added flags to interpolate daily albedo
!                                 and greenness fraction files
!=========================================================================

      module ldas_module 

      IMPLICIT NONE
      public ldasdec

      type ldasdec

!=== LDAS Parameters ======================================================
      INTEGER :: NCH          !actual number of tiles

!=== LDAS.crd Parameters: Model non-specific ==============================
      INTEGER :: DOMAIN       !Model domain, (1=NLDAS, 2=GLDAS)
      INTEGER :: LSM          !Land surface model (1=Mosaic,2=CLM1,3=CATCH,4=NOAH,5=CLM2)
      INTEGER :: FORCE        !Forcing data type (1=GDAS,2=GEOS,3=ETA,4=NCEP,5=NASA,6=CATCH,7=ECMWF,8=reanlECMWF)
      INTEGER :: SOIL	      !Soil parameter scheme (1=original veg-based, 2=Reynolds soils)
      INTEGER :: LAI          !LAI data source (1=original, 2=AVHRR satellite data, 3=MODIS satellite data)
      REAL*8  :: LAITIME      !Satellite LAI Time
      CHARACTER*40 :: AVHRRDIR  !AVHRR LAI Data directory
      CHARACTER*40 :: MODISDIR  !MODIS LAI Data directory

      INTEGER :: EXPCODE      !3 Digit experiment code 
      INTEGER :: NC           !Number of Columns in Grid
      INTEGER :: NR           !Number of Rows in Grid
      INTEGER :: NT           !Number of Vegetation Types 
      INTEGER :: RMOS         !Run Mosaic (0=no,1=yes)
      INTEGER :: RCAT         !Run Catchment Model (0=no,1=yes)
      INTEGER :: RCLM1        !Run CLM1 (0=no,1=yes)
      INTEGER :: RCLM2        !Run CLM2 (0=no,1=yes)
      INTEGER :: RNOAH        !Run NOAH LSM (0=no,1=yes)
      INTEGER :: WFOR         !Write Forcing (0=no,1=yes)
      INTEGER :: WTIL         !Write tile space data (0=no, 1=yes)
      INTEGER :: WHDF         !Write HDF output files (0=no,1=yes)
      INTEGER :: WGRB         !Write Grib output files (0=no,1=yes)
      INTEGER :: WBIN         !Write Binary output files (0=no,1=yes)
      INTEGER :: STARTCODE    !0=restart date, 1=card date
      INTEGER :: SSS          !Starting Second 
      INTEGER :: SDOY         !Starting Day of Year 
      INTEGER :: SMN          !Starting Minute 
      INTEGER :: SHR          !Starting Hour 
      INTEGER :: SDA          !Starting Day 
      INTEGER :: SMO          !Starting Month 
      INTEGER :: SYR          !Starting Year  
      INTEGER :: ENDCODE      !0=realtime, 1=specific date
      INTEGER :: ESS          !Ending Second
      INTEGER :: EMN          !Ending Minute
      INTEGER :: EDOY         !Ending Day of Year
      INTEGER :: EHR          !Ending Hour
      INTEGER :: EDA          !Ending Day
      INTEGER :: EMO          !Ending Month
      INTEGER :: EYR          !Ending Year
      INTEGER :: LDAS_KGDS(200)

      INTEGER :: TS           !Timestep (seconds) 
      INTEGER :: TSCOUNT      !Timestep Count
      INTEGER :: NUMOUTF      !Counts number of output times for forcing data
      INTEGER :: NUMOUTM      !Counts number of output times for Mosaic
      INTEGER :: NUMOUTC1     !Counts number of output times for CLM1
      INTEGER :: NUMOUTC2     !Counts number of output times for CLM2
      INTEGER :: NUMOUTCT     !Counts number of output times for CAT
      INTEGER :: NUMOUTNH     !Counts number of output times for NOAH
      REAL :: UDEF            !Undefined value
      REAL :: WRITEINTM       !Mosaic Output Interval (hours)
      REAL :: WRITEINTC1      !CLM1 Output Interval (hours)
      REAL :: WRITEINTC2      !CLM2 Output Interval (hours)      
      REAL :: WRITEINTCT      !CAT Output Interval (hours)
      REAL :: WRITEINTN       !NOAH Output Interval (hours)
      REAL :: WRITEINTF       !Forcing Output Interval (hours)
      CHARACTER*40 :: IDIR    !Input Data Base Directory
      CHARACTER*40 :: ODIR    !Output Data Base Directory
      INTEGER :: FETA,FNCEP   !Forcing Data Base Source Options 
      INTEGER :: FNASA        !Forcing Data Base Source Options
      INTEGER :: FGDAS,FGEOS  !Forcing Data Base Source Options
      INTEGER :: FCATCH       !Forcing data base 
      INTEGER :: FECMWF,FR_ECMWF !Forcing Options
      CHARACTER*40 :: CATCHDIR !Catchment data directory
      CHARACTER*40 :: NCEPDIR  !NCEP Forcing Directory
      CHARACTER*40 :: NASADIR  !NASA Forcing Directory
      CHARACTER*80 :: EDASDIR  !EDAS Forcing Directory
      CHARACTER*40 :: GDASDIR  !GDAS Forcing Directory
      CHARACTER*40 :: NRLDIR   !NRL Precip directory
      CHARACTER*40 :: PERSDIR  ! PERSIANN Precip directory
      CHARACTER*40 :: HUFFDIR  ! HUFFMAN Precip directory
      CHARACTER*40 :: CMAPDIR  ! CMAP Precip directory
      CHARACTER*40 :: GEOSDIR  !GEOS Forcing Directory
      CHARACTER*40 :: ECMWFDIR !ECMWF Forcing Directory
      CHARACTER*40 :: ETA3HRDIR  !ETA 3HR Forecast Directory
      CHARACTER*40 :: ETA6HRDIR  !ETA 6HR Forecast Directory
      CHARACTER*22 :: EVT1    !valid time of Before ETA data
      CHARACTER*22 :: EVT2    !valid time of After ETA data
      CHARACTER*40 :: PINKDIR !Pinker Hourly Radiation Forcing Directory
      CHARACTER*40 :: PINKDIROLD !Pinker Hourly Radiation Forcing Directory
      CHARACTER*80 :: RADBCDIR !Hourly Radiation Forcing Bias Correction Directory
      CHARACTER*40 :: NDDIR   !NESIDS Hourly Radiation Forcing Directory
      CHARACTER*40 :: AGRMDIR !AGRMET Hourly Radiation Forcing Directory
      CHARACTER*40 :: BRTTMPDIR !Pinker Hourly Radiation Forcing Directory
      CHARACTER*40 :: PARDIR  !PAR Hourly Forcing Directory
      CHARACTER*40 :: BRTTMPDIROLD !Pinker Hourly Radiation Forcing Directory
      CHARACTER*40 :: PARDIROLD  !PAR Hourly Forcing Directory
      CHARACTER*80 :: STAGEIV !STAGE IV Hourly Precipitation Forcing Directory
      CHARACTER*80 :: HIGGINS !Higgins 24hr Gauge Directory
      CHARACTER*80 :: HIGGINSUNI !Higgins Unified precip data Directory
      CHARACTER*80 :: HPD      ! HPD hourly precip data Directory
      CHARACTER*80 :: CMORPH  !cmorph  precip data Directory
      CHARACTER*80 :: MERGEH  !HIGGINS GAGE W/STAGEIV Wts merged product DIR
      CHARACTER*40 :: PRECIPWEIGHT !Name of precip weighting file  
      INTEGER :: SHORTFLAG    !Shortwave radiation source flag
                              !0=No radiation
                              !1=Instantaneous SW
                              !2=Time Averaged SW
      INTEGER :: LONGFLAG     !Longwave radiation source flag
                              !0=No radiation
                              !1=Instantaneous LW
                              !2=Time Averaged LW
      INTEGER :: REGENERATE   !0=DO NOT regenerate precip files if already exist
                              !1=REGENERATE precip files even if already exist
      INTEGER :: CPC          !0=DO NOT use CPC 24hr Gage Data in merged precip product
                              !1=Use CPC 24hr Gage Data 1st in merged precip product
                              !2=Use second in merged precip product
      INTEGER :: UNIFIED      !0=DO NOT use Higgins Unified data in merged precip product
                              !1=Use Higgins Unified 1st in merged precip product
                              !2=Use second in merged precip product
      INTEGER :: PCPCOUNT     !0=Hourly precip files have NOT recently
                              !been regenerated in a REALTIME run
                              !1=Hourly precip files HAVE recently
                              !been regenerated in a REALTIME run
      INTEGER :: NF           !Number of forcing variables
      INTEGER :: NMIF         !Number of forcing variables for model initialization option                        
      INTEGER :: PRECSOR(4)   !Precip Data Source Option
      INTEGER :: GPCPSRC(4)   !Global precipitation flags NRL,PERSIANN,HUFFMAN,CMAP
      INTEGER :: HIERFLAG     !Precipitation hierarchy flag (highest ranked precip data source)
      INTEGER :: MASTERPCP(4) !Program flow precipitation source flags
      INTEGER :: TEMPADJ      !Adjust temperature for elevation differences
      INTEGER :: PRESADJ      !Adjust pressure for elevation differences
      INTEGER :: HUMIDADJ     !Adjust humidity for elevation differences
      INTEGER :: LWRADADJ     !Adjust long wave radiation for elevation differences
      INTEGER :: PINKER       !Rank of Pinker SWdwn data source
			      !0=Don't Use, 1=Use First, 2=Use Second
      INTEGER :: NESDIS       !Rank of NESDIS SWdwn data source
		              !0=Don't Use, 1=Use First, 2=Use Second
      INTEGER :: AGRMETSW     !Rank of AGRMET SWdwn data source
			      !0=Don't Use, 1=Use
      INTEGER :: AGRMETLW     !Rank of AGRMET LWdwn data source
			      !0=Don't Use, 1=Use
      INTEGER :: PAR          !Flag for PAR data 
                              !0=Don't Use, 1=Use
      INTEGER :: BRTTMP       !1=Read in NESDIS BRT TEMP, 0=Don't read in
      INTEGER :: PRECIPMASK   !0=Don't apply weighting mask to merged precip product
                              !1=Apply weighting mask to merged precip product
      INTEGER :: MAXT         !Maximum tiles per grid  
      REAL :: MINA            !Min grid area for tile (%)
      INTEGER :: VCLASS       !Vegetation Classification (1=UMD)
      CHARACTER*40 :: MFILE   !land/water mask file for MODELLING (AVHRR)
      CHARACTER*40 :: M2FILE  !land/water mask file for MODELLING (MODIS)
      CHARACTER*40 :: FMFILE  !land/water mask file for FORCING DATA (AVHRR)
      CHARACTER*40 :: FM2FILE !land/water mask file for FORCING DATA (MODIS)
      CHARACTER*40 :: ELEVFILE !LDAS-EDAS elevation difference file
      CHARACTER*40 :: NARRELEV !NARR elevation file
      CHARACTER*40 :: CATELEVFILE !LDAS-CATCHMENT elevation difference file
      CHARACTER*40 :: INTERPA !Interpolation file awips212 to 1/8th, FOR MODELLING
      CHARACTER*40 :: INTERPAF!Interpolation file awips212 to 1/8th, FOR FORCING DATA
      CHARACTER*40 :: INTERPP !Interpolation file 1/2 to 1/8th, FOR MODELLING
      CHARACTER*40 :: INTERPPF!Interpolation file 1/2 to 1/8th, FOR FORCING DATA
      CHARACTER*40 :: VFILE   !vegetation classification file (AVHRR)
      CHARACTER*40 :: V2FILE  !vegetation classification file (MODIS)
      CHARACTER*40 :: SFILE   !soil classification file
      CHARACTER*40 :: TXFILE  !Name of LDAS 11 Layer Soil Class file from Yun Duan
      CHARACTER*40 :: SAFILE  !Sand fraction map file
      CHARACTER*40 :: SIFILE  !Silt fraction map file
      CHARACTER*40 :: CLFILE  !Clay fraction map file
      CHARACTER*40 :: POFILE  !Porosity map file
      CHARACTER*40 :: SLFILE  !Slope map file
      CHARACTER*40 :: ISCFILE !Soil color map file
      CHARACTER*40 :: DFILE   !runtime diagnostics file
      CHARACTER*40 :: FFILE   !runtime forcing source file
      CHARACTER*40 :: EVTFILE !ETA Valid Time File
      CHARACTER*40 :: KVFILE  ! Koster tilespace file
      CHARACTER*40 :: EMASKFILE !1/2deg Reanal-ECMWF Land-Sea Mask File
      INTEGER :: REMASK1D(720*360)  ! 1D Reanal-ECMWF Land-Sea Mask
      INTEGER :: KOSTER ! Flag to use Koster tilespace at 2.0x2.5 degrees
      INTEGER :: NKTYPE ! Number of Mosaic vegetation types

!=== LDAS.crd Parameters: Mosaic =======================================
      CHARACTER*40 :: MOS_RFILE  !Mosaic Active Restart File
      CHARACTER*40 :: MOS_MFILE  !Mosaic model init. restart file
      CHARACTER*40 :: MOS_VFILE  !Mosaic Static Vegetation Parameter File
      CHARACTER*40 :: MOS_MVFILE !Mosaic Monthly Vegetation Parameter File   
      CHARACTER*40 :: MOS_KVFILE !Mosaic Static Koster Vegetation Parameter File
      CHARACTER*40 :: MOS_KMVFILE !Mosaic Monthly Koster Vegetation Parameter File
      CHARACTER*40 :: MOS_SFILE  !Mosaic Soil Parameter File
      CHARACTER*40 :: MOS_KSFILE !Mosaic Koster Soil Parameter File
      INTEGER :: MOS_IC          !Mosaic Initial Condition Source
      REAL :: MOS_ISM            !Mosaic Initial Soil Moisture (m3/m3)
      REAL :: MOS_IT             !Mosaic Initial Soil Temperature (K)
      INTEGER :: MOS_SMDA        !Mosaic SM Assimilation Option
      INTEGER :: MOS_TDA         !Mosaic Temperature Assimilation Option
      INTEGER :: MOS_SDA         !Mosaic Snow Assimilation Option
      INTEGER :: MOS_NVEGP       !Number of static vegetation parameters
      INTEGER :: MOS_NMVEGP      !Number of monthly vegetation parameters
      INTEGER :: MOS_NSOILP      !Number of static soil parameters
      INTEGER :: MOS_NRET        !Number of Mosaic output variables

!=== LDAS.crd Parameters: Catchment ====================================
      CHARACTER*40:: CAT_RFILE   !Catchment Active Restart File
      CHARACTER*40 :: CAT_VFILE  !Catchment Static Vegetation Parameter File
      CHARACTER*40 :: CAT_MVFILE !Catchment Month Vegetation Parameter File   
      CHARACTER*40 :: CAT_SFILE  !Catchment Soil Parameter File
      INTEGER :: CAT_IC          !Catchment Initial Condition Source
      REAL :: CAT_ISM            !Catchment Initial Soil Moisture (m3/m3)
      REAL :: CAT_IT             !Catchment Initial Soil Temperature (K)
      INTEGER :: CAT_SMDA        !Catchment SM Assimilation Option
      INTEGER :: CAT_TDA         !Catchment Temperature Assimilation Option
      INTEGER :: CAT_SDA         !Catchment Snow Assimilation Option

!=== LDAS.crd Parameters: CLM ==========================================
      CHARACTER*40 :: CLM1_RFILE  !CLM Active Restart File
      CHARACTER*40 :: CLM1_VFILE  !CLM Vegetation Tile Specification File
      CHARACTER*40 :: CLM1_MVFILE !CLM Vegetation Type Parameter File 
      CHARACTER*40 :: CLM1_PFILE  !CLM 1D Parameter Output File
      CHARACTER*40 :: CLM1_SFILE  !CLM Soil Parameter File
      CHARACTER*40 :: CLM1_CFILE  !CLM Constant Parameter File
      CHARACTER*40 :: CLM1_OFILE  !CLM Output File
      INTEGER :: CLM1_IC          !CLM Initial Condition Source
      REAL :: CLM1_ISM            !CLM Initial Soil Moisture (m3/m3)
      REAL :: CLM1_IT             !CLM Initial Soil Temperature (K)
      REAL :: CLM1_ISCV           !CLM Initial Snow Mass (kg/m2)
      INTEGER :: CLM1_SMDA        !CLM SM Assimilation Option
      INTEGER :: CLM1_TDA         !CLM Temperature Assimilation Option
      INTEGER :: CLM1_SDA         !CLM Snow Assimilation Option

!=== LDAS.crd Parameters: NOAH ==========================================
      CHARACTER*40 :: NOAH_RFILE   !NOAH Active Restart File
      CHARACTER*40 :: NOAH_MFILE   !NOAH model init. restart file
      CHARACTER*40 :: NOAH_VFILE   !NOAH Static Vegetation Parameter File
      CHARACTER*40 :: NOAH_SFILE   !NOAH Soil Parameter File
      CHARACTER*40 :: NOAH_MGFILE  !NOAH Monthly Veg. Green Frac.
      CHARACTER*40 :: NOAH_ALBFILE !NOAH Quart. Snow-free albedo
      CHARACTER*40 :: NOAH_NMXSNAL !NOAH NLDAS 0.125 deg max snow albedo
      CHARACTER*40 :: NOAH_GMXSNAL !NOAH GLDAS 0.25 max snow albedo
      CHARACTER*40 :: NOAH_SMXSNAL !NOAH GLDAS 2x2.5 deg max snow albedo
      CHARACTER*40 :: NOAH_OMXSNAL !NOAH GLDAS 1.0 max snow albedo
      CHARACTER*40 :: NOAH_EMXSNAL !NOAH GLDAS 0.5 max snow albedo
      CHARACTER*40 :: NOAH_NTBOT   !NOAH NLDAS 0.125 Bottom Temp
      CHARACTER*40 :: NOAH_GTBOT   !NOAH GLDAS 0.25 Bottom Temp
      CHARACTER*40 :: NOAH_STBOT   !NOAH GLDAS 2x2.5Bottom Temp
      CHARACTER*40 :: NOAH_OTBOT   !NOAH GLDAS 1.0 Bottom Temp
      CHARACTER*40 :: NOAH_ETBOT   !NOAH GLDAS 0.5 Bottom Temp

      INTEGER :: NOAH_IC           !NOAH Initial Condition Source
      REAL :: NOAH_ISM             !NOAH Initial Soil Moisture (m3/m3)
      REAL :: NOAH_IT              !NOAH Initial Soil Temperature (K)
      INTEGER :: NOAH_SMDA         !NOAH SM Assimilation Option
      INTEGER :: NOAH_SDA          !NOAH Snow Assimilation Option
      INTEGER :: NOAH_TDA          !NOAH Temperature Assimilation Option
      INTEGER :: NOAH_NVEGP        !Number of static vegetation parameters
      INTEGER :: NOAH_NMVEGP       !Number of monthly vegetation parameters
      INTEGER :: NOAH_NSOILP       !Number of static soil parameters
      INTEGER :: NOAH_ZST          !Number of Zobler soil classes
      INTEGER :: NOAH_NRET         !Number of NOAH output variables
      INTEGER :: NOAH_ALBTIME      !Time flag to update albedo files
      REAL*8  :: NOAH_GFRACTIME    !Time flag to update gfrac files 
      INTEGER :: NOAH_ALBDCHK      !Flag to interpolate alb files
      INTEGER :: NOAH_GFRACDCHK    !Flag to interpolate gfrac values

!=== LDAS.crd Parameters: CLM2 ==========================================
      CHARACTER*40 :: CLM2_RFILE  !CLM2 Active Restart File
      CHARACTER*40 :: CLM2_VFILE  !CLM2 Vegetation Tile Specification File
      CHARACTER*40 :: CLM2_PFILE  !CLM2 1D Parameter Output File
      CHARACTER*40 :: CLM2_OFILE  !CLM2 Output File
      INTEGER :: CLM2_IC          !CLM2 Initial Condition Source
      REAL :: CLM2_ISM            !CLM2 Initial Soil Moisture (m3/m3)
      REAL :: CLM2_IT             !CLM2 Initial Soil Temperature (K)
      REAL :: CLM2_ISCV           !CLM2 Initial Snow Mass (kg/m2)
      INTEGER :: CLM2_SMDA        !CLM2 SM Assimilation Option
      INTEGER :: CLM2_TDA         !CLM2 Temperature Assimilation Option
      INTEGER :: CLM2_SDA         !CLM2 Snow Assimilation Option      

!=== END LDAS CARD PARAMETERS ===========================================

!=== LDAS Timing Variables (others above in ldas.crd section)  ==========
      REAL*8 :: TIME               !LDAS Current Model Time in Years
      REAL*8 :: ETIME              !LDAS End Time in Years
      INTEGER :: PDA               !LDAS Previous Timestep Day
      INTEGER :: DOY,YR,MO,DA,HR,MN,SS !LDAS Current Model Timing Variables   
      INTEGER :: ENDTIME           !LDAS Stop (0=continue time looping)
      REAL :: GMT,EGMT,SGMT

!=== MOSAIC variables ===================================================
      REAL :: MOS_STIME            !Mosaic Restart Time 
      INTEGER :: MYR,MMO,MDA       !LDAS Restart Model Timing Variables
      INTEGER :: MHR,MMN,MSS       !LDAS Restart Model Timing Variables

!=== CLM variables ======================================================
      REAL :: CLM_STIME            !CLM Restart Time 
      INTEGER :: CYR,CMO,CDA       !LDAS Restart Model Timing Variables
      INTEGER :: CHR,CMN,CSS       !LDAS Restart Model Timing Variables

!=== Catchment Model variables (K for Koster)  ===========================
      REAL :: KOS_STIME            !Mosaic Restart Time 
      INTEGER :: KYR,KMO,KDA       !LDAS Restart Model Timing Variables
      INTEGER :: KHR,KMN,KSS       !LDAS Restart Model Timing Variables
      INTEGER :: NCAT,NCATM        !Catchment Model Size Variables

!=== NOAH LSM variables ===================================================
      REAL :: NOAH_STIME           !NOAH Restart Time
      INTEGER :: NYR,NMO,NDA       !LDAS Restart Model Timing Variables
      INTEGER :: NHR,NMN,NSS       !LDAS Restart Model Timing Variables

!=== Forcing Data Variables ==============================================
      REAL*8 :: CATCHTIME1,CATCHTIME2
      REAL*8 :: NCEPTIME1,NCEPTIME2
      REAL*8 :: ETATIME1,ETATIME2
      REAL*8 :: GDASTIME1, GDASTIME2
      REAL*8 :: GEOSTIME1,GEOSTIME2
      REAL*8 :: FMODELTIME1,FMODELTIME2
      REAL*8 :: SW_TIME1,SW_TIME2
      REAL*8 :: PINKTIME1,PINKTIME2
      REAL*8 :: RADBCTIME1,RADBCTIME2
      REAL*8 :: NESTIME1,NESTIME2
      REAL*8 :: AGRMTIME1,AGRMTIME2
      REAL*8 :: PARTIME1,PARTIME2
      REAL*8 :: BRTTMPTIME1,BRTTMPTIME2
      REAL*8 :: PRECIPTIME1,PRECIPTIME2
      REAL*8 :: STAGEIVTIME,MERGETIME
      REAL*8 :: NRLTIME,PERSTIME,HUFFTIME,CMAPTIME
      INTEGER :: NCold,NRold        !AWIPS 212 dimensions
      INTEGER :: PSTAT1,PSTAT2
      INTEGER :: RADBCSTAT1,RADBCSTAT2
      INTEGER :: NESSTAT1,NESSTAT2
      INTEGER :: PARSTAT1,PARSTAT2
      INTEGER :: BRTTMPSTAT1,BRTTMPSTAT2
      INTEGER :: FSOURCE(16)        !Forcing Status Vector
      INTEGER :: SKIPINTP1,SKIPINTP2           
      INTEGER :: FTYPE            ! Raw forcing file type indicator
                                  ! 1=EDAS, 2=ETA3hr, 3=ETA6hr
      INTEGER :: latmax           !Per hemisphere, for AGRMET intepolation
      INTEGER :: NRGPCP           !Global precip interpolated grid row number
      INTEGER :: GRIDCHANGE       !Flag to check for forcing grid changes

!=== HDF Output Variables
      INTEGER :: FIDGF,RCGF,FIDTF,RCTF,FIDGC,RCGC,FIDTC,RCTC
      INTEGER :: FIDT3C,RCT3C,FIDGM,RCGM,FIDTM,RCTM
      INTEGER :: FIDGCT,RCGCT
      INTEGER :: YYYYMMDD,HHMMSS
      CHARACTER*200 :: RLATSGF,RLATSTF,RMHDFTF,RMHDFGF
      CHARACTER*200 :: RLATSGM,RLATSTM,RMHDFTM,RMHDFGM
      CHARACTER*200 :: RLATSGC,RLATSTC,RLATST3C,RMHDFTC,RMHDFGC,RMHDFT3C
      CHARACTER*200 :: RLATSGCT,RMHDFGCT
      CHARACTER*80  :: LATSDIAGF,LATSDIAGM,LATSDIAGC,LATSDIAGCT

!=== OUTPUT Variables:
      INTEGER :: MOSopen,FORopen,CLM1open,CLM2open,CATopen,NOAHopen ! Keeps track of opening files

!=== PSAS assimilation variable
      INTEGER :: RPSAS,RBIAS,RIBC,RDBC,RSDBC

!=== END LDAS parameters ============================

       end type
       end module ldas_module

