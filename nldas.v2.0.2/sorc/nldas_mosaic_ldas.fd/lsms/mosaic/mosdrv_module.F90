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
!BOP
!
! !MODULE: mosdrv_module.F90
!
! !DESCRIPTION:
!  Module for runtime specific Mosaic variables
!
! !REVISION HISTORY:
!
! 15 Oct 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
MODULE mosdrv_module
!EOP
  type mosdrvdec
     integer :: mosopen         !Keeps track of opening files
     integer :: numoutm         !Counts number of output times for mosaic
     INTEGER :: MOS_IC          !Mosaic Initial Condition Source
     INTEGER :: MOS_SMDA        !Mosaic SM Assimilation Option
     INTEGER :: MOS_TDA         !Mosaic Temperature Assimilation Option
     INTEGER :: MOS_SDA         !Mosaic Snow Assimilation Option
     INTEGER :: MOS_NVEGP       !Number of static vegetation parameters
     INTEGER :: MOS_NMVEGP      !Number of monthly vegetation parameters
     INTEGER :: MOS_NSOILP      !Number of static soil parameters
     
     CHARACTER*40 :: MOS_RFILE  !Mosaic Active Restart File
     CHARACTER*40 :: MOS_VFILE  !Mosaic Static Vegetation Parameter File
     CHARACTER*40 :: MOS_MVFILE !Mosaic Monthly Vegetation Parameter File
     CHARACTER*40 :: MOS_KVFILE !Mosaic Static Koster Vegetation Parameter File
     CHARACTER*40 :: MOS_KMVFILE !Mosaic Monthly Koster Vegetation Parameter File
     CHARACTER*40 :: MOS_SFILE  !Mosaic Soil Parameter File
     CHARACTER*40 :: MOS_KSFILE !Mosaic Koster Soil Parameter File
     REAL :: MOS_ISM            !Mosaic Initial Soil Moisture (m3/m3)
     REAL :: MOS_IT             !Mosaic Initial Soil Temperature (K)
     REAL :: WRITEINTM       !Mosaic Output Interval (hours)
  end type mosdrvdec
end MODULE mosdrv_module
