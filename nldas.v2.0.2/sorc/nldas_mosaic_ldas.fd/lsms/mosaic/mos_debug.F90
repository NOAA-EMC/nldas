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
! !ROUTINE: mos_totinit.F90
!
! !DESCRIPTION:
!  Initialize Mosaic output arrays
!
! !REVISION HISTORY:
!
!  1 Aug 2003  Sujay Kumar  Initial Specification
!
! !INTERFACE:
      subroutine mos_debug()
! !USES:
      use mos_varder      ! Mosaic LSM module
      use tile_spmdMod
      use lisdrv_module, only : lis
!EOP      
  IMPLICIT NONE

!=== End Variable List ===================================================
  integer t
!BOC 

     open(unit=31,file='mosold.bin',form='unformatted') 

     write(31) mos%VEGT !int field
     write(31) mos%VEGT !int field
     do t=1,24 
        write(31) mos%VEGP(t)    !Static vegetation parameter values, dim(MOS_NVEGP)
     enddo
     do t=1,6
        write(31) mos%VEGIP(t)   !Interpolated from monthly parameters (MOS_NMVEGP)
     enddo
     do t=1,10
        write(31) mos%SOILP(t)   !Static soil parameter values, dim(NOS_NSOILP)
     enddo

     write(31) mos%LAI           !AVHRR Derived LAI
     write(31) mos%GREEN         !AVHRR Derived GREENNESS
     write(31) mos%DSAI          !AVHRR Derived DSAI
     write(31) mos%DTCANAL                !MOSAIC Change in Temperature based on Analysis

!=== LDAS-MOSAIC States ===============================================
     write(31) mos%CT                     !MOSAIC Canopy/Soil Temperature 
     write(31) mos%QA                     !MOSAIC Canopy Humidity
     write(31) mos%ICS                    !MOSAIC Interception Canopy Storage
     write(31) mos%SNOW                   !MOSAIC Snow Depth
     write(31) mos%SoT    

     do t=1,3                !MOSAIC Deep Soil Temperaure
        write(31) mos%SoWet(t)      !MOSAIC Soil Wetness (3 layers)
     enddo
!=== Analysis and Bias Correction Variables ===========================

      close(31)

      STOP

!EOC
end subroutine mos_debug
