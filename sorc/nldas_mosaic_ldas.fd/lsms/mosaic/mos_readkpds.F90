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
! !ROUTINE: mos_readkpds.F90
! 
! !DESCRIPTION: 
!   Reads the ipds array from the grib table
!
! !REVISION HISTORY:
!
! 23 Oct 2003; Sujay Kumar; Initial Version
! 17 Oct 2013; Youlong Xia; Revised for grib2 0utput - read IPDS(13)
!
! !INTERFACE:
  subroutine mos_readkpds(ftn, ipds)
!
! !INTERFACE:
    use mos_varder, only : mosdrv
!EOP
    integer     :: ipds(13),ftn,i
    character*19 :: vname  
    read (ftn,*) vname,(ipds(i),i=1,13)
    end subroutine mos_readkpds
