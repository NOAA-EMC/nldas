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
! !ROUTINE: readkpds.F90
! 
! !DESCRIPTION: 
!   Reads the kpds array from the grib table
!
! !REVISION HISTORY:
!
! 23 Oct 2003; Sujay Kumar; Initial Version
!
! !INTERFACE:
  subroutine readkpds(ftn, kpds)
!
! !INTERFACE:
    use noah_varder, only : noahdrv
!EOP
    integer     :: kpds(25),ftn

    read (ftn, 15) kpds(5), kpds(6), kpds(7), kpds(14), &
         kpds(15), kpds(16), kpds(22)
    if(kpds(16).ne.0) kpds(15)=noahdrv%writeintn
    
15   format (29x, 7i6)     
  end subroutine readkpds
