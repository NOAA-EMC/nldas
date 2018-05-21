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
! !ROUTINE: hyssib_mapvegc.F90
!
! !DESCRIPTION: 
!  This subroutine converts the UMD classes to the SIB classes 
!  used by SSIB LSM (v 2.5).
!  (Originally from Dag Lohmann at NCEP)
!
! !REVISION HISTORY:
! 28 Apr 2002: K Arsenault, Added SSIB LSM to LDAS
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
      SUBROUTINE HYSSIB_MAPVEGC(VEGT)
!EOP
      implicit none

!  Local Variables
      INTEGER :: VEGT, SIBVEG

!  Convert UMD Classes to SIB Classes.
      IF (VEGT.EQ.1)  SIBVEG = 4
      IF (VEGT.EQ.2)  SIBVEG = 1
      IF (VEGT.EQ.3)  SIBVEG = 5
      IF (VEGT.EQ.4)  SIBVEG = 2
      IF (VEGT.EQ.5)  SIBVEG = 3
      IF (VEGT.EQ.6)  SIBVEG = 3
      IF (VEGT.EQ.7)  SIBVEG = 6
      IF (VEGT.EQ.8)  SIBVEG = 8
      IF (VEGT.EQ.9)  SIBVEG = 9
      IF (VEGT.EQ.10) SIBVEG = 7
      IF (VEGT.EQ.11) SIBVEG = 12
      IF (VEGT.EQ.12) SIBVEG = 11
      IF (VEGT.EQ.13) SIBVEG = 11
      IF (VEGT.GT.13) THEN
         SIBVEG = 7
      ENDIF

      VEGT = SIBVEG

      return
      end

