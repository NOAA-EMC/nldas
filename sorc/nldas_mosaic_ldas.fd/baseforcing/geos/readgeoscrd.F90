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
! !ROUTINE: readgeoscrd.F90
!
! !DESCRIPTION:
!  Routine to read GEOS specific parameters from the card file. 
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readgeoscrd(geosdrv,gridDesci)
! !USES:
  use geosdrv_module
#if ( defined OPENDAP )
    use geosopendap_module, only : opendap_geos_init, &
                                   def_gridDesc
#endif
!EOP
  implicit none
  integer :: lsm
  type(geosdrvdec) :: geosdrv
  namelist /geos/geosdrv
  real, intent(inout) :: gridDesci(50)
!BOC
  open(11,file='lis.crd',form='formatted',status='old')
  read(unit=11,NML=geos)
  print*,'Using GEOS forcing'
  print*, 'GEOS forcing directory :',geosdrv%GEOSDIR
  geosdrv%geostime1 = 3000.0
  geosdrv%geostime2 = 0.0
  geosdrv%gridchange = .true.
#if ( defined OPENDAP )
    call opendap_geos_init()
!    call def_gridDesc(gridDesci)
#endif
#if (! defined OPENDAP) || ( defined OPENDAP )
    gridDesci(1) = 0
    gridDesci(2) = geosdrv%ncold
    gridDesci(3) = geosdrv%nrold
    gridDesci(4) = -90.000
    gridDesci(5) = -180.000
    gridDesci(7) = 90.000
    gridDesci(8) = 179.000
    gridDesci(6) = 128
    gridDesci(9) = 1.000
    gridDesci(10) = 1.000
    gridDesci(20) = 255
#endif
  close(11)
!EOC
end subroutine readgeoscrd
