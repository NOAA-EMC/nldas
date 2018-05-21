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
! define_gds.F90
!
! DESCRIPTION:
! Assigns a grid definition section (GDS) array appropriate to the 
! G-LDAS resolution used.
!
!===========================================================================
! REVISION HISTORY:
!  20 Jul 2001: Urszula Jambor; Initial code
!  12 Feb 2002: Urszula Jambor; Added latmax variable assignment
!  06 Mar 2002: Urszula Jambor; Added 1 & 1/2 degree resolution GDS arrays
!  30 Jul 2002: Jon Gottschalck; Added nrgpcp to differentiate domains, used with global precip
!===========================================================================

subroutine define_gds ( ldas )

  use ldas_module      ! LDAS non-model-specific 1-D variables
  IMPLICIT NONE
  type (ldasdec)                   ldas              

!===  Begin declarations

  integer :: i, nc

!=== End declarations

!=== Initialize GDS array
  do i=1,200
     ldas%ldas_kgds(i) = 0
  end do

!=== NOTES from one of Matt's GDS arrays
!      kgds(1) = 4		!Input grid type (4=Gaussian)
!      kgds(2) = 128		!Number of points on a lat circle
!      kgds(3) = 64		!Number of points on a meridian
!      kgds(4) = -87864		!Latitude of origin x1000
!      kgds(5) = 0		!Longitude of origin x1000
!      kgds(6) = 128		!8 bits (1 byte) related to resolution
!				!(recall that 10000000 = 128), Table 7
!      kgds(7) = 87864		!Latitude of extreme point x1000
!      kgds(8) = -2812		!Longitude of extreme point x1000
!      kgds(9) = 2812		!N/S direction increment x1000
!      kgds(10) = 32		!(Gaussian) # lat circles pole-equator
!      kgds(11) = 64		!8 bit scanning mode flag (Table 8)

!=== Number of columns for this resolution
  nc = ldas%nc

  select case (nc)

     case ( 2880 )  ! 1/8 degree resolution

        ldas%ldas_kgds(1)  = 0
        ldas%ldas_kgds(2)  = ldas%nc  != 2880 
        ldas%ldas_kgds(3)  = ldas%nr  != 1200
        ldas%ldas_kgds(4)  = -59939
        ldas%ldas_kgds(5)  = -179938
        ldas%ldas_kgds(6)  = 128
        ldas%ldas_kgds(7)  = 89938
        ldas%ldas_kgds(8)  = 179938
        ldas%ldas_kgds(9)  = 125
        ldas%ldas_kgds(10) = 125
        ldas%ldas_kgds(11) = 64
        ldas%ldas_kgds(20) = 255

	ldas%latmax = 720
        ldas%nrgpcp = 960

     case ( 1440 )  ! 1/4 degree resolution

        ldas%ldas_kgds(1)  = 0
        ldas%ldas_kgds(2)  = ldas%nc  != 1440
        ldas%ldas_kgds(3)  = ldas%nr  !=  600
        ldas%ldas_kgds(4)  = -59875
        ldas%ldas_kgds(5)  = -179875
        ldas%ldas_kgds(6)  = 128
        ldas%ldas_kgds(7)  = 89875
        ldas%ldas_kgds(8)  = 179875
        ldas%ldas_kgds(9)  = 250
        ldas%ldas_kgds(10) = 250
        ldas%ldas_kgds(11) = 64
        ldas%ldas_kgds(20) = 255

	ldas%latmax = 360
        ldas%nrgpcp = 480

     case ( 720 )  ! 1/2 degree resolution

        ldas%ldas_kgds(1)  = 0
        ldas%ldas_kgds(2)  = ldas%nc  != 720
        ldas%ldas_kgds(3)  = ldas%nr  != 300
        ldas%ldas_kgds(4)  = -59750
        ldas%ldas_kgds(5)  = -179750
        ldas%ldas_kgds(6)  = 128
        ldas%ldas_kgds(7)  = 89750
        ldas%ldas_kgds(8)  = 179750
        ldas%ldas_kgds(9)  = 500
        ldas%ldas_kgds(10) = 500
        ldas%ldas_kgds(11) = 64
        ldas%ldas_kgds(20) = 255

	ldas%latmax = 180
        ldas%nrgpcp = 240

     case ( 360 )  ! 1 degree resolution

        ldas%ldas_kgds(1)  = 0
        ldas%ldas_kgds(2)  = ldas%nc  != 360
        ldas%ldas_kgds(3)  = ldas%nr  != 150
        ldas%ldas_kgds(4)  = -59500
        ldas%ldas_kgds(5)  = -179500
        ldas%ldas_kgds(6)  = 128
        ldas%ldas_kgds(7)  = 89500
        ldas%ldas_kgds(8)  = 179500
        ldas%ldas_kgds(9)  = 1000
        ldas%ldas_kgds(10) = 1000
        ldas%ldas_kgds(11) = 64
        ldas%ldas_kgds(20) = 255

	ldas%latmax = 90
        ldas%nrgpcp = 120

     case ( 144 )   ! 2 x 2.5 degree resolution

        ldas%ldas_kgds(1)  = 0
        ldas%ldas_kgds(2)  = ldas%nc  ! = 144
        ldas%ldas_kgds(3)  = ldas%nr  ! = 76
        ldas%ldas_kgds(4)  = -60000
        ldas%ldas_kgds(5)  = -180000
        ldas%ldas_kgds(6)  = 128
        ldas%ldas_kgds(7)  = 90000
        ldas%ldas_kgds(8)  = 177500
        ldas%ldas_kgds(9)  = 2500
        ldas%ldas_kgds(10) = 2000
        ldas%ldas_kgds(11) = 64
        ldas%ldas_kgds(20) = 255

	ldas%latmax = 46   !Actually, 45 + 45 + 1
        ldas%nrgpcp = 60

     case default
        
        print *, "No valid global grid defined for given resolution"
        print *, "columns: ", ldas%nc, " rows: ", ldas%nr
        print *, "Stopping..."
        stop
        
  end select

end subroutine define_gds

