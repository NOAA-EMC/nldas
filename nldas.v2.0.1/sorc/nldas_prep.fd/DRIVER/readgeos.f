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
! readgeos.f:
!
! DESCRIPTION:
!  Reads in GEOS data and performs interpolation to the LDAS domain.
!
! GEOS FORCING VARIABLES (unless noted, fields are 3-hr upstream averaged):
!  1. T 2m      Temperature interpolated to 2 metres [K]
!  2. q 2m      Instantaneous specific humidity interpolated to 2 metres[kg/kg]
!  3. radswg    Downward shortwave flux at the ground [W/m^2]
!  4. lwgdwn    Downward longwave radiation at the ground [W/m^2]
!  5. u 10m     Instantaneous zonal wind interpolated to 10 metres [m/s]
!  6. v 10m     Instantaneous meridional wind interpolated to 10 metres[m/s]
!  7. ps        Instantaneous Surface Pressure [Pa]
!  8. preacc    Total precipitation [mm/s]
!  9. precon    Convective precipatation [mm/s]
! 10. albedo    Surface albedo (0-1)
! 11. sfctyp    Surface types
! 12. snow      Snow depth (mm water equivalent)
! 13. gwet(top) Top soil layer wetness (fraction)
! 14. gwetroot  Root soil layer wetness (fraction)
! 14/15. T 10m  Temperature interpolated to 10 metres [K]
! 15/16. q 10m  Instantaneous specific humidity interpolated to 10 metres[kg/kg]
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  11 Apr 2000: Brian Cosgrove; Added read statements for forcing interpolation
!  17 Apr 2001: Jon Gottschalck; Added code to perform initialization of
!                                Mosaic with GEOS forcing and new intp. scheme
!  14 Aug 2001: Urszula Jambor; Added ferror flag as a routine argument
!  07 Dec 2001: Urszula Jambor; Began used of LDAS%LDAS_KGDS array
!  30 Oct 2002: Urszula Jambor; Modified to accomodate GEOS4 format
!               Switched order of snow & gwet (12 & 13), to allow for
!               GEOS4 gwettop & gwetroot.
!  25 Feb 2003: Urszula Jambor; Added file status check when opening
!               GEOS file and error message.
!  25 Mar 2003: Urszula Jambor; Modified argument list passed to GEOGFILL.
!=========================================================================

      subroutine readgeos(order,name,ldas,grid,ferror)

      use ldas_module      ! LDAS non-model-specific 1-D variables
      use grid_module      ! LDAS non-model-specific grid variables
      implicit none
      type (ldasdec) ldas
      type (griddec) grid(ldas%nc,ldas%nr)

!=== Local Variables =====================================================
      character*80 name
      integer nforce,i,j,v,order,ii
      integer :: ferror           ! set to zero if there's an error
      integer :: ios,ioerror      ! set to non-zero if there's an error
      real :: tempgeos(ldas%nc,ldas%nr,ldas%nmif) 
                                  ! ldas%nmif = Max.# parameters to retrieve
      real :: tempvar(ldas%ncold,ldas%nrold,ldas%nmif)
      integer :: gldas,ngeos                ! Size of I/O 1D fields
      real :: f(ldas%ncold*ldas%nrold),go(ldas%nc*ldas%nr) ! 1D I/O fields
      real :: gmask(ldas%ncold*ldas%nrold)  ! GEOS forcing mask
      real :: tgeos(ldas%nc,ldas%nr)        ! Interpolated 2D data field
      real :: rlat(ldas%nc*ldas%nr)
      real :: rlon(ldas%nc*ldas%nr)         ! Output lat & lon
      integer :: ibi,no,ibo,ipopt(20),count
      integer :: iret,km,ip,c
      integer :: kgds(200),kgdso(200)       ! Input,output grid info arrays
      logical*1 :: lb(ldas%ncold*ldas%nrold)
      logical*1 :: lo(ldas%nc*ldas%nr)      ! Input and output bitmaps
      logical*1 :: geogmask(ldas%nc,ldas%nr)! 2D output bitmap

      ngeos = ldas%ncold*ldas%nrold
      gldas = ldas%nc*ldas%nr

!=== End Defining Local Variables =========================================

      tempvar = 0.0
      ferror = 1                ! if a problem, ferror is set to zero
!=== Open GEOS forcing file
      print*,'open ', trim(name)
      open(40,file=name,form='unformatted',status="old",iostat=ios)
      if (ios /= 0) then
         print*, "Problem opening GEOS file"
         ferror = 0
      else
         read(40,iostat=ioerror)tempvar
         if (ioerror /= 0) then
            ferror = 0
            close(40)
         else
            close(40)

            !=== Finding number of forcing variables 
            !=== (13 if time step is 0, otherwise the normal 10)
            if (ldas%tscount .eq. 0) then
               nforce = ldas%nmif
            else
               nforce = ldas%nf
            endif

            do v=1,nforce

               !=== Transferring current data to 1-D array for interpolation
               c=0
               do i=1,ldas%nrold
                  do j=1,ldas%ncold
                     c = c + 1
                     f(c) = tempvar(j,i,v)
                     if (ldas%tscount .eq. 0 .and. order .eq. 1
     &                    .and. v .eq. 11) then
                        gmask(c) = f(c) ! Storing geos land mask for later use
                     endif            
                  enddo
               enddo

               !=== Initializing input and output grid arrays
               kgds  = 0
               kgdso = 0

               !=== Set input & output grid array values (GEOS to GLDAS)
               if (ldas%ncold==360) then !=== using GEOS3 data
                  kgds(1)  = 0
                  kgds(2)  = ldas%ncold
                  kgds(3)  = ldas%nrold
                  kgds(4)  = -90000
                  kgds(5)  = -180000
                  kgds(6)  = 128
                  kgds(7)  = 90000
                  kgds(8)  = 179000
                  kgds(9)  = 1000
                  kgds(10) = 1000
                  kgds(20) = 255
               else                      !=== using GEOS4 data
                  kgds(1)  = 0
                  kgds(2)  = ldas%ncold
                  kgds(3)  = ldas%nrold
                  kgds(4)  = -90000
                  kgds(5)  = -180000
                  kgds(6)  = 128
                  kgds(7)  = 90000
                  kgds(8)  = 178750
                  kgds(9)  = 1250
                  kgds(10) = 1000
                  kgds(20) = 255
               endif

               kgdso = ldas%ldas_kgds

               !=== Setting interpolation options 
               !=== (ip=0, bilinear),(iopt=0, no options), 
               !=== (km=1, one variable),(ibi=1, use bitmap)
               !=== Adjust to budget-bilinear for precip forcing fields
               if (v .eq. 8 .or. v .eq. 9) then
                 ip = 3
                 ipopt(1) = -1
                 ipopt(2) = -1
                 km = 1
                 ibi = 1
               else                 
                ip = 0
                do ii=1,20
                  ipopt(ii) = 0
                enddo
                km = 1
                ibi = 1
               endif

               !=== Defining input data bitmap
               do i=1,ngeos
                  lb(i)=.true.
               enddo

               !=== Alter default bitmap prescribed above for 
               !=== surface parameters (snow, soil wetness)
               if (v .eq. 12 .or. v .eq. 13) then
                  do i=1,ngeos

                     if (ldas%ncold==360) then !=== using GEOS3 data
                        if (gmask(i)==100.0 .or. gmask(i)==101.0) then
                           lb(i)=.false.
                        else
                           lb(i)=.true.
                        endif
                     else                      !=== using GEOS4 data
                        if (gmask(i)==0.0 .or. gmask(i)==2.0) then
                           lb(i)=.false.
                        else
                           lb(i)=.true.
                        endif
                     endif

                 enddo
               endif

               !=== Defining output data bitmap
               do i=1,gldas
                  lo(i)=.true.
               enddo
        
               !=== Interpolate data from GEOS grid to GLDAS grid
               CALL IPOLATES(ip,ipopt,kgds,kgdso,ngeos,gldas,km,ibi,lb,
     &              f,no,rlat,rlon,ibo,lo,go,iret)

               !=== Convert data to original 3D array & a 2D array to 
               !=== fill in of missing points due to geography difference  
               count = 0
               do j = 1, ldas%nr
                  do i = 1, ldas%nc
                     tempgeos(i,j,v) = go(i+count)
                     if (v .eq. 12 .or. v .eq. 13) then !snow or soil wetness
                        geogmask(i,j) = lo(i+count)
                        tgeos(i,j) = go(i+count)
                     endif
                  enddo
                  count = count + ldas%nc
               enddo

               !=== Fill in of missing points due to geography difference
               if (v .eq. 12 .or. v .eq. 13) then !snow or soil wetness only
                  CALL GEOGFILL(ldas%nc,ldas%nr,grid%fimask,
     &                 geogmask,tgeos,v)
                  do j = 1, ldas%nr
                     do i = 1, ldas%nc
                        if (tgeos(i,j) .eq. -1.0) then
                           if (v == 12) tgeos(i,j)=0.0 ! Snow
                           if (v == 13) tgeos(i,j)=0.5 ! mid-pt %-saturation
                        endif
                        tempgeos(i,j,v) = tgeos(i,j)
                     enddo
                  enddo
               endif

               !=== Fill in undefined and ocean points
               DO j = 1, ldas%nr
                  DO i = 1, ldas%nc
                     if ((grid(i,j)%mask < 1.0) .or. 
     &                    (tempgeos(i,j,v) == 9.9999999e+14)) then
                        tempgeos(i,j,v) = ldas%udef
                     endif

                     if(order.eq.1)then
                        grid(i,j)%glbdata1(v)=tempgeos(i,j,v)
                     else
                        grid(i,j)%glbdata2(v)=tempgeos(i,j,v)
                     endif

                  enddo !c
               enddo !r
            enddo !v

         end if
      end if
      return
      end

