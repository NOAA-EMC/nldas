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
#include "misc.h"
!BOP
!
! !ROUTINE: readgeos.F90
! 
! !DESCRIPTION: 
!
!  Reads in GEOS data and performs interpolation to the LDAS domain.
!
! \subsection{Core Functions of readgeos}
!  \begin{description}
!  \item[bilinear_interp]
!      Interpolates GEOS data to LDAS grid using bilinear interpolation
!  \end{description}
!
!
! GEOS FORCING VARIABLES (unless noted, fields are 3-hr upstream averaged): \\
!  1. T 2m    Temperature interpolated to 2 metres [$K$] \\
!  2. q 2m    Instantaneous specific humidity interpolated to 2 metres[$kg/kg$] \\
!  3. radswg  Downward shortwave flux at the ground [$W/m^2$] \\
!  4. lwgdwn  Downward longwave radiation at the ground [$W/m^2$]\\
!  5. u 10m   Instantaneous zonal wind interpolated to 10 metres [$m/s$] \\
!  6. v 10m   Instantaneous meridional wind interpolated to 10 metres[$m/s$]\\
!  7. ps      Instantaneous Surface Pressure [$Pa$] \\
!  8. preacc  Total precipitation [$mm/s$] \\
!  9. precon  Convective precipatation [$mm/s$] \\
! 10. albedo  Surface albedo (0-1) \\
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  11 Apr 2000: Brian Cosgrove; Added read statements for forcing interpolation
!  17 Apr 2001: Jon Gottschalck; Added code to perform initialization of
!                                Mosaic with GEOS forcing and new intp. scheme
!  14 Aug 2001: Urszula Jambor; Added ferror flag as a routine argument
!  07 Dec 2001: Urszula Jambor; Began used of LDAS$%$LDAS_GRIDDESC array
!
! !INTERFACE:      
subroutine readgeos(order,name,tscount,ferror)
! !USES:
  use lisdrv_module, only : lis
  use spmdMod
  use baseforcing_module, only : glbdata1, glbdata2 
  use bilinear_interpMod, only : w110,w120,w210,w220, &
                                 n110,n120,n210,n220, & 
                                 rlat0,rlon0
  use conserv_interpMod, only : w113,w123,w213,w223, &
                                n113,n123,n213,n223, & 
                                rlat3,rlon3           
  use geosdomain_module, only : geosdrv,mi 

#if ( defined OPENDAP )
  use geosopendap_module, only : geos_nr,geos_nc
#endif
  use lisdrv_module, only : gindex ! LDAS non-model-specific 1-D variables
  use lis_indices_module, only : lis_nc_working, lis_nr_working
  use lis_openfileMod

!EOP
  implicit none
!=== Local Variables =====================================================
  character*80, intent(in) :: name
  integer, intent(in)      :: order, tscount
  integer, intent(out)     :: ferror          
  integer :: mo
  integer :: nforce,i,j,v
  integer :: ios,ioerror      ! set to non-zero if there's an error
  real, allocatable :: tempgeos(:,:) 
                                  ! lf%nmif = Max.# parameters to retrieve
  real, allocatable :: tempvar(:,:,:)
  integer :: glis,ngeos                ! Size of I/O 1D fields
  real, allocatable :: f(:),go(:) ! 1D I/O fields
  real, allocatable :: gmask(:)  ! GEOS forcing mask
  integer :: ibi,ibo,count
  integer :: iret,c
  real :: gridDesco(50)       ! Input,output grid info arrays
  logical*1, allocatable :: lb(:)
  logical*1, allocatable :: lo(:)      ! Input and output bitmaps
!BOC

  integer :: nr_index, nc_index
  character :: corder
#if ( defined OPENDAP )
  integer :: geos_index
  nr_index = geos_nr
  nc_index = geos_nc
#else
  nr_index = geosdrv%nrold
  nc_index = geosdrv%ncold
#endif
  write(corder, '(i1)') order
  
  allocate(tempvar(nc_index,nr_index,geosdrv%nmif), stat=ios)
  call check_error(ios,'Error allocating tempvar.',iam)
  
  allocate(tempgeos(lis%d%ngrid,geosdrv%nmif), stat=ios)
  call check_error(ios,'Error allocating tempgeos.',iam)
  
  allocate(f(nc_index*nr_index), stat=ios)
  call check_error(ios,'Error allocating f.',iam)
  
  allocate(go(lis_nc_working*lis_nr_working), stat=ios)
  call check_error(ios,'Error allocating go.',iam)
  
  allocate(gmask(nc_index*nr_index), stat=ios)
  call check_error(ios,'Error allocating gmask.',iam)
  
  allocate(lb(nc_index*nr_index), stat=ios)
  call check_error(ios,'Error allocating lb.',iam)
  
  allocate(lo(lis_nc_working*lis_nr_working), stat=ios)
  call check_error(ios,'Error allocating lo.',iam)
  
  gmask = 0.0
  ngeos = nc_index*nr_index
  glis = lis_nc_working*lis_nr_working
  ferror = 1                
!------------------------------------------------------------------------
! Open GEOS forcing file
!------------------------------------------------------------------------
#if ( defined OPENDAP )
  call lis_open_file(40,file=name,status='old',form='unformatted',       & 
                     recl=geosdrv%ncold*geosdrv%nrold*4,access='direct', &
                     script='getgeos.pl',time_offset=corder)
#else
  call lis_open_file(40,file=name,form='unformatted')
!  call lis_open_file(40,file=name,status='old',form='unformatted',       &
!                     recl=geosdrv%ncold*geosdrv%nrold*4,access='direct', &
!                     script='none',time_offset=corder)
#endif


  call lis_read_geos_forcing(40,name,tempvar,nc_index,nr_index,geosdrv%nmif)

!------------------------------------------------------------------------
! Finding number of forcing variables 
! (13 if time step is 0, otherwise the normal 10)
!------------------------------------------------------------------------
  if (tscount .eq. 0) then
     nforce = geosdrv%nmif
  else
     nforce = 10 
  endif
  do v=1,nforce
!------------------------------------------------------------------------
! Transferring current data to 1-D array for interpolation
!------------------------------------------------------------------------
     c=0
     do i=1,nr_index    !i=1,lf%nrold
        do j=1,nc_index !j=1,lf%ncold
           c = c + 1
           f(c) = tempvar(j,i,v)
           if (tscount .eq. 0 .and. order .eq. 1 & 
                .and. v .eq. 11) then
              gmask(c) = f(c) ! Storing geos land mask for later use
           endif
        enddo
     enddo
!------------------------------------------------------------------------     
! Initializing input and output grid arrays
!------------------------------------------------------------------------
     gridDesco = 0
     gridDesco = lis%d%gridDesc
     mo = lis_nc_working * lis_nr_working
     if (v .eq. 8 .or. v .eq. 9) then
        ibi = 1
     else                 
        ibi = 1
     endif
!------------------------------------------------------------------------
! Defining input data bitmap
!------------------------------------------------------------------------
     do i=1,ngeos
        lb(i)=.true.
     enddo
!------------------------------------------------------------------------     
! Alter default bitmap prescribed above for 
! surface parameters (soil wetness, snow)
!------------------------------------------------------------------------
     if (v .eq. 12 .or. v .eq. 13) then
        do i=1,ngeos
           if (gmask(i)==100.0 .or. gmask(i)==101.0) then
              lb(i)=.false.
           else
              lb(i)=.true.
           endif
        enddo
     endif
!------------------------------------------------------------------------
! Defining output data bitmap
!------------------------------------------------------------------------
     do i=1,glis
        lo(i)=.true.
     enddo
!------------------------------------------------------------------------
! Interpolate data from GEOS grid to GLDAS grid
!------------------------------------------------------------------------
     if(lis%f%interp.eq.1) then 
        call bilinear_interp(gridDesco,ibi,lb,f,ibo,lo,go,mi,mo, & 
             rlat0,rlon0,w110,w120,w210,w220,n110,n120,n210,n220,iret)
     elseif(lis%f%interp.eq.2) then 
        if(v.eq.8 .or. v.eq. 9) then            
           call conserv_interp(gridDesco,ibi,lb,f,ibo,lo,go,mi,mo,& 
                rlat3,rlon3,w113,w123,w213,w223,n113,n123,n213,n223,iret)
        else 
           call bilinear_interp(gridDesco,ibi,lb,f,ibo,lo,go,mi,mo, & 
                rlat0,rlon0,w110,w120,w210,w220,n110,n120,n210,n220,iret)
        endif
     endif
!------------------------------------------------------------------------
! Convert data to original 3D array & a 2D array to 
! fill in of missing points due to geography difference  
!------------------------------------------------------------------------
     count = 0
     do j = 1, lis_nr_working
        do i = 1, lis_nc_working
           if(gindex(i,j) .ne. -1) then
              tempgeos(gindex(i,j),v) = go(i+count)
           endif
        enddo
        count = count + lis%d%lnc
     enddo
!------------------------------------------------------------------------
! Fill in undefined and ocean points
!------------------------------------------------------------------------
     do i = 1, lis%d%ngrid
        if (tempgeos(i,v) >= 9.9e+14) then
           tempgeos(i,v) = lis%d%udef
        endif
        if(order.eq.1)then
           glbdata1(v,i)=tempgeos(i,v)
        else
           glbdata2(v,i)=tempgeos(i,v)
        endif
     enddo               !i
  enddo                  !v
  
  print*,'DBG: readgeos -- Deallocating arrays'
  
  deallocate(tempvar, stat=ios)
  call check_error(ios,'Error deallocating tempvar.',iam)
  
  deallocate(tempgeos, stat=ios)
  call check_error(ios,'Error deallocating tempgeos.',iam)
  
  deallocate(f, stat=ios)
  call check_error(ios,'Error deallocating f.',iam)
  
  deallocate(go, stat=ios)
  call check_error(ios,'Error deallocating go.',iam)
  
  deallocate(gmask, stat=ios)
  call check_error(ios,'Error deallocating gmask.',iam)
  
  deallocate(lb, stat=ios)
  call check_error(ios,'Error deallocating lb.',iam)
  
  deallocate(lo, stat=ios)
  call check_error(ios,'Error deallocating lo.',iam)

  print*, 'MSG: readgeos -- Closing GEOS forcing file -', & 
       trim(name), ' (',iam,')'
  close(40, iostat=ios)
  if ( ios /= 0 ) then
     print*,'ERR: readgeos -- Error closing ', trim(name), & 
          '. Stopping.', ' (', iam, ')'
     call endrun
  endif
  
  print*,'DBG: readgeos -- leaving', ' (', iam, ')'
  
  return
!EOC

end subroutine readgeos

subroutine lis_read_geos_forcing(unit, name, tempvar,nc, nr, nf)

   use geosdomain_module, only : geosdrv
   use spmdMod,           only : iam

   implicit none

   integer, intent(in)                    :: unit, nc, nr, nf
   character(len=*), intent(in)           :: name
   real, dimension(nc,nr,nf), intent(out) :: tempvar

   integer :: i, ioerror

   call lis_log_msg('MSG: readgeos -- Reading GEOS forcing file: '//trim(name)) 

#if ( defined OPENDAP )
   do i = 1,geosdrv%nmif
!      read(unit=unit,iostat=ioerror)tempvar(:,:,i)         ! sequential access
      read(unit=unit,rec=i,iostat=ioerror)tempvar(:,:,i)  ! direct access -- normal
   enddo
#else
  read(unit=unit,iostat=ioerror)tempvar  ! sequential access -- normal
!   do i = 1,geosdrv%nmif                !direct access -- for debugging
!      read(unit=unit,rec=i,iostat=ioerror)tempvar(:,:,i)
!   enddo
#endif

   call check_error(ioerror,'readgeos -- Error reading '//trim(name),iam)
   call lis_log_msg('MSG: readgeos -- Read GEOS forcing file: '//trim(name)) 

end subroutine lis_read_geos_forcing
