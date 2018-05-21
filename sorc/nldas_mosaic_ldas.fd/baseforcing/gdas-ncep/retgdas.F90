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
! !ROUTINE: retgdas.F90
!
! !DESCRIPTION:
!  Defines forcing parameters, retrieves the fields using calls to
!  getgb, and interpolates the fields to LDAS specifications
! 
! !REVISION HISTORY:
!  14 Dec 2000: Urszula Jambor; Rewrote geteta.f to use GDAS forcing in GLDAS
!  15 Mar 2001: Jon Gottschalck; Added additional parameters and octets in 
!               which to search in GRIB files
!  01 Jun 2001: Urszula Jambor; Added option to get forcing from different 
!               files (F00 instantaneous and F06 six hour means)
!  29 Jan 2003: Urszula Jambor; Rewrote code, uses GETGB call to replace
!               ungribgdas.  Interpolation now occurs in interp_gdas.  
!               Using GETGB avoids problems with the Oct2002 GDAS 
!               grid update.
!  12 Nov 2003: Matt Rodell; Check to make sure input file exists before
!		opening and thereby creating a new, empty file.
!  14 Nov 2003: Matt Rodell; Ensure lugb varies in call to baopen
!  05 Feb 2004: James Geiger; Added GrADS-DODS Server functionality
!
! !INTERFACE:
#include "misc.h"
subroutine retgdas( order, name, nameF06, F00flag,ferror,try )
! !USES:  
  use lisdrv_module, only : lis, gindex
  use time_manager
  use baseforcing_module, only: glbdata1, glbdata2
  use gdasdomain_module, only : gdasdrv
  use spmdMod
#if ( defined OPENDAP )
  use opendap_module, only : ciam
  use gdasopendap_module
#endif
  use lis_indices_module, only : lis_nc_working, lis_nr_working
  implicit none
! !ARGUMENTS:
  integer, intent(in)           :: order    
  character(len=80), intent(in) :: name, nameF06
  integer, intent(in)           :: F00flag
  integer, intent(out)          :: ferror 
  integer, intent(inout)        :: try
!EOP
!==== Local Variables=======================
  
  character(len=80) :: fname
  integer :: iv, c, r
  integer :: errorflag
  integer :: endloop, nforce
  integer :: j, lugb,iret,gbret,jret,jpds(200),jgds(200)
  integer :: lubi,kf,k,kpds(200),gridDesc(200)
  integer :: ngdas
  integer, dimension(gdasdrv%nmif) :: pds5, pds7
  logical*1, allocatable :: lb(:)
  logical :: file_exists
  real, allocatable :: f(:)
  real, dimension(lis_nc_working, lis_nr_working) :: varfield
  character(len=22) :: validtime
  integer :: nstep

#if ( defined OPENDAP )
  real :: dummy
  real, dimension(gdasdrv%ncold,gdasdrv%nrold) :: f2d
  integer, dimension(gdasdrv%nmif) :: forcing_index
#endif

!=== End Variable Definition =======================
!BOC
  varfield = 0.0
  ngdas = (gdasdrv%ncold*gdasdrv%nrold)
!--------------------------------------------------------------------------
! Set the GRIB parameter specifiers
!--------------------------------------------------------------------------
  if ( masterproc ) then
     nstep = get_nstep(lis%t)
  endif
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
  call MPI_BCAST(nstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,errorflag)
#endif
  if ( nstep == 0 ) then
     pds5 = (/011,051,204,205,033,034,001,059,214,084,144,144, 011,011, 065/) !parameter
     pds7 = (/002,002,000,000,010,010,000,000,000,000,010,2760,010,2760,000/) !htlev2
     nforce = gdasdrv%nmif
#if ( defined OPENDAP )
     forcing_index = (/15,13,4,3,18,21,9,8,2,1,11,12,16,17,23/)
#endif
  else
     pds5 = (/ 011,051,204,205,033,034,001,059,214,084 /) !parameter
     pds7 = (/ 002,002,000,000,010,010,000,000,000,000 /) !htlev2
     nforce = 10 ! for now
#if ( defined OPENDAP )
     forcing_index = (/15,13,4,3,18,21,9,8,2,1/)
#endif

  endif

  ferror = 1  
!--------------------------------------------------------------------------
! if there's a problem then ferror is set to zero
!--------------------------------------------------------------------------
  iv = 0
  errorflag = 0 
  endloop = 0
  allocate(lb(gdasdrv%ncold*gdasdrv%nrold))
  allocate(f(gdasdrv%ncold*gdasdrv%nrold))
  do
     if ( endloop == 1 ) exit
     iv = iv+1
     if ( (F00flag > 0)  .and. &
          (iv==3.or.iv==4.or.iv==8.or.iv==9.or.iv==10) ) then
        fname = nameF06
     else
        fname = name
     end if
#if ( defined OPENDAP )
!     print*,'DBG: retgdas -- fname = ',fname
!     inquire (file=fname, exist=file_exists)
!     if ( .not. file_exists ) then
        print*, 'MSG: retgdas -- Retrieving GDAS forcing file ', & 
                trim(fname), ' (',iam,')'
        call system("opendap_scripts/getgdas.pl "//ciam//" "// & 
                     trim(fname)//" "//                        & 
                     cgdas_slat//" "//cgdas_nlat//" "//        &
                     cgdas_wlon//" "//cgdas_elon)
!        call perror("ERR: retgdas --")
!     endif
     open(unit=10,file=fname,access='sequential',form='unformatted')
     do j = 1, forcing_index(iv)-1
        read(10) dummy
     enddo
     read(10) f2d
     close(10)
     gbret=0
     ! Now make things look like grib data
     call twotoone(f2d,f,gdasdrv%ncold,gdasdrv%nrold)
     call get_kpds(kpds,iv)
     !call get_kgds(kgds)
     call get_gridDesc(gridDesc)
     if ( iv == 1 ) then
        call set_lb(f,lb,gdasdrv%ncold,gdasdrv%nrold)
     endif
#else
     inquire (file=fname, exist=file_exists)
     if (file_exists) then      
!--------------------------------------------------------------------------
! Set up to open file and retrieve specified field 
!--------------------------------------------------------------------------
        lugb = iv + try
        j = 0
        jpds = -1
        jpds(5) = pds5(iv)
        jpds(7) = pds7(iv)

        call baopen(lugb,fname,iret)
        if(iret==0) then 
           call getgb(lugb,lubi,ngdas,j,jpds,jgds,kf,k,kpds,gridDesc,lb,f,gbret)
        else 
           gbret = 99
        endif
        call baclose(lugb,jret)
     else
        ferror = 0
        deallocate(f)
        deallocate(lb)
        return
     endif
#endif
!--------------------------------------------------------------------------
! If field successfully retrieved, interplate to LIS domain
!--------------------------------------------------------------------------
     if (gbret==0) then
        call interp_gdas(kpds,ngdas,f,lb,lis%d%gridDesc, &
             lis_nc_working,lis_nr_working,varfield)
     else
        errorflag = 1
     endif
     
     if ( errorflag == 1 ) then 
        endloop = 1
        ferror = 0
     else
        do c =1, lis_nc_working
           do r = 1, lis_nr_working
             if (gindex(c,r).ne. -1) then 
              if ( order == 1 ) then 
                 glbdata1(iv,gindex(c,r)) = varfield(c,r)
              else
                 glbdata2(iv,gindex(c,r)) = varfield(c,r)
              end if
             endif
           end do 
        end do 
     end if 

     if ( errorflag == 1 ) then
        print *, 'Could not find correct forcing parameter in file',name
     end if 
     if ( iv == nforce ) endloop = 1
  end do
  deallocate(lb)
  deallocate(f)
  return
!EOC
end subroutine retgdas


!BOP
! !ROUTINE: interp_gdas
!
! !DESCRIPTION:
!   This subroutine interpolates a given GDAS field 
!   to the LIS domain.  Special treatment for some
!   initialization fields.
!   Code based on old ungribgdas.f
!
! !INTERFACE:
subroutine interp_gdas(kpds,ngdas,f,lb,lis_gds,nc,nr, &
     varfield)
! !USES:
  use lisdrv_module, only : lis
  use bilinear_interpMod, only : w110,w120,w210,w220,n110,n120,n210,n220,& 
       rlat0,rlon0
  use conserv_interpMod, only : w113,w123,w213,w223,n113,n123,n213,n223,rlat3,rlon3
  use gdasdomain_module, only : mi
!EOP
  implicit none

!=== Begin variable declarations
  integer, intent(in) :: nc, nr, ngdas
  integer, intent(in) :: kpds(200)
  real :: lis_gds(200)
  real, intent(out)   :: f(ngdas)
  real, intent(out)   :: varfield(nc,nr)
  integer :: ip, ipopt(20),ibi,km,iret
  integer :: no, ibo
  integer :: count,i,j,v,nglis, mo

  real :: ism, udef
  real, dimension(nc*nr) :: lis1d
  logical*1 :: lb(ngdas)
  logical*1 :: lo(nc*nr)

!=== End variable declarations
!BOC
!-----------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!-----------------------------------------------------------------------
  mo = nc*nr
  if (kpds(5)==59 .or. kpds(5)==214) then
     ip=3
     ipopt(1)=-1
     ipopt(2)=-1
     km=1
     ibi=1          
  else
     ip=0
     do i=1,20
        ipopt(i)=0
     enddo
     km=1
     ibi=1
  endif
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.
!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------  
  if(lis%f%interp.eq.1) then 
     call bilinear_interp(lis_gds,ibi,lb,f,ibo,lo,lis1d,mi,mo,&
          rlat0, rlon0,w110,w120,w210,w220,n110,n120,n210,n220,iret)
  elseif(lis%f%interp.eq.2) then 
     if (kpds(5)==59 .or. kpds(5)==214)then     
        call conserv_interp(lis_gds,ibi,lb,f,ibo,lo,lis1d,mi,mo, & 
             rlat3,rlon3,w113,w123,w213,w223,n113,n123,n213,n223,iret)
     else 
        call bilinear_interp(lis_gds,ibi,lb,f,ibo,lo,lis1d,mi,mo,&
             rlat0, rlon0,w110,w120,w210,w220,n110,n120,n210,n220,iret)
     endif
  endif
!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GDAS & LDAS. For LDAS land 
! points not included in GDAS geography dataset only.
!-----------------------------------------------------------------------    
  count = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count)
     enddo
     count = count + nc
  enddo
!  if(kpds(5).eq.1) then 
!     open(232, file='forcing1.bin',form='unformatted')
!     write(232) varfield
!     close(232)
!  endif

!-----------------------------------------------------------------------    
! Save air tempertaure interpolated field for later use in
! initialization of soil temp where geography differs 
! between GDAS and LDAS
!-----------------------------------------------------------------------    
!  if (kpds(5) .eq. 11 .and. kpds(6) .eq. 105) then
!     do i = 1, nc
!        do j = 1, nr
!           geogtemp(i,j) = varfield(i,j)
!        enddo
!     enddo
!  endif
!EOC
end subroutine interp_gdas
