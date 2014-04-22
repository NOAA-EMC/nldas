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
! !ROUTINE: retnldas2.F90
!
! !DESCRIPTION:
!  Retrieves the name of a file from getncep.f and then gets the NCEP-LDAS
!  forcing data using the zterp.f subroutine.  
!
!    TIME1 = most recent past data
!    TIME2 = most recent future data 
!
!  The strategy for missing data is to backwards up to 10 days to get
!  forcing at the same time of day.
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  11 Apr 2000: Brian Cosgrove; changed code to use Forcing Mask (With inland
!               water filled in).  Deleteted unused variables.
!  27 Apr 2000: Brian Cosgrove; changed code to use the original 
!               mask again since that is the
!               mask which  NCEP has already applied to the forcing data
!               by the time NASA gets it......not possible to use the 
!               expanded NASA forcing mask
!  1  May 2000: Brian Cosgrove; changed code so that if parameter 11 (sw)
!               is not found in hourly ncep data, it will just use
!               edas-based shortwave from the hourly ncep files
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!                not a hard-wired undefined value of -999.9 and -999.0
!  18 Aug 2000: Brian Cosgrove; changed code so that FMASK and not MASK
!               is used when ungribbing.  NCEP data already has a mask applied
!               to it and so may not be able to supply forcing data to
!               all LDAS land forcing points.  In areas where LDAS
!               forcing mask states that land exists, but where NCEP forcing
!               data is non-existant, assign undefined value to forcing data.
!  22 Aug 2000: Brian Cosgrove; Altered code for US/Mexico/Canada Mask
!  05 Sep 2001: Brian Cosgrove; Removed dirnom and infile variables, changed
!               call to ungribncep to match removal.  Added code to make use
!               of precip weighting mask
!  02 Feb 2004: Sujay Kumar; Initial Specification in LIS
! 
! !INTERFACE:
subroutine retnldas2(order,name,ferror,ftype,dataflag)
! !USES:
  use lisdrv_module, only : lis, grid, gindex
  use nldas2domain_module, only : nldas2drv
  use baseforcing_module, only : glbdata1, glbdata2,precipweight
!EOP
  implicit none
  integer, intent(in)      :: order,ftype,dataflag
  character*80, intent(in) :: name
  integer, intent(out)     :: ferror

  integer :: errorflag,count
  integer :: endloop,iv,iret,lugb
  logical :: file_exists
  integer :: gbret,jret,jpds(200),jgds(200),gridDesc(200),lubi
!  real :: kf,k,kpds(200)
  integer :: kf,k,kpds(200)
  integer :: nldas2,c,r,j,col,row
  integer, dimension(10) :: pds5, pds7,pds2
  logical*1, allocatable :: lb(:)
  real, allocatable :: f(:)
  real, allocatable :: oldprecip(:)
  real, dimension(lis%d%lnc, lis%d%lnr) :: varfield
  ferror = 1
  endloop=0
  errorflag = 0
  iv = 0
  nldas2 = (nldas2drv%ncold*nldas2drv%nrold)
  pds5 = (/ 011,051,204,205,033,034,001,061,153,  157/) !parameter
  pds7 = (/ 002,002,000,000,010,010,000,000,000,46080/) !htlev2
  pds2 = (/ 84, 84, 84, 84, 84, 84, 84,  84, 84,   84/)
  allocate(lb(nldas2drv%ncold*nldas2drv%nrold))
  allocate(f(nldas2drv%ncold*nldas2drv%nrold))
  allocate(oldprecip(nldas2drv%ncold*nldas2drv%nrold))

  do
     if ( endloop == 1 ) exit
     iv = iv+1
     lugb = iv+10
     inquire (file=name, exist=file_exists)
     if (file_exists) then      
!--------------------------------------------------------------------------
! Set up to open file and retrieve specified field 
!--------------------------------------------------------------------------
        lugb = iv+10
	lubi=0
        j = 0 
        jpds = -1
        jpds(5) = pds5(iv)
        jpds(7) = pds7(iv)
        jpds(2) = pds2(iv)
        call baopenr(lugb,name,iret)
	!print *,'iret=',iret
        if(iret==0) then 
	!print *,'in loop after iret'
           call getgb(lugb,lubi,nldas2,j,jpds,jgds,kf,k,kpds, &
                gridDesc,lb,f,gbret)
        !print *,'zero setting gbret to',gbret,'from getgb call'
	if(iv.eq.8) then

	c=1
!apply precip mask to precip data, blending between EDAS and Stage2
!around border of CONUS
	do row=1,lis%d%lnr
	 do col=1,lis%d%lnc
           oldprecip(c)=f(c)
           c=c+1
         enddo
        enddo
        endif !iv.eq.8
!compute ratio of EDAS convective to EDAS total precip
	if(iv.eq.9) then
	do c=1,nldas2
         if (oldprecip(c).gt.0.0) then
          f(c)=f(c)*oldprecip(c)
	 endif
	enddo
	endif !iv.eq.9
        else
	!print *,'now in else'
        !print *,'setting gbret now to 99'
           gbret = 99
        endif
        call baclose(lugb,jret)
	!print *,'here, gbret=',gbret
     else
        ferror = 0
	!print *,'now, just set ferror to 0'
        deallocate(f)
        deallocate(lb)
        deallocate(oldprecip)

     endif
	!print *,'gbret=',gbret
     if(gbret == 0) then 
! no need to interpolate since nldas2 forcing is on nldas2 domain	
!        call interp_nldas2(kpds,nldas2,f,lb,lis%d%gridDesc,& 
!             lis%d%lnc,lis%d%lnr,varfield)
	count=1
        do r =1, lis%d%lnr
           do c = 1, lis%d%lnc
              varfield(c,r)=f(count)
              count=count+1
	   enddo
	enddo

     else
        errorflag = 1
     endif
!move data into global arrays and adjust convective precip 
!to match merged total precip product
     if(errorflag == 1) then 
        endloop = 1
        ferror = 0
	!print *,'just set ferror to 0,errorflag=',errorflag
     else
        do r =1, lis%d%lnr
           do c = 1, lis%d%lnc
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
     if(errorflag ==9) then 
        print*, 'Could not find forcing parameter in file, ',name
     endif
     if(iv == 9) endloop = 1
  enddo
  deallocate(lb)
  deallocate(f)
  deallocate(oldprecip)
  return
end subroutine retnldas2


!BOP
! !ROUTINE: interp_nldas2
!
! !DESCRIPTION:
!   This subroutine interpolates a given NLDAS2 field 
!   to the LIS domain.  Special treatment for some
!   initialization fields.
!   Code based on old ungribgdas.f
!
! !INTERFACE:
subroutine interp_nldas2(kpds, nldas2,f,lb,lis_gds,nc,nr, &
     varfield)
! !USES:
  use lisdrv_module, only : lis
  use bilinear_interpMod, only : w110,w120,w210,w220,n110,n120,n210,n220,& 
       rlat0,rlon0
  use conserv_interpMod, only : w113,w123,w213,w223,n113,n123,n213,n223,&
       rlat3,rlon3
  use nldas2domain_module, only : mi
!EOP
  implicit none

!=== Begin variable declarations
  integer, intent(in)   :: nc, nr, nldas2,kpds(200)
  real, intent(in) :: lis_gds(200)
  logical*1, intent(in) :: lb(nldas2)
  real, intent(out)     :: f(nldas2),varfield(nc,nr)

  integer :: ip, ipopt(20),ibi,km,iret
  integer :: no, ibo,mo
  integer :: count,i,j,v


  real :: ism, udef
  real, dimension(nc*nr) :: lis1d
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
  if (kpds(5)==61 .or. kpds(5)==214) then
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
     if (kpds(5)==61 .or. kpds(5)==214) then     
        call conserv_interp(lis_gds,ibi,lb,f,ibo,lo,lis1d,mi,mo,& 
             rlat3,rlon3,w113,w123,w213,w223,n113,n123,n213,n223,iret)
     else 
        call bilinear_interp(lis_gds,ibi,lb,f,ibo,lo,lis1d,mi,mo,&
             rlat0, rlon0,w110,w120,w210,w220,n110,n120,n210,n220,iret)
     endif
  endif
!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between LDAS & LDAS. For LDAS land 
! points not included in LDAS geography dataset only.
!-----------------------------------------------------------------------    
  count = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count)
     enddo
     count = count + nc
  enddo
!-----------------------------------------------------------------------    
! Save air tempertaure interpolated field for later use in
! initialization of soil temp where geography differs 
! between LDAS and LDAS
!-----------------------------------------------------------------------    
!  if (kpds(5) .eq. 11 .and. kpds(6) .eq. 105) then
!     do i = 1, nc
!        do j = 1, nr
!           geogtemp(i,j) = varfield(i,j)
!        enddo
!     enddo
!  endif
!EOC
end subroutine interp_nldas2
