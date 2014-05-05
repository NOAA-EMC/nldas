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
! !ROUTINE: interp_cmap.F90
!
! !DESCRIPTION:
!
!  Interpolates CMAP observed precipitation forcing
!
! !INTERFACE: 
subroutine interp_cmap(kpds,ngdas,f,lb,lis_gds,nc,nr, &
     varfield)
! !USES:
  use cmapdomain_module, only : mi,w11,w12,w21,w22,&
       n11,n12,n21,n22,rlat,rlon

  implicit none
! !ARGUMENTS:
  integer :: nc, nr, ngdas
  integer :: kpds(200)
  real :: lis_gds(50)
  real :: f(ngdas)
  logical*1 :: lb(ngdas)
  real, dimension(nc,nr) :: varfield
!EOP

  integer :: ip, ipopt(20),ibi,km,iret
  integer :: no, ibo, mo
  integer :: count,i,j,v

  real :: ism, udef
  real, dimension(nc,nr) :: geogtemp
  real, dimension(nc*nr) :: lis1d

  logical*1 :: geogmask(nc,nr)
  logical*1 :: lo(nc*nr)

!=== End variable declarations
!BOC
!--------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!--------------------------------------------------------------------
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
!--------------------------------------------------------------------  
! Initialize output bitmap. Important for soil moisture and temp.
!--------------------------------------------------------------------  
  lo = .true.
  

!  call ipolates (ip,ipopt,gridDesc,lis_gds,ngdas,nglis, &
!       km,ibi,lb,f,no,rlat,rlon,ibo,lo,lis1d,iret)
  mi = ngdas
!  call polates0 (lis_gds,ibi,lb,f,ibo,lo,lis1d,mi,&
!       rlat, rlon,w11,w12,w21,w22,n11,n12,n21,n22,iret)
  call conserv_interp(lis_gds,ibi,lb,f,ibo,lo,lis1d,mi,mo,&
       rlat, rlon,w11,w12,w21,w22,n11,n12,n21,n22,iret)
!--------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GDAS & LDAS. For LDAS land 
! points not included in GDAS geography dataset only.
!--------------------------------------------------------------------    
  count = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count)
        geogmask(i,j) = lo(i+count)
     enddo
     count = count + nc
  enddo
!--------------------------------------------------------------------    
! Save air tempertaure interpolated field for later use in
! initialization of soil temp where geography differs 
! between GDAS and LDAS
!--------------------------------------------------------------------    
  if (kpds(5) .eq. 11 .and. kpds(6) .eq. 105) then
     do i = 1, nc
        do j = 1, nr
           geogtemp(i,j) = varfield(i,j)
        enddo
     enddo
  endif
!EOC  
end subroutine interp_cmap

