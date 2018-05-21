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
! retgdas.F90: 
!
! DESCRIPTION:
!  Defines forcing parameters, retrieves the fields using calls to
!  getgb, and interpolates the fields to LDAS specifications (interp_gdas)
!
! REVISION HISTORY:
!  14 Dec 2000: Urszula Jambor; Rewrote geteta.f to use GDAS forcing in GLDAS
!  15 Mar 2001: Jon Gottschalck; Added additional parameters and octets in 
!               which to search in GRIB files
!  01 Jun 2001: Urszula Jambor; Added option to get forcing from different 
!               files (F00 instantaneous and F06 six hour means)
!  29 Jan 2003: Urszula Jambor; Rewrote code, uses GETGB call to replace
!               ungribgdas.  Interpolation now occurs in interp_gdas.  
!               Using GETGB avoids problems with the Oct2002 GDAS 
!               grid update.
!=========================================================================

subroutine retgdas( order, ldas, grid, name, nameF06, F00flag, ferror )
	  
  use ldas_module      ! LDAS non-model-specific 1-D variables
  use grid_module      ! LDAS non-model-specific grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc, ldas%nr)


!==== Local Variables=======================

  integer :: order    ! 1 indicates lesser interp. bdry, 2 indicates greater
  integer :: F00flag  ! if 1, need for data from 2 files (name, nameF06)
  integer :: ferror   ! set to zero if there's an error
  integer :: iv, c, r
  integer :: errorflag, try
  integer :: endloop, nforce

  integer :: j, lugb,iret,gbret,jret,jpds(200),jgds(200)
  integer :: lubi,kf,k,kpds(200),kgds(200)
  integer :: ngdas
  integer, dimension(ldas%nmif) :: pds5, pds7

  logical*1 :: lb(ldas%ncold*ldas%nrold)

  real, dimension(ldas%ncold*ldas%nrold) :: f
  real, dimension(ldas%nc, ldas%nr) :: varfield
  real :: ism

  character(len=80) :: name, nameF06, fname
  character(len=22) :: validtime

!=== End Variable Definition =======================

  try = ferror
  varfield = 0.0
  ngdas = (ldas%ncold*ldas%nrold)

!=== Set the GRIB parameter specifiers

  if (ldas%tscount .eq. 0) then
     pds5 = (/011,051,204,205,033,034,001,059,214,084,144,144, 011,011, 065/) !parameter
     pds7 = (/002,002,000,000,010,010,000,000,000,000,010,2760,010,2760,000/) !htlev2
     nforce = ldas%nmif
  else
     pds5 = (/ 011,051,204,205,033,034,001,059,214,084 /) !parameter
     pds7 = (/ 002,002,000,000,010,010,000,000,000,000 /) !htlev2
     nforce = ldas%nf
  endif

  ferror = 1  !if there's a problem then ferror is set to zero
  iv = 0
  errorflag = 0  !set to >0 when problem with ungribbing a variable
  endloop = 0

  do
     if ( endloop == 1 ) exit
     iv = iv+1

     if ( (F00flag > 0)  .and. &
          (iv==3.or.iv==4.or.iv==8.or.iv==9.or.iv==10) ) then
        !=== parameter averages 
        fname = nameF06
     else
        fname = name
     end if
      
     !=== Set up to open file and retrieve specified field 
     lugb = iv+try  ! must alternate file unit # 
     j = 0
     jpds = -1
     jpds(5) = pds5(iv)
     jpds(7) = pds7(iv)
     call baopen(lugb,fname,iret)
     if (iret==0) then
        call getgb(lugb,lubi,ngdas,j,jpds,jgds,kf,k,kpds,kgds,lb,f,gbret)
     else
        gbret = 99
     endif
     call baclose(lugb,jret)

     !=== If field successfully retrieved, interplate to LDAS domain
     if (gbret==0) then
        ism = 0.3
        select case (ldas%lsm) ! specify default initial soil moisture
           case(1)
              ism = ldas%mos_ism
           case(2)
              ism = ldas%clm1_ism
           case(3)
              ism = ldas%noah_ism
           case(5)
              ism = ldas%clm2_ism
        end select
        call interp_gdas(kpds,kgds,ngdas,f,lb,ldas%ldas_kgds, &
             ldas%nc,ldas%nr,varfield,grid%fimask,ism,ldas%udef)
     else
        ! problem with retrieving field
        errorflag = 1
     endif
     
     if ( errorflag == 1 ) then 
        endloop = 1
        ferror = 0
     else
        do c =1, ldas%nc
           do r = 1, ldas%nr
              if ( order == 1 ) then 
                 grid(c,r)%glbdata1(iv) = varfield(c,r)
              else
                 grid(c,r)%glbdata2(iv) = varfield(c,r)
              end if
           end do !r
        end do !c
     end if !errorflag equals 1 or zero

     if ( errorflag == 1 ) then
        print *, 'Could not find correct forcing parameter in file',name
     else
        do c =1, ldas%nc
           do r = 1, ldas%nr
              if ( (varfield(c,r) == ldas%udef)  &
                   .and. (grid(c,r)%fmask > 1.0)  ) then
                 print *, 'retgdas error: ',c,r,' IV=',iv,varfield(c,r)
              end if !parameter is undefined and specified as land
           end do !r
        end do !c
     end if !errorflag equals 1 or zero

     if ( iv == nforce ) endloop = 1
  end do !loop while endloop is not 1

  return
end subroutine retgdas


!!!!!SSSSS  SUBROUTINES    SUBROUTINES    SUBROUTINES   SSSSS
!
!
!=======================================================
!  interp_gdas
!
!  DESCRIPTION:
!   This subroutine interpolates a given GDAS field 
!   to the LDAS domain.  Special treatment for some
!   initialization fields.
!   Code based on old ungribgdas.f
!
!=======================================================

subroutine interp_gdas(kpds,kgds,ngdas,f,lb,ldas_gds,nc,nr, &
     varfield,fimask,ism,udef)

  implicit none

!=== Begin variable declarations
  integer :: nc, nr, ngdas, ngldas
  integer :: kpds(200),kgds(200), ldas_gds(200)
  integer :: ip, ipopt(20),ibi,km,iret
  integer :: no, ibo
  integer :: count,i,j,v
  integer :: fimask(nc,nr)

  real :: f(ngdas)
  real :: ism, udef
  real, dimension(nc,nr) :: varfield, geogtemp
  real, dimension(nc*nr) :: ldas1d,rlat,rlon

  logical*1 :: geogmask(nc,nr)
  logical*1 :: lb(ngdas)
  logical*1 :: lo(nc*nr)

!=== End variable declarations

  !===> Setting interpolation options (ip=0,bilinear)
  !===> (km=1, one parameter, ibi=1,use undefined bitmap
  !===> (needed for soil moisture and temperature only)
  !===> Use budget bilinear (ip=3) for precip forcing fields
  ngldas = nc*nr
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
  
  !===> Initialize output bitmap. Important for soil moisture and temp.
  
  lo = .true.
  
  !===> Interpolate to LDAS grid
  
  call ipolates (ip,ipopt,kgds,ldas_gds,ngdas,ngldas, &
       km,ibi,lb,f,no,rlat,rlon,ibo,lo,ldas1d,iret)
  
  !===> Create 2D array for main program. Also define a "soil" mask
  !===> due to different geography between GDAS & LDAS. For LDAS land 
  !===> points not included in GDAS geography dataset only.
  
  count = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = ldas1d(i+count)
        geogmask(i,j) = lo(i+count)
     enddo
     count = count + nc
  enddo
  
  !===> Save air tempertaure interpolated field for later use in
  !===> initialization of soil temp where geography differs 
  !===> between GDAS and LDAS
  
  if (kpds(5) .eq. 11 .and. kpds(6) .eq. 105) then
     do i = 1, nc
        do j = 1, nr
           geogtemp(i,j) = varfield(i,j)
        enddo
     enddo
  endif
  
  !===> Due to difference in land/ocean datasets between GDAS and LDAS,
  !===> fill in points along coasts, small islands etc. Performed only 
  !===> for soil moisture, temp, and snow        
  
  if ((kpds(5) .eq. 11 .and. kpds(6) .eq. 112) .or.   &
       (kpds(5) .eq. 144 .and. kpds(6) .eq. 112) .or. &
       (kpds(5) .eq. 65)) then
     if (kpds(5) .eq. 11 .or. kpds(5) .eq. 144) v=13
     if (kpds(5) .eq. 65) v=12
     call geogfill(nc,nr,fimask,geogmask,varfield,v)
     
     !===> For remaining points not taken care of by above procedure 
     !===> set to user defined soil moisture and 2 m air temperature.
     
     do j=1,nr
        do i=1,nc
           if (varfield(i,j) .eq. -1) then
              if (kpds(5) .eq. 144) varfield(i,j) = ism
              if (kpds(5) .eq. 11 .and. kpds(6) .eq. 112) &
                   varfield(i,j) = geogtemp(i,j)
              if (kpds(5) .eq. 65) varfield(i,j) = 0
           endif
        enddo
     enddo
     
  endif
  
  do j=1,nr
     do i=1,nc
        if ((kpds(5).eq.204).or.(kpds(5).eq.205)) then
           if(fimask(i,j).lt.1) then
              varfield(i,j) = -9999.9
           endif
        else
           if(fimask(i,j).lt.1) then
              varfield(i,j) = udef
           endif
        endif
     enddo
  enddo
  
end subroutine interp_gdas
