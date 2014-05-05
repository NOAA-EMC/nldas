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
! obsprecip.F90: 
!
! DESCRIPTION:
!   Includes reading routines for global observation precip products. 
!   Used instead of GDAS/GEOS precipitation forcing
!
! REVISION HISTORY:
!  17 Jul 2001: Jon Gottschalck; Initial code
!  04 Feb 2002: Jon Gottschalck; Added necessary code to use global precip
!               observations with domain 3 (2x2.5)
!  30 Jul 2002: Jon Gottschalck; Added code to use Huffman and Persiann precip data
!=========================================================================

!============================================================================
! NRL precipitation read subroutine
!============================================================================

 subroutine glbprecip_nrl( ldas, grid, name, name_trmm, ferror_nrl )
	  
  use ldas_module      ! LDAS non-model-specific 1-D variables
  use grid_module      ! LDAS non-model-specific grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc, ldas%nr)

!==== Local Variables=======================

  integer :: c,r,ferror_nrl,ios,ferrorg,ferrors                ! Loop indicies and error flags
  integer :: i,j,xd,yd,gyd,ibad,mm,nn,xdc,ydc,offsetx,offsety  ! Domain specific parameters
  parameter(xd=1440,yd=480)                                    ! Dimensions of the original NRL data
  parameter(ibad=-32768)                                       ! Bad (missing data) value 
  integer*2 :: rr(xd,yd)                                       ! Original 2 byte integer NRL data
  real :: precip(xd,yd)                                        ! Orginal real precipitation array
  real, pointer  :: precip_regrid(:,:)                         ! Interpolated precipitation array
  character(len=80) :: name, fname, fname2, name_trmm          ! Filename variables for both the IR and microwave datafiles

!=== End Variable Definition =======================
   
   fname = name

!=== Determine domain specific parameters necessary for interpolation 

   if (ldas%domain .eq. 1) then
    xdc = 2880
    ydc = ldas%nrgpcp
    offsetx = 440
    offsety = 680
   else
    xdc = ldas%nc
    ydc = ldas%nrgpcp
   endif
   
   allocate (precip_regrid(xdc,ydc))
   
!=== Fill necessary arrays to assure not using old NRL data
   
    precip = -1.0
    precip_regrid = -1.0
    
!=== Find NRL Geostationary precip, read it in and assign to forcing precip array.
!=== **Must reverse grid in latitude dimension to be consistent with LDAS grid**

   open(unit=10,file=fname, status='old', &
&          access='direct',recl=xd*yd*2, &
&          form='unformatted',iostat=ios)
   if (ios .eq. 0) then   
     read(10,rec=1) rr    ! all rainrates 
     do i = 1,yd
      do j = 1,xd
        precip(j,i) = float(rr(j,yd+1-i)) / 100.0
	if ( rr(j,yd+1-i) .eq. ibad ) precip(j,i) = -1.0
      enddo
     enddo

!=== Interpolating to desired domain and resolution
!=== ****** Global precip datasets not used currently to force NLDAS ******
     
     select case (ldas%domain)
     
!     case (1)
     
!      call interp_glbpcp(-59875,59875,-179875,179875,-59938,59938,-179938,179938,&
!                         xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
!      do i = 1,ldas%nr
!       do j = 1,ldas%nc
!	if (precip_regrid(offsetx+j,offsety+i) .ne. -1.0) then
!	 grid(j,i)%obsprecip = precip_regrid(offsetx+j,offsety+i)
!	endif
!       enddo
!      enddo								    
     
     case (3)
     
      call interp_glbpcp(-59875,59875,-179875,179875,-60000,58000,-180000,177500,&
      xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
       do i = 1,ydc
        do j = 1,xdc
          if (precip_regrid(j,i) .ne. -1.0) then
            grid(j,i)%obsprecip = precip_regrid(j,i)
          endif
        enddo
       enddo
       
     case (2)
     
      do i = 1,yd
       do j = 1,xd
        if (precip(j,i) .ne. -1.0) then         
          grid(j,i)%obsprecip = precip(j,i)
        endif
       enddo
      enddo
      
     case(4)
     
      call interp_glbpcp(-59875,59875,-179875,179875,-59500,59500,-179500,179500,&
                          xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
      do i = 1,ydc
       do j = 1,xdc
	 if (precip_regrid(j,i) .ne. -1.0) then
	     grid(j,i)%obsprecip = precip_regrid(j,i)
	 endif
       enddo
      enddo
      
     case(5)
      
       call interp_glbpcp(-59875,59875,-179875,179875,-59750,59750,-179750,179750,&
                          xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
	do i = 1,ydc
	 do j = 1,xdc
	  if (precip_regrid(j,i) .ne. -1.0) then
            grid(j,i)%obsprecip = precip_regrid(j,i)
	  endif
	 enddo
	enddo				   
								    
     case default
     
      print *, "No valid global grid defined for given resolution"
      print *, "columns: ", ldas%nc, " rows: ", ldas%nr
      print *, "Stopping..."
      stop
				     
     end select
     
     ferrorg = 1  
     close(10)
    print*, "Obtained NRL geostationary precipitation data ", fname
   else
    print*, "Missing NRL geostationary precipitation data ", fname
    ferrorg = 0
   endif
   
!=== Preparing for SSMI/TRMM filename
    fname2 = name_trmm

!=== Find NRL Merged TRMM/SSMI precip, read it in and assign to forcing precip array.
!=== **Must reverse grid in latitude dimension to be consistent with LDAS grid**
!=== Also only replace Geostationary with valid SSMI/TRMM

    open(unit=10,file=fname2, status='old', &
&          access='direct',recl=xd*yd*2, &
&          form='unformatted',iostat=ios)

    if (ios .eq. 0) then
     read(10,rec=1) rr    ! all rainrates
     do i = 1,yd
      do j = 1,xd
        if (rr(j,yd+1-i) .ne. ibad) then
	   precip(j,i) = float(rr(j,yd+1-i)) / 100.0
	endif
      enddo
     enddo

!=== Interpolating to desired domain and resolution
!=== ****** Global precip datasets not used currently to force NLDAS ******indent: Command not found.
     
     select case (ldas%domain)
     
!     case (1)
     
!      call interp_glbpcp(-59875,59875,-179875,179875,-59938,59938,-179938,179938,&
!                         xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
!       do i = 1,ldas%nr
!	do j = 1,ldas%nc
!	 if (precip_regrid(offsetx+j,offsety+i) .ne. -1.0) then
!          grid(j,i)%obsprecip = precip_regrid(offsetx+j,offsety+i)
!	 endif
!	enddo
!       enddo
     
      case (3)
     
      call interp_glbpcp(-59875,59875,-179875,179875,-60000,58000,-180000,177500,&
      xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
       do i = 1,ydc
        do j = 1,xdc
	 if (precip_regrid(j,i) .ne. -1.0) then
         grid(j,i)%obsprecip = precip_regrid(j,i)
	 endif
        enddo
       enddo

     case(2)
     
      do i = 1,yd
       do j = 1,xd
        if (precip(j,i) .ne. -1.0) then
         grid(j,i)%obsprecip = precip(j,i)
	endif	
       enddo
      enddo
      
     case(4)
     
      call interp_glbpcp(-59875,59875,-179875,179875,-59500,59500,-179500,179500,&
                          xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
      do i = 1,ydc
       do j = 1,xdc
	if (precip_regrid(j,i) .ne. -1.0) then
         grid(j,i)%obsprecip = precip_regrid(j,i)
        endif
       enddo
      enddo
				   
     case(5)
					
      call interp_glbpcp(-59875,59875,-179875,179875,-59750,59750,-179750,179750,&
                          xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
      do i = 1,ydc
       do j = 1,xdc
	if (precip_regrid(j,i) .ne. -1.0) then
          grid(j,i)%obsprecip = precip_regrid(j,i)
	endif
       enddo
      enddo

     case default
     
      print *, "No valid global grid defined for given resolution"
      print *, "columns: ", ldas%nc, " rows: ", ldas%nr
      print *, "Stopping..."
      stop
     
     end select

     ferrors = 1
     close(10)
     print*, "Obtained NRL TRMM/SSMI precipitation data ", fname2
   else
     print*, "Missing NRL TRMM/SSMI precipitation data ", fname2
     ferrors = 0
   endif
   
   if (ferrorg .eq. 1 .or. ferrors .eq. 1) ferror_nrl = 1
   
   deallocate (precip_regrid)
   
end subroutine glbprecip_nrl

!============================================================================
! HUFFMAN precipitation read subroutine
!============================================================================

subroutine glbprecip_huff ( ldas, grid, name_huff, ferror_huff )

 use ldas_module      ! LDAS non-model-specific 1-D variables
 use grid_module      ! LDAS non-model-specific grid variables
 implicit none
 type (ldasdec) ldas
 type (griddec) grid(ldas%nc, ldas%nr)
	 
!==== Local Variables=======================
	 
 integer :: c,r,ferror_huff,ios,lrec,offsetx,offsety  ! Loop indicies and error flags
 integer :: i,j,xd,yd,gyd,ibad,mm,nn,xdc,ydc          ! Domain specific parameters
 parameter(xd=1440,yd=480)                            ! Dimension of original HUFFMAN data
 parameter(ibad=-31999)                               ! Bad (missing data) value
 integer*2 :: rr(xd,yd)                               ! Original 2 byte integer HUFFMAN data 
 real :: precip(xd,yd)                                ! Original real precipitation array
 real, pointer :: precip_regrid(:,:)                  ! Interpolated precipitation array
 character(len=80) :: fname, name_huff                ! Filename variables
 character(len=2880) ::  head                         ! Header variable
		       
!=== End Variable Definition =======================
		       
 fname = name_huff

!=== Determine domain specific parameters necessary for interpolation

 if (ldas%domain .eq. 1) then
   xdc = 2880
   ydc = ldas%nrgpcp
   offsetx = 440
   offsety = 680
 else
   xdc = ldas%nc
   ydc = ldas%nrgpcp
 endif

 allocate (precip_regrid(xdc,ydc))
 
!=== Fill necessary arrays to assure not using old HUFFMAN data

 precip = -1.0
 precip_regrid = -1.0
	    
!=== Find HUFFMAN precip data, read it in and assign to forcing precip array.
!=== **Must reverse grid in latitude dimension to be consistent with LDAS grid**

  open(unit=10,file=fname, status='old', &
   &          access='direct',recl=xd*yd*4, &
   &          form='unformatted',iostat=ios)

  if (ios .eq. 0) then
     read (10,rec=1) head (1:2880),&
          ( ( rr (i, j), i = 1, xd ), j = 1, yd)
   do i = 1,yd
    do j = 1,xd
     if (j .lt. 721) then
      precip(j,i) = float(rr(j+720,yd+1-i)) / 100.0
      if ( rr(j+720,yd+1-i) .eq. ibad ) precip(j,i) = -1.0 
      if ( rr(j+720,yd+1-i) .lt. 0.0 )  precip(j,i) = -1.0
     else
      precip(j,i) = float(rr(j-720,yd+1-i)) / 100.0
      if ( rr(j-720,yd+1-i) .eq. ibad ) precip(j,i) = -1.0
      if ( rr(j-720,yd+1-i) .lt. 0.0 )  precip(j,i) = -1.0
     endif
    enddo
  enddo

!=== Interpolating to desired domain and resolution
!=== ****** Global precip datasets not used currently to force NLDAS ******

  select case (ldas%domain)
  
!  case (1)
  
!   call interp_glbpcp(-59875,59875,-179875,179875,-59938,59938,-179938,179938,&
!                         xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
!   do i = 1,ldas%nr
!    do j = 1,ldas%nc
!     if (precip_regrid(offsetx+j,offsety+i) .ne. -1.0) then
!      grid(j,i)%obsprecip = precip_regrid(offsetx+j,offsety+i)
!     endif
!    enddo
!   enddo
			      
  case (3)

   call interp_glbpcp(-59875,59875,-179875,179875,-60000,58000,-180000,177500,&
                       xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
  do i = 1,ydc
   do j = 1,xdc
    if (precip_regrid(j,i) .ne. -1.0) then
     grid(j,i)%obsprecip = precip_regrid(j,i)
    endif
   enddo
  enddo

  case(2)

  do i = 1,yd
   do j = 1,xd
     if (precip(j,i) .ne. -1.0) then
       grid(j,i)%obsprecip = precip(j,i)
     endif
   enddo
  enddo
  
  case(4)
  
   call interp_glbpcp(-59875,59875,-179875,179875,-59500,59500,-179500,179500,&
                          xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
   do i = 1,ydc
    do j = 1,xdc
     if (precip_regrid(j,i) .ne. -1.0) then
      grid(j,i)%obsprecip = precip_regrid(j,i)
     endif
    enddo
   enddo
			      
  case(5)
				   
   call interp_glbpcp(-59875,59875,-179875,179875,-59750,59750,-179750,179750,&
                          xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
   do i = 1,ydc
    do j = 1,xdc
     if (precip_regrid(j,i) .ne. -1.0) then
      grid(j,i)%obsprecip = precip_regrid(j,i)
     endif
    enddo
   enddo

  case default
  
   print *, "No valid global grid defined for given resolution"
   print *, "columns: ", ldas%nc, " rows: ", ldas%nr
   print *, "Stopping..."
   stop
  
  end select

  ferror_huff = 1
  close(10)
  print*, "Obtained HUFFMAN precipitation data ", fname
 else
   print*, "Missing HUFFMAN precipitation data ", fname
   ferror_huff = 0
 endif
 
 deallocate (precip_regrid)
 
end subroutine glbprecip_huff

!============================================================================
! PERSIANN precipitation read subroutine
!============================================================================

subroutine glbprecip_pers ( ldas, grid, name_pers, ferror_pers )

  use ldas_module      ! LDAS non-model-specific 1-D variables
  use grid_module      ! LDAS non-model-specific grid variables
  implicit none
  type (ldasdec) ldas
  type (griddec) grid(ldas%nc, ldas%nr)
      
!==== Local Variables=======================
  
  integer :: c,r,ferror_pers,ios,lrec,offsetx,offsety      ! Loop indicies and error flags
  integer :: i,j,xd,yd,gyd,ibad,mm,nn,xdc,ydc,offset       ! Domain specific parameters
  parameter(xd=1440,yd=400)                                ! Dimension of original PERSIANN data
  parameter(ibad=-9999.0)                                  ! Bad (missing data) value
  real :: rr(xd,yd)                                        ! Original precip array 
  real :: precip(xd,yd)                                    ! Reconfigured original precip array
  real, pointer :: precip_regrid(:,:)                      ! Interpolated precip array
  character(len=80) :: fname, name_pers                    ! Filename variables
 
!=== End Variable Definition =======================
 
  fname = name_pers

!=== Determine domain specific parameters necessary for interpolation

 if (ldas%domain .eq. 1) then
   xdc = 2880
   ydc = ldas%nrgpcp - (20.0 / (ldas%ldas_kgds(9)/1000.0))
   offsetx = 440
   offsety = 600
 else
   xdc = ldas%nc
   ydc = ldas%nrgpcp - (20.0 / (ldas%ldas_kgds(9)/1000.0))
 endif

  allocate (precip_regrid(xdc,ydc))
  
!=== Fill necessary arrays to assure not using old PERSIANN data
  
  precip = -1.0
  precip_regrid = -1.0

!=== Determine offset in number of rows from 60 S since PERSIANN starts at 50 S
  
  select case (ldas%domain)
  
  case(3)
   offset = 6
  case(2)
   offset = (ldas%nr / 150.0) * 10.0
  case(4)
   offset = (ldas%nr / 150.0) * 10.0
  case(5) 
   offset = (ldas%nr / 150.0) * 10.0 
  case default
 
    print *, "No valid global grid defined for given resolution"
    print *, "columns: ", ldas%nc, " rows: ", ldas%nr
    print *, "Stopping in global PERSIANN precip routine ... "
    stop
 
   end select
  
!=== Find PERSIANN precip data, read it in and assign to forcing precip array.
!=== **Must reverse grid in latitude dimension to be consistent with LDAS grid**

  open(unit=10,file=fname, status='old',access='direct', &
      form='unformatted',recl=xd*yd*4,iostat=ios)
      
  if (ios .eq. 0) then
    read (10,rec=1) rr

    do i = 1,yd
      do j = 1,xd
        if (j .lt. 721) then
	 precip(j,i) = rr(j+720,yd+1-i)
	 if ( rr(j+720,yd+1-i) .eq. ibad .or. rr(j+720,yd+1-i) .lt. 0.0) precip(j,i) = -1.0
	else
	 precip(j,i) = rr(j-720,yd+1-i)
	 if ( rr(j-720,yd+1-i) .eq. ibad .or. rr(j-720,yd+1-i) .lt. 0.0) precip(j,i) = -1.0
	endif
      enddo
    enddo

!=== Interpolating to desired domain and resolution
!=== ****** Global precip datasets not used currently to force NLDAS ******

  select case (ldas%domain)
  
!    case (1)
    
!    call interp_glbpcp(-49750,50000,-180000,179750,-49938,49938,-179938,179938,&
!                         xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
!     do i = 1,ldas%nr
!      do j = 1,ldas%nc
!       if (precip_regrid(offsetx+j,offsety+i) .ne. -1.0) then
!	grid(j,i)%obsprecip = precip_regrid(offsetx+j,offsety+i)
!       endif
!      enddo
!     enddo
  
    case (3)

    call interp_glbpcp(-49750,50000,-180000,179750,-48000,50000,-180000,177500,& 
                        xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
    do i = 1,ydc
      do j = 1,xdc
        if (precip_regrid(j,i) .ne. -1.0) then
          grid(j,i+offset)%obsprecip = precip_regrid(j,i)
	endif
      enddo
    enddo

   case(2)

   call interp_glbpcp(-49750,50000,-180000,179750,-49875,49875,-179875,179875,&
                       xd,yd,xd,yd,precip,precip,ldas%domain)
    do i = 1,yd
      do j = 1,xd
	if (precip(j,i) .ne. -1.0) then
	  grid(j,i+offset)%obsprecip = precip(j,i)
	endif
      enddo
    enddo
    
    case(4)
    
     call interp_glbpcp(-49750,50000,-180000,179750,-49500,49500,-179500,179500,&
                          xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
     do i = 1,ydc
      do j = 1,xdc
       if (precip_regrid(j,i) .ne. -1.0) then
	grid(j,i+offset)%obsprecip = precip_regrid(j,i)
       endif
      enddo
     enddo
							
    case(5)
							  
     call interp_glbpcp(-49750,50000,-180000,179750,-49750,49750,-179750,179750,&
                          xd,yd,xdc,ydc,precip,precip_regrid,ldas%domain)
     do i = 1,ydc
      do j = 1,xdc
       if (precip_regrid(j,i) .ne. -1.0) then
	grid(j,i+offset)%obsprecip = precip_regrid(j,i)
       endif
      enddo
     enddo

    case default
    
     print *, "No valid global grid defined for given resolution"
     print *, "columns: ", ldas%nc, " rows: ", ldas%nr
     print *, "Stopping..."
     stop
    
    end select

     ferror_pers = 1
     close(10)
     print*, "Obtained PERSIANN precipitation data ", fname
   else
     print*, "Missing PERSIANN precipitation data ", fname
     ferror_pers = 0
   endif
	
   deallocate (precip_regrid)
	
  end subroutine glbprecip_pers

!============================================================================
! CMAP precipitation read subroutine
!============================================================================

 subroutine glbprecip_cmap( ldas, grid, fname, ferror_cmap, filehr )

  use ldas_module      ! LDAS non-model-specific 1-D variables
  use grid_module      ! LDAS non-model-specific grid variables

  implicit none

  type (ldasdec) ldas
  type (griddec) grid(ldas%nc, ldas%nr)

!==== Local Variables=======================

  integer            :: i,j,ferror_cmap,ios,filehr,iret,jret              ! Loop indicies and error flags
  real, pointer      :: precip_regrid(:,:)                                ! Interpolated precipitation array
  character(len=80)  :: fname                                             ! Filename variable for datafile
  integer, parameter :: ncmap=131072
  integer            :: jj,lugb,lugi,kf,kpds(25),k,kgdscmap(200),jpds(25),jgds(22)
  real               :: cmapin(ncmap),ism,udef
  logical*1          :: lb(ncmap)
  
!=== End Variable Definition =======================

!=== Determine domain specific parameters necessary for interpolation

  allocate (precip_regrid(ldas%nc,ldas%nr))
  
!=== Fill necessary arrays to assure not using old CMAP data

  precip_regrid = -1.0    
    
!=== Set necessary parameters for call to interp_gdas    
  ism     = 0
  udef    = ldas%udef	    
  jj      = 0
  if (mod((filehr),12).eq.0) then
   lugb=1
  else
   lugb=2
  endif				      
  lugi    = 0
  jpds    = -1
  jpds(5) = 59
  jpds(6) = 1
  jpds(7) = 0
  jgds    = 0

  call baopen (lugb,fname,iret)

  if (iret == 0 ) then
  
  call getgb (lugb,lugi,ncmap,jj,jpds,jgds,kf,k,kpds,kgdscmap,lb,cmapin,iret)
  
!=== Interpolating to desired domain and resolution

        call interp_gdas (kpds,kgdscmap,ncmap,cmapin,lb,ldas%ldas_kgds, &
      	     ldas%nc,ldas%nr,precip_regrid,grid%fimask,ism,ldas%udef)
          
       do i = 1,ldas%nr
        do j = 1,ldas%nc
          if (precip_regrid(j,i) .ne. -1.0) then
!=== Precipitation rate is multiplied here by 3600 (since it is in kg/m2/s). It is
!=== done because in glbobspcp.F90, the returned array is divided by 3600.	  
            grid(j,i)%obsprecip = precip_regrid(j,i)*3600.0
          endif
        enddo
       enddo

     call baclose (lugb,jret)

     ferror_cmap = 1
     close(10)
    print*, "Obtained CMAP CPC precipitation data ", fname
   else
    print*, "Missing CMAP CPC precipitation data ", fname
    ferror_cmap = 0
   endif

   deallocate (precip_regrid)

  end subroutine glbprecip_cmap

      
!==============================================================================

subroutine interp_glbpcp(in_edges,in_edgen,in_edgew,in_edgee,out_edges,out_edgen,out_edgew,out_edgee,&
                         xd,yd,xdc,ydc,precip,precip_regrid,d)

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
! interp_glbpcp.f90:
!
! DESCRIPTION:
!   Interpolates 1/4 global observed precip products to coarser resolutions
!
! REVISION HISTORY:
!  04 Feb 2002: Jon Gottschalck; Initial code
!  30 Jul 2002: Jon Gottschalck; Added additional domains
!=========================================================================

  implicit none

  integer :: xd,yd,i,j,c,nquart,ncourse,iret,xdc,ydc,d,dx,dy     ! Size of orginal and interpolated precipitation arrays
  integer :: kgdsin(200),kgdspcp(200),ip,ibi,km,ipopt(20),no,ibo ! Variables used by IPOLATES
  real :: precip(xd,yd),precip_regrid(xdc,ydc)                   ! Precipitation arrays for orginal and interpolated resolutions
  real :: f(xd*yd),g(xdc*ydc),rlat(xdc*ydc),rlon(xdc*ydc)        ! 1-D arrays needed by IPOLATES for input data, interpolated data, latitude, and longitude
  integer :: in_edges,in_edgen,in_edgew,in_edgee                 ! Edges of input grids, used as input to IPOLATES
  integer :: out_edgen,out_edges,out_edgew,out_edgee             ! Edges of output grids, used as input to IPOLATES
  logical*1 :: li(xd*yd),lo(xdc*ydc)                             ! 1-D valid data bitmap arrays needed by IPOLATES
  logical*1 :: ligrid(xd,yd),logrid(xdc,ydc)                     ! 2-D valid data bitmap arrays 

!=== Determine interpolated grid spacing

  select case (d)

!=== Not used for NLDAS currently  
!  case(1)
  
!   dx = 125
!   dy = 125
  
  case(3)
  
   dx = 2500
   dy = 2000
  
  case(2)
  
   dx = 250
   dy = 250
   
  case(4)
   
   dx = 1000
   dy = 1000
   
  case(5)
  
   dx = 500
   dy = 500
   
  case default
  
   print *, "No valid global grid defined for given resolution"
   print *, "Stopping in global precip interpolation ... "
   stop
  
  end select

!=== Setting up parameters needed by IPOLATES

  nquart = xd*yd
  ncourse = xdc*ydc
  
  do i=1,200
   kgdsin(i) = 0
   kgdspcp(i) = 0
  enddo
  
  kgdsin(1) = 0
  kgdsin(2) = xd
  kgdsin(3) = yd
  kgdsin(4) = in_edges
  kgdsin(5) = in_edgew
  kgdsin(6) = 128
  kgdsin(7) = in_edgen
  kgdsin(8) = in_edgee
  kgdsin(9) = dy
  kgdsin(10) = dx
  kgdsin(11) = 64
  kgdsin(20) = 255
  
  c=0
  do i=1,yd
   do j=1,xd
     c = c + 1
     f(c) = precip(j,i)
   enddo
  enddo

  kgdspcp(1)=0
  kgdspcp(2)=xdc
  kgdspcp(3)=ydc
  kgdspcp(4)=out_edges
  kgdspcp(5)=out_edgew
  kgdspcp(6)=128
  kgdspcp(7)=out_edgen
  kgdspcp(8)=out_edgee
  kgdspcp(9)=dy
  kgdspcp(10)=dx
  kgdspcp(11)=64
  kgdspcp(20)=255

!=== Setting interpolation options (ip=3,budget bilinear)
!=== (km=1, one parameter, ibi=1,use undefined bitmap

  ip=3
  ipopt(1) = -1
  ipopt(2) = -1
  km=1
  ibi=1
  
  do c=1,nquart
   if (f(c) .eq. -1.0) then 
     li(c) = .FALSE.
   else
     li(c) = .TRUE.
   endif
  enddo
  
  c=0
  do j = 1,yd
   do i = 1,xd
     ligrid(i,j) = li(i+c)
   enddo
   c = c + xd
  enddo
  
!=== Initializing output bitmap.

  do i=1,ncourse
   lo(i) = .TRUE.
  enddo

!=== Interpolation to new grid

  call ipolates (ip,ipopt,kgdsin,kgdspcp,nquart,ncourse, &
&                             km,ibi,li,f,no,rlat,       &
&                             rlon,ibo,lo,g,iret)

!=== Reformatting to 2-D grid assigning -1.0 to undefined data 
  c=0
  do j = 1,ydc
   do i = 1,xdc
     precip_regrid(i,j) = g(i+c)
     logrid(i,j) = lo(i+c)
     if (.NOT. logrid(i,j)) precip_regrid(i,j) = -1.0
   enddo
   c = c + xdc
  enddo

end subroutine interp_glbpcp


