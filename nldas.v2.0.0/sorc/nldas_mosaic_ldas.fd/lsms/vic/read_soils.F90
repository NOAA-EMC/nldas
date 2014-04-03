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
! !ROUTINE: read_soils.F90
!
! !DESCRIPTION:
! Reads the soil files
!
! !REVISION HISTORY:
! 
! 29 Mar 2004; Sujay Kumar : Initial Specification
! 
! !INTERFACE:
subroutine read_soils()
!EOP

  use lisdrv_module, only : lis
  use vic_varder, only : vicdrv
  use lis_openfileMod
  use lis_indices_module

  implicit none

  real    :: lat, lon
  integer :: r, c, cindex, rindex
  integer :: mlat, mlon,line,line1,line2,glnc,glnr
  integer :: ios1
  real    :: temp1, temp2
  real, allocatable :: sand1(:,:),clay1(:,:), type(:)


  line1 = (lis%d%gridDesc(4)-lis%d%gridDesc(44))/lis%d%gridDesc(49)+1
  line2 = (lis%d%gridDesc(5)-lis%d%gridDesc(45))/lis%d%gridDesc(50)+1

  allocate(type(lis%d%nch))

!  if(lis%d%domain.le.6) then 
     allocate(sand1(lis_nc_data,lis_nr_data)) !gnc,gnr
     allocate(clay1(lis_nc_data,lis_nr_data)) !gnc,gnr
     
     call readsand(lis%d%soil, sand1)
     call readclay(lis%d%soil, clay1)

!     read(15) sand1
!     read(16) clay1
     call get_soil_type(lis_nc_data, lis_nr_data, lis%d%nch, sand1, clay1, type)
     deallocate(sand1)
     deallocate(clay1)
!     close(15)
!     close(16)
#if 0 
  else if(lis%d%domain.eq.8) then 
     allocate(sand1(lis_nc_data,lis_nr_data)) !lnc,lnr
     allocate(clay1(lis_nc_data,lis_nr_data)) !lnc,lnr
     !open(11,file=lis%p%safile,status='old',form='unformatted',&
     !     access ='direct',recl=4,iostat=ios1)
     !open(12,file=lis%p%clfile,status='old',form='unformatted',&
     !     access ='direct',recl=4,iostat=ios1)
     call lis_open_file(11,file=lis%p%safile,status='old',form='unformatted',&
                        access='sequential',script='getsand.pl')
     call lis_open_file(12,file=lis%p%clfile,status='old',form='unformatted',&
                        access='sequential',script='getclay.pl')
     call get_soil_type(lis_nc_data, lis_nr_data, lis%d%nch, sand1, clay1, type)
     deallocate(sand1)
     deallocate(clay1)
     close(11)
     close(12)
  endif
#endif  
  call set_soilparam(trim(vicdrv%vic_sfile), len(trim(vicdrv%vic_sfile)), & 
       type, lis%d%nch, vicdrv%vic_nlayer)

  deallocate(type)

end subroutine read_soils

subroutine get_soil_type(nc, nr, nch, sand, clay, tex)
  use lisdrv_module, only : tile
  use lis_indices_module

  implicit none
  
  integer :: nc,nr,nch
  real    :: sand(nc,nr),clay(nc,nr)
  integer :: tex(nch)
  
  integer :: i,j
  real    :: frac,sa,si,cl

  TEX = 0

  do i=1,nch
!	    Test for ocean points.
     if(clay(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.0) then 
        tex(i) = -99
     else if(clay(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.23) then 
        if(sand(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.50) then 
           tex(i) = 8 !Loam
        elseif(sand(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.75) then 
           tex(i) = 4  ! sandy loam
        else
           tex(i) = 1  !loamy sand
        endif
     elseif(clay(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.28) then 
        if(sand(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.45) then 
           tex(i) = 4 !sandy loam
        else
           tex(i) = 7  !sandy clay loam
        endif
     elseif(clay(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.37) then 
        if(sand(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.2) then 
           tex(i) = 2 ! silty clay loam
        elseif(sand(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.43) then 
           tex(i) = 6 !clay loam
        else
           tex(i) = 7 ! sandy clay loam
        endif
     elseif(clay(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.41) then 
        if(sand(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.2) then 
           tex(i) = 2 ! silty clay loam
        elseif(sand(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.43) then 
           tex(i) = 6 !clay loam
        else
           tex(i) = 5 ! sandy clay
        endif
     else
        if(sand(tile(i)%col, tile(i)%row-lis_tnroffset).lt.0.43) then 
           tex(i) = 3 ! light clay
        else
           tex(i) = 5 ! sandy clay
        endif

     endif

  enddo
  return
end subroutine get_soil_type
