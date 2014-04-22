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
! !ROUTINE: maketiles_gswp.F90
!
! !DESCRIPTION:
!  This primary goal of this routine is to determine tile space for 
!  GSWP data sets
!
! !REVISION HISTORY:
! 23Feb04, Sujay Kumar : Intial Specification
!
! !INTERFACE:
subroutine maketiles_gswp()
! !USES:
  use lisdrv_module, only: lis, grid, glbgindex, tile
  use grid_module
  use spmdMod
#if ( defined USE_NETCDF )
  use netcdf
#endif
!EOP
  IMPLICIT NONE
  real, allocatable :: mask(:,:)
  real, allocatable :: vegclass(:)
  real,allocatable :: fgrd(:)
  real, allocatable :: nav_lat(:,:),nav_lon(:,:)
! for writing dominant veg types..
  real, allocatable :: domveg(:,:)
!=== Local Variables =====================================================
  integer :: c,r,count  ! Loop counters
  integer :: ierr
  integer :: gnc, gnr
  integer :: cindex, rindex  
  integer :: ncid, status
  integer :: mvarid, lcvarid, latid, lonid
!=== End Variable Definition =============================================
!BOC
  if ( masterproc ) then
     if(lis%d%gridDesc(42) > lis%d%lnc .or. &
          lis%d%gridDesc(43) > lis%d%lnr)  then !using a subdomain
        gnc = lis%d%gridDesc(42)
        gnr = lis%d%gridDesc(43)
     else
        gnc = lis%d%lnc
        gnr = lis%d%lnr
     endif
     lis%d%gnc = gnc
     lis%d%gnr = gnr
     
     allocate(mask(lis%d%lnc, lis%d%lnr), stat=ierr)
     allocate(nav_lat(lis%d%lnc, lis%d%lnr), stat=ierr)
     allocate(nav_lon(lis%d%lnc, lis%d%lnr), stat=ierr)
     call check_error(ierr,'Error allocating mask.',iam)
          
#if ( defined USE_NETCDF )
     print*,'MSG: maketiles_gswp -- Reading landmask (',iam,')'
     status = nf90_open(path=lis%p%mfile,mode= nf90_nowrite,ncid = ncid)
     status = nf90_inq_varid(ncid, "landmask",mvarid)
     status = nf90_inq_varid(ncid, "nav_lat",latid)
     status = nf90_inq_varid(ncid, "nav_lon",lonid)
     status = nf90_get_var(ncid, mvarid,mask)
     status = nf90_get_var(ncid, latid,nav_lat)
     status = nf90_get_var(ncid, lonid,nav_lon)
     status = nf90_close(ncid)
     print*,'MSG: maketiles_gswp -- Done reading ',trim(lis%p%mfile), & 
          ' (',iam,')'
#endif

!----------------------------------------------------------------------
!  Make Tile Space
!----------------------------------------------------------------------
      lis%d%glbnch=0
      do r=1,lis%d%lnr      
         do c=1,lis%d%lnc   
            if(mask(c,r).gt.0.99.and. & 
                 mask(c,r).lt.3.01)then !we have land
               lis%d%glbnch=lis%d%glbnch+1 
            endif
         enddo
      enddo
      print*, 'DBG: maketiles_gswp -- glbnch',lis%d%glbnch,' (',iam,')'
      allocate(tile(lis%d%glbnch))

      lis%d%glbngrid=0
      do r=1,lis%d%lnr
         do c=1,lis%d%lnc
            if(mask(c,r).gt.0.99 .and. & 
                 mask(c,r).lt.3.01) then
               lis%d%glbngrid=lis%d%glbngrid+1
            endif
         enddo
      enddo
      count = 1
      print*, 'DBG: maketiles_gswp -- glbnch',lis%d%glbnch,' (',iam,')'
      allocate(grid(lis%d%glbngrid))
      allocate(glbgindex(lis%d%lnc, lis%d%lnr))
      print*, 'DBG: maketiles_gswp -- glbnch',lis%d%glbnch,' (',iam,')'

      allocate(vegclass(lis%d%glbngrid), stat=ierr)
      call check_error(ierr,'Error allocating vegclass.',iam)
      
      allocate(fgrd(lis%d%glbngrid), stat=ierr)
      call check_error(ierr,'Error allocating vegclass.',iam)

#if ( defined USE_NETCDF )
     print*,'MSG: maketiles_gswp -- Reading vegclass (',iam,')'
     status = nf90_open(path=lis%p%vfile,mode= nf90_nowrite,ncid = ncid)
     status = nf90_inq_varid(ncid, "VegClass",lcvarid)
     status = nf90_inq_varid(ncid, "nav_lat",latid)
     status = nf90_inq_varid(ncid, "nav_lon",lonid)
     status = nf90_get_var(ncid, lcvarid,vegclass)
     status = nf90_get_var(ncid, latid,nav_lat)
     status = nf90_get_var(ncid, lonid,nav_lon)

     status = nf90_close(ncid)
#endif
     print*,'MSG: maketiles_gswp -- Done reading ',trim(lis%p%vfile), & 
          ' (',iam,')'
     do c=1,lis%d%lnc
        do r=1,lis%d%lnr
           rindex = (nav_lat(c,r)+59.5) + 1
           cindex = (nav_lon(c,r)+179.5) + 1
           glbgindex(cindex,rindex) = -1
           if(mask(cindex,rindex).gt.0.99 .and. & 
                mask(cindex,rindex).lt.3.01) then
              grid(count)%lat = nav_lat(c,r)
              grid(count)%lon = nav_lon(c,r)
              grid(count)%fgrd = 1
              glbgindex(cindex,rindex) = count
              count = count+1
           endif
        enddo
     enddo
     print*, 'DBG: maketiles_gswp -- glbnch',lis%d%glbnch,' (',iam,')'
!--------------------------------------------------------------------
!   For writing dominant Vegetation types
!--------------------------------------------------------------------
     if(lis%o%wparam .eq.1) then 
        allocate(domveg(lis%d%lnc,lis%d%lnr))
        domveg = -9999.0
     endif
     count = 0
     do c=1,lis%d%lnc            
        do r=1,lis%d%lnr
           rindex = (nav_lat(c,r)+59.5) + 1
           cindex = (nav_lon(c,r)+179.5) + 1
            if(mask(cindex,rindex).gt.0.99.and. & 
                 mask(cindex,rindex).lt.3.01)then 
!               if(fgrd(c,r).gt.0.0)then
                  count = count+1
                  tile(count)%row=rindex
                  tile(count)%col=cindex    
                  tile(count)%index = glbgindex(cindex,rindex)
!                  print*, count, grid(count)%lat, grid(count)%lon
!                  print*, glbgindex(c,r),count,vegclass(glbgindex(c,r))
                  tile(count)%vegt= vegclass(count)
                  if(lis%o%wparam.eq.1) then 
                     domveg(cindex,rindex) = vegclass(count)
                  endif
                  tile(count)%fgrd=1
!               endif
            endif
         enddo
      enddo

      if(lis%o%wparam.eq.1) then 
         open(32,file="domvegtype.bin",form='unformatted')
         write(32) domveg
         close(32)
         deallocate(domveg)
      endif
      deallocate(mask,stat=ierr)
      deallocate(vegclass,stat=ierr)
      deallocate(fgrd, stat=ierr)
      deallocate(nav_lat)
      deallocate(nav_lon)
      call check_error(ierr,'Error allocating glbfgrd',iam)
      
      write(*,*) 'MSG: maketiles_gswp -- Actual Number of Tiles:', & 
           lis%d%glbnch,' (',iam,')'
      write(*,*)
      
      write(*,*) 'MSG: maketiles_gswp -- Size of Grid Dimension:', & 
           lis%d%glbngrid,' (',iam,')'
      write(*,*)
      
   endif
   print*,'MSG: maketiles_gswp -- done',' (',iam,')'   
   return
!EOC
 end subroutine maketiles_gswp
 
