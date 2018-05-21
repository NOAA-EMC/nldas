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
! !ROUTINE: glbprecip_cmap.F90
!
! !DESCRIPTION:
!  Includes reading routines for global CMAP precipitation product
!  Used instead of GDAS/GEOS precipitation forcing
!
! !REVISION HISTORY:
!  17 Jul 2001: Jon Gottschalck; Initial code
!  04 Feb 2002: Jon Gottschalck; Added necessary code to use global precip
!               observations with domain 3 (2x2.5)
!  30 Jul 2002: Jon Gottschalck; Added code to use Huffman and Persiann precip data
! !INTERFACE:
subroutine glbprecip_cmap( fname, ferror_cmap, filehr)
! !USES:
  use lisdrv_module, only : lis, gindex      
  use obsprecipforcing_module, only: obsprecip
  use cmapdomain_module, only : cmapdrv
  implicit none
! !ARGUMENTS:
  character(len=80)  :: fname           ! Filename variable for datafile
  integer            :: ferror_cmap
  integer            :: filehr
!EOP

  integer            :: i,j,ios,iret,jret  ! Loop indicies and error flags
  real, pointer      :: precip_regrid(:,:)                                ! Interpolated precipitation array
  integer            :: ncmap
  integer            :: jj,lugb,lugi,kf,kpds(200),k,gridDesccmap(200),jpds(200),jgds(200)
  real               :: ism,udef
  real, allocatable  :: cmapin(:)
  logical*1,allocatable  :: lb(:)
  integer            :: index  
!=== End Variable Definition =======================
!BOC
  allocate (precip_regrid(lis%d%lnc,lis%d%lnr))
  obsprecip     = -1.0
  precip_regrid = -1.0    
!------------------------------------------------------------------------    
! Set necessary parameters for call to interp_gdas    
!------------------------------------------------------------------------    
  ism     = 0
  udef    = lis%d%udef
  jj      = 0
  if (mod((filehr),12).eq.0) then
     lugb=134
  else
     lugb=138
  endif
  ncmap = cmapdrv%ncold*cmapdrv%nrold
  allocate(cmapin(ncmap))
  allocate(lb(ncmap)) 
  lugi    = 0
  jpds    = -1
  jpds(5) = 59
  jpds(6) = 1
  jpds(7) = 0
  jgds    = 0
  cmapin = 0.0
  call baopen (lugb,fname,iret)
  if (iret == 0 ) then
     call getgb (lugb,lugi,ncmap,jj,jpds,jgds,kf,k,kpds,&
          gridDesccmap,lb,cmapin,iret)
!     do j=1,ncmap
!        if(cmapin(j).ne.0) print*, cmapin
!     enddo
     print*, iret
     call interp_cmap(kpds,ncmap,cmapin,lb,lis%d%gridDesc, &
          lis%d%lnc,lis%d%lnr,precip_regrid)
     do j = 1,lis%d%lnr
        do i = 1,lis%d%lnc
!           if(precip_regrid(i,j) .eq. -1) then 
!              print*, j,i,precip_regrid(i,j)
!           endif
           if (precip_regrid(i,j) .ne. -1.0) then
              index = gindex(i,j)
              if(index .ne. -1) then 
                 obsprecip(index) = precip_regrid(i,j)*3600.0
              endif
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
  deallocate(lb)
  deallocate(cmapin)
  !EOC 
end subroutine glbprecip_cmap





