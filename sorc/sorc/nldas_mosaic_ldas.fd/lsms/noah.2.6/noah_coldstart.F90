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
! !ROUTINE: noah_coldstart.F90
!
! !DESCRIPTION:
!  
! Routine for noah initialization from cold start
! 
! !INTERFACE:
subroutine noah_coldstart()
! !USES:
   use lisdrv_module, only: lis, tile, gindex
   use noah_varder
   use time_manager
   use spmdMod, only: iam
#if ( defined USE_NETCDF )
   use netcdf
#endif
!EOP
   implicit none

   integer :: t,l
   real :: soiltemp(lis%d%lnc, lis%d%lnr)
   real :: soiltemp1(lis%d%glbnch)
   integer :: ncid, status
   integer :: stid
   integer :: cindex,rindex,c,r

!BOC
   if ( lis%o%startcode == 2 ) then
      call lis_log_msg('MSG: noah_coldstart -- cold-starting noah: '// &
                       'using ics from card file')
   elseif ( lis%o%startcode == 3 ) then
      call lis_log_msg('MSG: noah_coldstart -- cold-starting noah: '//         &
                       'using ics from card file, reading intial soil temp '// &
                       'from GSWP-2 file')
   endif
      
   if ( lis%o%startcode == 2 ) then
      do t=1,lis%d%nch
         do l=1,4
            noah(t)%stc(l)=noahdrv%noah_it
         enddo
      enddo
   endif

   if ( lis%o%startcode == 3 ) then
#if ( defined USE_NETCDF )
      call lis_log_msg('MSG: noah_coldstart -- Reading initial soil temp: ' &
                       //trim(lis%p%soiltemp_init))
      status = nf90_open(path=trim(lis%p%soiltemp_init),mode= nf90_nowrite,&
                         ncid = ncid)
      status = nf90_inq_varid(ncid, "SoilTemp_init",stid)
      status = nf90_get_var(ncid, stid,soiltemp1)
      status = nf90_close(ncid)
      do c=1,lis%d%lnc
         do r=1,lis%d%lnr
            rindex = 150-r+1
            cindex = c
            if(gindex(cindex,rindex).ne.-1) then 
               soiltemp(cindex,rindex) = soiltemp1(gindex(cindex,rindex))
            endif
         enddo
      enddo
      

      do t=1,lis%d%nch
         if((soiltemp(tile(t)%col, tile(t)%row).ne.-9999.000)) then 
!            testtemp(t) = soiltemp(tile(t)%col, tile(t)%row)
            do l=1,4
               noah(t)%stc(l)= soiltemp(tile(t)%col, tile(t)%row)
            enddo
         endif
      enddo
!      ftn = 110
!      ier = nf90_create("soiltemp.nc",cmode=nf90_clobber,ncid=ftn)
!      ier = nf90_def_dim(ftn, 'land',lis%d%nch, dimID(1))
!      ier = nf90_def_dim(ftn, 'time',1, dimID(2))
!      ier = nf90_def_var(ftn,"soiltemp",nf90_float,dimids=dimID,varid=varid)
!      ier = nf90_put_att(ftn,varid,"units","K")
!      ier = nf90_enddef(ftn)
!      ier = nf90_put_var(ftn,varid, testtemp,(/1,1/),(/1,lis%d%nch/))
!      call drv_writevar_netcdf(ftn,testtemp,1,varid)
!      ier = nf90_close(ftn)
!      stop
#endif
   endif

   if ( lis%o%startcode == 2 .or. lis%o%startcode == 3 ) then
      print*,'DBG: noah_coldstart -- nch',lis%d%nch, &
           ' (', iam, ')'
      do t=1,lis%d%nch
         noah(t)%t1=noahdrv%noah_it
         noah(t)%t1=280.0         
         noah(t)%cmc=0.0004
         noah(t)%snowh=0.0
         noah(t)%sneqv=0.0
         noah(t)%ch=0.0150022404
         noah(t)%cm=0.0205970779
         noah(t)%smc(1)=0.3252287
         noah(t)%smc(2)=0.3194746
         noah(t)%smc(3)=0.3172167
         noah(t)%smc(4)=0.3078052
         noah(t)%sh2o(1)=0.1660042
         noah(t)%sh2o(2)=0.2828006
         noah(t)%sh2o(3)=0.3172163
         noah(t)%sh2o(4)=0.3078025
      enddo  
      lis%t%yr=lis%t%syr
      lis%t%mo=lis%t%smo 
      lis%t%da=lis%t%sda
      lis%t%hr=lis%t%shr
      lis%t%mn=lis%t%smn
      lis%t%ss=lis%t%sss

      call date2time(lis%t%time,lis%t%doy,lis%t%gmt,lis%t%yr,&
                     lis%t%mo,lis%t%da,lis%t%hr,lis%t%mn,lis%t%ss) 
      write(*,*)'MSG: noah_coldstart -- Using lis.crd start time ',&
                lis%t%time, ' (', iam, ')'
   endif
!EOC
end subroutine noah_coldstart
