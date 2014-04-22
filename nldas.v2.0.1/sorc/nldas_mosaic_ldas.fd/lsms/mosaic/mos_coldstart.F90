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
! !ROUTINE: mos_coldstart.F90
!
! !DESCRIPTION:
!  
! Routine for mosaic initialization from cold start
! 
! !INTERFACE:
subroutine mos_coldstart()
! !USES:
   use lisdrv_module, only : lis
   use mos_varder
   use time_manager
   use spmdMod, only : iam
!EOP
   implicit none

   integer :: t,l
!BOC
   if ( lis%o%startcode == 2 ) then
      print*,'MSG: mos_coldstart -- cold-starting mosaic', &
             '...using ics from card file',' (', iam, ')'
      
      print*,'DBG: mos_coldstart -- nch',lis%d%nch, &
           ' (', iam, ')'
      do t=1,lis%d%nch
         mos(t)%ct=mosdrv%mos_it
         mos(t)%qa=0.0
         mos(t)%ics=0.0
         mos(t)%snow=0.0
         mos(t)%SoT=mosdrv%mos_it
         do l=1,3
            mos(t)%SoWET(l)=mosdrv%mos_ism
         enddo
      enddo  
      lis%t%yr=lis%t%syr
      lis%t%mo=lis%t%smo 
      lis%t%da=lis%t%sda
      lis%t%hr=lis%t%shr
      lis%t%mn=lis%t%smn
      lis%t%ss=lis%t%sss

      call date2time(lis%t%time,lis%t%doy,lis%t%gmt,lis%t%yr,&
                     lis%t%mo,lis%t%da,lis%t%hr,lis%t%mn,lis%t%ss) 
      write(*,*)'MSG: mos_coldstart -- Using lis.crd start time ',&
                lis%t%time, ' (', iam, ')'
	print *,'mos(1)%sot=',mos(1)%sot
   endif
!EOC
end subroutine mos_coldstart
