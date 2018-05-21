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
! !ROUTINE: template_writevars.F90
!
! !DESCRIPTION:  
!  LIS data writer: Writes template output
!
! !REVISION HISTORY:
! 21 Jul 2004; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine template_writestats(ld,ftn_stats)
! !USES:
  use lisdrv_module, only : gindex
  use lis_module
  use template_varder
 
  implicit none
  
  type(lisdec) :: ld
  integer :: ftn_stats,t
!EOP
  real :: vmean,vstdev,vmin,vmax
  real :: rainf(ld%d%glbnch)
  real :: snowf(ld%d%glbnch)
!BOC
  do t=1,ld%d%glbnch
     if(template(t)%forcing(1) < 273.15) then
        rainf(t) = 0.0
        snowf(t) = template(t)%forcing(8)
     else
        rainf(t) = template(t)%forcing(8)
        snowf(t) = 0.0
     endif
  enddo
   if(ld%o%wfor.eq.1) then
      call stats(sqrt(template%forcing(5)*template%forcing(5)+ & 
           template%forcing(6)*template%forcing(6)), & 
           ld%d%udef,ld%d%glbnch,vmean,vstdev,vmin, vmax)
      write(ftn_stats,999) 'Wind(m/s)', & 
           vmean,vstdev,vmin,vmax
      call stats(rainf, & 
           ld%d%udef,ld%d%glbnch,vmean,vstdev,vmin, vmax)
      write(ftn_stats,998) 'Rainf(kg/m2s)', & 
           vmean,vstdev,vmin,vmax
      call stats(snowf, & 
           ld%d%udef,ld%d%glbnch,vmean,vstdev,vmin, vmax)
      write(ftn_stats,998) 'Snowf(kg/m2s)', & 
           vmean,vstdev,vmin,vmax
      call stats(template%forcing(1),ld%d%udef,ld%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'Tair(K)', & 
           vmean,vstdev,vmin,vmax
      call stats(template%forcing(2),ld%d%udef,ld%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'Qair(kg/kg)', & 
           vmean,vstdev,vmin,vmax
      call stats(template%forcing(7),ld%d%udef,ld%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'PSurf(Pa)', & 
           vmean,vstdev,vmin,vmax
      call stats(template%forcing(3),ld%d%udef,ld%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'SWdown (W/m2)', & 
           vmean,vstdev,vmin,vmax
      call stats(template%forcing(4),ld%d%udef,ld%d%glbnch,vmean,vstdev, & 
           vmin, vmax)
      write(ftn_stats,999) 'LWdown(W/m2)', & 
           vmean,vstdev,vmin,vmax
   endif
998    FORMAT(1X,A18,4E14.3)
999    FORMAT(1X,A18,4F14.3)
!EOC
 end subroutine template_writestats
 
