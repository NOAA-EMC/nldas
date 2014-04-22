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
! !ROUTINE: vic_writevars.F90
!
! !DESCRIPTION:  
!  LIS VIC data writer: Writes vic output
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine vic_writestats(ftn_stats)
! !USES:
  use lisdrv_module, only : gindex, lis
  use vic_varder
 
  implicit none
  
  real :: tmp(lis%d%glbnch)
  integer :: ftn_stats,t, nvars
!EOP
  real :: vmean,vstdev,vmin,vmax
!BOC
  tmp = 0
  call get_vicvar(1, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'SWnet(W/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(2, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'LWnet(W/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(3, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'Qle(W/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(4, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'Qh(W/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(5, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'Qg(W/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(6, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'Rainf(kg/m2s)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(7, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'Snowf(kg/m2s)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(8, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'Evap(kg/m2s)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(9, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'Qs(kg/m2s)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(10, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'Qsb(kg/m2s)', & 
       vmean,vstdev,vmin,vmax
  
  call get_vicvar(11, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'Qfz(kg/m2s)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(12, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'SnowT(K)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(13, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'AvgSurfT(K)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(14, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'RadT(K)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(15, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'Albedo(-)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(16, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'SoilTemp1(K)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(17, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'SoilTemp2(K)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(18, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'SoilTemp3(K)', & 
       vmean,vstdev,vmin,vmax


  call get_vicvar(19, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'SoilMoist1(kg/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(20, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'SoilMoist2(kg/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(21, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'SoilMoist3(kg/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(22, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'TVeg', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(23, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'ESoil', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(24, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'SoilWet(-)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(25, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'RootMoist(kg/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(26, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,998) 'SWE(kg/m2)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(27, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'Qsm(kg/m2s)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(28, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'DelSoilMoist(kg/m2s)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(29, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'DelSWE(kg/m2s)', & 
       vmean,vstdev,vmin,vmax

  call get_vicvar(30, lis%d%glbnch, tmp)
  call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
       vmin, vmax)
  write(ftn_stats,999) 'ACond(m/s)', & 
       vmean,vstdev,vmin,vmax

  if(lis%o%wfor.eq.1) then 
     call get_vicvar(31, lis%d%glbnch, tmp)
     call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
          vmin, vmax)
     write(ftn_stats,999) 'Wind(m/s)', & 
          vmean,vstdev,vmin,vmax

     call get_vicvar(32, lis%d%glbnch, tmp)
     call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
          vmin, vmax)
     write(ftn_stats,999) 'Rainfall(kg/m2/s)', & 
          vmean,vstdev,vmin,vmax

     call get_vicvar(33, lis%d%glbnch, tmp)
     call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
          vmin, vmax)
     write(ftn_stats,999) 'Snowfall(kg/m2/s)', & 
          vmean,vstdev,vmin,vmax

     call get_vicvar(34, lis%d%glbnch, tmp)
     call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
          vmin, vmax)
     write(ftn_stats,999) 'Tair(K)', & 
          vmean,vstdev,vmin,vmax

     call get_vicvar(35, lis%d%glbnch, tmp)
     call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
          vmin, vmax)
     write(ftn_stats,999) 'Qair(kg/kg)', & 
          vmean,vstdev,vmin,vmax

     call get_vicvar(36, lis%d%glbnch, tmp)
     call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
          vmin, vmax)
     write(ftn_stats,999) 'Psurf(Pa)', & 
          vmean,vstdev,vmin,vmax

     call get_vicvar(37, lis%d%glbnch, tmp)
     call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
          vmin, vmax)
     write(ftn_stats,999) 'SWdown(W/m2)', & 
          vmean,vstdev,vmin,vmax

     call get_vicvar(38, lis%d%glbnch, tmp)
     call stats(tmp,lis%d%udef,lis%d%glbnch,vmean,vstdev, & 
          vmin, vmax)
     write(ftn_stats,999) 'LWdown(W/m2)', & 
          vmean,vstdev,vmin,vmax

  endif

998    FORMAT(1X,A18,4E14.3)
999    FORMAT(1X,A18,4F14.3)
!EOC
end subroutine vic_writestats
 
