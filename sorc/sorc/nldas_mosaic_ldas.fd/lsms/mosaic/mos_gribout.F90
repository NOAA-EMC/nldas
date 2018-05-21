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
! !ROUTINE: mos_binout.F90
!
! !DESCRIPTION:  
!  LIS NOAH data writer: Writes mos output in binary format
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 17 Oct 2013; Youlong Xia; Revised for grib2 output 
! !INTERFACE:
subroutine mos_gribout(ftn)
! !USES:
  use lisdrv_module, only : lis, gindex
  use drv_output_mod, only : drv_writevar_grib2
  use mos_varder
  use time_manager, only : tick
 
  implicit none
  
  integer :: ftn
!EOP
  real :: vmean,vstdev,vmin,vmax
  real :: rainf(lis%d%glbnch)
  real :: snowf(lis%d%glbnch)
  integer :: t,c,r,i,k
  logical*1 :: lismask(lis%d%lnc,lis%d%lnr)
  character*8 :: today, yesterday
  character*1 :: tod(8), yes(8)
  character(len=100) :: temp1
  real*8 :: dummytime   
  real  :: dummygmt
  integer:: ss1,ts,mn1,hr1,da1,mo1,yr1,ts1,doy1
  integer :: kpds(13)
  real :: interval
!BOC
  interval = mosdrv%writeintm
  hr1=lis%t%hr
  da1=lis%t%da
  mo1=lis%t%mo
  yr1=lis%t%yr
  mn1=lis%t%mn
  ss1=0
  ts1=-3600*24
  dummygmt=1.0
  dummytime=1.0
  write(unit=temp1,fmt='(i4,i2,i2)')yr1,mo1,da1
  read(unit=temp1,fmt='(8a1)')tod
  do i=1,8
     if(tod(i).eq.(' '))tod(i)='0'
  enddo
  today=tod(1)//tod(2)//tod(3)//tod(4)//tod(5) &
       //tod(6)//tod(7)//tod(8)
  
  call tick(dummytime,doy1,dummygmt,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  write(unit=temp1,fmt='(i4,i2,i2)')yr1,mo1,da1
  read(unit=temp1,fmt='(8a1)')yes
  do i=1,8
     if(yes(i).eq.(' '))yes(i)='0'
  enddo
  yesterday=yes(1)//yes(2)//yes(3)//yes(4)//yes(5) &
       //yes(6)//yes(7)//yes(8)
  do i=1,13
     kpds(i)=0
  enddo

  open (unit = 69, file = 'KPDS_completemos.tbl')
  do k = 1, 42
     read(69,*)
  end do
  
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        if(gindex(c,r).gt.0) then
           lismask(c,r)=.true.
        else
           lismask(c,r)=.false.
        endif
     enddo
  enddo

  do t=1,lis%d%glbnch
     if(mos(t)%forcing(1) < 273.15) then
        rainf(t) = 0.0
        snowf(t) = mos(t)%forcing(8)
     else
        rainf(t) = mos(t)%forcing(8)
        snowf(t) = 0.0
     endif
  enddo
!---------------------------------------------------------------------------
! General Energy Balance Components
!---------------------------------------------------------------------------
   call mos_readkpds(69,kpds)
   mos%swnet = mos%swnet/float(mos%count)
   call drv_writevar_grib2(ftn,mos%swnet,kpds,lismask,interval,today,yesterday) 

   call mos_readkpds(69,kpds)   
   mos%lwnet = (-1)*mos%lwnet/float(mos%count)
   call drv_writevar_grib2(ftn,mos%lwnet,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%qle = mos%qle/float(mos%count)
   call drv_writevar_grib2(ftn,mos%qle,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%qh = mos%qh/float(mos%count)
   call drv_writevar_grib2(ftn,mos%qh,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%qg = mos%qg/float(mos%count)
   call drv_writevar_grib2(ftn,mos%qg,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
   mos%snohf = mos%snohf/float(mos%count)
   call drv_writevar_grib2(ftn,mos%snohf,kpds,lismask,interval,today,yesterday)

!---------------------------------------------------------------------------
! General Water Balance Components
!---------------------------------------------------------------------------
   call mos_readkpds(69,kpds)   
!   mos%snowf = mos%snowf
   snowf=snowf*lis%t%ts*mos%count
   call drv_writevar_grib2(ftn,snowf,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
!   mos%rainf = mos%rainf
   rainf=rainf*lis%t%ts*mos%count
   call drv_writevar_grib2(ftn,rainf,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%evap = mos%evap
   call drv_writevar_grib2(ftn,mos%evap,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%qs = mos%qs
   call drv_writevar_grib2(ftn,mos%qs,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%qsb = mos%qsb
   call drv_writevar_grib2(ftn,mos%qsb,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%qsm = mos%qsm
   call drv_writevar_grib2(ftn,mos%qsm,kpds,lismask,interval,today,yesterday)

!---------------------------------------------------------------------------
! Surface State Variables
!---------------------------------------------------------------------------
   call mos_readkpds(69,kpds)   
   call drv_writevar_grib2(ftn,mos%avgsurft,kpds,lismask,interval,today,yesterday)
   call mos_readkpds(69,kpds)   
   call drv_writevar_grib2(ftn,mos%albedo,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%swe,kpds,lismask,interval,today,yesterday)
!---------------------------------------------------------------------------
! Subsurface State Variables
!---------------------------------------------------------------------------
   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%sot,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds) 
   call drv_writevar_grib2(ftn,mos%soilmtot,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%soilmr,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds) 
   call drv_writevar_grib2(ftn,mos%soilm1,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%soilmoist1,kpds,lismask,interval,today,yesterday)
   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%soilmoist2,kpds,lismask,interval,today,yesterday)
   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%soilmoist3,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%soilwet=mos%soilwet*100.0
   call drv_writevar_grib2(ftn,mos%soilwet,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
   mos%mstavr=mos%mstavr*100.0
   call drv_writevar_grib2(ftn,mos%mstavr,kpds,lismask,interval,today,yesterday)

!---------------------------------------------------------------------------
! Evaporation Components
!---------------------------------------------------------------------------
   call mos_readkpds(69,kpds)   
   mos%ecanop = mos%ecanop/float(mos%count)
   call drv_writevar_grib2(ftn,mos%ecanop,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%tveg= mos%tveg/float(mos%count)
   call drv_writevar_grib2(ftn,mos%tveg,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   mos%esoil= mos%esoil/float(mos%count)
   call drv_writevar_grib2(ftn,mos%esoil,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
   mos%sbsno= mos%sbsno/float(mos%count)
   call drv_writevar_grib2(ftn,mos%sbsno,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%acond,kpds,lismask,interval,today,yesterday)
   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%ccond,kpds,lismask,interval,today,yesterday)
   call mos_readkpds(69,kpds)
   mos%green=mos%green*100.0
   call drv_writevar_grib2(ftn,mos%green,kpds,lismask,interval,today,yesterday)
   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn,mos%laiout,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)   
   call drv_writevar_grib2(ftn, mos%canopint,kpds,lismask,interval,today,yesterday)

   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn, mos%snod,kpds,lismask,interval,today,yesterday)
   call mos_readkpds(69,kpds)
   call drv_writevar_grib2(ftn, mos%snwfrcout,kpds,lismask,interval,today,yesterday)


!---------------------------------------------------------------------------
! Forcing Components
!---------------------------------------------------------------------------
   if(lis%o%wfor.eq.1) then
      call mos_readkpds(69,kpds)   
      call drv_writevar_grib2(ftn, mos%forcing(5), & 
           kpds,lismask,interval,today,yesterday)

      call mos_readkpds(69,kpds)
      call drv_writevar_grib2(ftn, mos%forcing(6), &                    
           kpds,lismask,interval,today,yesterday)
      
      call mos_readkpds(69,kpds)
      call drv_writevar_grib2(ftn,mos%forcing(1),kpds,lismask,interval,today,yesterday)

      call mos_readkpds(69,kpds)
      call drv_writevar_grib2(ftn,mos%forcing(2),kpds,lismask,interval,today,yesterday)

      call mos_readkpds(69,kpds)
      call drv_writevar_grib2(ftn,mos%forcing(7)/100.0,kpds,lismask,interval,today,yesterday)

      call mos_readkpds(69,kpds)   
      mos%swrad=mos%swrad/float(mos%count)
      call drv_writevar_grib2(ftn,mos%swrad,kpds,lismask,interval,today,yesterday)

      call mos_readkpds(69,kpds)   
      mos%lwrad=mos%lwrad/float(mos%count)
      call drv_writevar_grib2(ftn,mos%lwrad,kpds,lismask,interval,today,yesterday)

      call mos_readkpds(69,kpds)
      call drv_writevar_grib2(ftn,mos%forcing(9)*lis%t%ts*mos%count,kpds,lismask,interval,today,yesterday)

   endif
   close(69)

!EOC
 end subroutine mos_gribout
 
