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
! !ROUTINE: noah_binout.F90
!
! !DESCRIPTION:  
!  LIS NOAH data writer: Writes noah output in binary format
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine noah_gribout(ftn)
! !USES:
  use lisdrv_module, only : lis, gindex
  use drv_output_mod, only : drv_writevar_grib
  use noah_varder
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
  integer :: kpds(25)
  real :: interval
!BOC
  interval = noahdrv%writeintn
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
  do i=1,25
     kpds(i)=0
  enddo
  kpds(1)=221               !id for gsfc products
  kpds(2)=221               !id for noah model (change value for other models)
  kpds(4)=192               !bms flag... don't worry about this.
  kpds(12)=0                !assume output time minute always = 0
  kpds(13)=1                !forecast time unit (hours)
  kpds(17)=int((noahdrv%writeintn*3600.0)/lis%t%ts) !number of time steps in
  !averaged/accum variables
  kpds(18)=0                !grib version -- left as 0 in ncep products
  kpds(19)=1                !version number of kpds.tbl for lisas.  
  kpds(20)=0                !none missing from averages/accumulations (always4)
  kpds(23)=221              !gsfc id#
  kpds(24)=0                !does not apply to lisas output
  kpds(25)=0   

  open (unit = 69, file = './src/tables/KPDS_completenoah.tbl')
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
     if(noah(t)%forcing(1) < 273.15) then
        rainf(t) = 0.0
        snowf(t) = noah(t)%forcing(8)
     else
        rainf(t) = noah(t)%forcing(8)
        snowf(t) = 0.0
     endif
  enddo
!---------------------------------------------------------------------------
! General Energy Balance Components
!---------------------------------------------------------------------------
   call readkpds(69,kpds)
   noah%swnet = noah%swnet/float(noah%count)
   call drv_writevar_grib(ftn,noah%swnet,kpds,lismask,interval,today,yesterday) 

   call readkpds(69,kpds)   
   noah%lwnet = (-1)*noah%lwnet/float(noah%count)
   call drv_writevar_grib(ftn,noah%lwnet,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%qle = noah%qle/float(noah%count)
   call drv_writevar_grib(ftn,noah%qle,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%qh = noah%qh/float(noah%count)
   call drv_writevar_grib(ftn,noah%qh,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%qg = noah%qg/float(noah%count)
   call drv_writevar_grib(ftn,noah%qg,kpds,lismask,interval,today,yesterday)
!---------------------------------------------------------------------------
! General Water Balance Components
!---------------------------------------------------------------------------
   call readkpds(69,kpds)   
   noah%snowf = noah%snowf/float(noah%count)
   call drv_writevar_grib(ftn,noah%snowf,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%rainf = noah%rainf/float(noah%count)
   call drv_writevar_grib(ftn,noah%rainf,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%evap = noah%evap/float(noah%count)
   call drv_writevar_grib(ftn,noah%evap,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%qs = noah%qs/float(noah%count)
   call drv_writevar_grib(ftn,noah%qs,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%qsb = noah%qsb/float(noah%count)
   call drv_writevar_grib(ftn,noah%qsb,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%qsm = noah%qsm/float(noah%count)
   call drv_writevar_grib(ftn,noah%qsm,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,(noah%smc(1)*1000.0*0.1+ &
        noah%smc(2)*1000.0*0.3 + & 
        noah%smc(3)*1000.0*0.6 + & 
        noah%smc(4)*1000.0*1.0 -noah%soilm_prev)/float(noah%count),kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,(noah%sneqv*1000.0-noah%swe_prev)/float(noah%count),kpds,lismask,interval,today,yesterday)
!---------------------------------------------------------------------------
! Surface State Variables
!---------------------------------------------------------------------------
   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%avgsurft,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%albedo,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%swe,kpds,lismask,interval,today,yesterday)
!---------------------------------------------------------------------------
! Subsurface State Variables
!---------------------------------------------------------------------------
   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%stc(1),kpds,lismask,interval,today,yesterday)
   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%stc(2),kpds,lismask,interval,today,yesterday)
   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%stc(3),kpds,lismask,interval,today,yesterday)
   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%stc(4),kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%soilmoist1,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%soilmoist2,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%soilmoist3,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%soilmoist4,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn,noah%soilwet,kpds,lismask,interval,today,yesterday)
!---------------------------------------------------------------------------
! Evaporation Components
!---------------------------------------------------------------------------
   call readkpds(69,kpds)   
   noah%ecanop = noah%ecanop/float(noah%count)
   call drv_writevar_grib(ftn,noah%ecanop,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%tveg= noah%tveg/float(noah%count)
   call drv_writevar_grib(ftn,noah%tveg,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   noah%esoil= noah%esoil/float(noah%count)
   call drv_writevar_grib(ftn,noah%esoil,kpds,lismask,interval,today,yesterday)

   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn, noah%rootmoist,kpds,lismask,interval,today,yesterday)
   call readkpds(69,kpds)   
   call drv_writevar_grib(ftn, noah%canopint,kpds,lismask,interval,today,yesterday)

!---------------------------------------------------------------------------
! Forcing Components
!---------------------------------------------------------------------------
   if(lis%o%wfor.eq.1) then
      call readkpds(69,kpds)   
      call drv_writevar_grib(ftn, sqrt(noah%forcing(5)*noah%forcing(5)+ & 
           noah%forcing(6)*noah%forcing(6)),kpds,lismask,interval,today,yesterday)

      call readkpds(69,kpds)   
      call drv_writevar_grib(ftn,rainf,kpds,lismask,interval,today,yesterday)

      call readkpds(69,kpds)   
      call drv_writevar_grib(ftn,snowf,kpds,lismask,interval,today,yesterday)

      call readkpds(69,kpds)   
      call drv_writevar_grib(ftn,noah%forcing(1),kpds,lismask,interval,today,yesterday)

      call readkpds(69,kpds)   
      call drv_writevar_grib(ftn,noah%forcing(2),kpds,lismask,interval,today,yesterday)

      call readkpds(69,kpds)   
      call drv_writevar_grib(ftn,noah%forcing(7),kpds,lismask,interval,today,yesterday)

      call readkpds(69,kpds)   
      call drv_writevar_grib(ftn,noah%forcing(3),kpds,lismask,interval,today,yesterday)

      call readkpds(69,kpds)   
      call drv_writevar_grib(ftn,noah%forcing(4),kpds,lismask,interval,today,yesterday)
   endif
   close(69)

!EOC
 end subroutine noah_gribout
 
