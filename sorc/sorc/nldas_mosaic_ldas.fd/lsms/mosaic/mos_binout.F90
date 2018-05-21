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
!  LIS MOS data writer: Writes mos output in binary format
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine mos_binout(ftn)
! !USES:
  use lisdrv_module, only : lis
  use drv_output_mod, only : drv_writevar_bin
  use mos_varder
 
  implicit none
  
  integer :: ftn
!EOP
  real :: rainf(lis%d%glbnch)
  real :: snowf(lis%d%glbnch)
  integer :: t,c,r
!BOC
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
   mos%swnet = mos%swnet/float(mos%count)
   call drv_writevar_bin(ftn,mos%swnet) !Net shortwave radiation (surface) (W/m2)
   mos%lwnet = (-1)*mos%lwnet/float(mos%count)
   call drv_writevar_bin(ftn,mos%lwnet)!Net longwave radiation (surface) (W/m2)
   mos%qle = mos%qle/float(mos%count)
   call drv_writevar_bin(ftn,mos%qle) !Latent Heat Flux (W/m2)
   mos%qh = mos%qh/float(mos%count)
   call drv_writevar_bin(ftn,mos%qh) !Sensible Heat Flux (W/m2)
   mos%qg = mos%qg/float(mos%count)
   call drv_writevar_bin(ftn,mos%qg)
!---------------------------------------------------------------------------
! General Water Balance Components
!---------------------------------------------------------------------------
   mos%snowf = mos%snowf/float(mos%count)
   call drv_writevar_bin(ftn,mos%snowf)
   mos%rainf = mos%rainf/float(mos%count)
   call drv_writevar_bin(ftn,mos%rainf)
   mos%evap = mos%evap/float(mos%count)
   call drv_writevar_bin(ftn,mos%evap)
   mos%qs = mos%qs/float(mos%count)
   call drv_writevar_bin(ftn,mos%qs)
   mos%qsb = mos%qsb/float(mos%count)
   call drv_writevar_bin(ftn,mos%qsb)
   mos%qsm = mos%qsm/float(mos%count)
   call drv_writevar_bin(ftn,mos%qsm)
   call drv_writevar_bin(ftn,(mos%water1 + & 
        mos%water2 + & 
        mos%water3  & 
        -mos%soilm_prev)/float(mos%count))
   call drv_writevar_bin(ftn,(mos%snow-mos%swe_prev)/float(mos%count))
!---------------------------------------------------------------------------
! Surface State Variables
!---------------------------------------------------------------------------
   call drv_writevar_bin(ftn,mos%avgsurft)
   call drv_writevar_bin(ftn,mos%soT)
   call drv_writevar_bin(ftn,mos%albedo)
   mos%swe= mos%swe/float(mos%count)
   call drv_writevar_bin(ftn,mos%swe)
!---------------------------------------------------------------------------
! Subsurface State Variables
!---------------------------------------------------------------------------
   mos%soilmoist1= mos%soilmoist1/float(mos%count)
   call drv_writevar_bin(ftn,mos%soilmoist1)
   mos%soilmoist2= mos%soilmoist2/float(mos%count)
   call drv_writevar_bin(ftn,mos%soilmoist2)
   mos%soilmoist3= mos%soilmoist3/float(mos%count)
   call drv_writevar_bin(ftn,mos%soilmoist3)
   mos%soilwet= mos%soilwet/float(mos%count)
   call drv_writevar_bin(ftn,mos%soilwet)
!---------------------------------------------------------------------------
! Evaporation Components
!---------------------------------------------------------------------------
   mos%ecanop = mos%ecanop/float(mos%count)
   call drv_writevar_bin(ftn,mos%ecanop)
   mos%tveg= mos%tveg/float(mos%count)
   call drv_writevar_bin(ftn,mos%tveg)
   mos%esoil= mos%esoil/float(mos%count)
   call drv_writevar_bin(ftn,mos%esoil)
   mos%rootmoist = mos%rootmoist/float(mos%count)
   call drv_writevar_bin(ftn, mos%rootmoist)
   mos%canopint = mos%canopint/float(mos%count)
   call drv_writevar_bin(ftn, mos%canopint)
   mos%acond = mos%acond
   call drv_writevar_bin(ftn, mos%acond)
   if(lis%o%wfor.eq.1) then
      call drv_writevar_bin(ftn, sqrt(mos%forcing(5)*mos%forcing(5)+ & 
           mos%forcing(6)*mos%forcing(6)))
      call drv_writevar_bin(ftn,rainf)
      call drv_writevar_bin(ftn,snowf)
      call drv_writevar_bin(ftn,mos%forcing(1))
      call drv_writevar_bin(ftn,mos%forcing(2))
      call drv_writevar_bin(ftn,mos%forcing(7))
      call drv_writevar_bin(ftn,mos%forcing(3))
      call drv_writevar_bin(ftn,mos%forcing(4))
   endif
!EOC
 end subroutine mos_binout
 
