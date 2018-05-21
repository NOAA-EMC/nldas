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
subroutine noah_binout(ftn)
! !USES:
  use lisdrv_module, only : lis
  use drv_output_mod, only : drv_writevar_bin
  use noah_varder
 
  implicit none
  integer :: ftn
!EOP
  real :: rainf(lis%d%glbnch)
  real :: snowf(lis%d%glbnch)
  integer :: t,c,r
!BOC
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
   noah%swnet = noah%swnet/float(noah%count)
   call drv_writevar_bin(ftn,noah%swnet) !Net shortwave radiation (surface) (W/m2)
   noah%lwnet = (-1)*noah%lwnet/float(noah%count)
   call drv_writevar_bin(ftn,noah%lwnet)!Net longwave radiation (surface) (W/m2)   
   noah%qle = noah%qle/float(noah%count)
   call drv_writevar_bin(ftn,noah%qle) !Latent Heat Flux (W/m2)
   noah%qh = noah%qh/float(noah%count)
   call drv_writevar_bin(ftn,noah%qh) !Sensible Heat Flux (W/m2)
   noah%qg = noah%qg/float(noah%count)
   call drv_writevar_bin(ftn,noah%qg)
!---------------------------------------------------------------------------
! General Water Balance Components
!---------------------------------------------------------------------------
   noah%snowf = noah%snowf/float(noah%count)
   call drv_writevar_bin(ftn,noah%snowf)
   noah%rainf = noah%rainf/float(noah%count)
   call drv_writevar_bin(ftn,noah%rainf)
   noah%evap = noah%evap/float(noah%count)
   call drv_writevar_bin(ftn,noah%evap)
   noah%qs = noah%qs/float(noah%count)
   call drv_writevar_bin(ftn,noah%qs)
   noah%qsb = noah%qsb/float(noah%count)
   call drv_writevar_bin(ftn,noah%qsb)
   noah%qsm = noah%qsm/float(noah%count)
   call drv_writevar_bin(ftn,noah%qsm)
   call drv_writevar_bin(ftn,(noah%smc(1)*1000.0*0.1+ &
        noah%smc(2)*1000.0*0.3 + & 
        noah%smc(3)*1000.0*0.6 + & 
        noah%smc(4)*1000.0*1.0 -noah%soilm_prev)/float(noah%count))
   call drv_writevar_bin(ftn,(noah%sneqv*1000.0-noah%swe_prev)/float(noah%count))
!---------------------------------------------------------------------------
! Surface State Variables
!---------------------------------------------------------------------------
   call drv_writevar_bin(ftn,noah%avgsurft)
   call drv_writevar_bin(ftn,noah%albedo)
   call drv_writevar_bin(ftn,noah%swe)
!---------------------------------------------------------------------------
! Subsurface State Variables
!---------------------------------------------------------------------------
   call drv_writevar_bin(ftn,noah%stc(1))
   call drv_writevar_bin(ftn,noah%stc(2))
   call drv_writevar_bin(ftn,noah%stc(3))
   call drv_writevar_bin(ftn,noah%stc(4))
   call drv_writevar_bin(ftn,noah%soilmoist1)
   call drv_writevar_bin(ftn,noah%soilmoist2)
   call drv_writevar_bin(ftn,noah%soilmoist3)
   call drv_writevar_bin(ftn,noah%soilmoist4)
   call drv_writevar_bin(ftn,noah%soilwet)
!---------------------------------------------------------------------------
! Evaporation Components
!---------------------------------------------------------------------------
   noah%ecanop = noah%ecanop/float(noah%count)
   call drv_writevar_bin(ftn, noah%ecanop)
   noah%tveg= noah%tveg/float(noah%count)
   call drv_writevar_bin(ftn,noah%tveg)
   noah%esoil= noah%esoil/float(noah%count)
   call drv_writevar_bin(ftn,noah%esoil)
   call drv_writevar_bin(ftn, noah%rootmoist)
   call drv_writevar_bin(ftn, noah%canopint)
!---------------------------------------------------------------------------
! Forcing Components
!---------------------------------------------------------------------------
   if(lis%o%wfor.eq.1) then
      call drv_writevar_bin(ftn, sqrt(noah%forcing(5)*noah%forcing(5)+ & 
           noah%forcing(6)*noah%forcing(6)))
      call drv_writevar_bin(ftn,rainf)
      call drv_writevar_bin(ftn,snowf)
      call drv_writevar_bin(ftn,noah%forcing(1))
      call drv_writevar_bin(ftn,noah%forcing(2))
      call drv_writevar_bin(ftn,noah%forcing(7))
      call drv_writevar_bin(ftn,noah%forcing(3))
      call drv_writevar_bin(ftn,noah%forcing(4))
   endif

!EOC
 end subroutine noah_binout
 
