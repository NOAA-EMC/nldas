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
! !ROUTINE: noah_binout.F90
!
! !DESCRIPTION:  
!  LIS NOAH data writer: Writes noah output in binary format
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine noah_netcdfout(check, ftn)
! !USES:
#if ( defined USE_NETCDF )
  use netcdf
#endif
  use lisdrv_module, only : lis
  use drv_output_mod, only : drv_writevar_netcdf, drv_writevar_netcdf3d
  use noah_varder
 
  implicit none
  
  integer :: ftn, dim1
  logical :: check 
!EOP
  real :: rainf(lis%d%glbnch)
  real :: snowf(lis%d%glbnch)
  integer :: t,c,r,dimID(2), days1(12),days2(12)
  integer :: dimID1(3)
  data days1 /31,28,31,30,31,30,31,31,30,31,30,31/
  data days2 /31,29,31,30,31,30,31,31,30,31,30,31/
  logical :: leap
!BOC
#if ( defined USE_NETCDF )
  if(.not. check) then 
     if((mod(lis%t%yr,4).eq.0.and.mod(lis%t%yr,100).ne.0) &    
          .or.(mod(lis%t%yr,400).eq.0)) then            
        leap = .true.                   
     else
        leap = .false.
     endif
     if(lis%t%mo .eq. lis%t%emo .and. &
          lis%t%yr .eq. lis%t%eyr) then
        dim1 = ((lis%t%eda-lis%t%da -1)*24+&
             lis%t%ehr+(24-lis%t%hr))/noahdrv%writeintn
     else
        if(leap) then 
           dim1 = ((24-lis%t%hr)+&
                (days2(lis%t%mo)-lis%t%da)*24)/noahdrv%writeintn
        else
           dim1 = ((24-lis%t%hr)+&
                (days1(lis%t%mo)-lis%t%da)*24)/noahdrv%writeintn
        endif
     endif

     call check_nc(nf90_def_dim(ftn, 'land',lis%d%glbnch, dimID(1)))
     call check_nc(nf90_def_dim(ftn, 'time',dim1, dimID(2)))
     
     call check_nc(nf90_def_dim(ftn, 'land1',lis%d%glbnch, dimID1(1)))
     call check_nc(nf90_def_dim(ftn, 'soil',4,dimID1(2)))
     call check_nc(nf90_def_dim(ftn, 'time1',dim1, dimID1(3)))

     call check_nc(nf90_def_var(ftn, "Swnet", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(1)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(1),"units","W/m2"))
     
     call check_nc(nf90_def_var(ftn, "Lwnet", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(2)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(2),"units","W/m2"))
     
     call check_nc(nf90_def_var(ftn, "Qle", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(3)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(3),"units","W/m2"))
     
     call check_nc(nf90_def_var(ftn, "Qh", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(4)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(4),"units","W/m2"))
     
     call check_nc(nf90_def_var(ftn, "Qg", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(5)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(5),"units","W/m2"))
     
     call check_nc(nf90_def_var(ftn, "Snowf", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(6)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(6),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "Rainf", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(7)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(7),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "Evap", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(8)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(8),"units","kg/m2s"))
     
     
     call check_nc(nf90_def_var(ftn, "Qs", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(9)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(9),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "Qsb", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(10)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(10),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "Qsm", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(11)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(11),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "DelSoilMoist", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(12)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(12),"units","kg/m2s"))
     
     
     call check_nc(nf90_def_var(ftn, "DelSWE", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(13)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(13),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "AvgSurfT", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(14)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(14),"units","K"))
     
     call check_nc(nf90_def_var(ftn, "Albedo", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(15)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(15),"units","-"))
     
     call check_nc(nf90_def_var(ftn, "SWE", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(16)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(16),"units","kg/m2"))
     
     call check_nc(nf90_def_var(ftn, "SoilTemp", nf90_float, &
          dimids = dimID1, varID = noahdrv%varid(17)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(17),"units","K"))

     call check_nc(nf90_def_var(ftn, "SoilMoist", nf90_float, &
          dimids = dimID1, varID = noahdrv%varid(18)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(18),"units","kg/m2"))

       
     call check_nc(nf90_def_var(ftn, "SoilWet", nf90_float, &
       dimids = dimID, varID = noahdrv%varid(19)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(19),"units","-"))

     call check_nc(nf90_def_var(ftn, "ECanop", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(20)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(20),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "TVeg", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(21)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(21),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "ESoil", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(22)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(22),"units","kg/m2s"))
     
     call check_nc(nf90_def_var(ftn, "RootMoist", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(23)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(23),"units","kg/m2"))

     call check_nc(nf90_def_var(ftn, "CanopInt", nf90_float, &
          dimids = dimID, varID = noahdrv%varid(24)))
     call check_nc(nf90_put_att(ftn,noahdrv%varid(24),"units","kg/m2"))
     call check_nc(nf90_enddef(ftn))
  endif
  dim1 = noahdrv%numout

!---------------------------------------------------------------------------
! General Energy Balance Components
!---------------------------------------------------------------------------
  noah%swnet = noah%swnet/float(noah%count)
  call drv_writevar_netcdf(ftn, noah%swnet, dim1, noahdrv%varid(1))
  noah%lwnet = (-1)*noah%lwnet/float(noah%count)
  call drv_writevar_netcdf(ftn, noah%lwnet, dim1, noahdrv%varid(2))
  noah%qle = noah%qle/float(noah%count)
  call drv_writevar_netcdf(ftn, noah%qle, dim1, noahdrv%varid(3))
  noah%qh = noah%qh/float(noah%count)
  call drv_writevar_netcdf(ftn, noah%qh, dim1, noahdrv%varid(4))
  noah%qg = noah%qg/float(noah%count)
  call drv_writevar_netcdf(ftn, noah%qg, dim1, noahdrv%varid(5))

!---------------------------------------------------------------------------
! General Water Balance Components
!---------------------------------------------------------------------------
   noah%snowf = noah%snowf/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%snowf, dim1, noahdrv%varid(6))
   noah%rainf = noah%rainf/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%rainf, dim1, noahdrv%varid(7))
   noah%evap = noah%evap/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%evap, dim1, noahdrv%varid(8))
   noah%qs = noah%qs/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%qs, dim1, noahdrv%varid(9))
   noah%qsb = noah%qsb/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%qsb, dim1, noahdrv%varid(10))
   noah%qsm = noah%qsm/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%qsm, dim1, noahdrv%varid(11))
   call drv_writevar_netcdf(ftn,  noah%smc(1)*1000.0*0.1 + &
                                  noah%smc(2)*1000.0*0.3 + & 
                                  noah%smc(3)*1000.0*0.6 + & 
                                  noah%smc(4)*1000.0 -     &
                                  noah%soilm_prev,dim1, noahdrv%varid(12))
   call drv_writevar_netcdf(ftn, noah%sneqv*1000.0 - noah%swe_prev,&
        dim1, noahdrv%varid(13))

!---------------------------------------------------------------------------
! Surface State Variables
!---------------------------------------------------------------------------
   call drv_writevar_netcdf(ftn, noah%avgsurft ,dim1, noahdrv%varid(14))
   call drv_writevar_netcdf(ftn, noah%albedo ,dim1, noahdrv%varid(15))
   call drv_writevar_netcdf(ftn, noah%swe ,dim1, noahdrv%varid(16))

!---------------------------------------------------------------------------
! Subsurface State Variables
!---------------------------------------------------------------------------
   call drv_writevar_netcdf3d(ftn, noah%stc(1), dim1, 1, noahdrv%varid(17))
   call drv_writevar_netcdf3d(ftn, noah%stc(2), dim1, 2, noahdrv%varid(17))
   call drv_writevar_netcdf3d(ftn, noah%stc(3), dim1, 3, noahdrv%varid(17))
   call drv_writevar_netcdf3d(ftn, noah%stc(4), dim1, 4, noahdrv%varid(17))

   call drv_writevar_netcdf3d(ftn, noah%soilmoist1, dim1, 1, noahdrv%varid(18))
   call drv_writevar_netcdf3d(ftn, noah%soilmoist2, dim1, 2, noahdrv%varid(18))
   call drv_writevar_netcdf3d(ftn, noah%soilmoist3, dim1, 3, noahdrv%varid(18))
   call drv_writevar_netcdf3d(ftn, noah%soilmoist4, dim1, 4, noahdrv%varid(18))

   call drv_writevar_netcdf(ftn, noah%soilwet ,dim1, noahdrv%varid(19))

!---------------------------------------------------------------------------
! Evaporation Components
!---------------------------------------------------------------------------
   noah%ecanop = noah%ecanop/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%ecanop, dim1, noahdrv%varid(20))
   noah%tveg= noah%tveg/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%tveg ,dim1, noahdrv%varid(21))
   noah%esoil= noah%esoil/float(noah%count)
   call drv_writevar_netcdf(ftn, noah%esoil, dim1, noahdrv%varid(22))
   call drv_writevar_netcdf(ftn, noah%rootmoist, dim1, noahdrv%varid(23))
   call drv_writevar_netcdf(ftn, noah%canopint, dim1, noahdrv%varid(24))
!EOC
#endif
 end subroutine noah_netcdfout

