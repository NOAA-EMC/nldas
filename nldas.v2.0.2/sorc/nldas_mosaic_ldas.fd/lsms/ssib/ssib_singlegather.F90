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
! !ROUTINE: ssib_singlegather.F90
!
! !DESCRIPTION:
!  Gather single variable for output 
!
! !REVISION HISTORY:
! 22 May 2004: David Mocko, Conversion from HY-SSiB to SSiB
!
! !INTERFACE:
      subroutine ssib_singlegather(index, var)
! !USES:
      use lisdrv_module, only : lis
      use tile_spmdMod
      use ssib_varder
      use ssibpardef_module

      implicit none
! !ARGUMENTS:
      integer :: index             ! Index of SSiB variable
      real    :: var(lis%d%glbnch) ! SSiB variable being gathered
!EOP
      real :: var_temp(di_array(iam))
      integer :: t
      integer ierr
!BOC
      do t = 1,di_array(iam)
         select case (index)
         case(1)
            var_temp(t) = ssib(t)%swnet/float(ssib(t)%count)
         case(2)
            var_temp(t) = ssib(t)%lwnet/float(ssib(t)%count)
         case(3)
            var_temp(t) = ssib(t)%qle/float(ssib(t)%count)
         case(4)
            var_temp(t) = ssib(t)%qh/float(ssib(t)%count)
         case(5)
            var_temp(t) = ssib(t)%qg/float(ssib(t)%count)
         case(6)
            var_temp(t) = ssib(t)%qf/float(ssib(t)%count)
         case(7)
            var_temp(t) = ssib(t)%qtau/float(ssib(t)%count)
         case(8)
            var_temp(t) = ssib(t)%delsurfheat
         case(9)
            var_temp(t) = ssib(t)%snowf/float(ssib(t)%count)
         case(10)
            var_temp(t) = ssib(t)%rainf/float(ssib(t)%count)
         case(11)
            var_temp(t) = ssib(t)%evap/float(ssib(t)%count)
         case(12)
            var_temp(t) = ssib(t)%qs/float(ssib(t)%count)
         case(13)
            var_temp(t) = ssib(t)%qsb/float(ssib(t)%count)
         case(14)
            var_temp(t) = ssib(t)%qsm/float(ssib(t)%count)
         case(15)
            var_temp(t) = ssib(t)%delsoilmoist
         case(16)
            var_temp(t) = ssib(t)%delswe
         case(17)
            var_temp(t) = ssib(t)%delintercept
         case(18)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%vegtc/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%vegtc
            endif
         case(19)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%baresoilt/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%baresoilt
            endif
         case(20)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%avgsurft/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%avgsurft
            endif
         case(21)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%radteff/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%radteff
            endif
         case(22)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               if (ssib(t)%albedocount.gt.0) then
                  var_temp(t) = ssib(t)%albedo/float(ssib(t)%albedocount)
               else
                  var_temp(t) = lis%d%udef
               endif
            else
               var_temp(t) = ssib(t)%albedo
            endif
         case(23)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%swe/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%swe
            endif
         case(24)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%sweveg/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%sweveg
            endif
         case(25)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%soilmoist1/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%soilmoist1
            endif
         case(26)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%soilmoist2/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%soilmoist2
            endif
         case(27)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%soilmoist3/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%soilmoist3
            endif
         case(28)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%soiltemp/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%soiltemp
            endif
         case(29)
            if (ssibdrv%STATEVAR_AVG.eq.1) then
               var_temp(t) = ssib(t)%soilwet/float(ssib(t)%count)
            else
               var_temp(t) = ssib(t)%soilwet
            endif
         case(30)
            var_temp(t) = ssib(t)%ecanop/float(ssib(t)%count)
         case(31)
            var_temp(t) = ssib(t)%tveg/float(ssib(t)%count)
         case(32)
            var_temp(t) = ssib(t)%esoil/float(ssib(t)%count)
         case(33)
            var_temp(t) = ssib(t)%rootmoist/float(ssib(t)%count)
         case(34)
            var_temp(t) = ssib(t)%canopint/float(ssib(t)%count)
         case(35)
            var_temp(t) = ssib(t)%acond/float(ssib(t)%count)
         case(36)
            var_temp(t) = ssib(t)%snowfrac/float(ssib(t)%count)
         case(37)
            if (lis%o%wfor.eq.1) then
               var_temp(t) = sqrt(ssib(t)%forcing(5)*ssib(t)%forcing(5)+ & 
               ssib(t)%forcing(6)*ssib(t)%forcing(6))
            endif
         case(38)
            if (lis%o%wfor.eq.1) then
               if (ssib(t)%forcing(1).lt.273.15) then
                  var_temp(t) = 0.0
               else 
                  var_temp(t) = ssib(t)%forcing(8)
               endif
            endif
         case(39)
            if (lis%o%wfor.eq.1) then
               if (ssib(t)%forcing(1).lt.273.15) then
                  var_temp(t) = ssib(t)%forcing(8)
               else
                  var_temp(t) = 0.0
               endif
            endif
         case(40)
            if (lis%o%wfor.eq.1) then
               var_temp(t) = ssib(t)%forcing(1)
            endif
         case(41)
            if (lis%o%wfor.eq.1) then
               var_temp(t) = ssib(t)%forcing(2)
            endif
         case(42)
            if (lis%o%wfor.eq.1) then
               var_temp(t) = ssib(t)%forcing(7)
            endif
         case(43)
            if (lis%o%wfor.eq.1) then
               var_temp(t) = ssib(t)%forcing(3)
            endif
         case(44)
            if (lis%o%wfor.eq.1) then
               var_temp(t) = ssib(t)%forcing(4)
            endif
         end select
      enddo
#if (defined SPMD)
      call MPI_GATHERV(var_temp(1:di_array(iam)),di_array(iam), & 
           MPI_REAL,var,di_array,displs,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
      return
!EOC
      end subroutine ssib_singlegather

