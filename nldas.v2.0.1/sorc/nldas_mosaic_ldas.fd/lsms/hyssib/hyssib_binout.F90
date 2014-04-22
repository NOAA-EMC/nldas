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
! !ROUTINE: hyssib_gridout.F90
!
! !DESCRIPTION:
!  LIS HY-SSiB data writer: Writes HY-SSiB output in grid space
!
! !REVISION HISTORY:
! 02 Dec 2003, Sujay Kumar, Initial Version
!    Feb 2004, David Mocko, Conversion from SSiB to HY-SSiB
!
! !INTERFACE:
subroutine hyssib_binout(ftn)
! !USES:
  use lisdrv_module, only : lis
  use drv_output_mod, only : drv_writevar_bin
  use hyssib_varder

  implicit none
! !ARGUMENTS:
  integer :: ftn
!EOP
  real :: rainf(lis%d%glbnch)
  real :: snowf(lis%d%glbnch)
  integer :: t
!BOC
  do t=1,lis%d%glbnch
     if (hyssib(t)%forcing(1).lt.273.15) then
        rainf(t) = 0.0
        snowf(t) = hyssib(t)%forcing(8)
     else
        rainf(t) = hyssib(t)%forcing(8)
        snowf(t) = 0.0
     endif
  enddo
  
!---------------------------------------------------------------------------
! General Energy Balance Components
!---------------------------------------------------------------------------
  hyssib%swnet = hyssib%swnet/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%swnet)
  
  hyssib%lwnet = (-1)*hyssib%lwnet/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%lwnet)

  hyssib%qle = hyssib%qle/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%qle)

  hyssib%qh = hyssib%qh/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%qh)

  hyssib%qg = hyssib%qg/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%qg)

  hyssib%qf = hyssib%qf/float(hyssib%count)
  call drv_writevar_bin(ftn, hyssib%qf)

  hyssib%qv = hyssib%qv/float(hyssib%count)
  call drv_writevar_bin(ftn, hyssib%qv)
  
  hyssib%qtau = hyssib%qtau/float(hyssib%count)
  call drv_writevar_bin(ftn, hyssib%qtau)
  
  hyssib%qa = hyssib%qa/float(hyssib%count)
  call drv_writevar_bin(ftn, hyssib%qa)
  
  call drv_writevar_bin(ftn, hyssib%delsurfheat)
  
  call drv_writevar_bin(ftn, hyssib%delcoldcont)

  !---------------------------------------------------------------------------
  ! General Water Balance Components
  !---------------------------------------------------------------------------
  hyssib%snowf = hyssib%snowf/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%snowf)

  hyssib%rainf = hyssib%rainf/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%rainf)

  hyssib%evap = hyssib%evap/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%evap)

  hyssib%qs = hyssib%qs/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%qs)

  hyssib%qrec = hyssib%qrec/float(hyssib%count)
  call drv_writevar_bin(ftn, hyssib%qrec)

  hyssib%qsb = hyssib%qsb/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%qsb)

  hyssib%qsm = hyssib%qsm/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%qsm)

  hyssib%qfz = hyssib%qfz/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%qfz)
  
  hyssib%qst = hyssib%qst/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%qst)
  
  call drv_writevar_bin(ftn,hyssib%delsoilmoist)
  
  call drv_writevar_bin(ftn,hyssib%delswe)
  
  call drv_writevar_bin(ftn,hyssib%delintercept)
  
  !---------------------------------------------------------------------------
  ! Surface State Variables
  !---------------------------------------------------------------------------
  if (hyssibdrv%STATEVAR_AVG.eq.1) then
     do t = 1,lis%d%glbnch
        if (hyssib(t)%snowtcount.gt.0) then
           hyssib(t)%snowt = hyssib(t)%snowt/float(hyssib(t)%snowtcount)
        else
           hyssib(t)%snowt = lis%d%udef
        endif
        if (hyssib(t)%albedocount.gt.0) then
           hyssib(t)%albedo = hyssib(t)%albedo/float(hyssib(t)%albedocount)
        else
           hyssib(t)%albedo = lis%d%udef
        endif
        if (hyssib(t)%sliqfraccount.gt.0) then
           hyssib(t)%sliqfrac = hyssib(t)%sliqfrac/float(hyssib(t)%sliqfraccount)
        else
           hyssib(t)%sliqfrac = lis%d%udef
        endif
     enddo
     hyssib%vegtc = hyssib%vegtc/float(hyssib%count)
     hyssib%baresoilt = hyssib%baresoilt/float(hyssib%count)
     hyssib%avgsurft = hyssib%avgsurft/float(hyssib%count)
     hyssib%radteff = hyssib%radteff/float(hyssib%count)
     hyssib%swe = hyssib%swe/float(hyssib%count)
     hyssib%sweveg = hyssib%sweveg/float(hyssib%count)
  endif
  
  call drv_writevar_bin(ftn,hyssib%snowt)
  
  call drv_writevar_bin(ftn,hyssib%vegtc)
  
  call drv_writevar_bin(ftn,hyssib%baresoilt)

  call drv_writevar_bin(ftn,hyssib%avgsurft)

  call drv_writevar_bin(ftn, hyssib%radteff)
  
  call drv_writevar_bin(ftn,hyssib%albedo)

  call drv_writevar_bin(ftn,hyssib%swe)

  call drv_writevar_bin(ftn, hyssib%sweveg)
      
  !---------------------------------------------------------------------------
  ! Subsurface State Variables
  !---------------------------------------------------------------------------
  if (hyssibdrv%STATEVAR_AVG.eq.1) then
     hyssib%soilmoist1 = hyssib%soilmoist1/float(hyssib%count)
     hyssib%soilmoist2 = hyssib%soilmoist2/float(hyssib%count)
     hyssib%soilmoist3 = hyssib%soilmoist3/float(hyssib%count)
     hyssib%soiltemp = hyssib%soiltemp/float(hyssib%count)
     hyssib%soilwet = hyssib%soilwet/float(hyssib%count)
  endif
  
  call drv_writevar_bin(ftn,hyssib%soilmoist1)

  call drv_writevar_bin(ftn,hyssib%soilmoist2)

  call drv_writevar_bin(ftn,hyssib%soilmoist3)

  call drv_writevar_bin(ftn,hyssib%soiltemp)  

  call drv_writevar_bin(ftn,hyssib%soilwet)

  !---------------------------------------------------------------------------
  ! Evaporation Components
  !---------------------------------------------------------------------------
  hyssib%potevap = hyssib%potevap/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%potevap)
  
  hyssib%ecanop = hyssib%ecanop/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%ecanop)
  
  hyssib%tveg = hyssib%tveg/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%tveg)

  hyssib%esoil = hyssib%esoil/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%esoil)

  hyssib%rootmoist = hyssib%rootmoist/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%rootmoist)

  hyssib%subsnow = hyssib%subsnow/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%subsnow)
  
  hyssib%subsurf = hyssib%subsurf/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%subsurf)

  hyssib%acond = hyssib%acond/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%acond)

  hyssib%snowfrac = hyssib%snowfrac/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%snowfrac)
!-----------------------------------------------------------------------
! Cold Season Processes
!-----------------------------------------------------------------------
  hyssib%snowdepth = hyssib%snowdepth/float(hyssib%count)
  call drv_writevar_bin(ftn,hyssib%snowdepth)

  call drv_writevar_bin(ftn,hyssib%sliqfrac)

  if (lis%o%wfor.eq.1) then
     call drv_writevar_bin(ftn,sqrt(hyssib%forcing(5)*hyssib%forcing(5)+ & 
          hyssib%forcing(6)*hyssib%forcing(6)))
     
     call drv_writevar_bin(ftn,rainf)
     
     call drv_writevar_bin(ftn,snowf)
     
     call drv_writevar_bin(ftn,hyssib%forcing(1))
     
     call drv_writevar_bin(ftn,hyssib%forcing(2))
     
     call drv_writevar_bin(ftn,hyssib%forcing(7))
     
     call drv_writevar_bin(ftn,hyssib%forcing(3))
         
     call drv_writevar_bin(ftn,hyssib%forcing(4))
     
  endif
      
!EOC
end subroutine hyssib_binout
    
