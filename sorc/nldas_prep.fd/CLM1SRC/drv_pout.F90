#include <misc.h>

subroutine drv_pout (drv, tile, clm1)

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  	
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:		  
!  Write out area-average CLM parameters
!
! REVISION HISTORY:
!  15 December 1999:  Paul Houser and Jon Radakovich; Initial Code 
!  15 November 2000:  Mariana Vertenstein
!=========================================================================
! $Id: drv_pout.F90,v 1.1.1.1 2003/02/06 16:10:44 jgottsch Exp $
!=========================================================================

  use precision
  use drv_module          ! 1-D Land Model Driver variables
  use drv_tilemodule      ! Tile-space variables
  use clm1type             ! 1-D CLM variables
  use clm1_varpar, only : nlevsoi, nlevsno
  use clm1_varcon, only : tcrit
  implicit none

!=== Arguments ===========================================================

  type (drvdec)      :: drv              
  type (clm_tiledec) :: tile(drv%nch)
  type (clm11d)       :: clm1 (drv%nch)

!=== Local Variables =====================================================

  integer  :: t,n,l         !Tile and grid space counters
  integer  :: nch           !Temporary counters
  integer  :: mask(drv%nch) !Water mask
  real(r8) :: drv_gridave   !Spatial Averaging Function
  real(r8) :: temp(drv%nch),avg !Temporary for converting integer->real arrays

!=== End Variable List ===================================================

  nch=drv%nch

! determine grid average weights - only soil points will be averaged for print out

  do t=1,drv%nch 
!     if (clm1(t)%lakpoi) then
!        mask(t) = 0
!     else
!        mask(t) = 1
!     endif
      mask(t)=1
  end do

! 
  n=5
  open(n,file=drv%poutf1d,form='formatted')

! Write out 1-D parameters

  write(n,1)'Start Time:                ', &
       drv%smo,'/',drv%sda,'/',drv%syr,drv%shr,':',drv%smn,':',drv%sss 
  write(n,1)'End Time:                  ', &
       drv%emo,'/',drv%eda,'/',drv%eyr,drv%ehr,':',drv%emn,':',drv%ess 

  write(n,3)'Timestep (sec):            ',drv%ts 
  write(n,3)'CLM startcode:             ',drv%startcode 
  write(n,3)'CLM IC Source:             ',drv%clm_ic
  write(n,3)'Number of Columns:         ',drv%nc 
  write(n,3)'Number of Rows:            ',drv%nr 
  write(n,3)'Veg Class Scheme:          ',drv%vclass
  write(n,3)'Number of Veg Classes:     ',drv%nt
  write(n,3)'Maximum tile per grid:     ',drv%maxt
  write(n,2)'Min area for tile (%):     ',drv%mina
  write(n,2)'Undefined Value:           ',drv%udef
  write(n,3)'Total Model Tiles:         ',drv%nch
  write(n,4)'Veg Tile Spec File:        ',drv%vegtf 
  write(n,4)'Veg Type Parameter Values: ',drv%vegpf 
  write(n,4)'CLM 1D Met input file:     ',drv%metf1d 
  write(n,4)'CLM 1D Output File:        ',drv%outf1d 
  write(n,4)'CLM 1D Para Output File:   ',drv%poutf1d 
  write(n,4)'CLM active restart file:   ',drv%rstf 

  write(n,*)
  write(n,3)'Number of Soil Layers      ',nlevsoi
  write(n,3)'Max number of snow layers  ',nlevsno
  write(n,2)'Critical Rain/Snow Temp:   ',tcrit
  write(n,*)

! Write out 2-D parameters  
  write(n,2)'Average Latitude:          ',drv_gridave(nch,mask,tile%fgrd,clm1%latdeg)
  write(n,2)'Average Longitude:         ',drv_gridave(nch,mask,tile%fgrd,clm1%londeg)
  write(n,2)'Average Wind Obs. Height   ',drv_gridave(nch,mask,tile%fgrd,clm1%forc_hgt_u)
  write(n,2)'Average Temp Obs. Height:  ',drv_gridave(nch,mask,tile%fgrd,clm1%forc_hgt_t)
  write(n,2)'Average Humid Obs. Height: ',drv_gridave(nch,mask,tile%fgrd,clm1%forc_hgt_q)
  write(n,2)'Average Max Allowed Dew    ',drv_gridave(nch,mask,tile%fgrd,clm1%dewmx)
  temp = real(tile%vegt)
  write(n,2)'Average Veg Index          ',drv_gridave(nch,mask,tile%fgrd,temp)
  write(n,2)'Average Aero Rough Length  ',drv_gridave(nch,mask,tile%fgrd,clm1%z0m)
  write(n,2)'Average Displacement Height',drv_gridave(nch,mask,tile%fgrd,clm1%displa)
  write(n,2)'Average Leaf Dimension     ',drv_gridave(nch,mask,tile%fgrd,clm1%dleaf)
  write(n,2)'Average Stem Area Index    ',drv_gridave(nch,mask,tile%fgrd,clm1%tsai)
  write(n,2)'Average High WTable Fract: ',drv_gridave(nch,mask,tile%fgrd,clm1%wtfact)
  write(n,2)'Average Wilting Point:     ',drv_gridave(nch,mask,tile%fgrd,clm1%smpmax)
  write(n,2)'Average Irr Snow Sat:      ',drv_gridave(nch,mask,tile%fgrd,clm1%ssi)
  write(n,2)'Average Imp Porosity:      ',drv_gridave(nch,mask,tile%fgrd,clm1%wimp)
  write(n,2)'Average Ponding Depth:     ',drv_gridave(nch,mask,tile%fgrd,clm1%pondmx)
  write(n,2)'Average Vert Scale Factor: ',drv_gridave(nch,mask,tile%fgrd,tile%scalez)
  write(n,2)'Average KSAT Scale decreas:',drv_gridave(nch,mask,tile%fgrd,tile%hkdepth)
  temp = real(clm1%isoicol)
  write(n,2)'Average Soil Color:        ',drv_gridave(nch,mask,tile%fgrd,temp)
  write(n,2)'Average Soil Rough Length: ',drv_gridave(nch,mask,tile%fgrd,clm1%zlnd)
  write(n,2)'Average Snow Rough Length: ',drv_gridave(nch,mask,tile%fgrd,clm1%zsno)
  write(n,2)'Average Under Canopy DCoef:',drv_gridave(nch,mask,tile%fgrd,clm1%csoilc)
  write(n,*)   
  write(n,*)'Note: IC values below are after first model timestep'
  write(n,*)
  write(n,2)'Average Initial Surf Temp: ',drv_gridave(nch,mask,tile%fgrd,clm1%t_grnd)
  write(n,2)'Average Initial Leaf Temp: ',drv_gridave(nch,mask,tile%fgrd,clm1%t_veg)
  write(n,2)'Average Initial Snow Water:',drv_gridave(nch,mask,tile%fgrd,clm1%h2osno)
  write(n,2)'Average Initial Snow Age:  ',drv_gridave(nch,mask,tile%fgrd,clm1%snowage)
  write(n,2)'Average Initial Snow Depth:',drv_gridave(nch,mask,tile%fgrd,clm1%snowdp)
  write(n,2)'Average Initial Intercept: ',drv_gridave(nch,mask,tile%fgrd,clm1%h2ocan)
  write(n,2)'Average LayT to SurfT Fact:',drv_gridave(nch,mask,tile%fgrd,clm1%capr)
  write(n,2)'Average Crank Nichol Fact: ',drv_gridave(nch,mask,tile%fgrd,clm1%cnfac)
  temp = real(clm1%snl)
  write(n,2)'Average Num of Snow Layers:',drv_gridave(nch,mask,tile%fgrd,temp)
  write(n,2)'Maximum leaf area index    ',drv_gridave(nch,mask,tile%fgrd,clm1%maxlai)
  write(n,2)'Minimum leaf area index    ',drv_gridave(nch,mask,tile%fgrd,clm1%minlai)
  write(n,*)'Restrict min of soil poten ',drv_gridave(nch,mask,tile%fgrd,clm1%smpmin)
  write(n,2)'ROOTA                      ',drv_gridave(nch,mask,tile%fgrd,tile%roota)
  write(n,2)'ROOTB                      ',drv_gridave(nch,mask,tile%fgrd,tile%rootb)

! Write out 3-D parameters

  do l = 1,nlevsoi        
     write(n,*)
     write(n,5)'Average Hydro Cond at Sat   lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%hksat(l))
     write(n,5)'Average Soil Void Fract     lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%watsat(l))
     write(n,5)'Average Min Soil Suction    lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%sucsat(l))
     write(n,5)'Percent SAND                lay: ', l, drv_gridave(nch,mask,tile%fgrd,tile%sand(l))
     write(n,5)'Percent CLAY                lay: ', l, drv_gridave(nch,mask,tile%fgrd,tile%clay(l))
     write(n,5)'Clapp and Horn "b" para     lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%bsw(l))
     write(n,5)'Heat cap of soil soilds     lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%csol(l))
     write(n,5)'Thermal conduct of soil min lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%tkmg(l))
     write(n,5)'Thermal conduct sat soil    lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%tksatu(l))
     write(n,5)'Thermal conduct dry soil    lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%tkdry(l))
     write(n,5)'Average Depth Soil          lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%z(l))
     write(n,5)'Average Thick Soil          lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%dz(l))
     write(n,5)'Average Init Temp           lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%t_soisno(l))
     write(n,5)'Average Init Water          lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%h2osoi_liq(l))
     write(n,5)'Average Init Ice            lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%h2osoi_ice(l))
     write(n,5)'Root Fraction               lay: ', l, drv_gridave(nch,mask,tile%fgrd,clm1%rootfr(l))
  enddo

  close(n)

1 format(a28,1(1x,i2,a1,i2,a1,i4),1(1x,i2,a1,i2,a1,i2))
2 format(a28,1x,f20.10)
3 format(a28,1x,i10)
4 format(a28,1x,a20)
5 format(A32,1x,i3,1x,f20.10)

  return
end subroutine drv_pout







