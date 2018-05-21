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
! !ROUTINE: mos_output.F90
! 
! !DESCRIPTION: This subroutines sets up methods to write mosaic output 
! 
! !INTERFACE:
subroutine mos_output
! !USES:
  use lisdrv_module, only : lis, tile, glbgindex
  use mos_varder, only : mosdrv
  use spmdMod, only : masterproc,npes
!EOP
  integer :: i 
  real :: var(lis%d%glbnch)
!BOC
  if(lis%o%wsingle ==1) then
    print*, "Mosaic not currently able to write output for each variable in a separate file"
    print*, "Please set lis%o%wsingle = 2"
    stop
!------------------------------------------------------------------
! Writes each output variable to a separate file
!------------------------------------------------------------------
!!!     if(mod(lis%t%gmt,mosdrv%writeintm).eq.0)then
!!!        do i=1,32
!!!           call noah_singlegather(i,var)
!!!           if(masterproc) then 
!!!              call noah_singleout(lis, tile, glbgindex, var, i)
!!!           endif
!!!        enddo
!!!        call noah_totinit()
!!!     endif
  else 
!------------------------------------------------------------------
! Writes bundled output
!------------------------------------------------------------------
     if(mod(lis%t%gmt,mosdrv%writeintm).eq.0)then
        if(npes > 1 ) then 
           call mos_gather()
        endif
        if(masterproc) then 
           call mos_almaout()
        endif
        call mos_totinit()
     endif
  endif
!EOC
end subroutine mos_output
