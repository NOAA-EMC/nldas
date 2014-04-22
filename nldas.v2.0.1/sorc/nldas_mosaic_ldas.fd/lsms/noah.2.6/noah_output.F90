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
! !ROUTINE: noah_output.F90
! 
! !DESCRIPTION: This subroutines sets up methods to write noah
! output 
! 
! !INTERFACE:
subroutine noah_output
! !USES:
  use lisdrv_module, only : lis, tile, glbgindex
  use noah_varder, only : noahdrv
  use spmdMod, only : masterproc,npes
!EOP
  integer :: i 
  real :: var(lis%d%glbnch)
!BOC
!------------------------------------------------------------------
! Number of variables to be outputted from Noah LSM 
!------------------------------------------------------------------
  if ( lis%o%wsingle == 1 ) then 
     if ( lis%o%wfor == 0 ) then
        num_vars = 30
     else
        num_vars = 38
     endif
!------------------------------------------------------------------
! Writes each output variable to a separate file
!------------------------------------------------------------------
     if ( mod(lis%t%gmt,noahdrv%writeintn) == 0 ) then
        do i=1,num_vars
           call noah_singlegather(i,var)
           if ( masterproc ) then 
              call noah_singleout(lis, tile, var, i)
           endif
        enddo
        call noah_totinit()
     endif
  else 
!------------------------------------------------------------------
! Writes bundled output
!------------------------------------------------------------------
     if(mod(lis%t%gmt,noahdrv%writeintn).eq.0)then
        if(npes > 1 ) then 
           call noah_gather()
        endif
        if(masterproc) then 
           call noah_almaout()
        endif
        call noah_totinit()
     endif
  endif
!EOC
end subroutine noah_output

