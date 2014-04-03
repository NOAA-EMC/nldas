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
! !ROUTINE: clm2_output.F90
! 
! !DESCRIPTION: This subroutines sets up methods to write CLM2
! output 
! 
! !INTERFACE:
subroutine clm2_output
! !USES:
  use lisdrv_module, only : lis, tile
  use clm_varctl, only : clmdrv
  use spmdMod, only : masterproc, npes
!EOP
  integer :: ier
  real :: var(lis%d%glbnch)
!BOC  
  if(lis%o%wsingle==1) then 
!------------------------------------------------------------------
! Write output with each variable in a single file
!------------------------------------------------------------------
     if(mod(lis%t%gmt,clmdrv%writeintc2).eq.0)then
        do i=1,11
           call clm2_singlegather(i,var)
           if(masterproc) then 
              call clm2_singleout(var, i)
           endif
        enddo
        do i=14,47
           call clm2_singlegather(i, var)
           if(masterproc) then 
              call clm2_singleout(var,i)
           endif
        enddo
        if ( lis%o%wfor == 1 ) then
           do i=48,55
              call clm2_singlegather(i, var)
              if(masterproc) then 
                 call clm2_singleout(var,i)
              endif
           enddo
        endif
        call clm2_totinit()
     endif
  else 
!------------------------------------------------------------------
! Write bundled output 
!------------------------------------------------------------------
     if(mod(lis%t%gmt,clmdrv%writeintc2).eq.0)then
        if(npes > 1 ) then 
           call clm2_gather()
        endif
        if(masterproc) then 
           call clm2_almaout()
        endif
        call clm2_totinit()
     endif
  endif
!EOC
end subroutine clm2_output
