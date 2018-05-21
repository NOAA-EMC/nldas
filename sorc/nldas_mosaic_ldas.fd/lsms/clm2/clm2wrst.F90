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
! !ROUTINE:clm2wrst.F90
!
! !DESCRIPTION:
!  This program writes the restart files for CLM
!
! !REVISION HISTORY:
! 20 Jan 2003; Sujay Kumar Initial Specification
! 26 Aug 2004; James Geiger, Added support for GrADS-DODS based and MPI based
!                            parallel simulations.
! 
! !INTERFACE:
subroutine clm2wrst()
! !USES:
  use spmdMod, only : masterproc, npes
  use restFileMod, only : restwrt
  use clm_varctl, only : clmdrv
  use lisdrv_module, only : lis
!EOP
  implicit none
!BOC


  if ( ( lis%t%gmt == (24-clmdrv%writeintc2) )  .or. &
       lis%t%endtime == 1 ) then

#if ( ! defined OPENDAP )
     if ( lis%o%wsingle == 1 .and. npes > 1 ) then 
        call clm2_gather()
     endif

     if ( masterproc ) then 
#endif
        call lis_log_msg('MSG: clm2wrst -- Writing CLM restart')
        call restwrt()
#if ( ! defined OPENDAP )
     endif 
#endif

  endif  
!EOC
end subroutine clm2wrst
