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
! !ROUTINE: vic_output.F90
!
! !DESCRIPTION:
! This subroutines sets up methods to write VIC output
!
! !REVISION HISTORY:
!   14 Apr 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:
subroutine vic_output
! !USES:
  use lisdrv_module, only : lis
  use vic_varder,    only : vicdrv
  use spmdMod,       only : masterproc, npes
!EOP
  implicit none
  integer :: i
!  character(len=( len(lis%o%odir_array(1)) + 1)*lis%o%odirn + 1) :: dir_list

  real, dimension(:), allocatable :: tmp
!BOC
  if ( lis%o%wsingle == 1 ) then 
     if ( mod(lis%t%gmt,vicdrv%writeintvic) == 0 ) then
!-------------------------------------------------------------------------
!  Writes a separate output file for each variable
!-------------------------------------------------------------------------
!        if ( masterproc ) then
!           call pack_string_f2c(lis%o%odirn, lis%o%odir_array, dir_list)
!        endif

        if ( masterproc ) then
           allocate(tmp(lis%d%glbnch))
        else
           allocate(tmp(1))
        endif

        do i=1,27
           call vic_singlegather(i,tmp)
           if ( masterproc ) then
              call vic_singleout(i,tmp)
              !call vic_singleout(lis%o%startcode, vicdrv%vicopen,&
              !     lis%d%glbnch,lis%o%wfor,lis%o%wout, &
              !     lis%o%expcode,vicdrv%vic_snowband,vicdrv%vic_nlayer,&
              !     i,lis%d%domain,lis%d%lnc,lis%d%lnr,lis%o%odirn,dir_list)
           endif
        enddo
        do i=30,30
           call vic_singlegather(i,tmp)
           if ( masterproc ) then
              call vic_singleout(i,tmp)
              !call vic_singleout(lis%o%startcode, vicdrv%vicopen,&
              !     lis%d%glbnch,lis%o%wfor,lis%o%wout, &
              !     lis%o%expcode,vicdrv%vic_snowband,vicdrv%vic_nlayer,&
              !     i,lis%d%domain,lis%d%lnc,lis%d%lnr,lis%o%odirn,dir_list)
           endif
        enddo
        if ( lis%o%wfor == 1 ) then
           do i=31,38
              call vic_singlegather(i,tmp)
              if ( masterproc ) then
                 call vic_singleout(i,tmp)
              endif
           enddo
        endif

        deallocate(tmp)

        call vic_totinit()

     endif
  else
!-------------------------------------------------------------------------
!  Write bundled output
!-------------------------------------------------------------------------
     if ( mod(lis%t%gmt,vicdrv%writeintvic) == 0 ) then
        if ( npes > 1 ) then 
           call vic_gather()
        endif
        if ( masterproc ) then 
           call vic_almaout()
!           call vic_almaout(lis%o%startcode, vicdrv%vicopen, &
!                lis%d%glbnch,lis%o%wfor, lis%o%wout,  &
!                lis%o%expcode,vicdrv%vic_snowband,vicdrv%vic_nlayer)
        endif
        call vic_totinit()
     endif
  endif
!EOC
end subroutine vic_output
