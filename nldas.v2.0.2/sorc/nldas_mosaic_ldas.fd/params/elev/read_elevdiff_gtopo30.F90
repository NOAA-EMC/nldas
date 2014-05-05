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
! !ROUTINE: read_elevdiff_gtopo30
!
! !DESCRIPTION:
!  This subroutine retrieves UMD-AVHRR landcover data

! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  12 May 2005: James Geiger; Added opendap support
!
! !INTERFACE:
subroutine read_elevdiff_gtopo30(elevdiff)
! !USES:
  use lisdrv_module,      only: lis
  use lis_indices_module, only: lis_nc_data, lis_nr_data
  use lis_openfileMod

  implicit none

  integer :: line1, line2
  integer :: nc_dom, ierr
  real    :: elevdiff(lis_nc_data, lis_nr_data)
  !real    :: elevdiff(lis%d%lnc, lis%d%lnr)
  integer :: c, r, glnc, glnr, line
  !logical :: first_read = .TRUE.
  character(len=40) :: file_name
!EOP      

  if ( lis%f%ecor .eq. 1 ) then 

     elevdiff = 0.0

     call lis_open_file(21,file=lis%p%elevfile,form='unformatted', &
                        access='direct',recl=4,status='old',       &
                        script='getelev.pl')

     call lis_log_msg('MSG: read_elevdiff_gtopo30 -- Reading elevation file: '//trim(lis%p%elevfile))

#if 0
     ! Note: This routine is called during initialization of the domain
     ! (making the tiles) before the opendap module is initialized.
     ! Therefore the first time the elevation data is read, I assume that
     ! only the master process calling this routine and processing the data.
     if ( first_read ) then
        line1 = nint((lis%d%gridDesc(4)-lis%d%elev_gridDesc(1))/lis%d%gridDesc(9)) + 1
        line2 = nint((lis%d%gridDesc(5)-lis%d%elev_gridDesc(2))/lis%d%gridDesc(10)) + 1
     
        nc_dom = nint((lis%d%elev_gridDesc(4)-lis%d%elev_gridDesc(2))/&
                      lis%d%elev_gridDesc(5))+1        
     
        do r=1,lis_nr_data
           do c=1,lis_nc_data
              glnc = line2+c-1
              glnr = line1+r-1
              line = (glnr-1)*nc_dom+glnc
              read(21,rec=line) elevdiff(c,r)
           enddo
        enddo
        first_read = .FALSE.
     else
        line = 1
        do r=1,lis_nr_data
           do c=1,lis_nc_data
              read(21,rec=line) elevdiff(c,r)
              line = line + 1
           enddo
        enddo
     endif
#else
     line1 = nint((lis%d%gridDesc(4)-lis%d%elev_gridDesc(1))/lis%d%gridDesc(9)) + 1
     line2 = nint((lis%d%gridDesc(5)-lis%d%elev_gridDesc(2))/lis%d%gridDesc(10)) + 1
     
     nc_dom = nint((lis%d%elev_gridDesc(4)-lis%d%elev_gridDesc(2))/&
          lis%d%elev_gridDesc(6))+1
     
     do r=1,lis_nr_data!lis%d%lnr
        do c=1,lis_nc_data!lis%d%lnc
           glnc = line2+c-1
           glnr = line1+r-1
           line = (glnr-1)*nc_dom+glnc
           read(21,rec=line) elevdiff(c,r)
        enddo
     enddo
#endif
     close(21)
     call lis_log_msg('MSG: read_elevdiff_gtopo30 -- Done reading elevation difference file')
  endif
  
!EOC
end subroutine read_elevdiff_gtopo30

