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
! !ROUTINE: read_nldas_soilclass
!
! !DESCRIPTION:
!  This subroutine retrieves STATSGO soils data

! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_nldas_soilclass(array)
! !USES:
  use lisdrv_module, only : lis, tile
  use lis_openfileMod
  use lis_indices_module
!EOP      
  implicit none

  real,          intent(inout) :: array(lis_nc_data, lis_nc_data, 11)

  integer :: line1, line2, line
  integer :: c,r, glnc, glnr,i,j,ii,jj,k
  integer :: nc_dom
 
  print*, 'MSG: Reading NLDAS texture file'
	print *,'now reading in nldas files ADSFADFASDFADFASDFASDF'
        print *,'now reading in nldas files ADSFADFASDFADFASDFASDF'
        print *,'now reading in nldas files ADSFADFASDFADFASDFASDF'
        print *,'now reading in nldas files ADSFADFASDFADFASDFASDF'
        print *,'now reading in nldas files ADSFADFASDFADFASDFASDF'
        print *,'now reading in nldas files ADSFADFASDFADFASDFASDF'
        print *,'now reading in nldas files ADSFADFASDFADFASDFASDF'
        print *,'now reading in nldas files ADSFADFASDFADFASDFASDF'

!    N-LDAS COSBY/RAWLS soil parameterization.
!+++   Open soil texture file
        OPEN(14,FILE=LIS%p%soilclass_file,STATUS='OLD')
        do k=1,11
        do j=1,lis%d%lnr
        do i=1,lis%d%lnc
        read (14,'(I3,1X,I3,1X,3X,I2)') &
       ii,jj,array(k,i,j)
        if (array(k,i,j).eq.13) array(k,i,j)=12
        enddo
        enddo
        enddo
        print *,'   '
        print *,'SETTING SOIL CLASS 13 to 12...since have no params for organic class'
        print *,'WARNING!!!!!!!!!!!!!!!!!!'
        print *,'   '
        CLOSE(14)

  print*, 'MSG: read soil texture file'

!EOC
end subroutine read_nldas_soilclass
