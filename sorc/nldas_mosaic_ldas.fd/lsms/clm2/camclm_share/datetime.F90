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
   subroutine datetime(cdate, ctime) 
!-----------------------------------------------------------------------
!
! Purpose:
!
!  A generic Date and Time routine
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
!
! $Id: datetime.F90,v 1.5 2004/05/07 22:18:33 jim Exp $
!
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
!-----------------------------Arguments---------------------------------
   character , intent(out) :: cdate*8 
   character , intent(out) :: ctime*8 
!-----------------------------------------------------------------------
!
!---------------------------Local Variables------------------------------
   integer, dimension(8) :: values 
   character :: date*8, time*10, zone*5 
!-----------------------------------------------------------------------
 
   call date_and_time (date, time, zone, values) 
   cdate(1:2) = date(5:6) 
   cdate(3:3) = '/' 
   cdate(4:5) = date(7:8) 
   cdate(6:6) = '/' 
   cdate(7:8) = date(3:4) 
   ctime(1:2) = time(1:2) 
   ctime(3:3) = ':' 
   ctime(4:5) = time(3:4) 
   ctime(6:6) = ':' 
   ctime(7:8) = time(5:6) 
 
   return  
   end subroutine datetime 
 
