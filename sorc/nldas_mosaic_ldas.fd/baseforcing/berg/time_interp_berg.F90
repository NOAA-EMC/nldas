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
! !MODULE: time_interp_berg
!
! !DESCRIPTION:
!  Opens, reads, and interpolates BERG forcing.  
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data 
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
! 
! !REVISION HISTORY: 
! 24 Nov 2003: Sujay Kumar; Initial version of code adapted from the GLDAS 
!                           version
! 
! !INTERFACE:
subroutine time_interp_berg()
! !USES:
  use lisdrv_module, only : lis, grid
  use baseforcing_module, only : glbdata1,glbdata2
  use time_manager
  use bergdomain_module, only : bergdrv
  use grid_spmdMod
  use spmdMod
!EOP
  implicit none
  real :: wt1,wt2,zw1,zw2,czb,cze,czm,gmt1,gmt2
  integer :: zdoy, c,f,ier
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8 :: btime
  integer, parameter :: nforce=9  ! # forcing variables
  
  lis%f%nforce = nforce
#if (defined SPMD)
  call MPI_BCAST(bergdrv%fmodeltime1,1,MPI_REAL8,0,& 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(bergdrv%fmodeltime2,1,MPI_REAL8,0,& 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%nforce,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%t%time,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%t%gmt,1,MPI_REAL,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%rstflag,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier) 
#endif
  btime=bergdrv%fmodeltime1
  call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=bergdrv%fmodeltime2
  call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  
  !===  Interpolate Data in Time
  wt1=(bergdrv%fmodeltime2-lis%t%time)/&
       (bergdrv%fmodeltime2-bergdrv%fmodeltime1)
  wt2=1.0-wt1
  do f=1,nforce
     
     if(f.eq.3) then     ! Time Averaged Shortwave       
        do c=1,gdi(iam)
           zdoy=lis%t%doy
           call zterp(0,grid(c)%lat,grid(c)%lon, &
                gmt1,gmt2,lis%t%gmt,zdoy, &
                zw1,zw2,czb,cze,czm,lis)
           grid(c)%forcing(f)=glbdata1(f,c)*zw1
           if ((grid(c)%forcing(f).ne.lis%d%udef).and.  &
                (grid(c)%forcing(f).lt.0)                ) then
              if (grid(c)%forcing(f) > -0.00001) then !arbitrary!!
                 grid(c)%forcing(f) = 0.0             !threshold!!
              else
                 print*, 'Stopping because forcing not udef but lt 0, '
                 print*, f,c,grid(c)%forcing(f)
                 stop
              end if
           endif
           
           if (grid(c)%forcing(f).gt.1367) then
              grid(c)%forcing(f)=glbdata1(f,c)
           endif
        enddo
     else if (f.eq.4) then     ! Time Averaged Longwave, Block Interpolation
        do c=1,gdi(iam)
           grid(c)%forcing(f)=glbdata1(f,c)
        enddo
     else if (f.eq.6) then     ! Set to Zero, f=5 is magnitude of wind
        do c=1,gdi(iam)
           grid(c)%forcing(f)=0.0
        enddo
     else if(f.eq.8.or.f.eq.9) then    ! precip variable Block Interpolation
        do c=1,gdi(iam)
           grid(c)%forcing(f)=glbdata1(f,c)
        enddo
     else     !Linearly interpolate everything else	
        do c=1,gdi(iam)
           grid(c)%forcing(f)=glbdata1(f,c)*wt1+ &
                glbdata2(f,c)*wt2
        enddo
     endif
  enddo   !the f loop

end subroutine time_interp_berg
