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
! !ROUTINE: time_interp_ecmwf.F90
!
! !DESCRIPTION:
!  Opens, reads, and interpolates ECMWF forcing.  
!
!    TIME1 = most recent past data\\
!    TIME2 = nearest future data \\
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! \subsection{Core Functions of time$_-$interp$_-$geos}
!  \begin{description}
!  \item[zterp]  
!   Performs zenith angle-based temporal interpolation
!  \end{description}
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  25 Oct 1999: Jared Entin; Significant F90 Revision
!  11 Apr 2000: Brian Cosgrove; Fixed name construction error 
!               in Subroutine ETA6HRFILE 
!  27 Apr 2000: Brian Cosgrove; Added correction for use of old shortwave
!               data with opposite sign convention from recent shortwave data.
!               Added capability to use time averaged shortwave & longwave data
!               Altered times which are passed into ZTERP--used to be GMT1 
!               and GMT2, now they are LDAS%ETATIME1 and LDAS%ETATIME2
!  30 Nov 2000: Jon Radakovich; Initial code based on geteta.f
!  17 Apr 2001: Jon Gottschalck; A few changes to allow model init.  
!  13 Aug 2001: Urszula Jambor; Introduced missing data replacement.     
!   5 Nov 2001: Urszula Jambor; Reset tiny negative SW values to zero. 
! !INTERFACE:
subroutine time_interp_ecmwf()
  ! !USES:
  use lisdrv_module, only : lis, grid    
  use baseforcing_module, only : glbdata1, glbdata2
  use time_manager
  use grid_spmdMod
  use spmdMod
  use ecmwfdomain_module, only : ecmwfdrv
  !EOP
  implicit none

  integer :: ier
  integer :: c,r,f,zdoy,idoy,iyr,imo,ida,ihr,imn,its,iss
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8  :: btime,inittime
  real    :: wt1,wt2,czb,cze,czm,gmt1,gmt2,igmt
  real    :: zw1,zw2
  !BOC
  if(masterproc) then 
     if (get_nstep(lis%t) .eq. 0) then
        lis%f%nforce = lis%f%nmif
     else
        lis%f%nforce = lis%f%nf
     endif
  endif
#if(defined SPMD)
  call MPI_BCAST(ecmwfdrv%ecmwftime1,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(ecmwfdrv%ecmwftime2,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%t%time,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%t%gmt,1,MPI_REAL,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%nforce,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%rstflag,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier) 
#endif

  !=== Reset GMT times
  btime=ecmwfdrv%ecmwftime1
  call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=ecmwfdrv%ecmwftime2
  call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  !=== Need lower SW time boundary for call to zterp based on time2
  inittime=btime
  call time2date( inittime, idoy, igmt, iyr, imo, ida, ihr, imn )
  select case (ihr)
  case(00,12,24)
     its = -12*60*60
  case(03,15)
     its =  -3*60*60
  case(06,18)
     its =  -6*60*60
  case(09,21)
     its =  -9*60*60
  end select
  call tick( inittime, idoy, igmt, iyr, imo, ida, ihr, imn, iss, its )

  !===  Interpolate Data in Time
  wt1=(ecmwfdrv%ecmwftime2-lis%t%time)/ & 
       (ecmwfdrv%ecmwftime2-ecmwfdrv%ecmwftime1)
  wt2=1.0-wt1
  do f=1,lis%f%nforce
     if (f == 3) then     ! Time Averaged Shortwave
        do c = 1, gdi(iam)
           zdoy = lis%t%doy
           call zterp(0,grid(c)%lat,grid(c)%lon, &
                igmt,gmt2,lis%t%gmt,zdoy, &
                zw1,zw2,czb,cze,czm,lis)
           grid(c)%forcing(f) = zw1 * glbdata2(f,c)
           if ((grid(c)%forcing(f).ne.lis%d%udef).and.  &
                (grid(c)%forcing(f).lt.0)                ) then
              if (grid(c)%forcing(f) > -0.00001) then !arbitrary!!
                 grid(c)%forcing(f) = 0.0             !threshold!!
              else
                 print*, 'Stopping because SW forcing not udef but lt 0, '
                 print*, f,c,r,grid(c)%forcing(f)
                 stop
              endif
           endif
           if (grid(c)%forcing(f).gt.1367) then
              print*, '!!! ',grid(c)%forcing(f)
              print*, '!!! zterp produced d-SW-flux > solar constant!!!'
              stop
              grid(c)%forcing(f) = glbdata2(f,c)
           endif
        enddo ! c

     else if (f == 4) then     
        do c=1,gdi(iam)	      
           grid(c)%forcing(f)=glbdata2(f,c)
        enddo

     else if (f == 8 .or. f == 9) then    ! precip variable Block Interpolation
        do c=1,gdi(iam)
           grid(c)%forcing(f)=glbdata2(f,c)
        enddo
     else     !Linearly interpolate everything else
        do c=1,gdi(iam)
           grid(c)%forcing(f)=glbdata1(f,c)*wt1+glbdata2(f,c)*wt2
        enddo
     endif
  enddo   !the f loop

end subroutine time_interp_ecmwf
