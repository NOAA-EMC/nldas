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
! !ROUTINE: time_interp_gdas
!
! !DESCRIPTION:
!  Opens, reads, and interpolates GDAS forcing.  
!
!    TIME1 = most recent past data\\
!    TIME2 = nearest future data \\
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! \subsection{Core Functions of time$_-$interp$_-$gdas}
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
subroutine time_interp_gdas()
! !USES:
  use lisdrv_module, only :lis, grid  
  use baseforcing_module, only: glbdata1, glbdata2
  use time_manager
  use grid_spmdMod
  use spmdMod
  use gdasdomain_module, only : gdasdrv
!EOP
  implicit none
!==== Local Variables=======================
  integer :: ier
  integer :: c,r,f,zdoy,idoy,iyr,imo,ida,ihr,imn,its,iss
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8 :: btime,inittime
  real :: wt1,wt2,czb,cze,czm,gmt1,gmt2,igmt
  real :: zw1,zw2
!BOC
  if(masterproc) then 
     if (get_nstep(lis%t) .eq. 0) then
        lis%f%nforce = gdasdrv%nmif
     else
        lis%f%nforce = lis%f%nf
     endif
  endif
#if(defined SPMD)
  call MPI_BCAST(gdasdrv%gdastime1,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(gdasdrv%gdastime2,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%t%time,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%t%gmt,1,MPI_REAL,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%nforce,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%F00_flag,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%F06_flag,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%rstflag,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier) 
#endif
  btime=gdasdrv%gdastime1
  call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=gdasdrv%gdastime2
  call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  if ( (lis%f%F06_flag==0) .and. (lis%f%F00_flag==1) ) then
     inittime = gdasdrv%gdastime2
     call time2date( inittime, idoy, igmt, iyr, imo, ida, ihr, imn )
     its = -6*60*60
     call tick( inittime, idoy, igmt,iyr,imo,ida,ihr,imn,iss,its)
  end if
  wt1 = (gdasdrv%gdastime2-lis%t%time) / & 
       (gdasdrv%gdastime2-gdasdrv%gdastime1)
  wt2 = 1.0 - wt1
  do f=1,lis%f%nforce
     if ( f == 3 ) then  !shortwave         
        do c = 1, gdi(iam)
           zdoy = lis%t%doy
!           print*, 'f0,f6',lis%f%F00_flag,lis%f%F06_flag
           if ( (lis%f%F06_flag==0) .and. (lis%f%F00_flag==1) ) then
              call zterp( 0, grid(c)%lat, grid(c)%lon, igmt, gmt2, & 
                   lis%t%gmt,zdoy,zw1,zw2,czb,cze,czm,lis)
!              if(c==200) then 
!                 print*, grid(c)%lat, grid(c)%lon
!                 print*, igmt
!                 print*, gmt2
!                 print*, lis%t%gmt
!                 print*, zdoy
!              endif
           else
              call zterp( 0, grid(c)%lat, grid(c)%lon, gmt1, gmt2, & 
                   lis%t%gmt,zdoy,zw1,zw2,czb,cze,czm,lis)
!                 print*, 'blk2 ',c,zw1,zw2
           end if

           grid(c)%forcing(f) = zw1 * glbdata2(f,c)
!           if(c==200) print*, c,grid(c)%forcing(f),glbdata2(f,c),zw1
           if (grid(c)%forcing(f) < 0) then
              print *, '2 warning!!!  SW radiation is negative!!'
              print *, 'sw=', grid(c)%forcing(f), '... negative'
              print *, 'gdas2=', glbdata2(f,c)
              call endrun
           end if
           
           if (grid(c)%forcing(f).gt.1367) then
!              print*, 'in this high block..',c,f,grid(c)%forcing(f)
              grid(c)%forcing(f)=glbdata2(f,c)
           endif
        end do
     else if ( (f==4) .or. (f==10) ) then         
        do c = 1, gdi(iam)
           if ( lis%f%F00_flag==1 ) then 
              grid(c)%forcing(f)=2*glbdata2(f,c) -glbdata1(f,c)
           else
              grid(c)%forcing(f) = glbdata2(f,c)
           end if
        end do
        
     else if ( (f==8) .or. (f==9 ) ) then
        do c = 1, gdi(iam)
           if ( lis%f%F00_flag == 1) then 
              if (2*glbdata2(f,c) >= glbdata1(f,c)) then
                 grid(c)%forcing(f)=2*glbdata2(f,c)-glbdata1(f,c)
              else
                 grid(c)%forcing(f) = 0.0
              end if
           else
              grid(c)%forcing(f) = glbdata2(f,c)
           end if
        end do
     else 
        do c = 1, gdi(iam)
           grid(c)%forcing(f) = wt1 * glbdata1(f,c) +  & 
                wt2 *glbdata2(f,c)
        end do
     end if
  end do
  return
!EOC    
end subroutine time_interp_gdas
