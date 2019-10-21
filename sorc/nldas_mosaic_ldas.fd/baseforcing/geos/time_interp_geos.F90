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
! !ROUTINE: time_interp_geos.F90
!
! !DESCRIPTION:
!  Opens, reads, and interpolates GEOS forcing.  
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
subroutine time_interp_geos()
! !USES:
  use lisdrv_module, only : lis, grid    
  use baseforcing_module, only : glbdata1, glbdata2
  use time_manager
  use grid_spmdMod
  use spmdMod
  use geosdomain_module, only :geosdrv
!EOP
  implicit none
!==== Local Variables=======================
  integer :: ier
  integer :: c,f,zdoy
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8 :: btime
  real :: wt1,wt2,czb,cze,czm,gmt1,gmt2
  real :: zw1,zw2

!BOC
  if(masterproc) then 
     if (get_nstep(lis%t) .eq. 0) then
        lis%f%nforce = lis%f%nmif
     else
        lis%f%nforce = lis%f%nf
     endif
  endif

#if(defined SPMD)
  call MPI_BCAST(geosdrv%geostime1,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(geosdrv%geostime2,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%nforce,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%t%time,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%t%gmt,1,MPI_REAL,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%shortflag,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%longflag,1,MPI_INTEGER,0, & 
         MPI_COMM_WORLD,ier)
  call MPI_BCAST(lis%f%rstflag,1,MPI_INTEGER,0, & 
       MPI_COMM_WORLD,ier) 
#endif
  btime=geosdrv%geostime1
  call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=geosdrv%geostime2
  call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
!-----------------------------------------------------------------------
!  Interpolate Data in Time
!-----------------------------------------------------------------------
  wt1=(geosdrv%geostime2-lis%t%time)/ & 
       (geosdrv%geostime2-geosdrv%geostime1)
  wt2=1.0-wt1
  do f=1,lis%f%nforce
     if(f.eq.3) then     
        if (lis%f%shortflag.eq.2) then
!-----------------------------------------------------------------------
! Got Time Averaged SW
!-----------------------------------------------------------------------
           do c=1,gdi(iam)
              zdoy=lis%t%doy
!              print*, c, zdoy, grid(c)%lat, grid(c)%lon, &
!                   gmt1, gmt2, lis%t%gmt
              call zterp(0,grid(c)%lat,grid(c)%lon, & 
                   gmt1,gmt2,lis%t%gmt,zdoy, & 
                   zw1,zw2,czb,cze,czm,lis)
              grid(c)%forcing(f)=glbdata2(f,c)*zw1
              if ((grid(c)%forcing(f).ne.lis%d%udef).and. & 
                   (grid(c)%forcing(f).lt.0) ) then
                 if (grid(c)%forcing(f) > -0.00001) then 
                    grid(c)%forcing(f) = 0.0 
                 else
                    print*,'ERR: time_interp_geos -- Stopping because ', & 
                         'forcing not udef but lt0,'
                    print*,'ERR: time_interp_geos -- ', & 
                         'f,c,grid(c)%forcing(f),glbdata2(f,c)',  & 
                         f,c,grid(c)%forcing(f),glbdata2(f,c), & 
                         ' (',iam,')'
                    call endrun
                 end if
              endif
              
              if (grid(c)%forcing(f).gt.1367) then
                 grid(c)%forcing(f)=glbdata2(f,c)
              endif
           enddo
        endif
        
     else if(f.eq.8.or.f.eq.9) then 
!-----------------------------------------------------------------------
! precip variable Block Interpolation
!-----------------------------------------------------------------------
        do c=1,gdi(iam)
           grid(c)%forcing(f)=glbdata2(f,c)
        enddo
        
     else if (f.eq.4) then     
        if (lis%f%longflag.eq.1) then 
!-----------------------------------------------------------------------
!    Got Instantaneous LW
!-----------------------------------------------------------------------
           do c=1,gdi(iam)
              grid(c)%forcing(f)=glbdata1(f,c)*wt1+ & 
                   glbdata2(f,c)*wt2 
           enddo
        endif
        
        if (lis%f%longflag.eq.2) then 
!-----------------------------------------------------------------------
!    Got Time Averaged LW
!-----------------------------------------------------------------------
           do c=1,gdi(iam)
              grid(c)%forcing(f)=glbdata2(f,c)
              
           enddo
          endif
          
       else
!-----------------------------------------------------------------------
!     Linearly interpolate everything else
!-----------------------------------------------------------------------
          do c=1,gdi(iam)
             grid(c)%forcing(f)=glbdata1(f,c)*wt1+ & 
                  glbdata2(f,c)*wt2
          enddo
       endif
    enddo               
84  format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)
    return 
!EOC
  end subroutine time_interp_geos
