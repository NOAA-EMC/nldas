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
!
!	MODIFIED TO READ HELIN'S GDAS T62 FORCING
!	JESSE 20040404
!
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
  integer, parameter :: c1=2960 !jesse
  
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
#endif
  btime=gdasdrv%gdastime1
  call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=gdasdrv%gdastime2
  call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)

  wt1 = (gdasdrv%gdastime2-lis%t%time) / & 
       (gdasdrv%gdastime2-gdasdrv%gdastime1)
  wt2 = 1.0 - wt1

  WRITE(*,'(1X,A,3(A,F8.3))') "J---TIME_INTERP_GDAS",&
       " BTIME=",GMT1," ETIME=",GMT2," MTIME=",lis%t%gmt

  do f=1,lis%f%nforce
     if ( f == 3 ) then  !shortwave         
        do c = 1, gdi(iam)
           zdoy = lis%t%doy
!           print*, 'f0,f6',lis%f%F00_flag,lis%f%F06_flag
!           call zterp( 0, grid(c)%lat, grid(c)%lon, gmt1, gmt2, &	!SARAH 
           call zterp( 1, grid(c)%lat, grid(c)%lon, gmt1, gmt2, &	!HELIN
                lis%t%gmt,zdoy,zw1,zw2,czb,cze,czm,lis)
           if( c == c1 ) &
              WRITE(*,'(1X,A,3(A,F8.3))') "J---TIME_INTERP_GDAS",&
                   " CZB  =",czb," CZE  =",cze," CZM  =",czm
           
!           grid(c)%forcing(f) = zw1 * glbdata2(f,c)				!SARAH
           grid(c)%forcing(f) = zw1 * glbdata1(f,c) + zw2 * glbdata2(f,c)	!HELIN
 
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
!           print*, c, grid(c)%forcing(f),glbdata2(f,c)

        end do
!      else if ( (f==4) .or. (f==10) ) then 	!SARAH'S LWDN IS AVERAGE 
      else if (  f==10 )              then 	!HELIN'S LWDN IS INSTANTANEOUS 
        do c = 1, gdi(iam)
           if ( lis%f%F00_flag==1 ) then 
              grid(c)%forcing(f)=2*glbdata2(f,c) -glbdata1(f,c)
           else
              grid(c)%forcing(f) = glbdata2(f,c)
           end if
        end do
        
     else if ( (f==8) .or. (f==9 ) ) then
        do c = 1, gdi(iam)
!           if ( lis%f%F00_flag == 1) then 
!              if (2*glbdata2(f,c) >= glbdata1(f,c)) then
!                 grid(c)%forcing(f)=2*glbdata2(f,c)-glbdata1(f,c)
!              else
!                 grid(c)%forcing(f) = 0.0
!              end if
!           else
!              grid(c)%forcing(f) = glbdata2(f,c)
!           end if
            grid(c)%forcing(f) = glbdata1(f,c) !jesse 20040404 for Helin's forcing
        end do
     else 
        do c = 1, gdi(iam)
           grid(c)%forcing(f) = wt1 * glbdata1(f,c) +  & 
                wt2 *glbdata2(f,c)
        end do
 
    end if

!J- ON SCREEN OUTPUT, JESSE 20040505

        IF ( F==1 ) then
        WRITE(*,'(1X,A,A,I8,2(A,F8.3))') "J---TIME_INTERP_GDAS",&
             " GDI  =",gdi(iam),&
             " LAT1 =",grid(1)%lat,  " LAT2 =",grid(gdi(iam))%lat
        WRITE(*,'(1X,A,A,I8,2(A,F8.3))') "J---TIME_INTERP_GDAS",&
             " GDI  =",gdi(iam),&
             " LON1 =",grid(1)%lon,  " LON2 =",grid(gdi(iam))%lon
        WRITE(*,'(1X,A,A,I8,2(A,F8.3))') "J---TIME_INTERP_GDAS",&
             " GRID =",c1,&
             " LAT  =",grid(c1)%lat, " LON  =",grid(c1)%lon
        ENDIF
        IF ( F==3 ) &
        WRITE(*,'(1X,A,3(A,F8.3))') "J---TIME_INTERP_GDAS",&
             " SWDN1=",glbdata1(f,c1)," SWDN2=",glbdata2(f,c1),&
             " SWDN =",grid(c1)%forcing(f)
        IF ( F==4 ) &
        WRITE(*,'(1X,A,3(A,F8.3))') "J---TIME_INTERP_GDAS",&
             " LWDN1=",glbdata1(f,c1)," LWDN2=",glbdata2(f,c1),&
             " LWDN =",grid(c1)%forcing(f) 
        IF ( F==8 ) &
        WRITE(*,'(1X,A,3(A,F8.3))') "J---TIME_INTERP_GDAS",&
             " PRAT1=",glbdata1(f,c1)*1000000.,&
             " PRAT2=",glbdata2(f,c1)*1000000.,&
             " PRAT =",grid(c1)%forcing(f)*1000000.
!J-

  end do
  return
!EOC    
end subroutine time_interp_gdas
