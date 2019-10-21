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
! !ROUTINE: time_interp_nldas2.F90
!
! !DESCRIPTION:
!  Opens, reads, and interpolates NLDAS2 forcing.  
!
!    TIME1 = most recent past data\\
!    TIME2 = nearest future data \\
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! \subsection{Core Functions of time$_-$interp$_-$nldas2}
!  \begin{description}
!  \item[zterp]  
!   Performs zenith angle-based temporal interpolation
!  \end{description}
!
! !REVISION HISTORY:
! 02Feb2004: Sujay Kumar; Initial Specification
! !INTERFACE:
subroutine time_interp_nldas2()
! !USES:
  use lisdrv_module, only : lis, grid
  use baseforcing_module, only : glbdata1, glbdata2
  use nldas2domain_module, only : nldas2drv
  use grid_spmdMod
  use time_manager, only : tick, time2date
  implicit none
!EOP
  integer :: zdoy,ier
  real :: zw1, zw2
  real :: czm, cze, czb
  real :: wt1, wt2,swt1,swt2
  real :: gmt1, gmt2
  integer :: c,f
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: time1,time2,btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss, tempbts

  lis%f%nforce = lis%f%nf
#if(defined SPMD)
  call MPI_BCAST(nldas2drv%nldas2time1,1,MPI_REAL8,0, & 
       MPI_COMM_WORLD,ier)
  call MPI_BCAST(nldas2drv%nldas2time2,1,MPI_REAL8,0, & 
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
  btime=nldas2drv%nldas2time1
  call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  tempbdoy=bdoy
  tempgmt1=gmt1

  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
   tempbts=000

  call tick(newtime1,tempbdoy,tempgmt1,& 
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn, & 
       tempbss,tempbts)
   
  btime=nldas2drv%nldas2time2
  call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
   tempbdoy=bdoy
   tempgmt2=gmt2
   tempbyr=byr
   tempbmo=bmo
   tempbda=bda
   tempbhr=bhr
   if (tempbhr.eq.24) tempbhr=0
   tempbmn=bmn
   tempbss=0
   tempbts=000

   call tick(newtime2,tempbdoy,tempgmt2,&
        tempbyr,tempbmo,tempbda,tempbhr,tempbmn,&
        tempbss,tempbts)
   
!=== Interpolate Data in time      
   wt1=(nldas2drv%nldas2time2-lis%t%time)/ & 
        (nldas2drv%nldas2time2-nldas2drv%nldas2time1)
   wt2=1.0-wt1
   swt1=(newtime2-lis%t%time)/(newtime2-newtime1)
   swt2=1.0-swt1
            
   do f=1,lis%f%nf
      if(f.eq.3)then
         do c=1,gdi(iam)
            zdoy=lis%t%doy
!	compute and apply zenith angle weights
            call zterp(1,grid(c)%lat,grid(c)%lon,&
                 gmt1,gmt2,lis%t%gmt,zdoy,zw1,zw2,czb,cze,czm,lis)
!            print*, c, glbdata1(f,c), glbdata2(f,c)
!            print*, c, grid(c)%lat, grid(c)%lon, gmt1, gmt2, lis%t%gmt, zdoy, zw1, zw2
            grid(c)%forcing(f)=glbdata1(f,c)*zw1+&
                 glbdata2(f,c)*zw2
	!if (c.eq.24373) print *,'forcing,glb1,glb2,zw1,zw2',grid(c)%forcing(f), &
             !glbdata1(f,c),glbdata2(f,c),zw1,zw2
!	In cases of small cos(zenith) angles, use linear weighting
!       to avoid overly large weights

            if((grid(c)%forcing(f).gt.glbdata1(f,c).and. & 
                 grid(c)%forcing(f).gt.glbdata2(f,c)).and. & 
                 (czb.lt.0.1.or.cze.lt.0.1))then
               grid(c)%forcing(f)=glbdata1(f,c)*swt1+ & 
                    glbdata2(f,c)*swt2
!        if (c.eq.24373) print *,'using swt1,swt2,forcing',swt1,swt2,grid(c)%forcing(f)
!	print *,'using swt1 and swt2',swt1,swt2
            endif

            if (grid(c)%forcing(f).gt.1367) then
               print*,'warning, sw radiation too high!!'
               print*,'it is',grid(c)%forcing(f)
               print*,'ncepdata1=',glbdata1(f,c)
               print*,'ncepdata2=',glbdata2(f,c)
               print*,'zw1=',zw1,'zw2=',zw2
               print*,'swt1=',swt1,'swt2=',swt2
               grid(c)%forcing(f)=glbdata1(f,c)*swt1+ & 
                    glbdata2(f,c)*swt2
            endif
               
            if ( (glbdata1(f,c).eq.lis%d%udef).and. & 
                 (glbdata2(f,c).eq.lis%d%udef) ) then
               grid(c)%forcing(f)=lis%d%udef
            endif
         enddo
      else if(f.eq.8.or.f.eq.9)then	!do block precipitation interpolation
         do c=1,gdi(iam)
            grid(c)%forcing(f)=glbdata2(f,c)	
         enddo
      else !linearly interpolate everything else
         do c=1,gdi(iam)
            grid(c)%forcing(f)=glbdata1(f,c)*wt1+ & 
                 glbdata2(f,c)*wt2
         enddo
      endif
   enddo
   
!=== ADJUST PRECIP TO VALUE PER SEC DATA SHOULD COME IN AS VALUE PER 
!=== TIME PERIOD ONE HOUR FOR NCEP, THREE HOUR FOR ETA
   do c=1,gdi(iam)
      grid(c)%forcing(8)=grid(c)%forcing(8)/(60.0*60.0)
      grid(c)%forcing(9)=grid(c)%forcing(9)/(60.0*60.0)
    enddo
!	print *,'at end of timeinterp, lw24373forcing=',grid(24373)%forcing(4)

!EOC
 end subroutine time_interp_nldas2
 
