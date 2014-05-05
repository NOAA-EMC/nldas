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
! !ROUTINE: time_interp_gswp.F90
!
! !DESCRIPTION:
!  Opens, reads, and interpolates GSWP forcing.  
!
!    TIME1 = most recent past data\\
!    TIME2 = nearest future data \\
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! \subsection{Core Functions of time$_-$interp$_-$gswp}
!  \begin{description}
!  \item[zterp]  
!   Performs zenith angle-based temporal interpolation
!  \end{description}
!
! !REVISION HISTORY:
! 20Feb2004; Sujay Kumar : Initial Specification
! !INTERFACE:
subroutine time_interp_gswp()
! !USES:
  use lisdrv_module, only : lis, grid    
  use baseforcing_module, only : glbdata1, glbdata2, glbdata3
  use time_manager
  use grid_spmdMod
  use spmdMod
  use gswpdomain_module, only :gswpdrv
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
  integer :: madtt
  integer :: sub_dt = 1
  character(len=1), dimension(10) :: trp_flag = &
                                     (/'I','I','L','L','I','X','I','L','L','X'/)

!BOC
  if ( masterproc ) then 
     if ( get_nstep(lis%t) == 0 ) then
        lis%f%nforce = lis%f%nmif
     else
        lis%f%nforce = lis%f%nf
     endif
  endif

  zdoy=lis%t%doy
  btime=gswpdrv%gswptime1
  call time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=gswpdrv%gswptime2
  call time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)

  if ( lis%d%domain == 1 ) then
   !-----------------------------------------------------------------------
   !  Interpolate Data in Time
   !-----------------------------------------------------------------------
     wt1=(gswpdrv%gswptime2-lis%t%time)/ & 
          (gswpdrv%gswptime2-gswpdrv%gswptime1)
     wt2=1.0-wt1
     
     do f=1,lis%f%nforce
        if(f.eq.3) then     
           if (lis%f%shortflag.eq.2) then
   !-----------------------------------------------------------------------
   ! Got Time Averaged SW
   !-----------------------------------------------------------------------
              do c=1,gdi(iam)
                 zdoy=lis%t%doy
                 call zterp(0,grid(c)%lat,grid(c)%lon, & 
                      gmt1,gmt2,lis%t%gmt,zdoy, & 
                      zw1,zw2,czb,cze,czm,lis)
                 grid(c)%forcing(f)=glbdata2(f,c)*zw1
                 
                 if ((grid(c)%forcing(f).ne.lis%d%udef).and. & 
                      (grid(c)%forcing(f).lt.0) ) then
                    if (grid(c)%forcing(f) > -0.00001) then 
                       grid(c)%forcing(f) = 0.0 
                    else
                       print*,'ERR: time_interp_gswp -- Stopping because ', & 
                            'forcing not udef but lt0,'
                       print*,'ERR: time_interp_gswp -- ', & 
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
  else
     ! Note:  the GSWP temporal interpolation routine, finterp, returns
     ! an array of length madtt of interpolated values.
     ! Helin Wei's version only returns the value corresponding to the
     ! current time-step, sub_dt, w.r.t. the given 3-hour forcing interval.

     madtt = 3600 / lis%t%ts * 3 ! number of time-steps in a 
                                 ! 3-hourly forcing interval
     do f = 1, lis%f%nforce
        do c = 1, gdi(iam)
#if 0
           if ( f == 3 ) then ! shortwave
              call zterp(0,grid(c)%lat,grid(c)%lon, & 
                         gmt1,gmt2,lis%t%gmt,zdoy, & 
                         zw1,zw2,czb,cze,czm,lis)
              grid(c)%forcing(f)=glbdata2(f,c)*zw1

              if ( czm <= 0.1 .and. grid(c)%forcing(f) >= 400 ) then
                 print*,'ERR: time_interp_gswp -- ', &
                        'Triggered Helins arbitrary cut off', &
                        c,f,grid(c)%lat,grid(c)%lon,czm,grid(c)%forcing(f)
                 call finterp(0, trp_flag(f), &
                              0.0, glbdata1(f,c), glbdata2(f,c), glbdata3(f,c),&
                              madtt, sub_dt,         &
                              grid(c)%forcing(f))
                 print*,'ERR: time_interp_gswp -- ', &
                        'Replacing swnet with', &
                        c,f,grid(c)%lat,grid(c)%lon,grid(c)%forcing(f)
              endif
           else
              call finterp(0, trp_flag(f), &
                           0.0, glbdata1(f,c), glbdata2(f,c), glbdata3(f,c), &
                           madtt, sub_dt,         &
                           grid(c)%forcing(f))
           endif
#else
           call finterp(0, trp_flag(f), &
                        0.0, glbdata1(f,c), glbdata2(f,c), glbdata3(f,c), &
                        madtt, sub_dt,         &
                        grid(c)%forcing(f))
#endif
        enddo
     enddo

     sub_dt = sub_dt + 1
     if ( sub_dt > madtt ) then
        sub_dt = 1
     endif
  endif
84  format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)
    return 
!EOC
  end subroutine time_interp_gswp

