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
!BOP
!
! !ROUTINE: noah_almaout.F90
!
! !DESCRIPTION:  
!  LIS NOAH data writer:  Binary and stat files in ALMA convention
!
! !REVISION HISTORY:
!  4 Nov. 1999: Jon Radakovich; Initial Code
! 28 Apr. 2002: Kristi Arsenault; Added NOAH LSM to LDAS
! 15 Jun  2003: Sujay Kumar; ALMA version 
! 
! !INTERFACE:
#include "misc.h"
subroutine noah_almaout
! !USES:
#if ( defined USE_NETCDF )
  use netcdf
#endif
  use lisdrv_module, only : lis
  use noah_varder      ! NOAH-specific variables
  use lis_openfileMod, only : create_output_directory, &
                              create_output_filename,  &
                              create_stats_filename
  
  implicit none 
! !ARGUMENTS:
!EOP
  integer           :: iret, ftn
  character(len=80) :: filengb
  character(len=80) :: name
  logical           :: file_exists
!BOC
!-------------------------------------------------------------------------
! Test to see if output writing interval has been reached
!-------------------------------------------------------------------------
  if(mod(lis%t%gmt,noahdrv%writeintn).eq.0)then
     noahdrv%numout=noahdrv%numout+1    
     call create_output_directory()
     call create_output_filename(filengb, model_name='Noah', &
                                 writeint=noahdrv%writeintn)
!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
     if(noahdrv%noahopen.eq.0)then
        call create_stats_filename(name,'Noahstats.dat')
        if(lis%o%startcode.eq.1)then
           open(65,file=name,form='formatted',status='unknown', & 
                position='append')
        else
           open(65,file=name,form='formatted',status='replace')       
        endif
        noahdrv%noahopen=1
     endif
          
       write(65,996)'       Statistical Summary of Noah output for:  ', & 
            lis%t%mo,'/',lis%t%da,'/',lis%t%yr,lis%t%hr,':',lis%t%mn,':',lis%t%ss
996    format(a47,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
       write(65,*)
       write(65,997)
997    format(t27,'Mean',t41,'Stdev',t56,'Min',t70,'Max')
       ftn = 58

!-----------------------------------------------------------------------
! Write Binary Output 
!-----------------------------------------------------------------------
       if(lis%o%wout.eq.1) then
          open(ftn,file=filengb,form='unformatted')
          call noah_binout(ftn)
          close(58)
!-----------------------------------------------------------------------
! Write Grib Output 
!-----------------------------------------------------------------------
       elseif(lis%o%wout.eq.2) then 
          call baopen (ftn,filengb, iret)
          call noah_gribout(ftn)
          call baclose(ftn,iret)
       elseif ( lis%o%wout == 3 ) then !netcdf
#if ( defined USE_NETCDF )
!-----------------------------------------------------------------------
! Write NetCDF Output 
!-----------------------------------------------------------------------
          iret = nf90_create(path=trim(filengb),cmode=nf90_clobber,ncid=ftn)
          call noah_netcdfout(.false.,ftn)
          iret = nf90_close(ftn)
#endif
       elseif ( lis%o%wout == 4 ) then !netcdf
#if ( defined USE_NETCDF )
          inquire (file=trim(filengb), exist=file_exists)
          if(file_exists) then
             iret = nf90_open(path=trim(filengb), mode=nf90_write, ncid=ftn)
             call noah_netcdfout(file_exists, ftn)
             iret = nf90_close(ftn)
          else
             noahdrv%numout = 1
             iret = nf90_create(path=trim(filengb), cmode=nf90_clobber, &
                                ncid=ftn)
             call noah_netcdfout(file_exists, ftn)
             iret = nf90_close(ftn)
          endif
#endif
       endif
       call noah_writestats(65)
       noah%count=0  !reset counters
       write(65,*)
       write(65,*)
    endif
!EOC
  end subroutine noah_almaout
  
