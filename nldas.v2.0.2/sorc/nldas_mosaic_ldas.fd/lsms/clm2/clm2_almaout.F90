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
! !ROUTINE: clm2_almaout.F90
! 
! !DESCRIPTION:
!  CLM output writer.
!
! !REVISION HISTORY:
! 29 Oct. 1999: Jon Radakovich; Initial code
! 27 Sep. 2000: Brian Cosgrove; Major revisions to enable CLM to 
!               output ALMA/LDAS variables
! 27 Sep. 2000: Added abitrary root zone cutoff value of .05 so 
!               that root zone is considered only those levels
!               with a rooting value >= .05
! 19 Mar  2000: Cosgrove; Changed code to allow LDAS GRIB output.  Removed
!               calls to LATS4D that had been used for grib output and 
!               added a call to griboutclm which outputs GRIB CLM data.
! 05 Apr 2002: Jon Gottschalck; Modified code to work with CLM2 
! 15 Jun 2003; Sujay Kumar; ALMA version of the output routine
! 
! !INTERFACE: 
subroutine clm2_almaout ()
! !USES:
  use lisdrv_module, only : lis, tile
  use clm_varctl, only : clmdrv
  use lis_openfileMod, only : create_output_directory, &
                              create_output_filename,  &
                              create_stats_filename
#if ( defined USE_NETCDF )
  use netcdf
#endif
!EOP
  implicit none
!=== Local variables =====================================================
  integer           :: ftn, iret
  character(len=80) :: filengb
  character(len=80) :: name
  logical           :: file_exists
      
!=== End Variable List =========================================================

!=== Total arrays hold a running total of each output variable for time
!=== averaging, between output writes
!BOC
!-----------------------------------------------------------------------
! Test to see if output writing interval has been reached
!-----------------------------------------------------------------------
  if(mod(lis%t%gmt,clmdrv%writeintc2).eq.0)then
!-----------------------------------------------------------------------
! Generate directory structure and file names for CLM Output 
!-----------------------------------------------------------------------
     clmdrv%numout=clmdrv%numout+1    !Counts number of output times
     call create_output_directory()
     call create_output_filename(filengb, model_name='CLM', &
                                 writeint=clmdrv%writeintc2)
!-----------------------------------------------------------------------
! Write statistical output
!-----------------------------------------------------------------------
     if(clmdrv%clm2open.eq.0)then
        call create_stats_filename(name,'CLMstats.dat')
        if(lis%o%startcode.eq.1)then
           open(60,file=name,form='formatted',status='unknown', & 
                position='append')
        else
           open(60,file=name,form='formatted',status='replace')
        endif
        clmdrv%clm2open=1
     endif
     
     write(60,996)'       Statistical Summary of CLM Output for:  ', & 
          lis%t%mo,'/',lis%t%da,'/',lis%t%yr,lis%t%hr,':', & 
          lis%t%mn,':',lis%t%ss
996  format(a47,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
     write(60,*)
     write(60,997)
997  format(t26,'Mean',t40,'StDev',t54,'Min',t68,'Max')

     ftn = 57
!-----------------------------------------------------------------------
! Write Binary Output 
!-----------------------------------------------------------------------
     if(lis%o%wout.eq.1)then
        open(ftn,file=filengb,form='unformatted')
        call clm2_binout(ftn)
        close(ftn)
     elseif(lis%o%wout.eq.2)then
!-----------------------------------------------------------------------
! Write Grib Output 
!-----------------------------------------------------------------------
          call baopen (ftn,filengb, iret)
          call clm2_gribout(ftn)
          call baclose(ftn,iret)
     elseif(lis%o%wout.eq.3) then !netcdf
!-----------------------------------------------------------------------
! Write NetCDF Output 
!-----------------------------------------------------------------------
#if ( defined USE_NETCDF )
          iret = nf90_create(path=trim(filengb),cmode=nf90_clobber,ncid=ftn)
          call clm2_netcdfout(.false., ftn)
          iret = nf90_close(ftn)
#endif
       elseif ( lis%o%wout == 4 ) then !netcdf
!-----------------------------------------------------------------------
! Write GSWP-2 styly NetCDF Output 
!-----------------------------------------------------------------------
#if ( defined USE_NETCDF )
          inquire (file=trim(filengb), exist=file_exists)
          if(file_exists) then
             iret = nf90_open(path=trim(filengb), mode=nf90_write, ncid=ftn)
             call clm2_netcdfout(file_exists, ftn)
             iret = nf90_close(ftn)
          else
             clmdrv%numout = 1
             iret = nf90_create(path=trim(filengb), cmode=nf90_clobber, &
                                ncid=ftn)
             call clm2_netcdfout(file_exists, ftn)
             iret = nf90_close(ftn)
          endif
#endif
     endif
     call clm2_writestats(60)
  end if
!EOC
end subroutine clm2_almaout
