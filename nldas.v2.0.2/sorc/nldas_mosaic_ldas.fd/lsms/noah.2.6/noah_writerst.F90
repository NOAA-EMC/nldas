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
! !ROUTINE: noah_writerst.F90
!
! !DESCRIPTION:
!  This program writes restart files for NOAH.  This
!   includes all relevant water/energy storages, tile information,
!   and time information.  It also rectifies changes in the tile space.  
!
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  05 Sep 2001: Brian Cosgrove; Modified code to use Dag Lohmann's NOAA
!               initial conditions if necessary.  This is controlled with
!               local variable NOAAIC.  Normally set to 0 in this subroutine
!               but set to 1 if want to use Dag's NOAA IC's.  Changed output
!               directory structure, and commented out if-then check so that
!               directory is always made.
!  28 Apr 2002: Kristi Arsenault; Added NOAH LSM into LDAS
!  28 May 2002: Kristi Arsenault; For STARTCODE=4, corrected SNEQV values  
!                and put SMC, SH2O, STC limit for GDAS and GEOS forcing.
!  14 Jun 2003: Sujay Kumar , Separated the write restart from the original 
!                code
!  24 Aug 2004: James Geiger, Added support for GrADS-DODS server and
!                             corrected MPI base parallel restart writing
! RESTART FILE FORMAT(fortran sequential binary):
!  YR,MO,DA,HR,MN,SS,VCLASS,NCH !Restart time,Veg class,no.tiles, no.soil lay 
!  TILE(NCH)%COL        !Grid Col of Tile   
!  TILE(NCH)%ROW        !Grid Row of Tile
!  TILE(NCH)%FGRD       !Fraction of Grid covered by tile
!  TILE(NCH)%VEGT       !Vegetation Type of Tile
!  NOAH(NCH)%STATES     !Model States in Tile Space
! 
#include "misc.h"
! !INTERFACE:
subroutine noah_writerst()
! !uses:
  use lisdrv_module, only : lis,tile
  use time_manager
  use noah_varder      ! NOAH tile variables
  use lis_openfileMod, only : create_output_directory, &
                              create_restart_filename
!EOP
  implicit none       

!=== Local Variables =====================================================
  character(len=80) :: filen
!=== End Variable Definition =============================================
!BOC

   if ( ( lis%t%gmt == (24-noahdrv%writeintn) ) .or. &
        lis%t%endtime == 1 ) then
!
! All processes of a multi-process run call this subroutine
!
! If you are not writing bundled output, then noah_gather has not been
! called.  Gather the Noah data structure, unless you are using a GrADS-DODS
! server (OPENDAP).  In this case, the variables will be gathered one at a time.
!
#if ( ! defined OPENDAP )
     if ( lis%o%wsingle == 1 .and. npes > 1 ) then 
        call noah_gather()
     endif

! If you are using a GrADS-DODS server then every process 
! should continue, else only the master process.
     if ( masterproc ) then 
#endif

!-------------------------------------------------------------------------
! Restart Writing (2 files are written = active and archive)
!-------------------------------------------------------------------------


! Eventhough all processes may continue to this point, only the master
! process should perform any writing.

        if ( masterproc ) then
           call create_output_directory()
           call create_restart_filename(filen,'.Noahrst')
           open(40,file=filen,status='unknown',form='unformatted')
        endif

        call noah_dump_restart(40)

        if ( masterproc ) then
           close(40)   
           write(*,*)'MSG: noah_writerst -- archive restart written: ',filen
        endif
        
#if ( ! defined OPENDAP )
      endif 
#endif
   endif

   return
!EOC
end subroutine noah_writerst

!BOP
! 
! !ROUTINE: noah_dump_restart
!
! !DESCRIPTION:
!    This routine gathers the necessary restart variables and performs
!    the actual write statements to create the restart files.
!
! !REVISION HISTORY:
!  24 Aug 2004: James Geiger, Initial Specification
! !INTERFACE:
subroutine noah_dump_restart(ftn)

! !USES:
   use lisdrv_module, only : lis
   use noah_varder
   use tile_spmdMod
   use time_manager
!EOP

   implicit none

   integer, intent(in) :: ftn

   integer :: l,t,ierr
   integer :: alloc_size, local_size, tile_loop
   real, allocatable :: tmptilen(:)

! All processes may be calling this routine.
! The master process must allocate a temporary array that is large
! enough to hold all the land points.  The slave processes should only
! allocate a temporary array the same size as the number of land points
! they computed.
   if ( masterproc ) then
      alloc_size = lis%d%glbnch  ! all land points
   else
      alloc_size = di_array(iam) ! number of land points slave process processed
   endif
   local_size = di_array(iam) ! number of land points each process (master or
                              ! slave) processed

! For several do loops below, you must loop over the number of land points.
! If you are using a GrADS-DODS server, each process must loop over the number 
! of land the process processed.  Else only the master process is in this
! routine, and it must loop over all the land points.
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
   tile_loop = di_array(iam)
#else
   tile_loop = lis%d%glbnch
#endif

   allocate(tmptilen(alloc_size),stat=ierr)
   call check_error(ierr,'ERR: noah_dump_restart -- Error allocating tmptilen',iam)

   if ( masterproc ) then
      call timemgr_write_restart(ftn)

      write(ftn) lis%p%vclass,lis%d%lnc,lis%d%lnr,lis%d%glbnch  !Veg class, no tiles       
   endif

! If you are using a GrADS-DODS server, then each process must fill in
! a temporary array.  Then these temporary arrays must the gathered 
! into a master copy held by the master process.  Then the master process
! writes this array out.  Else the master process may simply write the noah
! data structure element.
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
   tmptilen(1:local_size) = noah(1:local_size)%t1
   call noah_gather_restart(tmptilen)
   if ( masterproc ) then
      write(ftn) tmptilen       !NOAH Skin Temperature (K)
   endif
#else
   write(ftn) noah%t1        !NOAH Skin Temperature (K)
#endif

#if ( ( defined OPENDAP ) && ( defined SPMD ) )
   tmptilen(1:local_size) = noah(1:local_size)%cmc
   call noah_gather_restart(tmptilen)
   if ( masterproc ) then
      write(ftn) tmptilen       !NOAH Canopy Water Content
   endif
#else
   write(ftn) noah%cmc       !NOAH Canopy Water Content
#endif

#if ( ( defined OPENDAP ) && ( defined SPMD ) )
   tmptilen(1:local_size) = noah(1:local_size)%snowh
   call noah_gather_restart(tmptilen)
   if ( masterproc ) then
      write(ftn) tmptilen       !NOAH Actual Snow Depth
   endif
#else
   write(ftn) noah%snowh     !NOAH Actual Snow Depth
#endif

#if ( ( defined OPENDAP ) && ( defined SPMD ) )
   tmptilen(1:local_size) = noah(1:local_size)%sneqv
   call noah_gather_restart(tmptilen)
   if ( masterproc ) then
      write(ftn) tmptilen       !NOAH Water Equivalent Snow Depth
   endif
#else
   write(ftn) noah%sneqv     !NOAH Water Equivalent Snow Depth
#endif

   do l=1,4
     do t=1,tile_loop
      tmptilen(t)=noah(t)%stc(l)
     enddo
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
      call noah_gather_restart(tmptilen)
#endif
     if ( masterproc ) then
        write(ftn) tmptilen  !NOAH Soil Temperature (4 layers)
     endif
   enddo

   do l=1,4
     do t=1,tile_loop
       tmptilen(t)=noah(t)%smc(l)
     enddo
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
      call noah_gather_restart(tmptilen)
#endif
     if ( masterproc ) then
        write(ftn) tmptilen  !NOAH Total Soil Moist. (4 layers)
     endif
   enddo

   do l=1,4
     do t=1,tile_loop
       tmptilen(t)=noah(t)%sh2o(l)
     enddo
#if ( ( defined OPENDAP ) && ( defined SPMD ) )
      call noah_gather_restart(tmptilen)
#endif
     if ( masterproc ) then
        write(ftn) tmptilen  !NOAH Liquid Soil Moist. (4 layers)
     endif
   enddo

#if ( ( defined OPENDAP ) && ( defined SPMD ) )
   tmptilen(1:local_size) = noah(1:local_size)%ch
   call noah_gather_restart(tmptilen)
   if ( masterproc ) then
      write(ftn) tmptilen     !NOAH Heat/Moisture Sfc Exchange Coef.
   endif
#else
   write(ftn) noah%ch      !NOAH Heat/Moisture Sfc Exchange Coef.
#endif

#if ( ( defined OPENDAP ) && ( defined SPMD ) )
   tmptilen(1:local_size) = noah(1:local_size)%cm
   call noah_gather_restart(tmptilen)
   if ( masterproc ) then
      write(ftn) tmptilen     !NOAH Momentum Sfc Exchange Coef.
   endif
#else
   write(ftn) noah%cm      !NOAH Momentum Sfc Exchange Coef.
#endif

   deallocate(tmptilen)

end subroutine noah_dump_restart

!BOP
! 
! !ROUTINE: noah_gather_restart
!
! !DESCRIPTION:
!    This routine is a wrapper for the mpi_gatherv subroutine.  It
!    performs any necessary data gathering for the noah_writerst routine.
!
! !REVISION HISTORY:
!  24 Aug 2004: James Geiger, Initial Specification
! !INTERFACE:
subroutine noah_gather_restart(tmptilen)

   use spmdMod
   use tile_spmdMod

   implicit none

   real, dimension(*), intent(inout) :: tmptilen
   integer :: ierr

#if ( defined SPMD )

   call MPI_GATHERV(tmptilen(1:di_array(iam)), di_array(iam),       & 
                    MPI_REAL, tmptilen, di_array, displs, MPI_REAL, & 
                    0, MPI_COMM_WORLD, ierr)
#endif

end subroutine noah_gather_restart
