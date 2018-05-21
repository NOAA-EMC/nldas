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
! !ROUTINE: mos_writerst.F90
!
! !DESCRIPTION:
!  This program writes restart files for Mosaic.  This
!   includes all relevant water/energy storages, tile information,
!   and time information.  It also rectifies changes in the tile space.  
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  19 Jan 2001: Brian Cosgrove; Added CLOSE statement
!  15 Mar 2001: Jon Gottschalck; Added option 4 (initialize Mosaic with GDAS forcing)
!  17 Apr 2001: Jon Gottschalck; Modified option 4 to work with GEOS forcing
!  05 Sep 2001: Brian Cosgrove; Modified code to use Dag Lohmann's NOAA
!               initial conditions if necessary.  This is controlled with
!               local variable NOAAIC.  Normally set to 0 in this subroutine
!               but set to 1 if want to use Dag's NOAA IC's.  Changed output
!               directory structure, and commented out if-then check so that
!               directory is always made.
!  18 Sep 2001: Brian Cosgrove; Changed use of Dag's ICs so that third mosaic
!               soil layer is set to layer 2's soil moisture since layer
!               3 doesn't have a wilting point
!  04 Feb 2002: Jon Gottschalck; Added section to set all Koster tiles of veg
!               type 9 (land ice) to have tons of snow so it never melts since the soil/veg
!               parameters are set for bare soil (Koster -- personal communication).
!  05 Feb 2002: Brian Cosgrove; Added a few diagnostic vars to the NOAAIC=1 option
!  20 Nov 2002: Jon Radakovich; Updated for PSAS temperature assimilation and BC so
!               the forecast bias is included in the restart file.  Added conditionals
!               based on startcode=5 (Using spun-up restart so model time comes
!               from card file) and startcode=6 (Restarting a bias correction run).
! 12 Dec. 2002: Brian Cosgrove; Fixed usage of Wiltpoint variable.  Before,
!               Wiltpoint1 and Wiltpoint2 were used in calculation of
!               root zone soil moisture availability...now, only Wiltpoint2
!               is used since wiltpoint1 is not the correct wilting point
!               needed for the calculation.
!  14 Jan 2003: Urszula Jambor; Added deallocation statements near end of routine
!               and changed pointer variables to allocatable.
!  23 Jan 2003: Urszula Jambor; Switch index order of GEOS forcing
!               array.  Snow is 12, soil wetness is 13.
! RESTART FILE FORMAT(fortran sequential binary):
!  VCLASS,NC,NR,NCH !Veg class,no. columns, no. rows, no. tiles
!  MOS(NCH)%STATES  !Model States in Tile Space
! 
! !INTERFACE:
subroutine mos_writerst()
! !uses:
  use lisdrv_module, only : lis,tile
  USE mos_varder      ! Mosaic tile variables
  use time_manager
  use tile_spmdMod
!EOP
  IMPLICIT NONE      

!=== Local Variables =====================================================
  integer :: c,r,t,i,j,l,n,f ! loop counters

  character*80 filen,mkfyrmo
  character*1  fname(80),fbase(80),fsubs(80),fmkdir(80)
  character*1  ftime(10),fyrmodir(80)
!=== Temporary tile space transfer files (different than in lis_module)      
  real, allocatable :: tmptilen(:)
  CHARACTER (LEN=100) :: temp
!=== End Variable Definition =============================================
!BOC
  if(masterproc) then 
!-------------------------------------------------------------------------
! Restart Writing (2 files are written = active and archive)
!-------------------------------------------------------------------------
     if((lis%t%gmt.eq.(24-mosdrv%writeintm)) & 
          .or.lis%t%endtime.eq.1)then
        allocate(tmptilen(lis%d%nch))
        open(40,file=mosdrv%mos_rfile,form='unformatted') !Active archive restart
        !call timemgr_write_restart(40)
        !write(40) lis%p%vclass,lis%d%lnc,lis%d%lnr,lis%d%nch  !Veg class, no tiles       
        write(40) lis%d%lnc,lis%d%lnr,lis%d%nch  !Veg class, no tiles       
        write(40) mos%ct         !MOSAIC Canopy/Soil Temperature
	write(40) mos%qa         !MOSAIC Canopy Humidity
	write(40) mos%ics        !MOSAIC Interception Canopy Storage
	write(40) mos%snow       !MOSAIC Snow Depth
	write(40) mos%SoT        !MOSAIC Deep Soil Temperaure
        do l=1,3
          do t=1,lis%d%nch
            tmptilen(t)=mos(t)%SoWET(l)
          enddo
          write(40) tmptilen  !Mosaic Soil Wetness (3 layers)
        enddo
        close(40)   
        write(*,*)'Mosaic Active restart written: ',mosdrv%mos_rfile
        write(unit=temp,fmt='(i4,i2,i2,i2)') lis%t%yr,lis%t%mo, & 
            lis%t%da,lis%t%hr
        read(unit=temp,fmt='(10a1)')ftime
        do i=1,10
          if(ftime(i).eq.(' '))ftime(i)='0'
        enddo
        write(unit=temp,fmt='(a4,a3,a5,i4,a1,i4,i2,i2,a6,a3,a1)') & 
            '/EXP',lis%o%expcode,'/MOS/',lis%t%yr, & 
            '/',lis%t%yr,lis%t%mo, & 
            lis%t%da,'/LIS.E',lis%o%expcode,'.'
        read(unit=temp,fmt='(80a1)') (fname(i),i=1,35)
        do i=1,35
          if(fname(i).eq.(' '))fname(i)='0'
        enddo
        write(unit=temp,fmt='(a9)')'mkdir -p '
        read(unit=temp,fmt='(80a1)')(fmkdir(i),i=1,9)
        write(unit=temp,fmt='(a4,a3,a5,i4,a1,i4,i2,i2)') & 
            '/EXP',lis%o%expcode,'/MOS/', & 
            lis%t%yr,'/',lis%t%yr,lis%t%mo,lis%t%da
        read(unit=temp,fmt='(80a1)') (fyrmodir(i),i=1,25)
        do i=1,25
          if(fyrmodir(i).eq.(' '))fyrmodir(i)='0'
        enddo

        write(unit=temp,fmt='(a7)')'.MOSrst'
        read(unit=temp,fmt='(80a1)') (fsubs(i),i=1,7)

        write(unit=temp,fmt='(a40)') lis%o%odir                       
        read(unit=temp,fmt='(80a1)') (fbase(i),i=1,80)
        c=0
        do i=1,80
          if(fbase(i).eq.(' ').and.c.eq.0)c=i-1
        enddo
        write(unit=temp,fmt='(80a1)')(fbase(i),i=1,c),(fname(i),i=1,35), &
                        (ftime(i),i=1,10),(fsubs(i),i=1,7) 
        read(unit=temp,fmt='(a80)')filen
 
        write(unit=temp,fmt='(80a1)')(fmkdir(i),i=1,9),(fbase(i),i=1,c), &
         (fyrmodir(i),i=1,25)
        read(unit=temp,fmt='(a80)')mkfyrmo

!-------------------------------------------------------------------------
! Archive File Name Generation Complete
! Make the directories for the NOAH restart file
!-------------------------------------------------------------------------
        CALL SYSTEM(MKFYRMO)
!-------------------------------------------------------------------------
! Archive File Name Generation Complete
!-------------------------------------------------------------------------
        open(40,file=filen,status='unknown',form='unformatted')
        write(40) lis%d%lnc,lis%d%lnr,lis%d%nch  !veg class, no tiles       
        !write(40) lis%p%vclass,lis%d%lnc,lis%d%lnr,lis%d%nch  !veg class, no tiles       
        call timemgr_write_restart(40)
        write(40) mos%ct         !MOSAIC Canopy/Soil Temperature
	write(40) mos%qa         !MOSAIC Canopy Humidity
	write(40) mos%ics        !MOSAIC Interception Canopy Storage
	write(40) mos%snow       !MOSAIC Snow Depth
	write(40) mos%SoT        !MOSAIC Deep Soil Temperaure
        do l=1,3
	  do t=1,lis%d%nch
	   tmptilen(t)=mos(t)%SoWET(l)
	  enddo
	  write(40) tmptilen !Mosaic Soil Wetness (3 layers)
        enddo
        close(40)
        
        write(*,*)'mosaic archive restart written: ',filen
        deallocate(tmptilen)
      endif 
   endif
   return
!EOC
 end subroutine mos_writerst








