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
!
! DESCRIPTION:
!  LDAS MOSAIC data writer.
!
! REVISION HISTORY:
! 4 Nov. 1999: Jon Radakovich; Initial Code
! 6 Nov. 1999: Paul Houser; Revision for Moasic Writing
!22 Aug. 2000: Brian Cosgrove; Modified code for output of
!              standard LDAS output variables.  Replaced
!              old output variables with new calls and
!              variables
!07 Sep. 2000: Brian Cosgrove; changed code so that downward lw
!              and sw values are output as averged values and not
!              instantaneous
!16 Nov. 2000: Brian Cosgrove; changed code to allow for output of
!              GRIB files using modified subroutine from NCEP.  Old
!              references to GRIB output based on LATS4D were
!              removed
!05 Sep. 2001: Brian Cosgrove; Added Close/Open statements, altered
!              output directory structure to match LDAS standards
!20 Sep. 2001: Urszula Jambor; Altered hdf and binary output directory 
!              structure to match NLDAS standards
!==========================================================================

subroutine mos_almaout()
 
! Declare modules and data structures
  use lisdrv_module, only : lis, tile ! LDAS non-model-specific 1-D variables
  use tile_module      ! LDAS non-model-specific tile variables
  use grid_module      ! LDAS non-model-specific grid variables
  use mos_varder
  IMPLICIT NONE

  
!=== Local variables =====================================================
  INTEGER :: T,C,R,M,I
  integer :: ftn,iret
  CHARACTER*80 MKFYRMO,FILENMT,FILENMG,CDIR,NAMET,NAMEG
  CHARACTER*90 FILENGB
  CHARACTER*80 MKFYRMO2
  CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
  CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13),FTIMEC(4)
  CHARACTER*1  FYRMODIR(25),FSUBFT(80)
  CHARACTER*1  FSUBFG(80),FTIMEB(10),FSUBGB(8)
  
  CHARACTER (LEN=100) :: FBINNAME
  CHARACTER (LEN=100) :: temp1

!=== Variables used for writing output in HDF format
    
  CHARACTER*40 FILE
  CHARACTER*80 NAME
!=== End Variable List ===================================================
!-------------------------------------------------------------------------
! Test to see if output writing interval has been reached
!-------------------------------------------------------------------------
  if(mod(lis%t%gmt,mosdrv%writeintm).eq.0)then
     mosdrv%numoutm=mosdrv%numoutm+1    
     write(unit=temp1,fmt='(i4,i2,i2)')lis%t%yr,lis%t%mo,lis%t%da
     read(unit=temp1,fmt='(8a1)') ftime
     do i=1,8
        if(ftime(i).eq.(' '))ftime(i)='0'
     enddo
     write(unit=temp1,fmt='(i4)')lis%t%yr
     read(unit=temp1,fmt='(8a1)')ftimec
     do i=1,4
        if(ftimec(i).eq.(' '))ftimec(i)='0'
     enddo
!     write(unit=temp1,fmt='(a6,i3,a1)') '/LIS.E',lis%o%expcode,'.'
!     read(unit=temp1,fmt='(80a1)') (fname(i),i=1,10)
!     do i=1,10
!        if(fname(i).eq.(' '))fname(i)='0'
!     enddo
     write(unit=temp1,fmt='(a40)') lis%o%odir
     read(unit=temp1,fmt='(40a1)') (fbase(i),i=1,40)
     c=0
     do i=1,40
        if(fbase(i).eq.(' ').and.c.eq.0)c=i-1
     enddo
!	write(*,'(a3)') lis%o%expcode
!	write (*,'(i4)') lis%t%yr
!	write (*,'(i2)') lis%t%mo
!	write (*,'(i2)') lis%t%da
     write(unit=temp1,fmt='(a4,a3,a5,i4,a1,i4,i2,i2)')'/EXP', & 
          lis%o%expcode,'/MOS/', & 
          lis%t%yr,'/',lis%t%yr,lis%t%mo,lis%t%da
     read(unit=temp1,fmt='(80a1)') (fyrmodir(i),i=1,25)
     do i=1,25
        if(fyrmodir(i).eq.(' '))fyrmodir(i)='0'
     enddo
     write(unit=temp1,fmt='(a9)')'mkdir -p '
     read(unit=temp1,fmt='(80a1)')(fmkdir(i),i=1,9)
     
     write(unit=temp1,fmt='(80a1)')(fmkdir(i),i=1,9),(fbase(i),i=1,c), & 
          (fyrmodir(i),i=1,25)
     read(unit=temp1,fmt='(a80)')mkfyrmo
     call system(mkfyrmo)
!----------------------------------------------------------------------
! Generate file name for BINARY output
!----------------------------------------------------------------------
     write(unit=fbinname, fmt='(i4,i2,i2,i2)') lis%t%yr,lis%t%mo, & 
          lis%t%da,lis%t%hr
     read(unit=fbinname,fmt='(10a1)') ftimeb
     do i=1,10
        if(ftimeb(i).eq.(' '))ftimeb(i)='0'
     enddo
     if(lis%o%wout.eq.1) then
        write(unit=fbinname,fmt='(a5)') '.gs4r'
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,5)
     elseif(lis%o%wout.eq.2) then 
        write(unit=fbinname,fmt='(a5)') '.grb '
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,5)
     elseif(lis%o%wout.eq.3) then         
        write(unit=fbinname,fmt='(a5)') '.nc  '
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,5)
     endif
        write(unit=fbinname,fmt='(69a1)')(fbase(i),i=1,c), & 
             (fyrmodir(i),i=1,25),'/',& 
             (ftimeb(i),i=1,10), & 
             (fsubgb(i),i=1,5)
        read(unit=fbinname,fmt='(a69)')filengb
!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
     if(mosdrv%mosopen.eq.0)then
        file='MOSstats.dat'
        call openfile(name,lis%o%odir,lis%o%expcode,file)
        if(lis%o%startcode.eq.1)then
           open(65,file=name,form='formatted',status='unknown', & 
                position='append')
        else
           open(65,file=name,form='formatted',status='replace')       
        endif
        mosdrv%mosopen=1
     endif
          
       write(65,996)'       Statistical Summary of MOS output for:  ', & 
            lis%t%mo,'/',lis%t%da,'/',lis%t%yr,lis%t%hr,':',lis%t%mn,':',lis%t%ss
996    format(a47,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
       write(65,*)
       write(65,997)
997    format(t27,'Mean',t41,'Stdev',t56,'Min',t70,'Max')
       ftn = 58
       if(lis%o%wout.eq.1) then
          open(ftn,file=filengb,form='unformatted')
          call mos_binout(ftn)
          close(58)
!-----------------------------------------------------------------------
! Write Grib Output
!-----------------------------------------------------------------------
       elseif(lis%o%wout.eq.2) then
          call baopen (ftn,filengb, iret)
          call mos_gribout(ftn)
          call baclose(ftn,iret)
       endif
       call mos_writestats(65)
       mos%count=0  !reset counters
       write(65,*)
       write(65,*)
    endif
!EOC
  end subroutine mos_almaout








