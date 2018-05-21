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
! !ROUTINE: vic_almaout.F90
!
! !DESCRIPTION:  
!  LIS VIC data writer:  Binary and stat files in ALMA convention
!
! !REVISION HISTORY:
!  4 Nov. 1999: Jon Radakovich; Initial Code
! 28 Apr. 2002: Kristi Arsenault; Added VIC LSM to LDAS
! 15 Jun  2003: Sujay Kumar; ALMA version 
! 
! !INTERFACE:
subroutine vic_almaout()
! !USES:
!  use netcdf
  use lisdrv_module, only : lis, grid, gindex, tile
  use vic_varder      ! VIC-specific variables
  use time_manager, only : get_nstep
  
  implicit none 
!EOP

  integer :: t,c,r,m,i,n,iret,ftn
  character*80 mkfyrmo,filenmt,filenmg,cdir,namet,nameg,filengb
  character*80 mkfyrmo2
  character*1  fname(80),fbase(40),fmkdir(80)
  character*1  ftime(8),fcd(3),frm(3),flats(13),ftimec(4)
  character*1  fyrmodir(26),fsubft(80)
  character*1  fsubfg(80),ftimeb(10),fsubgb(9)
  
  character (len=100) :: fbinname
  character (len=100) :: temp1

  integer,parameter :: nvarsg=29
  character*80 :: vname(nvarsg)
  
  character*40 file
  character*80 name
!BOC
!-------------------------------------------------------------------------
! Test to see if output writing interval has been reached
!-------------------------------------------------------------------------
  if(mod(lis%t%gmt,vicdrv%writeintvic).eq.0)then
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
     
     write(unit=temp1,fmt='(a4,i3,a5,i4,a1,i4,i2,i2)')'/EXP', & 
          lis%o%expcode,'/VIC/', & 
          lis%t%yr,'/',lis%t%yr,lis%t%mo,lis%t%da
     read(unit=temp1,fmt='(80a1)') (fyrmodir(i),i=1,25)
     do i=1,25
        if(fyrmodir(i).eq.(' '))fyrmodir(i)='0'
     enddo
     
     write(unit=temp1,fmt='(a9)')'mkdir -p '
     read(unit=temp1,fmt='(80a1)')(fmkdir(i),i=1,9)
     
     write(unit=temp1,fmt='(80a1)')(fmkdir(i),i=1,9),(fbase(i),i=1,c), & 
          (fyrmodir(i),i=1,26)
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
        if(lis%o%wtil.eq.1) then 
           write(unit=fbinname,fmt='(a8)') '.ls4r   '
        else
           write(unit=fbinname,fmt='(a8)') '.gs4r   '
        endif
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,8)
     elseif(lis%o%wout.eq.2) then 
        write(unit=fbinname,fmt='(a8)') '.VIC.grb'
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,8)
     elseif(lis%o%wout.eq.3) then         
        write(unit=fbinname,fmt='(a8)') '.VIC.nc '
        read(unit=fbinname,fmt='(80a1)') (fsubgb(i),i=1,8)
     endif
!     write(unit=fbinname,fmt='(80a1)')(fbase(i),i=1,c), & 
!          (fyrmodir(i),i=1,26), & 
!          (fname(i),i=1,10),(ftimeb(i),i=1,10), & 
!          (fsubgb(i),i=1,8)
     write(unit=fbinname,fmt='(80a1)')(fbase(i),i=1,c), & 
          (fyrmodir(i),i=1,25), '/',& 
          (ftimeb(i),i=1,10), & 
          (fsubgb(i),i=1,8)
     read(unit=fbinname,fmt='(a80)')filengb
!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
     if(vicdrv%vicopen.eq.0)then
        file='VICstats.dat'
        call openfile(name,lis%o%odir,lis%o%expcode,file)
        if(lis%o%startcode.eq.1)then
           open(65,file=name,form='formatted',status='unknown', & 
                position='append')
        else
           open(65,file=name,form='formatted',status='replace')       
        endif
        vicdrv%vicopen=1
     endif
          
     write(65,996)'       Statistical Summary of Vic output for:  ', & 
          lis%t%mo,'/',lis%t%da,'/',lis%t%yr,lis%t%hr,':',lis%t%mn,':',lis%t%ss
996    format(a47,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
       write(65,*)
       write(65,997)
997    format(t27,'Mean',t41,'Stdev',t56,'Min',t70,'Max')
       ftn = 73
       if(lis%o%wout.eq.1) then
          open(ftn,file=trim(filengb),form='unformatted')
          call vic_binout(ftn)
          close(73)
!       elseif(lis%o%wout.eq.3) then !netcdf
!          iret = nf90_create(path=trim(filengb),cmode=nf90_clobber,ncid=ftn)
!          call vic_netcdfout(ld,ftn)
!          iret = nf90_close(ftn)
       endif
       call vic_writestats(65)
       write(65,*)
       write(65,*)
    endif
!EOC
  end subroutine vic_almaout
  
