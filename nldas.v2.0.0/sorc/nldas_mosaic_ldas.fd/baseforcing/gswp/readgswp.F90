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
! !ROUTINE: readgswp.F90
!  
! !DESCRIPTION: 
!
!  Reads GSWP forcing.  
!
!    TIME1 = most recent past data\\
!    TIME2 = nearest future data \\
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! !REVISION HISTORY:
!
! 20Feb2004; Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine readgswp(order,yr,mo,da,hr,mn,ss)
! !USES:
#if ( defined USE_NETCDF )
  use netcdf
#endif
  use lisdrv_module, only : lis,gindex
  use baseforcing_module, only: glbdata1,glbdata2,glbdata3
  use gswp_module, only : getgswp_timeindex
  use gswpdomain_module, only : gswpdrv
!EOP
  implicit none
  integer :: yr,mo,da,hr,mn,ss
  integer :: tmp_yr, tmp_mo
  integer :: index,c,r,order
  character*8 :: cyear,cmo
  character*8 :: cslat,cnlat
  character*8 :: cwlon,celon
  character*8 :: c_index
  real :: tmp
  character*100 :: ffile
  real :: fvar(lis%d%lnc,lis%d%lnr)
  real :: fvar1(lis%d%glbnch)
  integer :: tid,qid,pid,wid,swid,lwid,rid,sid
  integer :: ncid,status
  integer :: cindex,rindex
  logical :: file_exis

  print*,'readgswp1',yr,mo,da,hr,mn,ss
  call getgswp_timeindex(yr,mo,da,hr,mn,ss,index)
  print*,'readgswp2',yr,mo,da,hr,mn,ss,index
#if ( defined GSWP_OPENDAP )
  tmp = lis%d%gridDesc(4)/1000.0
  write(cslat, '(f8.2)') tmp
  tmp = lis%d%gridDesc(7)/1000.0
  write(cnlat, '(f8.2)') tmp
  tmp = lis%d%gridDesc(5)/1000.0
  write(cwlon, '(f8.2)') tmp
  tmp = lis%d%gridDesc(8)/1000.0
  write(celon, '(f8.2)') tmp
  write(c_index,'(i8)') index
  
  print*,'MSG: GSWP forcing -- Reading tair '
  call system("./gswp_scripts/gettair.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp_scripts/tair.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ! Note:  GSWP forcing data are written into monthly files with
  ! a frequency of 3 hours.  The first entry corresponds to 03:00:00
  ! of the first day of the given month.  The last entry corresponds
  ! to 00:00:00 of the first day of the next month.
  !
  ! E.g.; for Tair_cru198207.nc, the data run from 1982-07-01T03:00:00
  ! through 1982-08-01T00:00:00, inclusive.
  !
  ! So, when you are at hour 0 on the first day of the month,
  ! you need to open the file corresponding to the previous month.
  if ( da == 1 .and. hr == 0 .and. mn == 0 .and. ss == 0 ) then
     if ( mo == 1 ) then
        tmp_mo = 12
        tmp_yr = yr - 1
     else
        tmp_mo = mo - 1
        tmp_yr = yr
     endif
  else
     tmp_mo = mo
     tmp_yr = yr
  endif
  print*, tmp_yr,tmp_mo,da,hr,mn,ss,index
  write(cyear, '(I4)') tmp_yr
  write(cmo, '(I2.2)') tmp_mo
  ffile = trim(gswpdrv%tair)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  print*,'MSG: GSWP forcing -- Reading tair ',ffile

#if ( defined USE_NETCDF )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Tair",tid)
  status = nf90_get_var(ncid, tid,fvar1,&
       start = (/1,index/),&
       count = (/lis%d%glbnch,1/))
  status = nf90_close(ncid)
#endif

  fvar = -9999.0
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if(gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = fvar1(gindex(cindex,rindex))
        endif
     enddo
  enddo
#endif
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then
           if ( order == 1 ) then 
              glbdata1(1,gindex(c,r)) = fvar(c,r)
           elseif ( order == 2 ) then 
              glbdata2(1,gindex(c,r)) = fvar(c,r)
           else 
              glbdata3(1,gindex(c,r)) = fvar(c,r)
           endif
        endif
     enddo
  enddo

#if ( defined GSWP_OPENDAP )
  print*,'MSG: GSWP forcing -- Reading qair '
  call system("./gswp_scripts/getqair.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp_scripts/qair.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswpdrv%qair)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  print*,'MSG: GSWP forcing -- Reading qair ',ffile

#if ( defined USE_NETCDF )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Qair",qid)
  status = nf90_get_var(ncid, qid,fvar1,&
       start = (/1,index/),&
       count = (/lis%d%glbnch,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if(gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = fvar1(gindex(cindex,rindex))
        endif
     enddo
  enddo
#endif
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then
           if ( order == 1 ) then 
              glbdata1(2,gindex(c,r)) = fvar(c,r)
           elseif ( order == 2 ) then 
              glbdata2(2,gindex(c,r)) = fvar(c,r)
           else 
              glbdata3(2,gindex(c,r)) = fvar(c,r)
           endif
        endif
     enddo
  enddo

#if ( defined GSWP_OPENDAP )
  print*,'MSG: GSWP forcing -- Reading swdown '
  call system("./gswp_scripts/getswdown.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp_scripts/swdown.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswpdrv%swdown)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  print*,'MSG: GSWP forcing -- Reading swdown ',ffile

#if ( defined USE_NETCDF )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "SWdown",swid)
  status = nf90_get_var(ncid, swid,fvar1,&
       start = (/1,index/),&
       count = (/lis%d%glbnch,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if(gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = fvar1(gindex(cindex,rindex))
        endif
     enddo
  enddo
#endif
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then
           if ( order == 1 ) then 
              glbdata1(3,gindex(c,r)) = fvar(c,r)
           elseif ( order == 2 ) then 
              glbdata2(3,gindex(c,r)) = fvar(c,r)
           else 
              glbdata3(3,gindex(c,r)) = fvar(c,r)
           endif
        endif
     enddo
  enddo

#if ( defined GSWP_OPENDAP )
  print*,'MSG: GSWP forcing -- Reading lwdown '
  call system("./gswp_scripts/getlwdown.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp_scripts/lwdown.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswpdrv%lwdown)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  print*,'MSG: GSWP forcing -- Reading LWdown ',ffile

#if ( defined USE_NETCDF )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "LWdown",lwid)
  status = nf90_get_var(ncid, lwid,fvar1,&
       start = (/1,index/),&
       count = (/lis%d%glbnch,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if(gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = fvar1(gindex(cindex,rindex))
        endif
     enddo
  enddo
#endif
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then
           if ( order == 1 ) then 
              glbdata1(4,gindex(c,r)) = fvar(c,r)
           elseif ( order == 2 ) then 
              glbdata2(4,gindex(c,r)) = fvar(c,r)
           else 
              glbdata3(4,gindex(c,r)) = fvar(c,r)
           endif
        endif
     enddo
  enddo

#if ( defined GSWP_OPENDAP )
  print*,'MSG: GSWP forcing -- Reading Wind '
  call system("./gswp_scripts/getwind.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp_scripts/wind.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswpdrv%wind)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  print*,'MSG: GSWP forcing -- Reading Wind ',ffile

#if ( defined USE_NETCDF )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Wind",wid)
  status = nf90_get_var(ncid, wid,fvar1,&
       start = (/1,index/),&
       count = (/lis%d%glbnch,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if(gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = fvar1(gindex(cindex,rindex))
        endif
     enddo
  enddo
#endif
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then
           if ( order == 1 ) then 
              glbdata1(5,gindex(c,r)) = fvar(c,r)
           elseif ( order == 2 ) then 
              glbdata2(5,gindex(c,r)) = fvar(c,r)
           else 
              glbdata3(5,gindex(c,r)) = fvar(c,r)
           endif
        endif
     enddo
  enddo


#if ( defined GSWP_OPENDAP )
  print*,'MSG: GSWP forcing -- Reading Pressure '
  call system("./gswp_scripts/getpsurf.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp_scripts/psurf.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswpdrv%psurf)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  print*,'MSG: GSWP forcing -- Reading psurf ',ffile

#if ( defined USE_NETCDF )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "PSurf",pid)
  status = nf90_get_var(ncid, pid,fvar1,&
       start = (/1,index/),&
       count = (/lis%d%glbnch,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if(gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = fvar1(gindex(cindex,rindex))
        endif
     enddo
  enddo
#endif
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then
           if ( order == 1 ) then 
              glbdata1(7,gindex(c,r)) = fvar(c,r)
           elseif ( order == 2 ) then 
              glbdata2(7,gindex(c,r)) = fvar(c,r)
           else 
              glbdata3(7,gindex(c,r)) = fvar(c,r)
           endif
        endif
     enddo
  enddo

#if ( defined GSWP_OPENDAP )
  print*,'MSG: GSWP forcing -- Reading Rainf '
  call system("./gswp_scripts/getrainf.sh "// & 
       cslat//" "//cnlat//" "//cwlon//" "//celon//" "//c_index)
  ffile = './gswp_scripts/rainf.bin'
  open(30,file=ffile,form='unformatted',status='old')
  read(30) fvar
  close(30)
#else
  ffile = trim(gswpdrv%rainf)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  print*,'MSG: GSWP forcing -- Reading Rainf ',ffile

#if ( defined USE_NETCDF )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Rainf",rid)
  status = nf90_get_var(ncid, rid,fvar1,&
       start = (/1,index/),&
       count = (/lis%d%glbnch,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if(gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = fvar1(gindex(cindex,rindex))
        endif
     enddo
  enddo
#endif
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then
           if ( order == 1 ) then 
              glbdata1(8,gindex(c,r)) = fvar(c,r)
           elseif ( order == 2 ) then 
              glbdata2(8,gindex(c,r)) = fvar(c,r)
           else 
              glbdata3(8,gindex(c,r)) = fvar(c,r)
           endif
        endif
     enddo
  enddo

#if ( defined GSWP_OPENDAP )
#else
  ffile = trim(gswpdrv%snowf)//trim(adjustl(cyear))//trim(adjustl(cmo))//".nc"
  print*,'MSG: GSWP forcing -- Reading Snowf ',ffile

#if ( defined USE_NETCDF )
  status = nf90_open(path=ffile,mode= nf90_nowrite,&
       ncid = ncid)
  status = nf90_inq_varid(ncid, "Snowf",sid)
  status = nf90_get_var(ncid, sid,fvar1,&
       start = (/1,index/),&
       count = (/lis%d%glbnch,1/))
  status = nf90_close(ncid)
#endif
  fvar = -9999.0
  do c=1,lis%d%lnc
     do r=1,lis%d%lnr
        rindex = 150-r+1
        cindex = c
        if(gindex(cindex,rindex).ne.-1) then 
           fvar(cindex,rindex) = fvar1(gindex(cindex,rindex))
        endif
     enddo
  enddo
#endif
  
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then
           if ( order == 1 ) then 
              glbdata1(9,gindex(c,r)) = fvar(c,r)
           elseif ( order == 2 ) then 
              glbdata2(9,gindex(c,r)) = fvar(c,r)
           else 
              glbdata3(9,gindex(c,r)) = fvar(c,r)
           endif
        endif
     enddo
  enddo
        
end subroutine readgswp
