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
! !ROUTINE: vic_setup.F90
!
! !DESCRIPTION:
! Completes the setup routines for VIC
!
! !REVISION HISTORY:
! 
! 29 Mar 2004; Sujay Kumar : Initial Specification
! 
! !INTERFACE:
subroutine read_parammap()
!EOP

  use lisdrv_module, only : lis,tile,gindex
  use vic_varder, only : vicdrv

  implicit none

  real :: lat, lon
  integer :: r, c, cindex, rindex
  integer :: line,line1,line2,glnc,glnr
  integer :: ios1

  real :: Dst(lis%d%lnc,lis%d%lnr)
  real :: Ds(lis%d%nch)

  real :: Dsmaxt(lis%d%lnc,lis%d%lnr)
  real :: Dsmax(lis%d%nch)

  real :: Wst(lis%d%lnc,lis%d%lnr)
  real :: Ws(lis%d%nch)

  real :: b_infiltt(lis%d%lnc,lis%d%lnr)
  real :: b_infilt(lis%d%nch)

  real :: depth1t(lis%d%lnc,lis%d%lnr)
  real :: depth2t(lis%d%lnc,lis%d%lnr)
  real :: depth3t(lis%d%lnc,lis%d%lnr)

  real :: depth1(lis%d%nch)
  real :: depth2(lis%d%nch)
  real :: depth3(lis%d%nch)

  line1 = (lis%d%gridDesc(4)-lis%d%gridDesc(44))/lis%d%gridDesc(49)+1
  line2 = (lis%d%gridDesc(5)-lis%d%gridDesc(45))/lis%d%gridDesc(50)+1

  print*, 'MSG: read_parammap -- Reading VIC parameter files'

  open(15,file=vicdrv%vic_dsmapfile,status='old',form='unformatted',&
          access='direct',recl=4,iostat=ios1)
  open(16,file=vicdrv%vic_dsmaxmapfile,status='old',form='unformatted',&
          access='direct',recl=4,iostat=ios1)
  open(17,file=vicdrv%vic_wsmapfile,status='old',form='unformatted',&
          access='direct',recl=4,iostat=ios1)
  open(18,file=vicdrv%vic_infiltmapfile,status='old',form='unformatted',&
          access='direct',recl=4,iostat=ios1)
  open(19,file=vicdrv%vic_depth1mapfile,status='old',form='unformatted',&
          access='direct',recl=4,iostat=ios1)
  open(20,file=vicdrv%vic_depth2mapfile,status='old',form='unformatted',&
          access='direct',recl=4,iostat=ios1)
  open(21,file=vicdrv%vic_depth3mapfile,status='old',form='unformatted',&
          access='direct',recl=4,iostat=ios1)

  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        glnc = line2+c-1
        glnr = line1+r-1
        line = (glnr-1)*lis%d%gridDesc(42)+glnc
!        print*,'DBG: read_parammap -- c,r,line1,line2,glnc,glnr,line', &
!               c,r,line1,line2,glnc,glnr,line
        read(15,rec=line) Dst(c,r)
        read(16,rec=line) Dsmaxt(c,r)
        read(17,rec=line) Wst(c,r)
        read(18,rec=line) b_infiltt(c,r)
        read(19,rec=line) depth1t(c,r)
        read(20,rec=line) depth2t(c,r)
        read(21,rec=line) depth3t(c,r)
     enddo
  enddo

  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)

! convert to tile array
  do r=1,lis%d%lnr
     do c=1,lis%d%lnc
        if(gindex(c,r).ne.-1) then 
           Ds(gindex(c,r)) = Dst(c,r)
           Dsmax(gindex(c,r)) = Dsmaxt(c,r)
           Ws(gindex(c,r)) = Wst(c,r)
           b_infilt(gindex(c,r)) = b_infiltt(c,r)
           depth1(gindex(c,r)) = depth1t(c,r)
           depth2(gindex(c,r)) = depth2t(c,r)
           depth3(gindex(c,r)) = depth3t(c,r) 
        endif
     enddo
  enddo
#if 0 
  do c=1,lis%d%nch
     if(nint(grid(c)%lat*1000).ge.lis%d%gridDesc(4).and. & 
          nint(grid(c)%lat*1000).le.lis%d%gridDesc(7).and. & 
          nint(grid(c)%lon*1000).ge.lis%d%gridDesc(5).and. & 
          nint(grid(c)%lon*1000).le.lis%d%gridDesc(8)) then
        rindex = tile(c)%row - (lis%d%gridDesc(4)-lis%d%gridDesc(44)) &
             /lis%d%gridDesc(9)
        cindex = tile(c)%col - (lis%d%gridDesc(5)-lis%d%gridDesc(45)) &
             /lis%d%gridDesc(10)
        if(Dst(cindex,rindex).ne.-9999.0) then
           Ds(c) = Dst(cindex, rindex)
        endif
        if(Dsmaxt(cindex,rindex).ne.-9999.0) then
           Dsmax(c) = Dsmaxt(cindex, rindex)
        endif
        if(Wst(cindex,rindex).ne.-9999.0) then
           Ws(c) = Wst(cindex, rindex)
        endif
        if(b_infiltt(cindex,rindex).ne.-9999.0) then
           b_infilt(c) = b_infiltt(cindex, rindex)
        endif

        if(depth1t(cindex,rindex).ne.-9999.0) then
           depth1(c) = depth1t(cindex, rindex)
        endif
        if(depth2t(cindex,rindex).ne.-9999.0) then
           depth2(c) = depth2t(cindex, rindex)
        endif
        if(depth3t(cindex,rindex).ne.-9999.0) then
           depth3(c) = depth3t(cindex, rindex)
        endif
     end if
  end do
#endif
  call set_parammap(Ds, Dsmax, Ws, b_infilt, depth1, depth2, depth3, & 
                    lis%d%nch, vicdrv%vic_nlayer)

end subroutine read_parammap
