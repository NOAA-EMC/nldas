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
! !ROUTINE: clmlairead.F90:
!
! !DESCRIPTION:
!  This program reads in AVHRR LAI data for CLM
!
! !REVISION HISTORY:
!  27 Nov 2001: Jon Gottschalck; Initial code
!  20 Feb 2002: Jon Gottschalck; Modified to use for 1/4 and 2x2.5 using 1/8 degree monthly data
!  01 Oct 2002: Jon Gottschalck; Modified to add MODIS LAI data
! 
! !INTERFACE: 
subroutine clm2lairead ()
! !USES:
  use lisdrv_module, only : lis
  use clm_varder
!EOP
  implicit none
  integer :: t
  real(r8) :: lai1(lis%d%nch),lai2(lis%d%nch)
  real(r8) :: sai1(lis%d%nch),sai2(lis%d%nch)
  real :: wt1, wt2

  call readlai(lis%p%lai, lai1, lai2, wt1, wt2)
  call readsai(lis%p%lai, sai1, sai2, wt1, wt2)

  if(lis%p%laiflag.eq.1) then 
     clm(1:lis%d%nch)%lai1 = lai1
     clm(1:lis%d%nch)%lai2 = lai2
  endif

  if(lis%p%saiflag.eq.1) then 
     clm(1:lis%d%nch)%sai1 = sai1
     clm(1:lis%d%nch)%sai2 = sai2
  endif

   do t=1,lis%d%nch
      clm(t)%tlai = wt1 * clm(t)%lai1 + wt2 * clm(t)%lai2
      clm(t)%tsai = wt1 * clm(t)%sai1 + wt2 * clm(t)%sai2

      if ( lis%p%lai == 2 ) then ! for UMD
         if (clm(t)%itypveg .eq. 12 ) then 
            clm(t)%tlai=0.0
            clm(t)%tsai=0.0
            clm(t)%htop=0.0
            clm(t)%hbot=0.0  
         endif
         if (clm(t)%itypveg .eq. 13 ) then 
            clm(t)%tlai=0.0
            clm(t)%tsai=0.0
         endif
      elseif ( lis%p%lai == 4 ) then ! for IGBP
         if (clm(t)%itypveg == 11 .or. clm(t)%itypveg >= 15 ) then
            clm(t)%tlai = 0.0
            clm(t)%htop = 0.0
            clm(t)%hbot = 0.0
         endif
      else
         call lis_log_msg('MSG: clm2lairead -- Unrecognized lis%p%lai value')
      endif
   enddo

end subroutine clm2lairead
