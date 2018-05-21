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
! !ROUTINE: ssib_output.F90
!
! !DESCRIPTION:
!  This subroutine sets up methods to write SSiB output 
!
! !REVISION HISTORY:
!  1 Mar 2004: Luis Gustavo G de Goncalves, SSiB in LIS
!  6 May 2004: David Mocko, made compatible with SiB-lings
!
! !INTERFACE:
      subroutine ssib_output
! !USES:
      use lisdrv_module, only : lis, tile, glbgindex
      use ssib_varder, only : ssibdrv
      use spmdMod, only : masterproc, npes
!EOP
      integer :: i 
      real :: var(lis%d%glbnch)
!BOC
      if (lis%o%wsingle.eq.1) then
!------------------------------------------------------------------
! Writes each output variables to a separate file
!------------------------------------------------------------------
         if (mod(lis%t%gmt,ssibdrv%writeintn).eq.0) then
            do i = 1,33
               call ssib_singlegather(i,var)
               if (masterproc) then
                  call ssib_singleout(lis, tile, glbgindex, var, i)
               endif
            enddo
            call ssib_totinit
         endif
      else
!------------------------------------------------------------------
! Writes bundled output
!------------------------------------------------------------------
         if (mod(lis%t%gmt,ssibdrv%writeintn).eq.0) then
            if (npes.gt.1) then
               call ssib_gather
            endif
            if (masterproc) then
               call ssib_almaout()
            endif
            call ssib_totinit
         endif
      endif
      return
!EOC
      end subroutine ssib_output

