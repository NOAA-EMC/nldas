!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! ldasdrv.f: 
!
! DESCRIPTION:
!  LDAS Driver: this is the driver for the LDAS system.  
!
! REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  4  Apr 2000: Jeffrey Walker; Added catchment model
!  11 Apr 2000: Brian Cosgrove; Added elevation correction code
!  6  Jun 2000: Jon Radakovich; Updated for new version of CLM
!  23 Feb 2001: Urszula Jambor; Added GEOS & GDAS forcing option in GLDAS
!  27 Feb 2001: Brian Cosgrove; Added calls to make use of catchment
!               forcing data (getcatch subroutine)
!   5 Mar 2001: Jon Radakovich; Update for CLM version 2
!  07 Mar 2001: Brian Cosgrove; Changed if then checks to allow for use
!               of NASA-LDAS data
!  15 Mar 2001: Jon Gottschalck; Allocated space for new forcing parameters
!               when option 4 is selected
!  23 Mar 2001: Jon Radakovich; Added PSAS temperature assimilation
!   6 Apr 2001: Matt Rodell; Include tile module in call to mosdynp
!  13 Apr 2001: Separated sibalb() into 1 module (sibalb), and 2 
!               subroutines (mapsib2umd, umd_sibalb).
!  17 Jul 2001: Jon Gottschalck; Updated to allow global observed NRL precipitation forcing      
!  05 Sep 2001: Brian Cosgrove; Removed call to pink2ldas and awip2ldas, 
!               uses GRIB interpolation package now
!  19 Sep 2001: Brian Cosgrove; Added call to getpar
!  27 Nov 2001: Jon Gottschalck; Added AVHRR LAI data
!  07 Dec 2001: Urszula Jambor, added call to define_gds subroutine      
!  28 Jan 2002: Jon Gottschalck; Add section to initialize CLM with GDAS model forcing
!  06 Feb 2002: Jon Gottschalck; Modifications to use Koster tile space
!  28 Apr 2002: Kristi Arsenault; Added NOAH LSM, v 2.5, to LDAS
!  30 Jul 2002: Jon Gottschalck; Added additional observed precipitation sources
!  01 Oct 2002: Jon Gottschalck; Added MODIS LAI data       
!  20 Nov 2002: Jon Radakovich; Added assimilation and bias correction for CLM
!               Updated assimilation and BC in MOS.  Added options for 
!               startcode=5 (Using spun-up restart so model time comes
!               from card file) and startcode=6 (Restarting a bias correction run).
!  05 Feb 2003: Jon Gottschalck; Added CLM LSM, v 2.0 to LDAS  
!  04 Mar 2003: Urszula Jambor; Added call to set_forcinggrid
!  03 Jun 2003: Jon Gottschalck; Added GEOS initialization for CLM1  
!  23 Jun 2003: Urszula Jambor; Added call to getecwmf     
!=========================================================================

!=== CLM2 pre-processor directive input files
#include <misc.h>
#include <preproc.h>
      program ldasdrv

      use precision 
      use ldas_module      ! LDAS non-model-specific 1-D variables
      use tile_module      ! LDAS non-model-specific tile variables
      use grid_module      ! LDAS non-model-specific grid variables
      use mos_module       ! Mosaic tile variables
      use sibalb_module    ! Mosaic SiB Albedo coefficients
      use drv_module       ! CLM1 driver variables
      use drv_tilemodule   ! CLM1 Tile-space variables
      use drv_gridmodule   ! CLM1 Grid-space variables
      use clm1type          ! CLM1 tile variables
      use cat_module       ! Catchment model variables
      use noah_module      ! NOAH model variables

      implicit none


      type (ldasdec)                   ldas              
      type (tiledec),       pointer :: tile(:)
      type (griddec),       pointer :: grid(:,:)   
      type (mosdec),        pointer :: mos(:)
      type (sibalbdec)                 sib
      type (drvdec)                 :: drv
      type (clm_tiledec),   pointer :: clm_tile(:)
      type (clm_griddec),   pointer :: clm_grid(:,:)
      type (clm11d),        pointer :: clm1(:)
      type (catdec),        pointer :: cat(:)
      type (noahdec),       pointer :: noah(:)

!=== Local Variables =====================================================

      integer :: c,r,t,i,j,iii,m,l    ! Loop counters
      integer :: ttr,ttc,nstep
      integer  :: yr                    !current year (0, ...)
      integer  :: mon                   !current month (1 -> 12)
      integer  :: day                   !current day (1 -> 31)
      integer  :: ncsec                 !current time of day [seconds]
      

      
!=== End Variable Definition =============================================

!===  Set the Timestep count to zero
      ldas%tscount=0


!=== Read in LDAS Driver Boundary Conditions
      call readcard(ldas)                 !Read in user secified options

!=== Define Grid Definition Section (GDS) array
      if ( LDAS%DOMAIN > 1 ) call define_gds(ldas)

!=== Determine TILE MODULE dimension, allocate memory, and make tile space
      
c     if (ldas%rmos.eq.1 .or. ldas%rclm1.eq.1) then
       if (ldas%domain==1) then
          ldas%nch = ldas%nc*ldas%nr*ldas%maxt
       else !Global, consider only land areas
          ldas%nch = 0.33*(ldas%nc*ldas%nr*ldas%maxt)
       endif
       allocate (tile(ldas%nch),grid(ldas%nc,ldas%nr))
       call maketiles(ldas,grid,tile)      !Determine actual NCH
       deallocate (tile)                   !save memory
       allocate (tile(ldas%nch)) 
       call maketiles(ldas,grid,tile)
c     endif

!=== Allocate memory for catchment model
      ldas%ncatm=2677
      if (ldas%rcat.eq.1) allocate (cat(ldas%ncatm))
c     if (ldas%rcat.eq.1) call system('lamboot') 
      if (ldas%rcat.eq.1) print*,'Catchment model does not use LDAS',
     &   ' soil or vegetation; 87 ISLSCP used!'

!=== Allocate memory for select variables
!=== *NOTE*: Allocation for CLM2 tile variables are performed in clm_varder.F90 module      
      if(ldas%rmos.eq.1) allocate (mos(ldas%nch))   !MOSAIC
      if(ldas%rclm1.eq.1)then                        !CLM
       allocate (clm_tile(ldas%nch),clm_grid(ldas%nc,ldas%nr),
     &   clm1(ldas%nch))
      endif
      if(ldas%rnoah.eq.1) then                      !NOAH
        allocate (noah(ldas%nch))
      endif

      do c=1,ldas%nc                                !FORCING
       do r=1,ldas%nr
        if((ldas%feta.eq.1).or.(ldas%precsor(1).gt.0))then 
         allocate (grid(c,r)%etadata1(ldas%nf))
         allocate (grid(c,r)%etadata2(ldas%nf))
         allocate (grid(c,r)%precipdata1(ldas%nf))
         allocate (grid(c,r)%precipdata2(ldas%nf))
         allocate (grid(c,r)%precip(ldas%nf))
         allocate (grid(c,r)%etaprecip(24))
        endif
        if((ldas%fncep.eq.1).or.(ldas%fnasa.eq.1))then
         allocate (grid(c,r)%ncepdata1(ldas%nf))
         allocate (grid(c,r)%ncepdata2(ldas%nf))
        endif
        if(ldas%fcatch.eq.1)then
         allocate (grid(c,r)%catchdata1(ldas%nf))
         allocate (grid(c,r)%catchdata2(ldas%nf))
        endif
        if(ldas%domain >= 2 ) then
          if (ldas%startcode == 4) then
            allocate (grid(c,r)%glbdata1(ldas%nmif))
            allocate (grid(c,r)%glbdata2(ldas%nmif))
          else
            allocate (grid(c,r)%glbdata1(ldas%nf))
            allocate (grid(c,r)%glbdata2(ldas%nf))
          endif
        end if
        if (ldas%startcode .eq. 4) then
          allocate (grid(c,r)%forcing(ldas%nmif))
        else
          allocate (grid(c,r)%forcing(ldas%nf))
        endif
        
       enddo
      enddo 

       if(ldas%rmos.eq.1 .or. ldas%rclm1.eq.1 .or. ldas%rclm2 .eq.1 .or.
     &             ldas%rnoah.eq.1) then
        if (ldas%startcode .eq. 4) then
          do t=1,ldas%nch
            allocate (tile(t)%forcing(ldas%nmif))
          enddo
        else
          do t=1,ldas%nch
            allocate (tile(t)%forcing(ldas%nf))
          enddo
        endif
      endif

!===  Read in precip mask if necessary
      if(ldas%precipmask.eq.1)call precipmask(ldas,grid)

!===  Read in Reanalysis ECMWF land-sea mask if necessary
      if (ldas%fr_ecmwf.eq.1)call recmwfmask(ldas)

!=== SETUP MODEL PARAMETERS 
!=== Mosaic LSM
      if(ldas%rmos.eq.1)call setmosp(ldas,tile,mos)
!=== Map MOSAIC Albedo Coefficients      
      if(ldas%rmos.eq.1)call mapsib2umd(sib, ldas)

!=== CLM
      if(ldas%rclm1.eq.1)then

!--- Read in the clm input file (drv_clmin.dat)
       call drv_readclmin(ldas,tile,grid,drv,clm_tile,clm_grid,clm1)

!--- Initialize clm derived type components
       call clm1_typini(ldas%nch, clm1)

!--- Read in vegetation data and set tile information accordingly
       call drv_readvegtf (drv, clm_grid, tile, clm_tile, clm1, ldas)

!--- Transfer grid variables to tile space      
       call drv_g2clm (drv%udef, drv, clm_grid, clm_tile, clm1)   

!--- Read vegetation parameter data file for IGBP classification

       call drv_readvegpf (ldas, drv, clm_grid, clm_tile, clm1)

!--- Read in forcing to use in initialization
       if (ldas%startcode == 4) then

!--- Determine model forcing
        if (ldas%fgdas == 1) call getgdas(ldas,grid)
        if (ldas%fgeos == 1) call getgeos(ldas,grid) 

!--- Convert model interpolated forcing data from grid to CLM1 tile space
         call force2tile(ldas,tile,grid)
         call drv_getforce(ldas,drv,grid,clm_tile,clm1)

       endif ! End of CLM1 initialization section

!--- Initialize CLM and DIAG variables

!       do t=1,drv%nch 
!        clm%kpatch = t
!        call drv_clmini (drv, clm_grid, clm_tile(t), clm(t))      !Initialize CLM Variables
!       enddo
       call drv_clmini (ldas, drv, clm_grid, clm_tile, clm1)        !Initialize CLM Variables

!   Adjust observation heights
       clm1%forc_hgt_u=10.+clm1%displa+clm1%z0m
       clm1%forc_hgt_t=10.+clm1%displa+clm1%z0m
       clm1%forc_hgt_q=10.+clm1%displa+clm1%z0m

      endif   !End of Setting up CLM Parameters

      
!=== Catchment LSM
      if(ldas%rcat.eq.1)call setcatp(ldas,cat)

!=== NOAH LSM
      if(ldas%rnoah.eq.1) call setnoahp(ldas,tile,noah)

!=== Read Restart Files (1 to read, 2 to write restart)
      if(ldas%rmos.eq.1)call mosrst(1,ldas,tile,grid,mos)
      if(ldas%rclm1.eq.1)call drv_restart(1,drv,ldas,grid,clm_tile,clm1)
      if(ldas%rcat.eq.1)call catrst(1,ldas,cat)
      if(ldas%rnoah.eq.1)call noahrst(1,ldas,tile,grid,noah) 

!=== Determine if forcing grid dimensions change during run
      if(ldas%domain>1) call set_forcinggrid(ldas)

!=== Initialize Output Arrays and Analysis Terms
      if(ldas%rmos.eq.1)then
       call mos_initout(ldas,mos)
       mos%dtcanal=0.    
       do r=1,ldas%nr
        do c=1,ldas%nc
         do m=1,5 
          if(ldas%startcode.ne.6)grid(c,r)%mosbetak(m)=0.
         enddo
         if(ldas%startcode.ne.6)grid(c,r)%mosdelt=0.
        enddo
       enddo        
      endif

      if(ldas%rclm1.eq.1)then
        call clm1_initout(drv,clm1)
       clm1%dtcanal=0.
       do r=1,ldas%nr
        do c=1,ldas%nc
         do m=1,5 
          if(ldas%startcode.ne.6)grid(c,r)%clmbetak(m)=0.
         enddo
         if(ldas%startcode.ne.6)grid(c,r)%clmdelt=0.
        enddo
       enddo        
      endif

!=== Initialize Output Arrays
      if(ldas%rcat.eq.1) call cat_initout(ldas,cat)
      if(ldas%rnoah.eq.1) then 
        call noah_initout(ldas,noah)
      endif
      if(ldas%wfor.eq.1) call for_initout(ldas,grid)

!=== Check to make sure all models have the same start (TBD)

!=== Must initially set north of 60 N data to -1 to ensure model data is left
!=== when using global observed precipitation
      grid%obsprecip = -1.0
!=== Time Loop
      ldas%endtime=0

      do while (ldas%endtime.eq.0) !Continue until Endtime is reached
        call ticktime(ldas)
        ldas%tscount=ldas%tscount+1
        
       if(ldas%rclm1.eq.1) call drv_tick(drv)
!=== Get LDAS Base Forcing
       if(ldas%feta.eq.1)  call geteta(ldas,grid)  !Get EDAS/ETA 3-6 Hr Data
       if(ldas%feta.eq.1)  call getradbc(ldas,grid)
       if(ldas%fcatch.eq.1)  call getcatch(ldas,grid)  !Get Aarons Catchment Data
       if((ldas%fncep.eq.1).or.(ldas%fnasa.eq.1)) 
     &     call getncep(ldas,grid)                 !Get NCEP-LDAS 1 Hr Data
       if(ldas%fgdas.eq.1) call getgdas(ldas,grid) !Get NCEP-GDAS 3 Hr Data
       if(ldas%fgeos.eq.1) call getgeos(ldas,grid) !Get GEOS 3 Hr Data 
       if(ldas%fr_ecmwf.eq.1) call getreanlecmwf(ldas,grid) !Get Aarons 6hr 
                                                            !reanal-ECMWF 
       if(ldas%fecmwf.eq.1) call getecmwf(ldas,grid)!Get ECMWF 3-hr data

!=== Get LDAS Observed Radiation Forcing
       if ( LDAS%DOMAIN == 1 ) then !NLDAS
        if((ldas%pinker.gt.0).or.(ldas%nesdis.gt.0).or.
     &    (ldas%brttmp.gt.0)) then
             call getrad(ldas,grid)
        if(ldas%wfor.gt.0) call getpar(ldas,grid)
        end if
       else if (LDAS%DOMAIN >= 2) then !GLDAS 
         if (ldas%agrmetsw.gt.0) then
            call getgrad(ldas,grid)
         end if
       end if
!=== Get LDAS Observed Precipitation Forcing 

        if(ldas%precsor(1).gt.0) then
c        call precdata(ldas,grid)   !call to stageiv routine, not needed now?  
         call makeprecip(ldas,grid)
        endif

!=== Get GLDAS Observed Precipitation Forcing

        if (ldas%domain .ge. 2 .and. (ldas%gpcpsrc(1) .gt. 0
     &    .or. ldas%gpcpsrc(2).gt. 0 .or. ldas%gpcpsrc(3) .gt. 0)
     &    .or. ldas%gpcpsrc(4).gt. 0)
     &   call glbobspcp(ldas,grid)       

!=== Apply Elevation Corrections To Forcing Data
        if ((ldas%tempadj.eq.1).or.(ldas%presadj.eq.1).or.
     &       (ldas%humidadj.eq.1).or.(ldas%lwradadj.eq.1)) then
           if(ldas%fsource(6).ne.1) then
             if (ldas%fnasa.ne.1) then
              call elevadjust(ldas,grid)
             endif
           endif
        endif
!=== Transfer LDAS Forcing to Model Tiles
        if(ldas%rclm1.eq.1 .or. ldas%rmos.eq.1 .or. 
     &  ldas%rnoah.eq.1)
     &    call force2tile(ldas,tile,grid)
        if(ldas%rclm1.eq.1)
     &    call drv_getforce(ldas,drv,grid,clm_tile,clm1)
!=== Set MOSAIC Dynamic Variables
       if(ldas%rmos.eq.1) call mosdynp(ldas,mos,tile)
!=== Read in AVHRR/MODIS satellite LAI data for CLM

       if(ldas%rclm1 .eq. 1 .and. ldas%lai .ge. 2)
     &   call clm1lairead(ldas,drv,clm_tile,clm1)

!=== Read and Interpolate NOAH Greenness Fraction (monthly) 
       if(ldas%rnoah.eq.1) call noah_gfrac(ldas,noah,tile) 

!=== Read and Interpolate NOAH Albedo (quarterly) 
       if(ldas%rnoah.eq.1) call noah_alb(ldas,noah,tile) 

!=== Read & interpolate Maps of time/space varying observed parameters (TBD)
!=== Tile Loop & Run Models

!=== Call Mosaic Main Subroutine

       if(ldas%rmos.eq.1)then
        if(ldas%rpsas.eq.1)then
        if(mod(ldas%gmt,6.).eq.0)then
!          call mos_psas(ldas,grid,tile,mos)
         else
	  grid%tovsts=ldas%udef
          mos%dtcanal=0.   !Set analysis increment to zero if not assimilating
         endif
        endif
        if(ldas%rpsas.eq.1.and.ldas%rbias.eq.1)then
          call mos_bias(ldas,grid,tile,mos)
        else
           grid%fbias=0.
        endif
        do t=1,ldas%nch
         call mos_main(t,ldas,tile(t),mos(t),sib)
        enddo	

       endif
!=== Call CLM Main Subroutine

       if(ldas%rclm1.eq.1)then
        if(ldas%rpsas.eq.1)then
        if(mod(ldas%gmt,6.0).eq.0)then
!          call clm1_psas(ldas,grid,tile,clm1)
         else
	  grid%tovsts=ldas%udef
          clm1%dtcanal=0.   !Set analysis increment to zero if not assimilating
         endif
        endif
        if(ldas%rpsas.eq.1.and.ldas%rbias.eq.1)then
          call clm1_bias(ldas,grid,tile,clm1)
        else
          grid%fbias=0.
        endif
        do t=1,ldas%nch
          call clm1_main(clm1(t),drv%day)
        enddo 

       endif
           
     
!=== Call Catchment Main Subroutine

       if(ldas%rcat.eq.1) call cat_main(ldas,grid,cat)

!=== Call NOAH LSM Main Subroutine

       if(ldas%rnoah.eq.1)then
        do t=1,ldas%nch
          call noah_main(t,ldas,tile(t),noah(t))
        enddo
       endif


!=== Perform Land-surface Data Assimilation(TBD)
!=== Write Output
       call fsource(ldas)
       if(ldas%rmos.eq.1) call mos_out(ldas,tile,grid,mos)
       if(ldas%rclm1.eq.1) call clm1_out(ldas,tile,grid,clm_tile,clm1)
       if(ldas%rcat.eq.1) call cat_out(ldas,grid,cat)
       if(ldas%rnoah.eq.1) call noah_out(ldas,tile,grid,noah)
       if(ldas%wfor.eq.1) call for_out(ldas,tile,grid)
!=== Write Daily Restarts
       if((ldas%gmt.eq.(24.-ldas%writeintm).or.ldas%endtime.eq.1).and.
     &  ldas%rmos.eq.1)call mosrst(2,ldas,tile,grid,mos)
       if((ldas%gmt.eq.(24.-ldas%writeintc1).or.ldas%endtime.eq.1).and.
     &   ldas%rclm1.eq.1)call drv_restart(2,drv,ldas,grid,clm_tile,clm1)
       if((ldas%gmt.eq.(24.-ldas%writeintct).or.ldas%endtime.eq.1).and.
     &  ldas%rcat.eq.1) call catrst(2,ldas,cat)
       if((ldas%gmt.eq.(24.-ldas%writeintn).or.ldas%endtime.eq.1).and.
     &  ldas%rnoah.eq.1) call noahrst(2,ldas,tile,grid,noah) 
       if(ldas%rclm1.eq.1) then
!=== Return required surface fields to atmospheric model (return to grid space)
         call drv_clm2g (drv, clm_grid, clm_tile, clm1)
!=== Write spatially-averaged BC's and IC's to file for user
         if (clm1(1)%istep==1) call drv_pout(drv,clm_tile,clm1)
       endif
      enddo  !End Time Loop

      
c     if (ldas%rcat.eq.1) call system('lamclean')

      end

     
