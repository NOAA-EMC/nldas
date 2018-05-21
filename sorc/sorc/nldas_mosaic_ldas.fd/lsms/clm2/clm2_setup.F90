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
! !ROUTINE: clm2_setup.F90
! 
! !DESCRIPTION: 
! 
! Completes the CLM2 setup routines. 
! 
! !REVISION HISTORY: 
! 
! 20 Jan 2003  Sujay Kumar Initial Specification
! 
! !INTERFACE:
subroutine clm2_setup()
! !USES:
  use lisdrv_module, only: lis, tile
  use spmdMod
  use time_manager
  use clm_varder
  use clm_varcon, only: eccen, obliqr, lambm0 , mvelpp
  use clm_varctl, only: nsrest
!EOP
  implicit none
  logical  :: readini               !true if read in initial data set
  integer  :: yr                    !current year (0, ...)
  integer  :: mon                   !current month (1 -> 12)
  integer  :: day                   !current day (1 -> 31)
  integer  :: k
  integer  :: i, ncsec
!BOC
#if ( ! defined OPENDAP )
  if ( masterproc ) then
#endif

! ----------------------------------------------------------------------
! Get current date
! ----------------------------------------------------------------------
     
     if ( masterproc ) then
!        call get_curr_date(yr, mon, day, ncsec)
        yr = lis%t%yr
        mon = lis%t%mo
        day = lis%t%da
        ncsec = lis%t%ss
     endif
! ----------------------------------------------------------------------
! If initial run: initialize time-varying data 
! If continuation run: end of initialization because time varying
! read in from restart file
! ----------------------------------------------------------------------

     if (nsrest == 0) then
        call canhtset()
        if (masterproc) then
           write (6,*) ('Attempting to initialize time variant variables .....')
        endif

        readini = .false.
        call iniTimeVar (readini, eccen, obliqr, lambm0 , mvelpp, lis, tile)
        
        if (masterproc) then
           write (6,*) ('Successfully initialized time variant variables')
           write (6,*)
        endif
        
     endif

! ----------------------------------------------------------------------
! End initialization
! ----------------------------------------------------------------------
     
     if (masterproc) then
        write (6,*) ('Successfully initialized the land model')
        if (nsrest == 0) then
           write (6,*) 'begin initial run at: '
        else
           write (6,*) 'begin continuation run at:'
        end if
        write (6,*) '   nstep= ',get_nstep(lis%t), &
             ' year= ',yr,' month= ',mon,' day= ',day,' seconds= ',ncsec
        write (6,*)
        write (6,'(72a1)') ("*",i=1,60)
        write (6,*)
     endif
     clm%totfsa=0.              ! solar absorbed solar radiation [W/m2]
     clm%toteflx_lwrad_net=0.   ! net longwave radiation [W/m2]
     clm%toteflx_lh_tot=0.      ! total latent heat flux [W/m2]
     clm%toteflx_sh_tot=0.      ! total sensible heat flux [W/m2]      
     clm%toteflx_soil_grnd=0.   ! ground heat flux [W/m2]
     clm%totqflx_snomelt=0.     ! snowmelt heat flux [W/m2]
     clm%totrain=0.             ! accumulation of rain [mm]
     clm%totsnow=0.             ! accumulation of snow [mm]
     clm%totqflx_evap=0.        ! total evaporation [mm]
     clm%totqflx_surf=0.        ! surface runoff [mm]
     clm%totqflx_drain=0.       ! subsurface runoff [mm]
     clm%totqflx_ecanop=0.      ! interception evaporation [W/m2]
     clm%totqflx_tran_veg=0.    
     clm%totqflx_evap_grnd=0.
     clm%totqflx_sub_snow=0.
     clm%count=0
     clm%canopint = 0
     clm%acond=0.
#if ( ! defined OPENDAP )
     endif    
     if ( npes > 1 ) then 
        call clm2_scatter
     endif
#endif
     return
!EOC
   end subroutine clm2_setup

