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
! BOP
!
! !ROUTINE: setnoahp.F90
!
! !DESCRIPTION:
!  This subroutine retrieves MOSAIC parameters - Significant F90 revisions
!  below this subroutine will be requiret in the future.  
!
! !REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Code
!  31 Jul 2001: Matt Rodell; Updated for new soil parameter definition.
!  14 Feb 2002: Jon Gottschalck; Added allocated space for AVHRR LAI/DSAI     
!  07 Mar 2002: Brian Cosgrove; Corrected declaration of TEX1 var from real to int
!  14 Jan 2003: Urszula Jambor; Added conditional to check if need exists 
!               to allocate for AVHRR LAI/DSAI variables.
!
! !INTERFACE:
subroutine setmosp()
! !USES:      

  use lisdrv_module, only : lis, tile, gindex
  use mos_varder      ! NOAH tile variables
  use lis_indices_module
#if ( defined OPENDAP )
  use opendap_module
#endif
! EOP
   IMPLICIT NONE
 
   integer  :: line1,line2
   integer :: glnr,glnc,line,ios1
   real    :: temp1, temp2
   integer :: nc_dom
!=== Local Parameters for new soil definition ==============================
!==== Layer thicknesss correspond to soil property map files - do not change!
!       Thicknesses of soil layers for any simulations
!       These are compatible with the FIRST set of Matt Rodell's 
!       soil files, i.e. sand_nlis3.5.bfsa

      REAL, PARAMETER :: D1 = 0.02	! thickness of soil layer 1, m
      REAL, PARAMETER :: D2 = 1.48	! thickness of soil layer 2, m
      REAL, PARAMETER :: D3 = 2.00	! thickness of soil layer 3, m

!	thicknesss of soil layers (semi-official) for NLDAS simulations
!	These are compatible with either Yun Duans's soil maps
!       or the SECOND set of Matt Rodell's soil files, i.e., sand_nlis2m.bfsa
      REAL, PARAMETER :: DNLDAS1 = 0.1    ! thickness of soil layer 1, m
      REAL, PARAMETER :: DNLDAS2 = 0.3    ! thickness of soil layer 2, m
      REAL, PARAMETER :: DNLDAS3 = 1.6    ! thickness of soil layer 3, m


!==== Minimum values of sin(theta) based on vegetation type.
      REAL, PARAMETER :: S8 = 0.57787	! Closed shrubland
      REAL, PARAMETER :: S9 = 0.95504	! Open shrubland
      REAL, PARAMETER :: S12 = 0.1736	! Bare soil
      REAL, PARAMETER :: S0 = 0.05	! All others 
!==== Maximum allowable porosity.
      REAL, PARAMETER :: PORMAX=0.70

!=== Local Variables =====================================================
      INTEGER :: C,R,T,I,J,K,JJ,II,count     ! Loop counters
      REAL :: VEGP(LIS%D%NCH,MOSDRV%MOS_NVEGP)   
!      REAL :: VEGMP(LIS%D%NCH,MOSDRV%MOS_NMVEGP,12)   
!      REAL :: SOILP(LIS%D%NCH,MOSDRV%MOS_NSOILP)   
      REAL :: VALUE(LIS%P%NT,MOSDRV%MOS_NVEGP)
!      REAL :: VALUEMON(LIS%P%NT,MOSDRV%MOS_NMVEGP,12)
      REAL :: BASICSET(LIS%P%NT,MOSDRV%MOS_NSOILP)
      integer, allocatable :: soiltextij(:,:,:)
      integer, allocatable :: soiltext(:,:)
      integer, allocatable :: tex1(:)
      real, allocatable :: sand1(:,:)
      real, allocatable :: silt1(:,:)
      real, allocatable :: clay1(:,:)
      real, allocatable :: por1(:,:)
      real, allocatable :: por2(:,:)
      real, allocatable :: por3(:,:)
      real :: por1a,por2a,por3a
      real, allocatable :: slope(:,:)
      real, allocatable :: ksat1(:)
      real, allocatable :: psi1(:)
      real, allocatable :: b1(:)
      real :: hycon(12)
      real :: soilpot(12)
      real :: poros(12)
      real :: bparam(12)     
!=== End Variable Definition =============================================

!=== Get Vegetation Parameters for Mosaic Model in Tile Space

!=== Read in the Mosaic Static and Monthly Vegetation Parameter Files
      open(unit=15,file=mosdrv%mos_vfile,status='old')
      do j=1,mosdrv%mos_nvegp
       read(15,*)(value(i,j),i=1,lis%p%nt)
      enddo 
      close(15)
   
!=== Assign STATIC vegetation parameters to each tile based on the
!=== type of vegetation present in that tile.
!=== These parameters will be stored in one long array--structured
!=== as follows: Tile 1, all the parameters (1 through numparam)
!=== then Tile 2, all the parameters. 
!=== Then Tile 3, all the parameters etc.
      do i=1,lis%d%nch
       do j=1,mosdrv%mos_nvegp
        mos(i)%vegp(j)=value(tile(i)%vegt,j)				
       enddo !j
      enddo !i

!=== Original vegetation-based soil parameterization.
      if (LIS%D%soil .eq. 1) then

!===   Get Soil Parameters (Based on Vegetation) for Mosaic Model in Tile Space
!===   Read in Soil Parameter Data
	    
        open(10,file=mosdrv%mos_sfile,status='old', &
        access='sequential')
     
        do i=1,mosdrv%mos_nsoilp
         read(10,*)(basicset(jj,i),jj=1,lis%p%nt)
        enddo
        close(10)
    
        do i=1,lis%d%nch
         k=tile(i)%vegt
          do j=1,mosdrv%mos_nsoilp
           mos(i)%soilp(j)=basicset(k,j)                  
          enddo !j 
        enddo !i	 

     end if	!soil=1

!    New soil parameterization.
     line1 = nint((lis%d%gridDesc(4)-lis%d%gridDesc(44))/lis%d%gridDesc(9))+1
     line2 = nint((lis%d%gridDesc(5)-lis%d%gridDesc(45))/lis%d%gridDesc(10))+1
     allocate(soiltextij(11,lis%d%lnc,lis%d%lnr))
     allocate(soiltext(11,lis%d%nch))
     allocate(tex1(lis%d%nch))
     allocate(ksat1(lis%d%nch))
     allocate(psi1(lis%d%nch))
     allocate(b1(lis%d%nch))
            
      if ((lis%d%soil.eq.2).or.(lis%d%soil.eq.3).or.(lis%d%soil.eq.6)) then
        OPEN(14,FILE=LIS%p%soilclass_file,STATUS='OLD')
        do k=1,11
        do j=1,lis%d%lnr
        do i=1,lis%d%lnc
        read (14,'(I3,1X,I3,1X,3X,I2)') &
       ii,jj,soiltextij(k,i,j)
        if (soiltextij(k,i,j).eq.13) soiltextij(k,i,j)=12
        enddo
        enddo
        enddo

	do k=1,11
	do i=1,lis%d%nch
	soiltext(k,i)=soiltextij(k,tile(i)%col,tile(i)%row)
	enddo
	enddo

        print *,'   '
        print *,'SETTING SOIL CLASS 13 to 12...since have no params for organic class'
        print *,'WARNING!!!!!!!!!!!!!!!!!!'
        print *,'   '
        CLOSE(14)


!===    Determine saturated hydraulic conductivity at the surface.
        DATA HYCON/3.75E-05,1.81E-05,6.36E-06,2.14E-06,2.14E-06, &
                  1.44E-06, &
                  9.72E-07, &
                  1.19E-06,3.33E-07,2.78E-07,4.44E-07,5.28E-07/
        do i=1,lis%d%nch
        KSAT1(i)=HYCON(soiltext(1,i))
        enddo


!===    Determine saturated soil potential based on texture class.
        DATA SOILPOT/-0.0692,-0.0363,-0.1413,-0.7586,-0.7586,-0.3548, &
                    -0.1349,-0.6166,-0.263,-0.0977,-0.3236,-0.4677/
        do i=1,lis%d%nch
        PSI1(i)=SOILPOT(soiltext(1,i))
        enddo

!===    Determine parameter b based on texture class.
        DATA BPARAM/2.79,4.26,4.74,5.33,5.33,5.25,6.77, &
                    8.72,8.17,10.73,10.39,11.55/

        do i=1,lis%d%nch
	B1(i)=BPARAM(soiltext(1,i))
        enddo


!====   assign soil parameters in tile space.
        write(*,*) 'Assigning soil parameters in tile space.'
        do i=1,lis%d%nch
          mos(i)%soilp(1) = b1(i)
          	  if (mos(i)%soilp(1) .le. 0.0) then 
	    write(*,*) 'col,row,tile,b1',i,& 
     	     mos(i)%soilp(1)
	    stop
         end if
         mos(i)%soilp(2) = psi1(i)
         if (mos(i)%soilp(2) .ge. 0.0) then 
	    write(*,*) 'col,row,tile,psi1', &
     	     i,mos(i)%soilp(2)
	    stop
         end if
         mos(i)%soilp(3) = ksat1(i)
         if (mos(i)%soilp(3) .le. 0.0) then 
	    write(*,*) 'col,row,tile,ksat1', & 
     	     i,mos(i)%soilp(3)
	    stop
         end if
!	set layer thicknesses based on soil scheme chosen in card file
         if (lis%d%soil.eq.2) then
            mos(i)%soilp(4) = d1
            mos(i)%soilp(5) = d2
            mos(i)%soilp(6) = d3
        endif
        if (lis%d%soil.eq.3) then
           mos(i)%soilp(4) = dnldas1
           mos(i)%soilp(5) = dnldas2
           mos(i)%soilp(6) = dnldas3
        endif
        if (lis%d%soil.eq.6) then
           mos(i)%soilp(4) = dnldas1
           mos(i)%soilp(5) = dnldas2
           mos(i)%soilp(6) = dnldas3
        endif
     enddo

     line1 = nint((lis%d%gridDesc(4)-lis%d%gridDesc(44))/lis%d%gridDesc(49))+1
     line2 = nint((lis%d%gridDesc(5)-lis%d%gridDesc(45))/lis%d%gridDesc(50))+1

     open(21,file=lis%p%slfile,status='old',form='unformatted')
     allocate(slope(lis%d%lnc,lis%d%lnr))
           read(21)slope
     close(21)

!===     Special for slope: use vegetation-based values as minima for
!===      desert-like vegetation classes.
!     if(lis%d%gridDesc(9) .ne. 0.01) then 
!        do r=1,lis%d%lnr
!           do c=1,lis%d%lnc
!              if(gindex(c,r).ne.-1) then 
!                 i = gindex(c,r)
!                 select case (tile(i)%vegt)
!                 case (8)	! closed shrubland
!                    mos(i)%soilp(7) = max(s8, sin(slope(c,r)))
!                 case (9)	! open shrubland
!                    mos(i)%soilp(7) = max(s9, sin(slope(c,r)))
!                 case (12)	! bare ground
!                    mos(i)%soilp(7) = max(s12, sin(slope(c,r)))
!                 case default	! all other vegetation classes
!                    mos(i)%soilp(7) = max(s0, sin(slope(c,r)))
!                 end select
!                 if (MOS(I)%SOILP(7) .le. 0.0) then 
!                    write(*,*) 'COL,ROW,TILE,SIN(THETA)',MOS(I)%SOILP(7)
!                 end if
!              endif
!           enddo
!        enddo
!     elseif(lis%d%gridDesc(9) .eq. 0.01) then 

        do i=1,lis%d%nch
           select case (tile(i)%vegt)
           case (8)	! closed shrubland
              mos(i)%soilp(7) = max(s8, sin(slope(tile(i)%col, &
                   tile(i)%row)))
           case (9)	! open shrubland
              mos(i)%soilp(7) = max(s9, sin(slope(tile(i)%col, &
                   tile(i)%row)))
           case (12)	! bare ground
              mos(i)%soilp(7) = max(s12, sin(slope(tile(i)%col, &
                   tile(i)%row)))
           case default	! all other vegetation classes
              mos(i)%soilp(7) = max(s0, sin(slope(tile(i)%col, &
                   tile(i)%row)))
           end select
           if (MOS(I)%SOILP(7) .le. 0.0) then 
              write(*,*) 'COL,ROW,TILE,SIN(THETA)',TILE(I)%COL,TILE(I)%ROW, &
                   I,MOS(I)%SOILP(7)
           end if
        enddo
!     endif
       deallocate(slope)
!===    Determine Porosity at 3 Levels. HARDWIRED FOR 0-10,10-40,40-200cm
!===    Porosity data for each of 12 soil textures
        DATA POROS/0.339,0.421,0.434,0.476,0.476,0.439,0.404,0.464, &
                  0.465,0.406,0.468,0.457/
          allocate(por1(lis%d%lnc,lis%d%lnr))
          allocate(por2(lis%d%lnc,lis%d%lnr))
          allocate(por3(lis%d%lnc,lis%d%lnr))

        do j=1,lis%d%lnr
        do i=1,lis%d%lnc
        POR1(i,j)=(0.5*poros(soiltextij(1,i,j)))+(0.5* &
                  poros(soiltextij(2,i,j)))
        POR2(i,j)=((1.0/3.0)*poros(soiltextij(3,i,j)))+ &
                 ((1.0/3.0)*poros(soiltextij(4,i,j)))+ &
                 ((1.0/3.0)*poros(soiltextij(5,i,j)))
        POR3(i,j)=(0.125*poros(soiltextij(6,i,j)))+(0.125* &
                  poros(soiltextij(7,i,j)))+ &
                 (0.125*poros(soiltextij(8,i,j)))+(0.3125* &
                 poros(soiltextij(9,i,j)))+ &
                 (0.3125*poros(soiltextij(10,i,j)))
        enddo
        enddo


          line1 = nint((lis%d%gridDesc(4)-lis%d%soil_gridDesc(1))/lis%d%gridDesc(9))+1
          line2 = nint((lis%d%gridDesc(5)-lis%d%soil_gridDesc(2))/lis%d%gridDesc(10))+1
          nc_dom = nint((lis%d%soil_gridDesc(4)-lis%d%soil_gridDesc(2))/&
               lis%d%soil_gridDesc(5))+1
          do r=1,lis%d%lnr
             do c=1,lis%d%lnc
                glnc = line2+c-1
                glnr = line1+r-1
                line = (glnr-1)*nc_dom+glnc
             enddo
          enddo
          close(15)
!       endif
       do i=1,lis%d%nch
          por1a=amin1(por1(tile(i)%col,tile(i)%row),pormax)
          por2a=amin1(por2(tile(i)%col,tile(i)%row),pormax)
          por3a=amin1(por3(tile(i)%col,tile(i)%row),pormax)
          if (lis%d%soil.eq.2) then
             mos(i)%soilp(8) = por1a * d1 * 1000.0
          endif
          if (lis%d%soil.eq.3) then
             mos(i)%soilp(8) = por1a * dnldas1 * 1000.0
          endif
          if (lis%d%soil.eq.6) then
             mos(i)%soilp(8) = por1a * dnldas1 * 1000.0
          endif
	  if (mos(i)%soilp(8) .le. 0.0) then 
             write(*,*) 'col,row,tile,wsat1',& 
                  i,mos(i)%soilp(8)
             stop
	  end if
          if (lis%d%soil.eq.2) then
             mos(i)%soilp(9) = por2a * d2 * 1000.0
          endif
          if (lis%d%soil.eq.3) then
             mos(i)%soilp(9) = por2a * dnldas2 * 1000.0
          endif
          if (lis%d%soil.eq.6) then
             mos(i)%soilp(9) = por2a * dnldas2 * 1000.0
          endif

          
	  if (mos(i)%soilp(9) .le. 0.0) then 
             write(*,*) 'col,row,tile,wsat2',& 
                  i,mos(i)%soilp(9)
             stop
	  end if

          if (lis%d%soil.eq.2) then
             mos(i)%soilp(10) = por3a * d3 * 1000.0
          endif

          if (lis%d%soil.eq.3) then
             mos(i)%soilp(10) = por3a * dnldas3 * 1000.0
          endif
          if (lis%d%soil.eq.6) then
             mos(i)%soilp(10) = por3a * dnldas3 * 1000.0
          endif

	  if (mos(i)%soilp(10) .le. 0.0) then 
             write(*,*) 'col,row,tile,wsat3',& 
                  i,mos(i)%soilp(10)
             stop
	  end if
       enddo !i	 
!	do i=1,lis%d%nch
!	if ((tile(i)%row.eq.117).and.(tile(i)%col.eq.77)) print *,'i=',i
!	enddo
!	stop
		
       deallocate(por1)
       deallocate(por2)
       deallocate(por3)
    end if	!soil=2 or 3 or 6
    return
    end subroutine setmosp

