!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. CONTINENTAL-SCALE
!  D                                      L  LAND MODELING AND DATA
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  ASSIMILATION PROJECT.
!  S                                      A  THIS IS THE GSFC-LDAS CODE.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  HTTP://LDAS.GSFC.NASA.GOV
!
!   GSFC - NCEP - OH - PRINCETON - WASHINGTON - RUTGERS
!
!=========================================================================
! MAKEPRECIP.F:
!
! DESCRIPTION:
!  Corrects Temperature, Pressure, Humidity and Longwave Radiation
!  values for differences in elevation between EDAS and LDAS grids.
!
! REVISION HISTORY:
!  11 Apr 2000: Brian Cosgrove; Initial Code
!  12 May 2000: Brian Cosgrove; Corrected for zero humidities
!  09 Aug 2000: Brian Cosgrove; Corrected program so that
!               it only performs calculations if both
!               the elevation difference file and the forcing
!               data file (use temperature data as check for all
!               fields) contain defined values
!  25 Jan 2001: Matt Rodell; Compute number of input and output
!		grid points, use to allocate local arrays
!  27 Feb 2001: Brian Cosgrove; Added statement to check for use of
!               catchment data so that correct elevation correction
!               files is used
!  15 Mar 2001: Jon Gottschalck; if-then to handle negative vapor
!		pressures in long wave correction
!  15 Mar 2001: Matt Rodell; merge NLDAS and GLDAS versions
!
!
!=========================================================================


	SUBROUTINE ELEVADJUST(LDAS,GRID)

      USE ldas_module      ! LDAS non-model-specific 1-D variables
      USE grid_module      ! LDAS non-model-specific grid variables
      IMPLICIT NONE
      type (ldasdec) LDAS
      type (griddec) GRID(LDAS%NC,LDAS%NR)

!        INTEGER NLDAS,NEDAS,P,X,Y,BB,COUNT,I
        INTEGER NLDAS,NEDAS,P,X,Y,COUNT,I
	integer, parameter :: BB=2016
        integer err !iostat error code
C
	real, allocatable, dimension (:) :: ELEVDIFF, TEMP, PRES,NARRELEV,
     &	 HUMID, LWRAD, TEDAS, TBAR, PEDAS, EEDAS, ELDAS,RHO,
     &   temphybrid,humidhybrid,preshybrid
	real ::  LAPSE, GRAV, RDRY, RATIO, ESATEDAS, QSATEDAS, RH, 
     &	 ESATLDAS, QSATLDAS, EMISSLDAS, EMISSEDAS
	real rv(LDAS%NC*LDAS%NR),rvedas(LDAS%NC*LDAS%NR)
	real tv(LDAS%NC*LDAS%NR),tvedas(LDAS%NC*LDAS%NR)
        real tvedashybrid(LDAS%NC*LDAS%NR),vtemp(LDAS%NC*LDAS%NR)
        real vtemphybrid(LDAS%NC*LDAS%NR)
C
C	write(*,*) 'Begin elevadjust'

        GRAV = 9.81
        RDRY = 287.
        LAPSE = -0.0065

!	Compute number of points in input and output grids.
	NEDAS = LDAS%NCold * LDAS%NRold
	NLDAS = LDAS%NC * LDAS%NR

!	Allocate the local arrays.
        allocate (NARRELEV(NLDAS))
	allocate (temphybrid(NLDAS))
        allocate (preshybrid(NLDAS))
        allocate (humidhybrid(NLDAS))
	allocate (ELEVDIFF(NLDAS))
	allocate (TEMP(NLDAS))
	allocate (PRES(NLDAS))
	allocate (HUMID(NLDAS))
	allocate (LWRAD(NLDAS))
	allocate (TEDAS(NLDAS))
	allocate (TBAR(NLDAS))
	allocate (PEDAS(NLDAS))
	allocate (EEDAS(NLDAS))
	allocate (ELDAS(NLDAS))
        allocate (RHO(NLDAS))

	IF (LDAS%FSOURCE(16).NE.1) THEN
!         OPEN (UNIT=21,FILE=LDAS%ELEVFILE,FORM='UNFORMATTED',
!     &	       STATUS='OLD',IOSTAT=err)
         OPEN (UNIT=21,FILE=LDAS%ELEVFILE,
     &         STATUS='OLD',IOSTAT=err)
         OPEN (UNIT=22,FILE=LDAS%NARRELEV,
     &         STATUS='OLD',IOSTAT=err)


	 IF (err /= 0) THEN
	    write(79,*) "STOP: problem opening elevation difference file" 
	    write(79,*) "Try running without elevation correction option."
	    print*, "STOP: problem opening elevation difference file" 
	    print*, "Try running without elevation correction option."
	    STOP
	 ELSE
	 do i=1,nldas
	    READ (21,*) ELEVDIFF(I)
	 enddo
	 ENDIF
         CLOSE(21)

         IF (err /= 0) THEN
            write(79,*) "STOP problem opening elevation difference file"
            write(79,*) "Try running without elevation correction 
     &  option"
            print*, "STOP: problem opening elevation difference file"
            print*, "Try running without elevation correction option."
            STOP
         ELSE
         do i=1,nldas
            READ (22,*) NARRELEV(I)
         enddo
         ENDIF
         CLOSE(22)

        ELSEIF (LDAS%FSOURCE(16).EQ.1) THEN
         OPEN (UNIT=21, FILE=LDAS%CATELEVFILE,FORM = 'UNFORMATTED')
         READ (21) (ELEVDIFF(I), I = 1, NLDAS)
         CLOSE(21)
	ENDIF

          COUNT = 0
          DO Y = 1, LDAS%NR
            DO X = 1, LDAS%NC
              GRID(X,Y)%FORCING(9)=GRID(X,Y)%FORCING(9)-
     &  NARRELEV(X+COUNT)
	      TEMP(X+COUNT) = GRID(X,Y)%FORCING(14)
	      PRES(X+COUNT) = GRID(X,Y)%FORCING(3)
	      HUMID(X+COUNT) = GRID(X,Y)%FORCING(15)
	      LWRAD(X+COUNT) = GRID(X,Y)%FORCING(21)

C       Correct for humidities of 0.0
              IF (HUMID(X+COUNT).EQ.0) HUMID(X+COUNT)=1E-4
            END DO
            COUNT = COUNT+LDAS%NC
          END DO

! compute virtual temperature
        do p=1,nldas
        IF (TEMP(P).GE.-9998.) THEN
        EEDAS(P)=(HUMID(P)*PRES(P))/0.622
        RvEDAS(p)=( (0.622*EEDAS(p))/(PRES(p)-EEDAS(p)) )
        TvEDAS(p) = (1 + 0.61* RvEDAS(p))*TEMP(p)
c        TvEDAS(p) = (1 + 0.61* RvEDAS(p))*TEMP(p)
c	print *,'rvedas,eedas,pres,temp,tvedas',rvedas(p),
c     &  eedas(p),pres(p),temp(p),tvedas(p)
        endif
        enddo
c	open (unit=16,file='out.bin',form='unformatted')
c	write(16) temp
c	write(16) tvedas
	
C	Correct Temperature
C	print *,'correcting temperature'
          DO P = 1, NLDAS
c            IF ( (ELEVDIFF(P).GE.-9998.) THEN
            IF ((ELEVDIFF(P).GE.-9998.).AND.
     &            (TEMP(P).GE.-9998.) ) THEN


                TEDAS(P)=TEMP(P)
! compute elevation corrected temperature
                TEMP(P)=TEMP(P)+(LAPSE*ELEVDIFF(P))
! compute elevation corrected virtual temperature
                VTEMP(P)=TvEDAS(P)+(LAPSE*ELEVDIFF(P))
                TBAR(P)=(TvEDAS(P)+VTEMP(P))/2.
            END IF
          END DO

c	write(16)tbarold	
c	write(16)tbar
c	close(16)

C	write(*,*) 'Temperature correction complete'

C	Correct Pressure
C	print*,'correcting pressure'
          DO P = 1, NLDAS
            IF ((ELEVDIFF(P).GE.-9998.).AND.
     &            (TEMP(P).GE.-9998.) ) THEN

c            IF (ELEVDIFF(P).GE.-9998.) THEN

                PEDAS(P)=PRES(P)
                PRES(P)=PRES(P)/(EXP((GRAV*ELEVDIFF(P))/
     +                                   (RDRY*TBAR(P))))
c		PRES(P)=PRES(P)/(EXP((GRAV*ELEVDIFF(P))/
c     +                                   (RDRY*(TBAR(P)+1))))

c            IF(P.EQ.41362) PRINT*,PRES(P),ELEVDIFF(P),
c     +        TBAR(P),PEDAS(P),TEMP(P),TEDAS(P),TvEDAS(P),
c     +        VTEMP(P)

            END IF
          END DO

C	write(*,*) 'Pressure correction complete'


C	Correct Humidity
C	print *,'correcting humidity'
          DO P = 1, NLDAS
c            IF (ELEVDIFF(P).GE.-9998.) THEN
            IF ((ELEVDIFF(P).GE.-9998.).AND.
     &            (TEMP(P).GE.-9998.) ) THEN

                EEDAS(P)=(HUMID(P)*PEDAS(P))/0.622
                ESATEDAS=611.2*(EXP((17.67*(TEDAS(P)-273.15))/
     +                              ((TEDAS(P)-273.15)+243.5)))
                QSATEDAS=(0.622*ESATEDAS)/(PEDAS(P)-(0.378*
     +                                                 ESATEDAS))
                RH=(HUMID(P)/QSATEDAS)*100.
                ESATLDAS=611.2*(EXP((17.67*(TEMP(P)-273.15))/
     +                              ((TEMP(P)-273.15)+243.5)))
                QSATLDAS=(0.622*ESATLDAS)/(PRES(P)-(0.378*
     +                                                 ESATLDAS))
                HUMID(P)=(RH*QSATLDAS)/100.
                ELDAS(P)=(HUMID(P)*PRES(P))/0.622
            END IF
          END DO

c	write(*,*) 'Humidity correction complete'

C	Correct Long Wave
c	print *,'correcting long wave'
          DO P = 1, NLDAS
            IF ((ELEVDIFF(P).GE.-9998.).AND.
     &            (TEMP(P).GE.-9998.) ) THEN
c            IF (ELEVDIFF(P).GE.-9998.) THEN
                EEDAS(P)=EEDAS(P)/100.
                ELDAS(P)=ELDAS(P)/100.
c      Correct for negative vapor pressure at very low temperatures at
c	high latitudes                
                IF (EEDAS(P) .LE. 0) THEN
                  EEDAS(P) = 1E-20
                ENDIF
                IF (ELDAS(P) .LE. 0) THEN
                  ELDAS(P) = 1E-20
                ENDIF
                EMISSEDAS=1.08*(1-EXP(-EEDAS(P)**(TEDAS(P)/BB)))
                EMISSLDAS=1.08*(1-EXP(-ELDAS(P)**(TEMP(P)/BB)))
                RATIO=(EMISSLDAS*(TEMP(P)**4))/
     +                (EMISSEDAS*(TEDAS(P)**4))
                LWRAD(P)=LWRAD(P)*RATIO
            END IF
          END DO

c	write(*,*) 'Longwave correction complete'

          COUNT = 0
          DO Y = 1, LDAS%NR
            DO X = 1, LDAS%NC
              IF ((ELEVDIFF(X+COUNT).GE.-9998.).AND.
     &            (TEMP(X+COUNT).GE.-9998.) ) THEN
               GRID(X,Y)%FORCING(14)=TEMP(X+COUNT)
               GRID(X,Y)%FORCING(3)=PRES(X+COUNT)
               GRID(X,Y)%FORCING(15)=HUMID(X+COUNT)
               GRID(X,Y)%FORCING(21)=LWRAD(X+COUNT)
              ENDIF
            END DO
            COUNT = COUNT+LDAS%NC
          END DO


c now correct fields on hybrid levels
          COUNT = 0
          DO Y = 1, LDAS%NR
            DO X = 1, LDAS%NC
              TEMPhybrid(X+COUNT) = GRID(X,Y)%FORCING(7)
              PREShybrid(X+COUNT) = GRID(X,Y)%FORCING(8)
              HUMIDhybrid(X+COUNT) = GRID(X,Y)%FORCING(10)

C       Correct for humidities of 0.0
              IF (HUMIDhybrid(X+COUNT).EQ.0) HUMIDhybrid(X+COUNT)=1E-4
            END DO
            COUNT = COUNT+LDAS%NC
          END DO


! compute virtual temperature
        do p=1,nldas
        IF (TEMPhybrid(P).GE.-9998.) then
        EEDAS(P)=(HUMIDhybrid(P)*PREShybrid(P))/0.622
        RvEDAS(p)=( (0.622*EEDAS(p))/(PREShybrid(p)-EEDAS(p)) )
        TvEDAShybrid(p) = (1 + 0.61* RvEDAS(p))*TEMPhybrid(p)
        endif
        enddo

c        write(16) tvedashybrid
c        close(16) 
c        stop

C       Correct Temperature
          DO P = 1, NLDAS
c            IF ( (ELEVDIFF(P).GE.-9998.) THEN
            IF ((ELEVDIFF(P).GE.-9998.).AND.
     &            (TEMPhybrid(P).GE.-9998.) ) THEN


                TEDAS(P)=TEMPhybrid(P)

! compute elevation corrected temperature
                TEMPhybrid(P)=TEMPhybrid(P)+(LAPSE*ELEVDIFF(P))
! compute elevation corrected virtual temperature
                VTEMPhybrid(P)=TvEDAShybrid(P)+(LAPSE*ELEVDIFF(P))
                TBAR(P)=(TvEDAShybrid(P)+VTEMPhybrid(P))/2.
            END IF
          END DO


C       Correct Pressure
          DO P = 1, NLDAS
            IF ((ELEVDIFF(P).GE.-9998.).AND.
     &            (TEMPhybrid(P).GE.-9998.) ) THEN

c            IF (ELEVDIFF(P).GE.-9998.) THEN

                PEDAS(P)=PREShybrid(P)
                PREShybrid(P)=PREShybrid(P)/(EXP((GRAV*ELEVDIFF(P))/
     +                                   (RDRY*TBAR(P))))
            END IF
          END DO


C       Correct Humidity
          DO P = 1, NLDAS
c            IF (ELEVDIFF(P).GE.-9998.) THEN
            IF ((ELEVDIFF(P).GE.-9998.).AND.
     &            (TEMPhybrid(P).GE.-9998.) ) THEN

                EEDAS(P)=(HUMIDhybrid(P)*PEDAS(P))/0.622
                ESATEDAS=611.2*(EXP((17.67*(TEDAS(P)-273.15))/
     +                              ((TEDAS(P)-273.15)+243.5)))
                QSATEDAS=(0.622*ESATEDAS)/(PEDAS(P)-(0.378*
     +                                                 ESATEDAS))
                RH=(HUMIDhybrid(P)/QSATEDAS)*100.
                ESATLDAS=611.2*(EXP((17.67*(TEMPhybrid(P)-273.15))/
     +                              ((TEMPhybrid(P)-273.15)+243.5)))
                QSATLDAS=(0.622*ESATLDAS)/(PREShybrid(P)-(0.378*
     +                                                 ESATLDAS))
                HUMIDhybrid(P)=(RH*QSATLDAS)/100.
                ELDAS(P)=(HUMIDhybrid(P)*PREShybrid(P))/0.622
            END IF
          END DO



          COUNT = 0
          DO Y = 1, LDAS%NR
            DO X = 1, LDAS%NC
              IF ((ELEVDIFF(X+COUNT).GE.-9998.).AND.
     &            (TEMPhybrid(X+COUNT).GE.-9998.) ) THEN
               GRID(X,Y)%FORCING(7)=TEMPhybrid(X+COUNT)
               GRID(X,Y)%FORCING(8)=PREShybrid(X+COUNT)
               GRID(X,Y)%FORCING(10)=HUMIDhybrid(X+COUNT)

c compute density and adjust surface exchange coefficient
!               RHO(X+COUNT)=(PRES(X+COUNT))
!     &  /(RDRY*VTEMP(X+COUNT))
!               GRID(X,Y)%FORCING(17)=GRID(X,Y)%FORCING(17)/RHO(X+COUNT)
              ENDIF
            END DO
            COUNT = COUNT+LDAS%NC
          END DO


        deallocate (RHO)
	deallocate (ELEVDIFF)
        deallocate (TEMPhybrid)
        deallocate (PREShybrid)
        deallocate (HUMIDhybrid)
	deallocate (NARRELEV)
	deallocate (TEMP)
	deallocate (PRES)
	deallocate (HUMID)
	deallocate (LWRAD)
	deallocate (TEDAS)
	deallocate (TBAR)
	deallocate (PEDAS)
	deallocate (EEDAS)
	deallocate (ELDAS)

c	write(*,*) 'End elevadjust'


        RETURN
        END

