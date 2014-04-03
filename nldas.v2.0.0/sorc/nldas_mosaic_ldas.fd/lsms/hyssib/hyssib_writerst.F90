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
! !ROUTINE: hyssib_writerst.F90
!
! !DESCRIPTION:
!  This program writes restart files for HY-SSiB.  This includes all relevant
!  water/energy storages, tile information, and time information.  It also
!  rectifies changes in the tile space.  
!
! !REVISION HISTORY:
! 1  Oct 1999: Jared Entin, Initial code
! 15 Oct 1999: Paul Houser, Significant F90 Revision
! 05 Sep 2001: Brian Cosgrove, Modified code to use Dag Lohmann's NOAA
!              initial conditions if necessary.  This is controlled with
!              local variable NOAAIC.  Normally set to 0 in this subroutine
!              but set to 1 if want to use Dag's NOAA IC's.  Changed output
!              directory structure, and commented out if-then check so that
!              directory is always made.
! 28 Apr 2002: Kristi Arsenault, Added SSIB LSM into LDAS
! 28 May 2002: Kristi Arsenault, For STARTCODE=4, corrected SNEQV values  
!               and put SMC, SH2O, STC limit for GDAS and GEOS forcing.
! 14 Jun 2003: Sujay Kumar, Separated the write restart from the original code
! 30 Oct 2003: Matt Rodell, Added back COL,ROW,FGRD,VEGT write statements
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
!
! RESTART FILE FORMAT(fortran sequential binary):
!  YR,MO,DA,HR,MN,SS,VCLASS,NCH !Restart time,Veg class,no.tiles, no.soil lay 
!  TILE(NCH)%COL        !Grid Col of Tile
!  TILE(NCH)%ROW        !Grid Row of Tile
!  TILE(NCH)%FGRD       !Fraction of Grid covered by tile
!  TILE(NCH)%VEGT       !Vegetation Type of Tile
!  HYSSIB(NCH)%STATES     !Model States in Tile Space
!
! !INTERFACE:
      SUBROUTINE HYSSIB_WRITERST()
! !USES:
      use lisdrv_module, only : lis,tile
      use hyssib_varder         ! HY-SSiB tile variables
      use time_manager
      use tile_spmdMod
!EOP
      implicit none      

!=== Local Variables =====================================================
      INTEGER :: C,R,T,I,J,L,N,F ! Loop counters
      CHARACTER*80 FILEN,MKFYRMO
      CHARACTER*1  FNAME(80),FBASE(80),FSUBS(80),FMKDIR(80)
      CHARACTER*1  FTIME(10),FYRMODIR(80)
!=== Temporary tile space transfer files (different than in lis_module)      
      real, allocatable :: tmptilen(:)
      CHARACTER (LEN=100) :: temp
!=== End Variable Definition =============================================

      if (masterproc) then 
!=== Restart Writing (2 files are written = active and archive)

         IF ((lis%t%gmt.eq.(24-hyssibdrv%writeint)) & 
            .or.lis%t%endtime.eq.1) THEN
            allocate(tmptilen(lis%d%nch))
 !Active archive restart
            OPEN(40,FILE=HYSSIBDRV%HYSSIB_RFILE,FORM='unformatted')
!   OPEN(33,FILE="timmgr.bin",FORM='unformatted') !Active archive restart
            call timemgr_write_restart(40)
!     WRITE(40) LDAS%T%YR,LDAS%T%MO,LDAS%T%DA,LDAS%T%HR,LDAS%T%MN, &
! LDAS%T%SS,LDAS%P%VCLASS,LDAS%D%lnc,LDAS%D%lnr,LDAS%D%NCH !Veg class,no.tiles

  !Veg class, no tiles
            WRITE(40) LIS%P%VCLASS,LIS%D%lnc,LIS%D%lnr,LIS%D%NCH

!== NOTE: Next four write statements originally commented out in LIS
            WRITE(40) TILE%COL  !Grid Col of Tile   
            WRITE(40) TILE%ROW  !Grid Row of Tile
            WRITE(40) TILE%FGRD !Fraction of Grid covered by tile
            WRITE(40) TILE%VEGT !Vegetation Type of Tile
            WRITE(40) HYSSIB%TC
            WRITE(40) HYSSIB%TG
            WRITE(40) HYSSIB%TSN
            WRITE(40) HYSSIB%TD
            DO L=1,3
               DO T=1,LIS%D%NCH
                  TMPTILEN(T)=HYSSIB(T)%WWW(L)
               ENDDO
               WRITE(40) TMPTILEN !HY-SSiB Soil Moisture (3 layers)
            ENDDO
            DO L=1,2
               DO T=1,LIS%D%NCH
                  TMPTILEN(T)=HYSSIB(T)%CAPAC(L)
               ENDDO
               WRITE(40) TMPTILEN !HY-SSiB Water on canopy/ground (2 layers)
            ENDDO
            DO L=1,2
               DO T=1,LIS%D%NCH
                  TMPTILEN(T)=HYSSIB(T)%SNOW(L)
               ENDDO
               WRITE(40) TMPTILEN !HY-SSiB Snow on canopy/ground (2 layers)
            ENDDO
            WRITE(40) HYSSIB%SGFG
            WRITE(40) HYSSIB%SDENS

            CLOSE(40)

            WRITE(*,*) 'HY-SSiB Active Restart Written: ', &
                       HYSSIBDRV%HYSSIB_RFILE

!        call ESMF_LogInfo(diaglog, 
!     &       "SSIB Active Restart Written: %40s\n",ssibdrv%ssib_rfile)
!        WRITE(79,*)'SSIB Active Restart Written: ',SSIBDRV%SSIB_RFILE

            WRITE(unit=temp,fmt='(i4,i2,i2,i2)') LIS%T%YR,LIS%T%MO, & 
                                                 LIS%T%DA,LIS%T%HR
            READ(UNIT=temp,fmt='(10A1)') FTIME
            DO I=1,10
               IF (FTIME(I).EQ.(' ')) FTIME(I)='0'
            ENDDO
            WRITE(UNIT=temp,fmt='(A4,I3,A8,I4,A1,I4,I2,I2,A8,I3,A1)') &
                  '/EXP',LIS%O%EXPCODE,'/HYSSIB/',LIS%T%YR,'/',LIS%T%YR, &
                  LIS%T%MO,LIS%T%DA,'/LIS.EXP',LIS%O%EXPCODE,'.'
            READ(UNIT=temp,fmt='(80a1)') (FNAME(I),I=1,40)
            DO I=1,75
               IF (FNAME(I).EQ.(' ')) FNAME(I)='0'
            ENDDO

            WRITE(UNIT=temp,fmt='(a9)') 'mkdir -p '
            READ(UNIT=temp,fmt='(80a1)') (FMKDIR(I),I=1,9)
            WRITE(UNIT=temp,fmt='(A4,I3,A8,I4,A1,I4,I2,I2)') & 
            '/EXP',LIS%O%EXPCODE,'/HYSSIB/', & 
            LIS%T%YR,'/',LIS%T%YR,LIS%T%MO,LIS%T%DA
            READ(UNIT=temp,fmt='(80a1)') (FYRMODIR(I),I=1,28)
            DO I=1,28
               IF (FYRMODIR(I).EQ.(' ')) FYRMODIR(I)='0'
            ENDDO

            WRITE(UNIT=temp,fmt='(a10)')'.HYSSIBrst'
            READ(UNIT=temp,FMT='(80A1)') (FSUBS(I),I=1,10)

            WRITE(UNIT=temp,fmt='(a40)') LIS%O%ODIR
            READ(UNIT=temp,FMT='(80A1)') (FBASE(I),I=1,80)
            C=0
            DO I=1,80
               IF (FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1
            ENDDO
            WRITE(UNIT=temp,FMT='(80A1)') (FBASE(I),I=1,C), &
                 (FNAME(I),I=1,40),(FTIME(I),I=1,10),(FSUBS(I),I=1,10)
            READ(UNIT=temp,fmt='(a80)') FILEN

            WRITE(UNIT=temp,FMT='(80A1)')(FMKDIR(I),I=1,9), &
                 (FBASE(I),I=1,C),(FYRMODIR(I),I=1,28)
            READ(UNIT=temp,fmt='(a80)') MKFYRMO

!== Archive File Name Generation Complete
!== Make the directories for the SSIB restart file
!       IF(LDAS%GMT.LE.LDAS%WRITEINTN)THEN             
            CALL SYSTEM(MKFYRMO)
!       ENDIF

!== Archive File Name Generation Complete
            print *,filen
            OPEN(40,FILE=FILEN,status='unknown',FORM='unformatted')
!        WRITE(40) LDAS%T%YR,LDAS%T%MO,LDAS%T%DA,LDAS%T%HR,LDAS%T%MN, &
! LDAS%T%SS, LDAS%P%VCLASS,LDAS%D%lnc,LDAS%D%lnr,LDAS%D%NCH !Veg class,no.tiles
            call timemgr_write_restart(40)
!Veg class, no. tiles
            WRITE(40) LIS%P%VCLASS,LIS%D%lnc,LIS%D%lnr,LIS%D%NCH
!== NOTE: Next four write statements originally commented out in LIS
            WRITE(40) TILE%COL  !Grid Col of Tile   
            WRITE(40) TILE%ROW  !Grid Row of Tile
            WRITE(40) TILE%FGRD !Fraction of Grid covered by tile
            WRITE(40) TILE%VEGT !Vegetation Type of Tile
            WRITE(40) HYSSIB%TC
            WRITE(40) HYSSIB%TG
            WRITE(40) HYSSIB%TSN
            WRITE(40) HYSSIB%TD
            DO L=1,3
               DO T=1,LIS%D%NCH
                  TMPTILEN(T)=HYSSIB(T)%WWW(L)
               ENDDO
               WRITE(40) TMPTILEN !HY-SSiB Soil Moisture (3 layers)
            ENDDO
            DO L=1,2
               DO T=1,LIS%D%NCH
                  TMPTILEN(T)=HYSSIB(T)%CAPAC(L)
               ENDDO
               WRITE(40) TMPTILEN !HY-SSiB Water on canopy/ground (2 layers)
            ENDDO
            DO L=1,2
               DO T=1,LIS%D%NCH
                  TMPTILEN(T)=HYSSIB(T)%SNOW(L)
               ENDDO
               WRITE(40) TMPTILEN !HY-SSiB Water on canopy/ground (2 layers)
            ENDDO
            WRITE(40) HYSSIB%SGFG
            WRITE(40) HYSSIB%SDENS

            CLOSE(40)

            WRITE(*,*)'HY-SSiB Archive Restart Written: ',FILEN
!        call ESMF_LogInfo(diaglog, 
!     &       "SSIB Archive Restart Written: %40s\n",filen)
!        WRITE(79,*)'SSIB Archive Restart Written: ',FILEN
            deallocate(tmptilen)
         ENDIF                  !End restart writing (active and archive files)
      endif                     ! if masterproc

!      endif

      return
      end SUBROUTINE HYSSIB_WRITERST

