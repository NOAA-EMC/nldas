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
! !ROUTINE: hyssib_almaout.F90
!
! !DESCRIPTION:
!  LIS HY-SSiB data writer:  Binary and stat files in ALMA convention
!
! !REVISION HISTORY:
!  4 Nov 1999: Jon Radakovich, Initial Code
! 28 Apr 2002: Kristi Arsenault, Added SSIB LSM to LDAS
! 15 Jun 2003: Sujay Kumar, ALMA version
!    Feb 2004: David Mocko, Conversion from SSiB to HY-SSiB
! 
! !INTERFACE:
      subroutine hyssib_almaout()
! !USES:
      use lisdrv_module, only : lis
      use hyssib_varder         ! HY-SSiB-specific variables

      implicit none 
!EOP

      INTEGER :: T,C,R,M,I,N           
      CHARACTER*80 MKFYRMO,FILENMT,FILENMG,CDIR,NAMET,NAMEG,FILENGB
      CHARACTER*80 MKFYRMO2
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13),FTIMEC(4)
      CHARACTER*1  FYRMODIR(26),FSUBFT(80)
      CHARACTER*1  FSUBFG(80),FTIMEB(10),FSUBGB(11)

      CHARACTER (LEN=100) :: FBINNAME
      CHARACTER (LEN=100) :: temp1

      INTEGER,PARAMETER :: NVARSG=29,NVARST=29,KMG=1
      CHARACTER*80 :: VNAME(NVARSG)
      INTEGER      :: PREC,KBEGT,KOUNTT

      DATA VNAME / "SWnet(W/m2)","LWnet(W/m2)", &
                   "Qle(W/m2)","Qh(W/m2)","Qg(W/m2)", &
                   "Snowf(kg/m2s)","Rainf(kg/m2s)","Evap(kg/m2s)", &
                   "Qs(kg/m2s)","Qsb(kg/m2s)","Qsm(kg/m2s)", &
                   "DelSoilMoist(kg/m2)","DelSWE(kg/m2)", &
                   "SnowT(K)","VegT(K)","BareSoilT(K)","AvgSurfT(K)", &
                   "RadT(K)","Albedo(-)","SWE(kg/m2)", &
                   "SoilMoist1(kg/m2)","SoilMoist2(kg/m2)", &
                   "SoilMoist3(kg/m2)","SoilMoist4(kg/m2)","SoilWet(-)", &
                   "TVeg(kg/m2s)","ESoil(kg/m2s)","RootMoist(kg/m2)", &
                   "ACond(m/s)" /

      CHARACTER*40 FILE
      CHARACTER*80 NAME
!BOC
!-------------------------------------------------------------------------
! Test to see if output writing interval has been reached
!-------------------------------------------------------------------------
      IF (MOD(LIS%T%GMT,hyssibdrv%WRITEINT).EQ.0) THEN
         hyssibdrv%NUMOUT = hyssibdrv%NUMOUT+1    
         WRITE(UNIT=temp1,FMT='(I4,I2,I2)') LIS%T%YR,LIS%T%MO,LIS%T%DA
         READ(UNIT=temp1,FMT='(8A1)') FTIME
         DO I=1,8
            IF (FTIME(I).EQ.(' ')) FTIME(I)='0'
         ENDDO
         WRITE(UNIT=temp1,FMT='(I4)') LIS%T%YR
         READ(UNIT=temp1,FMT='(8A1)') FTIMEC
         DO I=1,4
            IF (FTIMEC(I).EQ.(' ')) FTIMEC(I)='0'
         ENDDO

         WRITE(UNIT=temp1,FMT='(A8,I3,A1)') '/LIS.EXP',LIS%O%EXPCODE,'.'
         READ(UNIT=temp1,FMT='(80A1)') (FNAME(I),I=1,12)
         DO I=1,12
            IF (FNAME(I).EQ.(' ')) FNAME(I)='0'
         ENDDO

         WRITE(UNIT=temp1,FMT='(A40)') LIS%O%ODIR
         READ(UNIT=temp1,FMT='(40A1)') (FBASE(I),I=1,40)
         C=0
         DO I=1,40
            IF (FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1
         ENDDO

         WRITE(UNIT=temp1,FMT='(A4,I3,A8,I4,A1,I4,I2,I2)') '/EXP', &
               LIS%O%EXPCODE,'/HYSSIB/', & 
               LIS%T%YR,'/',LIS%T%YR,LIS%T%MO,LIS%T%DA
         READ(UNIT=temp1,FMT='(80A1)') (FYRMODIR(I),I=1,28)
         DO I=1,28
            IF (FYRMODIR(I).EQ.(' ')) FYRMODIR(I)='0'
         ENDDO

         WRITE(UNIT=temp1,FMT='(A9)')'mkdir -p '
         READ(UNIT=temp1,FMT='(80A1)')(FMKDIR(I),I=1,9)

         WRITE(UNIT=temp1,FMT='(80A1)')(FMKDIR(I),I=1,9),(FBASE(I),I=1,C), &
                                       (FYRMODIR(I),I=1,28)
         READ(UNIT=temp1,FMT='(A80)') MKFYRMO
         CALL SYSTEM(MKFYRMO)
!----------------------------------------------------------------------
! Generate file name for BINARY output
!----------------------------------------------------------------------
         IF (LIS%O%WOUT.EQ.1) THEN
            WRITE(UNIT=FBINNAME, FMT='(I4,I2,I2,I2)') LIS%T%YR,LIS%T%MO, &
                                                      LIS%T%DA,LIS%T%HR
            READ(UNIT=FBINNAME,FMT='(10A1)') FTIMEB
            DO I=1,10
               IF (FTIMEB(I).EQ.(' ')) FTIMEB(I)='0'
            ENDDO
            WRITE(UNIT=FBINNAME,FMT='(A11)') '.HYSSIBgbin'
            READ(UNIT=FBINNAME,FMT='(80A1)') (FSUBGB(I),I=1,11)

            WRITE(UNIT=FBINNAME,FMT='(80A1)')(FBASE(I),I=1,C), & 
                       (FYRMODIR(I),I=1,28), & 
                       (FNAME(I),I=1,12),(FTIMEB(I),I=1,10), & 
                       (FSUBGB(I),I=1,11)
            READ(UNIT=FBINNAME,FMT='(A80)') FILENGB
!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
            IF (hyssibdrv%HYSSIBopen.EQ.0) THEN
               FILE='HYSSIBstats.dat'
               CALL OPENFILE(NAME,LIS%O%ODIR,LIS%O%EXPCODE,FILE)
               IF (LIS%O%STARTCODE.EQ.1) THEN
                  OPEN(65,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN', &
                          POSITION='APPEND')
               ELSE
                  OPEN(65,FILE=NAME,FORM='FORMATTED',STATUS='REPLACE')
               ENDIF
               hyssibdrv%HYSSIBopen=1
            ENDIF

            WRITE(65,996)'     Statistical Summary of HY-SSiB Output for:  ', &
                  LIS%T%MO,'/',LIS%T%DA,'/',LIS%T%YR,LIS%T%HR,':',LIS%T%MN,':', &
                  LIS%T%SS
 996        FORMAT(A47,I2,A1,I2,A1,I4,1X,I2,A1,I2,A1,I2)
            WRITE(65,*) 
            WRITE(65,997)
 997        FORMAT(T27,'Mean',T41,'StDev',T56,'Min',T70,'Max')
         endif
         if (lis%o%wout.eq.1) then
            open(58,file=filengb,form='unformatted')
         endif
         if (lis%o%wout.eq.1) then
            call hyssib_binout(58) 
         endif
         call hyssib_writestats(65)
         hyssib%count=0           !reset counters

         write(65,*)
         write(65,*)
      endif
!EOC
      end subroutine hyssib_almaout

