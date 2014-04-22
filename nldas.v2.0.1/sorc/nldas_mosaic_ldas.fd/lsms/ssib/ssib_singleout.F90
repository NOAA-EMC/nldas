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
! !ROUTINE: ssib_singleout.F90
!
! !DESCRIPTION:
!  Write output file for a single HY-SSiB variable
!
! !REVISION HISTORY:
! 14 Jun 2002: Sujay Kumar, Initial Specification
! 22 May 2004: David Mocko, Conversion from HY-SSiB to SSiB
!
! !INTERFACE:
      subroutine ssib_singleout(ld,tile,gindex,var_array,index)
! !USES:
      use lis_module            ! LIS non-model-specific 1-D variables
      use tile_module           ! LIS non-model-specific tile variables
      use ssib_varder, only : ssibdrv
      use drv_output_mod, only : t2gr

      implicit none
! !ARGUMENTS:
      type (lisdec)  :: ld       !data structure for lis domain specific variables
      type (tiledec) :: tile(ld%d%glbnch) !tile array for the modeled domain
      integer        :: gindex(ld%d%lnc, ld%d%lnr) !2-d array for mapping from 2d to 1d
      real           :: var_array(ld%d%glbnch) !array of variable that is being output
      integer        :: index   !Index of the output varible in the ALMA list.
!EOP
!=== Local Variables ==================================================
      INTEGER :: T,C,R,M,I,N,length
      CHARACTER*80 MKFYRMO,FILENMT,FILENMG,CDIR,NAMET,NAMEG,FILENGB
      CHARACTER*80 MKFYRMO2
      CHARACTER*1  FNAME(80),FBASE(40),FMKDIR(80)
      CHARACTER*1  FTIME(8),FCD(3),FRM(3),FLATS(13),FTIMEC(4)
      CHARACTER*1  FYRMODIR(26),FSUBFT(80)
      CHARACTER*1  FVARNAME(10)
      CHARACTER*1  FSUBFG(80),FTIMEB(10),FSUBGB(9)

      CHARACTER (LEN=100) :: FBINNAME
      CHARACTER (LEN=100) :: temp1
!=== Variables used for writing output in HDF format
      INTEGER,PARAMETER :: NVARSG=36,NVARST=36,KMG=1
      CHARACTER*80 :: VNAME(NVARSG)
      CHARACTER*80 :: VNAME1(NVARSG)
      INTEGER      :: PREC,KBEGT,KOUNTT
      real :: gtmp (ld%d%glbngrid)
      REAL :: VMEAN,VSTDEV,VMIN,VMAX

      DATA VNAME / "SWnet(W/m2)","LWnet(W/m2)","Qle(W/m2)","Qh(W/m2)", &
                   "Qg(W/m2)","Qf(W/m2)","Qtau(N/m2)",                 &
                   "DelSurfHeat(J/m2)","Snowf(kg/m2s)",                &
                   "Rainf(kg/m2s)","Evap(kg/m2s)","Qs(kg/m2s)",        &
                   "Qsb(kg/m2s)","Qsm(kg/m2s)","DelSoilMoist(kg/m2)",  &
                   "DelSWE(kg/m2)","DelIntercept(kg/m2)","VegT(K)",    &
                   "BaresoilT(K)","AvgSurfT(K)","RadT(K)","Albedo(-)", &
                   "SWE(kg/m2)","SWEVeg(kg/m2)","SoilMoist1(kg/m2)",   &
                   "SoilMoist2(kg/m2)","SoilMoist3(kg/m2)",            &
                   "SoilTemp(K)","SoilWet(-)","ECanop(kg/m2s)",        &
                   "TVeg(kg/m2s)","ESoil(kg/m2s)","RootMoist(kg/m2)",  &
                   "CanopInt(kg/m2)","ACond(m/s)","SnowFrac(-)"        /

      DATA VNAME1 / "SWnet","LWnet","Qle","Qh","Qg","Qf","Qtau",       &
                    "DelSurfHeat","Snowf","Rainf","Evap","Qs","Qsb",   &
                    "Qsm","DelSoilMoist","DelSWE","DelIntercept",      &
                    "VegT","BaresoilT","AvgSurfT","RadT","Albedo",     &
                    "SWE","SWEVeg","SoilMoist1","SoilMoist2",          &
                    "SoilMoist3","SoilTemp","SoilWet","ECanop","TVeg", &
                    "ESoil","RootMoist","CanopInt","ACond","SnowFrac"  /

      CHARACTER*40 FILE
      CHARACTER*80 NAME
!BOC
!-------------------------------------------------------------------------
! Test to see if output writing interval has been reached
!-------------------------------------------------------------------------
      IF (MOD(LD%T%GMT,ssibdrv%WRITEINTN).EQ.0) THEN
         ssibdrv%NUMOUTNH=ssibdrv%NUMOUTNH+1
         length = len(trim(vname1(index)))
         WRITE(UNIT=temp1,FMT='(A10)') VNAME1(index)
         READ(UNIT=temp1,FMT='(10A1)') (FVARNAME(I),I=1,length)
         WRITE(UNIT=temp1,FMT='(I4,I2,I2)') LD%T%YR,LD%T%MO,LD%T%DA
         READ(UNIT=temp1,FMT='(8A1)') FTIME
         DO I=1,8
            IF (FTIME(I).EQ.(' ')) FTIME(I)='0'
         ENDDO
         WRITE(UNIT=temp1,FMT='(I4)') LD%T%YR
         READ(UNIT=temp1,FMT='(8A1)') FTIMEC
         DO I=1,4
            IF (FTIMEC(I).EQ.(' ')) FTIMEC(I)='0'
         ENDDO

         WRITE(UNIT=temp1,FMT='(A8,I3,A1)') '/LIS.EXP',LD%O%EXPCODE,'.'
         READ(UNIT=temp1,FMT='(80A1)') (FNAME(I),I=1,12)
         DO I=1,12
            IF (FNAME(I).EQ.(' ')) FNAME(I)='0'
         ENDDO

         WRITE(UNIT=temp1,FMT='(A40)') LD%O%ODIR
         READ(UNIT=temp1,FMT='(40A1)') (FBASE(I),I=1,40)
         C=0
         DO I=1,40
            IF (FBASE(I).EQ.(' ').AND.C.EQ.0) C=I-1
         ENDDO

         WRITE(UNIT=temp1,FMT='(A4,I3,A6,I4,A1,I4,I2,I2)')'/EXP', &
               LD%O%EXPCODE,'/SSIB/', & 
               LD%T%YR,'/',LD%T%YR,LD%T%MO,LD%T%DA
         READ(UNIT=temp1,FMT='(80A1)') (FYRMODIR(I),I=1,26)
         DO I=1,26
            IF (FYRMODIR(I).EQ.(' ')) FYRMODIR(I)='0'
         ENDDO

         WRITE(UNIT=temp1,FMT='(A9)')'mkdir -p '
         READ(UNIT=temp1,FMT='(80A1)')(FMKDIR(I),I=1,9)

         WRITE(UNIT=temp1,FMT='(80A1)')(FMKDIR(I),I=1,9),(FBASE(I),I=1,C), &
                                       (FYRMODIR(I),I=1,26)
         READ(UNIT=temp1,FMT='(A80)') MKFYRMO
         CALL SYSTEM(MKFYRMO)
!----------------------------------------------------------------------
! Generate file name for BINARY output
!----------------------------------------------------------------------
         IF (LD%O%WOUT.EQ.1) THEN
            WRITE(UNIT=FBINNAME,FMT='(I4,I2,I2,I2)') LD%T%YR,LD%T%MO, &
                                                     LD%T%DA,LD%T%HR
            READ(UNIT=FBINNAME,FMT='(10A1)') FTIMEB
            DO I=1,10
               IF (FTIMEB(I).EQ.(' ')) FTIMEB(I)='0'
            ENDDO
            WRITE(UNIT=FBINNAME,FMT='(A9)') '.SSIBgbin'
            READ(UNIT=FBINNAME,FMT='(80A1)') (FSUBGB(I),I=1,9)
            WRITE(UNIT=FBINNAME,FMT='(80A1)')(FBASE(I),I=1,C), & 
                       (FYRMODIR(I),I=1,26), & 
                       (FNAME(I),I=1,12),(FTIMEB(I),I=1,10), & 
                       (FVARNAME(I),I=1,length),(FSUBGB(I),I=1,9)
            READ(UNIT=FBINNAME,FMT='(A80)') FILENGB
!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
            IF (ssibdrv%SSIBopen.EQ.0) THEN
               FILE='SSIBstats.dat'
               CALL OPENFILE(NAME,LD%O%ODIR,LD%O%EXPCODE,FILE)
               IF (LD%O%STARTCODE.EQ.1) THEN
                  OPEN(65,FILE=NAME,FORM='FORMATTED',STATUS='UNKNOWN', &
                          POSITION='APPEND')
               ELSE
                  OPEN(65,FILE=NAME,FORM='FORMATTED',STATUS='REPLACE')
               ENDIF
               ssibdrv%SSIBopen=1
            ENDIF

            WRITE(65,996)'        Statistical Summary of SSiB Output for:  ', &
                  LD%T%MO,'/',LD%T%DA,'/',LD%T%YR,LD%T%HR,':',LD%T%MN,':', &
                  LD%T%SS
 996        FORMAT(A47,I2,A1,I2,A1,I4,1X,I2,A1,I2,A1,I2)
            WRITE(65,*)
            WRITE(65,997)
 997        FORMAT(T27,'Mean',T41,'StDev',T56,'Min',T70,'Max')
         ENDIF

         if (ld%o%wout.eq.1) then
            open(58,file=filengb,form='unformatted')
         endif
         if (ld%o%wout.eq.1) then
            call t2gr(var_array,gtmp,ld%d%glbngrid,ld%d%glbnch,tile); 
            write(58) gtmp
            call stats(var_array,ld%d%udef,ld%d%glbnch, vmean,vstdev, vmin, vmax); 
            write(65,999) vname(index),vmean,vstdev,vmin,vmax
         endif

 998     FORMAT(1X,A18,4E14.3)
 999     FORMAT(1X,A18,4F14.3)
      endif
      return
!EOC
      end subroutine ssib_singleout

