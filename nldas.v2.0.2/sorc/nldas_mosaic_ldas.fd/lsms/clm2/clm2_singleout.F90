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
! !ROUTINE: clm2_singleout.F90
!
! !DESCRIPTION:  
!  Write output file for a single CLM variable
!
! !REVISION HISTORY:
! 14 Jun 2002; Sujay Kumar  Initial Specification
!
! !INTERFACE:
subroutine clm2_singleout (var_array, index)
! !USES:
  use lisdrv_module, only : lis, tile
  use clm_varcon, ONLY : denh2o, denice, hvap, hsub, hfus, istwet 
  use clm_varpar, ONLY : nlevsoi
  use clm_varmap, ONLY : patchvec
  use clm_varctl, only : clmdrv
  use drv_output_mod, only : t2gr
!EOP
  implicit none
  real :: var_array(lis%d%glbnch)
  integer :: index,length
!=== Local variables =====================================================
  integer :: t,c,r,m,i,j,flag,tt,pp
  character(len=100)::temp
!=== Temporary transfer variables 
  !real :: tempvar(lis%d%glbnch)
  real,allocatable,dimension(:)  :: gtmp!(lis%d%glbngrid)
  real,allocatable,dimension(:,:) :: g2tmp!(lis%d%lnc,lis%d%lnr)
  integer :: idisk, nrec ! disk to use for output, direct access record number

!=== Variables used for statistical summary
  real vmean,vstdev,vmin,vmax

  character*80 fileng,filent,filent3,mkfyrmo,cdir,namet,nameg,namet3
  character*82 filengb
  character*1  fname(80),fbase(40),fsubst(80),fmkdir(80)
  character*1  ftime(8),fyrmodir(80),fsubsg(80),fsubst3(80)
  character*1  fcd(3),frm(3),flats(13),ftimeb(10),fsubgb(9)
  character*1 ftimec(4)
  character*80 :: vname1(55)
  character*80 :: vname(55)
  character*1  fvarname(12)
  integer,parameter :: nvarsg=46,nvarst=45,nvarst3=9
  character*12 :: levunits
  integer :: kbeg,timinc
  integer :: prec,kount
  
  DATA VNAME / "SWnet(W/m2)",          &
               "LWnet(W/m2)",          &
               "Qle(W/m2)",            &
               "Qh(W/m2)",             &
               "Qg(W/m2)",             &
               "Snowf(kg/m2s)",        &
               "Rainf(kg/m2s)",        &
               "Evap(kg/m2s)",         &
               "Qs(kg/m2s)",           &
               "Qsb(kg/m2s)",          &
               "Qsm(kg/m2s)",          &
               "DelSoilMoist(kg/m2)",  &
               "DelSWE(kg/m2)",        &
               "SnowT(K)",             &
               "VegT(K)",              &
               "BareSoilT(K)",         &
               "AvgSurfT(K)",          &
               "RadT(K)",              &
               "Albedo(-)",            &
               "SWE(kg/m2)",           &
               "SoilTemp1(K)",         &
               "SoilTemp2(K)",         &
               "SoilTemp3(K)",         &
               "SoilTemp4(K)",         &
               "SoilTemp5(K)",         &
               "SoilTemp6(K)",         &
               "SoilTemp7(K)",         &
               "SoilTemp8(K)",         &
               "SoilTemp9(K)",         &
               "SoilTemp10(K)",        &
               "SoilMoist1(kg/m2)",    &
               "SoilMoist2(kg/m2)",    &
               "SoilMoist3(kg/m2)",    &
               "SoilMoist4(kg/m2)",    &
               "SoilMoist5(kg/m2)",    &
               "SoilMoist6(kg/m2)",    &
               "SoilMoist7(kg/m2)",    &
               "SoilMoist8(kg/m2)",    &
               "SoilMoist9(kg/m2)",    &
               "SoilMoist10(kg/m2)",   &
               "RootMoist(kg/m2)",     &
               "SoilWet(-)",           &
               "TVeg(kg/m2s)",         &
               "ECanop(kg/m2s)",       &
               "ESoil(kg/m2s)",        &
               "Canopint(kg/m2)",      &
               "ACond(m/s)",           &
               "Wind(m/s)",            &
               "Rainfforc(kg/m2s)",    &
               "Snowfforc(kg/m2s)",    &
               "Tair(K)",              &
               "Qair(kg/kg)",          &
               "Psurf(Pa)",            &
               "SWdown(W/m2)",         &
               "LWdown(W/m2)" /

  DATA VNAME1 / "swnet",         &
                "lwnet",         &
                "qle",           &
                "qh",            &
                "qg",            &
                "snowf",         &
                "rainf",         &
                "evap",          &
                "qs",            &
                "qsb",           &
                "qsm",           &
                "delsoilmoist",  &
                "delswe",        &
                "snowt",         &
                "vegt",          &
                "baresoilt",     &
                "avgsurft",      &
                "radt",          &
                "albedo",        &
                "swe",           &
                "soiltemp1",     &
                "soiltemp2",     &
                "soiltemp3",     &
                "soiltemp4",     &
                "soiltemp5",     &
                "soiltemp6",     &
                "soiltemp7",     &
                "soiltemp8",     &
                "soiltemp9",     &
                "soiltemp10",    &
                "soilmoist1",    &
                "soilmoist2",    &
                "soilmoist3",    &
                "soilmoist4",    &
                "soilmoist5",    &
                "soilmoist6",    &
                "soilmoist7",    &
                "soilmoist8",    &
                "soilmoist9",    &
                "soilmoist10",   &
                "rootmoist",     &
                "soilwet",       &
                "tveg",          &
                "ecanop",        &
                "esoil",         &
                "canopint",      &
                "acond",         &
                "wind",          &
                "rainfforc",     &
                "snowfforc",     &
                "tair",          &
                "qair",          &
                "psurf",         &
                "swdown",        &
                "lwdown" /

  character*40 file
  character*80 name
  
!=== End Variable List =========================================================
!BOC
!----------------------------------------------------------------------
! Test to see if output writing interval has been reached
!----------------------------------------------------------------------
  if(mod(lis%t%gmt,clmdrv%writeintc2).eq.0)then
!----------------------------------------------------------------------
! Generate directory structure and file names for CLM Output 
!----------------------------------------------------------------------
     length = len(trim(vname1(index)))
     WRITE(UNIT=temp, FMT='(A12)') VNAME1(index)
     READ(UNIT=temp,FMT='(12A1)') (FVARNAME(I), I=1,length)
     WRITE(unit=temp,fmt='(I4,I2,I2)')lis%T%YR,lis%T%MO,lis%T%DA
     READ(unit=temp,fmt='(8a1)')FTIME
     DO I=1,8
        IF(FTIME(I).EQ.(' '))FTIME(I)='0'
     ENDDO
     
     WRITE(unit=temp,fmt='(I4)')lis%T%YR
     READ(unit=temp,fmt='(8a1)')FTIMEC
     DO I=1,4
        IF(FTIMEC(I).EQ.(' '))FTIMEC(I)='0'
     ENDDO
     
#if 0
     WRITE(unit=temp,fmt='(a6,i3,a1)')'/LIS.E',lis%O%EXPCODE,'.'
     READ(unit=temp,fmt='(80a1)') (FNAME(I),I=1,10)
     DO I=1,10
        IF(FNAME(I).EQ.(' '))FNAME(I)='0'
     ENDDO
#endif
     
!     idisk = mod(index, lis%o%odirn)
     idisk = 1
!     if ( idisk == 0 ) then
!        idisk = lis%o%odirn
!     endif
!     WRITE(unit=temp,fmt='(a40)') lis%O%ODIR_ARRAY(idisk)
     WRITE(unit=temp,fmt='(a40)') lis%O%ODIR
     READ(unit=temp,fmt='(40a1)') (FBASE(I),I=1,40)
     C=0
     DO I=1,40
        IF(FBASE(I).EQ.(' ').AND.C.EQ.0)C=I-1
     ENDDO
     
     WRITE(unit=temp,fmt='(A4,I3,A6,I4,A1,I4,I2,I2)')'/EXP', & 
          lis%O%EXPCODE,'/CLM2/', & 
          lis%t%YR,'/',lis%T%YR,lis%T%MO,lis%T%DA
     READ(unit=temp,fmt='(80A1)') (FYRMODIR(I),I=1,26)
     DO I=1,26
        IF(FYRMODIR(I).EQ.(' '))FYRMODIR(I)='0'
     ENDDO
     
     WRITE(unit=temp,fmt='(A9)')'mkdir -p '
     READ(unit=temp,fmt='(80A1)')(FMKDIR(I),I=1,9)
     
     WRITE(unit=temp,fmt='(80A1)')(FMKDIR(I),I=1,9),(FBASE(I),I=1,C), & 
          (FYRMODIR(I),I=1,26)
     READ(unit=temp,fmt='(A80)')MKFYRMO
     
!----------------------------------------------------------------------     
! Make the directories for the CLM2 output files              
!----------------------------------------------------------------------
     call system(mkfyrmo)

!----------------------------------------------------------------------
! Generate file name for binary output
!----------------------------------------------------------------------
     if(lis%o%wout.eq.1)then
        write(unit=temp,fmt='(I4,I2,I2,I2)')lis%t%yr, & 
             lis%t%mo,lis%t%da,lis%t%hr
        read(unit=temp,fmt='(10A1)')ftimeb
        do i=1,10
           if(ftimeb(i).eq.(' '))then
              ftimeb(i)='0'
           endif
        enddo
        
        !write(unit=temp,fmt='(A9)')'.CLM2gbin'
#if ( defined FARMER_DOG_BONES )
        write(unit=temp,fmt='(A5)')'.gd4r'
        read(unit=temp,fmt='(80A1)') (fsubgb(i),i=1,5)
#else
        write(unit=temp,fmt='(A5)')'.ls4r'
        read(unit=temp,fmt='(80A1)') (fsubgb(i),i=1,5)
#endif
        
#if 0
        write(unit=temp,fmt='(82A1)')(fbase(i),i=1,c), & 
             (FYRMODIR(I),I=1,26), & 
             (fname(i),i=1,10),(ftimeb(i),i=1,10), & 
             (fvarname(i),i=1,length),(fsubgb(i),i=1,9 )
        read(unit=temp,fmt='(A82)')filengb
#endif
        write(unit=temp,fmt='(69A1)')(fbase(i),i=1,c), & 
             (FYRMODIR(I),I=1,26), '/', & 
             (ftimeb(i),i=1,10), & 
             (fvarname(i),i=1,length),(fsubgb(i),i=1,5 )
        read(unit=temp,fmt='(A69)')filengb
     endif
     if(lis%o%wout.eq.1)then
       clmdrv%numout=clmdrv%numout+1 
#if ( defined FARMER_DOG_BONES )
       allocate(g2tmp(lis%d%lnc,lis%d%lnr))
       g2tmp = lis%d%UDEF
       do i = 1, lis%d%glbnch
          g2tmp(tile(i)%col, tile(i)%row) = var_array(i)
       enddo

       open(58,file=filengb,form='unformatted',access='direct', &
            recl=lis%d%lnc * lis%d%lnr * 4)
       write(58, rec=1) g2tmp
       deallocate(g2tmp)
#else
       open(57,file=filengb,form='unformatted')
       allocate(gtmp(lis%d%glbngrid))
       call t2gr(var_array,gtmp,lis%d%glbngrid,lis%d%glbnch,tile)
       write(57) gtmp
       deallocate(gtmp)
#endif
!----------------------------------------------------------------------        
! Write statistical output
!----------------------------------------------------------------------     
#if ( defined FARMER_DOG_BONES )
        call lis_log_msg("DBG: clm2_singleout -- farmer-dog-bones "// &
                         "running mode cannot write stats file")
#else
        if(clmdrv%clm2open.eq.0)then
           file='CLMstats.dat'
           call openfile(name,lis%o%odir,lis%o%expcode,file)
           if(lis%o%startcode.eq.1)then
              open(60,file=name,form='formatted',status='unknown', & 
                   position='append')
           else
              open(60,file=name,form='formatted',status='replace')
           endif
           clmdrv%clm2open=1
        endif
       
        write(60,996)'       Statistical Summary of CLM Output for:  ', & 
             lis%t%mo,'/',lis%t%da,'/',lis%t%yr,lis%t%hr,':', & 
             lis%t%mn,':',lis%t%ss 
        
996     format(a47,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
997     format(t26,'Mean',t40,'StDev',t54,'Min',t68,'Max')

        call stats(var_array,lis%d%udef,lis%d%glbnch,vmean,vstdev,vmin, & 
             vmax)
        write(60,999) vname(index),vmean,vstdev,vmin,vmax
#endif
     endif
  endif
995 format (1x,a10,I1,a9,4f14.3)
999 format (1x,a15,4f14.3)
998 format (1x,a15,4e14.3)
!EOC  
end subroutine clm2_singleout
