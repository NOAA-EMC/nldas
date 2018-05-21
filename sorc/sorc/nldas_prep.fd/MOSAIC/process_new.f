      subroutine process(nch,dtstep,sfrac,jday,sunang,alon,alat,ityp
     i,        tmp2m,dpt2m,sfcprs,cpcp,tpcp,alwdn,swdn,wnd
     i,        bf1,bf2,bf3,vgwmax,cdcr1,cdcr2
     i,        psis,bee,poros,wpwet,cond,gnu
     i,        ars1,ars2,ars3,ara1,ara2,ara3,ara4,arw1,arw2,arw3,arw4 
     I,        tsa1,tsa2,tsb1,tsb2
     i,        albed,  green,    zlt,     z0,      d
     i,        bug,bug1,bug2
     u,        tc1,tc2,tc4,qa1,qa2,qa4,catdef,rzex,srfex,capac
     u,        ght1,ght2,ght3,ght4,ght5,ght6,TSURF
     u,        wesn1,wesn2,wesn3, htsnn1,htsnn2,htsnn3
     u,        sndz1,sndz2,sndz3
     o,        evap,shflux,runoff,eint,esoi,eveg,esno,snowterm
     o,        hlwup, smelt, bflow, runsrf, ar1, ar2, rzeq, qinfil
     o,        ghflux, tmsnow, asnow0, hlatn, totalb
     O,        tp1,tp2,tp3,tp4,tp5,tp6
     o,        es,ev,et,qin,srflw,rzflw,bflw,adj,cor,exc1,exc2)
c****
c****                   driver for lsm_pilps.f
c****

      integer nch
      logical*1 bug,bug1,bug2
      parameter(ntyps=10,mxchp=3500)
      real sfrac,alon(nch),alat(nch)
      dimension tc1(nch),TC2(nch),TC4(nch),catdef(nch),
     &    rzex(nch),srfex(nch),capac(nch),
     &    qa1(nch),qa2(nch),qa4(nch),hlatn(nch),totalb(nch)
      dimension tmp2m(nch),dpt2m(nch),sfcprs(nch),cpcp(nch),
     &   tpcp(nch),alwdn(nch),swdn(nch),wnd(nch),sunang(nch)
      dimension bf1(nch),bf2(nch),bf3(nch),
     &      vgwmax(nch), cdcr1(nch), cdcr2(nch),
     &      ARS1(NCH), ARS2(NCH), ARS3(NCH), ARA1(NCH), ARA2(NCH), 
     &      ARA3(NCH), ARA4(NCH), ARW1(NCH), ARW2(NCH), ARW3(NCH),
     &      ARW4(NCH), tsa1(NCH), tsa2(NCH), tsb1(NCH), tsb2(NCH),
     &      psis(nch), bee(nch), poros(nch), wpwet(nch), 
     &      cond(nch), gnu(nch)
      dimension albed(nch),green(nch),zlt(nch),z0(nch), d(nch)
      integer ityp(nch)

      dimension          evap(nch), shflux(nch),  runoff(nch), 
     &      eint(nch),   esoi(nch),   eveg(nch), 
     &      esno(nch),  hlwup(nch),  smelt(nch),   bflow(nch),
     &    runsrf(nch), qinfil(nch),    ar1(nch),    rzeq(nch),
     &    ghflux(nch), tmsnow(nch), asnow0(nch),     ar2(nch),
     &  snowterm(nch)

      dimension ea1(mxchp),ea2(mxchp),ea4(mxchp),eas(mxchp),
     &       sqscat(mxchp),z2(mxchp),dzm(mxchp),rsoil1(mxchp),
     &  rsoil2(mxchp),satcap(mxchp),rdc(mxchp),u2fac(mxchp)
      dimension em(mxchp),tm(mxchp),um(mxchp),trainl(mxchp),
     &  trainc(mxchp),tsnow(mxchp),hlwdwn(mxchp),psur(mxchp)

      REAL RA1(MXCHP),RA2(MXCHP),RA4(MXCHP),RAS(MXCHP)
      REAL ETURB1(mxchp),DEDEA1(mxchp),HSTURB1(mxchp),DHSDTC1(mxchp),
     &     DHSDEA1(mxchp),DEDTC1(mxchp),QSAT1(mxchp),DQS1(mxchp)
      REAL ETURB2(mxchp),DEDEA2(mxchp),HSTURB2(mxchp),DHSDTC2(mxchp),
     &     DHSDEA2(mxchp),DEDTC2(mxchp),QSAT2(mxchp),DQS2(mxchp)
      REAL ETURB4(mxchp),DEDEA4(mxchp),HSTURB4(mxchp),DHSDTC4(mxchp),
     &     DHSDEA4(mxchp),DEDTC4(mxchp),QSAT4(mxchp),DQS4(mxchp)
      REAL ETURBS(mxchp),DEDEAS(mxchp),HSTURBS(mxchp),DHSDTCS(mxchp),
     &     DHSDEAS(mxchp),DEDTCS(mxchp),QSATS(mxchp),DQSS(mxchp)

      dimension avisdr(mxchp),anirdr(mxchp),avisdf(mxchp),anirdf(mxchp)
      dimension swnets(mxchp),swnetf(mxchp),par(mxchp),pdir(mxchp),
     &  ghtcnt(6,mxchp),wesn(3,mxchp),htsnn(3,mxchp),sndz(3,mxchp)
      dimension ght1(nch),ght2(nch),ght3(nch),ght4(nch),ght5(nch),
     &    ght6(nch),tp1(nch),tp2(nch),tp3(nch),tp4(nch),tp5(nch),
     &    tp6(nch),TSURF(NCH),wesn1(nch),wesn2(nch),wesn3(nch),
     &    htsnn1(nch),htsnn2(nch),htsnn3(nch),sndz1(nch),sndz2(nch),
     &    sndz3(nch)

      dimension snwmid(ntyps)
      parameter (epsilon = 18.01/28.97)
      parameter (alhe = 2.4548e6, alhs = 2.8368e6, alhm = alhs-alhe)
      data snwmid /50.,50.,50.,2.,50.,2.,2.,2.,2.,2./
c
cc jpw 8/4/99 ----------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
      integer ch
      parameter (ch=4810)
      real a11(ch),a12(ch),a21(ch),a22(ch),a23(ch),a32(ch)
      real a33(ch),q11(ch),q22(ch),q33(ch)
      logical*1 assim
      data assim /.false./
c
      real es(nch),ev(nch),et(nch),qin(nch),srflw(nch),rzflw(nch)
      real bflw(nch),adj(nch),cor(nch),exc1(nch),exc2(nch)
c
cc ---------------------------------------------------------------------
c
c**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c****       pmonth   updates the seasonally varying parameters
c****       turb     computes the turbulent fluxes and tendencies
c****       chip     is the modular si2b.


c     print *, 'process'
c     print *, nch,dtstep,sfrac,jday,sunang(1),alon(1),alat(1),ityp(1)
c     print *, tmp2m(1),dpt2m(1),sfcprs(1),cpcp(1),tpcp(1),alwdn(1)
c     print *, swdn(1),wnd(1),bf1(1),bf2(1),bf3(1),vgwmax(1),cdcr1(1)
c     print *, cdcr2(1),psis(1),bee(1),poros(1),wpwet(1),cond(1),gnu(1)
c     print *, ars1(1),ars2(1),ars3(1),ara1(1),ara2(1),ara3(1),ara4(1)
c     print *, arw1(1),arw2(1),arw3(1),arw4(1),tsa1(1),tsa2(1),tsb1(1)
c     print *, tsb2(1),albed(1),green(1),zlt(1),z0(1),d(1)
c     print *, bug,bug1,bug2
c     print *, tc1(1),tc2(1),tc4(1),qa1(1),qa2(1),qa4(1),catdef(1)
c     print *, rzex(1),srfex(1),capac(1)
c     print *, ght1(1),ght2(1),ght3(1),ght4(1),ght5(1),ght6(1),TSURF(1)
c     print *, wesn1(1),wesn2(1),wesn3(1),htsnn1(1),htsnn2(1),htsnn3(1)
c     print *, sndz1(1),sndz2(1),sndz3(1)
c     print *, evap(1),shflux(1),runoff(1),eint(1),esoi(1),eveg(1)
c     print *, esno(1),snowterm(1)
c     print *, hlwup(1), smelt(1), bflow(1), runsrf(1), ar1(1), ar2(1)
c     print *, rzeq(1), qinfil(1)
c     print *, ghflux(1), tmsnow(1), asnow0(1), hlatn(1), totalb(1)
c     print *, tp1(1),tp2(1),tp3(1),tp4(1),tp5(1),tp6(1)
c     print *, es(1),ev(1),et(1),qin(1),srflw(1),rzflw(1),bflw(1),adj(1)
c     print *, cor(1),exc1(1),exc2(1)


      bug1=.false.
      bug2=.false.
      if (bug) then
         bug1=.true.
         bug=.false.
      endif

      if (bug) then 
         WRITE(*,*) 'start of process OK'
      endif

c****
      call pmonth(
     i  nch, ityp, jday, alat, green, zlt, z0, d,.false.,
     o  sqscat, z2, dzm,
     o  rsoil1, rsoil2, satcap, rdc, u2fac
     &           )

      if (bug) then 
         WRITE(*,*) 'call pmonth OK'
      endif

c****
      do n=1,nch
         tm(n)=tmp2m(n)
         um(n)=wnd(n)
         tsnow(n)=0.
         trainc(n)=amin1(cpcp(n),tpcp(n))
         trainl(n)=tpcp(n)-trainc(n)

c air temperature criterion for snow different in GSWP and PILPS2c
         if(tmp2m(n) .lt. 275.16) then
            tsnow(n)=tpcp(n)
            trainl(n)=0.
            trainc(n)=0.
         endif

         snowterm(n)=((tm(n)-273.16)*(1.9e6/920.)-alhm)*
     &              tsnow(n)/dtstep

         hlwdwn(n)=alwdn(n)
         sunang(n)=amax1(sunang(n),0.01)

c pressure unit different in GSWP and PILPS2c
         psur(n)=sfcprs(n)/100.
         em(n)=esat(dpt2m(n))
         ea1(n)=qa1(n)*psur(n)/epsilon
         ea2(n)=qa2(n)*psur(n)/epsilon
         ea4(n)=qa4(n)*psur(n)/epsilon
         eas(n)=qsat(TMSNOW(N),psur(n),alhs)*psur(n)/epsilon

         ghtcnt(1,n)=ght1(n)
         ghtcnt(2,n)=ght2(n)
         ghtcnt(3,n)=ght3(n)
         ghtcnt(4,n)=ght4(n)
         ghtcnt(5,n)=ght5(n)
         ghtcnt(6,n)=ght6(n)

         wesn(1,n)=wesn1(n)
         wesn(2,n)=wesn2(n)
         wesn(3,n)=wesn3(n)

         htsnn(1,n)=htsnn1(n)
         htsnn(2,n)=htsnn2(n)
         htsnn(3,n)=htsnn3(n)

         sndz(1,n)=sndz1(n)
         sndz(2,n)=sndz2(n)
         sndz(3,n)=sndz3(n)

      enddo

      if (bug) then 
         WRITE(*,*) 'loop 10 OK'
      endif

c**** account for fact that forcing can lead the land:
      do n=1,nch
         
         esattc = qsat(TC1(N),psur(n),alhe)  * psur(n) / epsilon
         if (EA1(n).gt.esattc .and. ea1(n).gt.em(n) ) 
     &        ea1(n) = amax1( em(n), esattc )
         if (ea1(n).lt.esattc .and. ea1(n).lt.em(n) )
     &        ea1(n) = amin1( em(n), esattc )
         
         esattc = qsat(TC2(N),psur(n),alhe)  * psur(n) / epsilon
         if (ea2(n).gt.esattc .and. ea2(n).gt.em(n) ) 
     &        ea2(n) = amax1( em(n), esattc )
         if (ea2(n).lt.esattc .and. ea2(n).lt.em(n) )
     &        ea2(n) = amin1( em(n), esattc )
         
         esattc = qsat(TC4(N),psur(n),alhe)  * psur(n) / epsilon
         if (ea4(n).gt.esattc .and. ea4(n).gt.em(n) ) 
     &        ea4(n) = amax1( em(n), esattc )
         if (ea4(n).lt.esattc .and. ea4(n).lt.em(n) )
     &        ea4(n) = amin1( em(n), esattc )

      enddo

      if (bug) then 
         WRITE(*,*) 'loop 15 OK'
      endif

c****
c****
      call turb(
     i  nch, em,tm,um,TC1,ea1,psur,z0,bug,
     o  ra1,eturb1,dedea1,dedtc1,hsturb1,dhsdea1,dhsdtc1
     &         )

      call turb(
     i  nch, em,tm,um,TC2,ea2,psur,z0,bug,
     o  ra2,eturb2,dedea2,dedtc2,hsturb2,dhsdea2,dhsdtc2
     &         )

      call turb(
     i  nch, em,tm,um,TC4,ea4,psur,z0,bug,
     o  ra4,eturb4,dedea4,dedtc4,hsturb4,dhsdea4,dhsdtc4
     &         )

      call turb(
     i  nch, em,tm,um,TMSNOW,eaS,psur,z0,bug,
     o  raS,eturbS,dedeaS,dedtcS,hsturbS,dhsdeaS,dhsdtcS
     &         )

      if (bug) then 
         WRITE(*,*) 'call turb OK'
      endif                      

c****
c****               (process radiation quantities:)
c****
      do n=1,nch
         avisdr(n)=albed(n)
         anirdr(n)=albed(n)
         avisdf(n)=albed(n)
         anirdf(n)=albed(n)

c         snowx=(wesn(1,n)+wesn(2,n)+wesn(3,n))*1000.
c         if(snowx .gt. 0.) then
c            snwfrc=snowx/(snowx+snwmid(ityp(n)))
c            avisdr(n)=avisdr(n)*(1.-snwfrc) + 0.85*snwfrc
c            anirdr(n)=anirdr(n)*(1.-snwfrc) + 0.50*snwfrc  ! GSWP =/= PILPS2c
cc            anirdr(n)=anirdr(n)*(1.-snwfrc) + 0.65*snwfrc ! PILPS2c
c            avisdf(n)=avisdf(n)*(1.-snwfrc) + 0.85*snwfrc
c            anirdf(n)=anirdf(n)*(1.-snwfrc) + 0.50*snwfrc  ! GSWP =/= PILPS2c
cc            anirdf(n)=anirdf(n)*(1.-snwfrc) + 0.65*snwfrc ! PILPS2c
c         endif

         albave=0.25*(avisdr(n)+anirdr(n)+avisdf(n)+anirdf(n))
         swnetf(n)=(1.-albave)*swdn(n)
         par(n)=0.5*swdn(n)
         pdir(n)=0.5

         avisdrs=0.85
         avisdfs=0.85
         anirdfs=0.50    ! GSWP =/= PILPS2c
         anirdrs=0.50    ! GSWP =/= PILPS2c
c         anirdrs=0.65   ! PILPS2c
c         anirdfs=0.65   ! PILPS2c
         albsn=0.25*(avisdrs+anirdrs+avisdfs+anirdfs)
         swnets(n)=(1.-albsn)*swdn(n)
         totalb(n)=1.-((swnetf(n)*(1.-asnow0(n))
     &        +swnets(n)*asnow0(n))/(swdn(n)+1.e-20))

      enddo

c****
c****             (update prognostic variables with chip:)
c****

      if (bug) then
         WRITE(*,*) 'Before chip OK'
      endif
c
cc jpw 8/4/99 ----------------------------------------------------------
c  added variables assim, a11, a12, a21, a22, a23, a32, a33, q11, q22
c  and q33
c ----------------------------------------------------------------------
c
      CALL CHIP(       nch, dtstep, sfrac, ityp, trainl, trainc, tsnow, 
     i                 um, 
     I                 ETURB1, DEDEA1, DEDTC1, HSTURB1,DHSDEA1, DHSDTC1,
     I                 ETURB2, DEDEA2, DEDTC2, HSTURB2,DHSDEA2, DHSDTC2,
     I                 ETURB4, DEDEA4, DEDTC4, HSTURB4,DHSDEA4, DHSDTC4,
     I                 ETURBS, DEDEAS, DEDTCS, HSTURBS,DHSDEAS, DHSDTCS,
     i                 tm, em, ra1, ra2, ra4, ras, sunang, par, pdir,
     i                 swnetf,swnets,hlwdwn, psur, zlt, green, z2,
     i                 sqscat, rsoil1, rsoil2, satcap, rdc, u2fac,
     i                 bf1, bf2, bf3, vgwmax, cdcr1, cdcr2,
     i                 psis, bee, poros, wpwet, cond,gnu,ars1,ars2,ars3,
     i                 ara1,ara2,ara3,ara4,arw1,arw2,arw3,arw4,
     I                 tsa1,tsa2,tsb1,tsb2,bug1,assim,
     u                 tc1, TC2, TC4, ea1, EA2, EA4, 
     u                 capac,catdef,rzex,srfex,ghtcnt,tsurf,
     U                 WESN, HTSNN, SNDZ,
     o                 evap, shflux, runoff, runsrf, smelt,
     o                 hlwup, hlatn, qinfil, bflow, ar1, ar2, rzeq,
     o                 eint,esoi,eveg, esno, ghflux, 
     O                 tmsnow, asnow0, tp1,tp2,tp3,tp4,tp5,tp6,
     o                 a11,a12,a21,a22,a23,a32,a33,q11,q22,q33,
     o                 es,ev,et,qin,srflw,rzflw,bflw,adj,cor,exc1,exc2
     &         )
c
cc ---------------------------------------------------------------------
c
      if (bug) then
         WRITE(*,*) 'After chip OK'
      endif
      
c****
c****                    (accumulate diagnostics:)
c****
c**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c****
      do n=1,nch
         qa1(n)=ea1(n)*epsilon/psur(n)
         qa2(n)=ea2(n)*epsilon/psur(n)
         qa4(n)=ea4(n)*epsilon/psur(n)
c        cond(n)=1./ra(n)
      enddo

      do n=1,nch
        ght1(n)=ghtcnt(1,n)
        ght2(n)=ghtcnt(2,n)
        ght3(n)=ghtcnt(3,n)
        ght4(n)=ghtcnt(4,n)
        ght5(n)=ghtcnt(5,n)
        ght6(n)=ghtcnt(6,n)

        wesn1(n)=wesn(1,n)
        wesn2(n)=wesn(2,n)
        wesn3(n)=wesn(3,n)

        htsnn1(n)=htsnn(1,n)
        htsnn2(n)=htsnn(2,n)
        htsnn3(n)=htsnn(3,n)

        sndz1(n)=sndz(1,n)
        sndz2(n)=sndz(2,n)
        sndz3(n)=sndz(3,n)
      enddo

      if (bug) then
         WRITE(*,*) 'After loop 80'
      endif

      return
      end
c****
c**** -----------------------------------------------------------------
c**** /////////////////////////////////////////////////////////////////
c**** -----------------------------------------------------------------
c****
      subroutine turb(
     i      nch,em,tm,um,tc,ea,psur,z0,bug,
     o      ra,eturb,dedea,dedtc,hsturb,dhsdea,dhsdtc)
c****
      logical*1 bug
      parameter(mxchp=5022)
      dimension em(nch),tm(nch),um(nch),tc(nch),
     &      ea(nch),psur(nch),z0(nch),
     &      ra(nch),eturb(nch),dedea(nch),dedtc(nch),
     &      hsturb(nch),dhsdea(nch),dhsdtc(nch)

      data grav/9.81/,rgas/287./,cp/1010./
      data vkrmn/0.41/,epsi/0.611/,b/5./,c/5./,d/5./,
     &    zthick/50./
c****
      do 100 n=1,nch

         
         if (bug .and. n .eq. 134) then
            write(*,*) 'turb: n=',n
         endif

         rhoair=psur(n)*100./(rgas*tc(n))
         qm=em(n)*epsi/psur(n)
         qc=ea(n)*epsi/psur(n)
         dtv=tm(n)*(1.+epsi*qm)-tc(n)*(1.+epsi*qc)
         ri=grav*zthick*dtv/(tm(n)*um(n)*um(n))
         dridtc=-grav*zthick*(1.+epsi*qc)/(tm(n)*um(n)*um(n))
         dridqc=-grav*zthick*epsi*tc(n)/(tm(n)*um(n)*um(n))
         dummy=alog(zthick/z0(n)+1.)
         cn=vkrmn*vkrmn/(dummy*dummy)
c****
         if(ri.ge.0.) then
            t=sqrt(1.+d*ri)
            fu=1./(1.+2.*b*ri/t)
            ft=1./(1.+3.*b*ri*t)
            dftdri=-3.*b*ft*ft*( (d*ri)/(2.*t) + t )
         endif
         if(ri.lt.0.) then
            cs=cn*sqrt((zthick/z0(n))+1.)
            r=3.*b*cs*sqrt(-ri)
            s=1./(1.+c*r)
            t=b*ri*s
            ft=1.-3.*t
            dftdri=-1.5*b*s*s*(2.+c*r)
         endif
c****
         dftdqc=dftdri*dridqc
         dftdtc=dftdri*dridtc
         ct=cn*ft
         hsturb(n)=rhoair*cp*ct*um(n)*(tc(n)-tm(n))
         eturb(n)=rhoair*ct*um(n)*(qc-qm)
         dedqa=rhoair*cn*um(n)*( dftdqc*(qc-qm) + ft )
         dedtc(n)=rhoair*cn*um(n)*dftdtc*(qc-qm)
         dhsdqa=rhoair*cp*cn*um(n)*dftdqc*(tc(n)-tm(n))
         dhsdtc(n)=rhoair*cp*cn*um(n)*( dftdtc*(tc(n)-tm(n)) + ft )
c**** 
         dedea(n)=dedqa*epsi/psur(n)
         dhsdea(n)=dhsdqa*epsi/psur(n)
c****
         ra(n)=1/(ct*um(n))
c****
         if(dhsdtc(n).lt.0.) dhsdtc(n)=0.
         if(dhsdea(n).lt.0.) dhsdea(n)=0.
         if(dedtc(n).lt.0.) dedtc(n)=0.
         if(dedea(n).lt.0.) dedea(n)=0.
c****
 100  continue

      return
      end
c****
c**** -----------------------------------------------------------------
c**** /////////////////////////////////////////////////////////////////
c**** -----------------------------------------------------------------
*---<|>--1---------2---------3---------4---------5---------6---------7-<
c****
c**** [ begin chip ]
c**** 
      subroutine chip (
     i                 nch, dtstep, sfrac, ityp, trainl, trainc, tsnow, 
     i                 um, 
     I                 ETURB1, DEDEA1, DEDTC1, HSTURB1,DHSDEA1, DHSDTC1,
     I                 ETURB2, DEDEA2, DEDTC2, HSTURB2,DHSDEA2, DHSDTC2,
     I                 ETURB4, DEDEA4, DEDTC4, HSTURB4,DHSDEA4, DHSDTC4,
     I                 ETURBS, DEDEAS, DEDTCS, HSTURBS,DHSDEAS, DHSDTCS,
     i                 tm, em, ra1, ra2, ra4, ras, sunang, par, pdir,
     i                 swnetf,swnets,  hlwdwn, psur, zlai, green, z2,
     i                 sqscat, rsoil1, rsoil2, satcap, rdc, u2fac,
     i                 bf1, bf2, bf3, vgwmax, cdcr1, cdcr2,
     i                 psis, bee, poros, wpwet, cond,gnu,ars1,ars2,ars3,
     i                 ara1,ara2,ara3,ara4,arw1,arw2,arw3,arw4,
     I                 tsa1,tsa2,tsb1,tsb2,bug,assim,
     u                 tc1, TC2, TC4, ea1, EA2, EA4, 
     u                 capac,catdef,rzex,srfex,ghtcnt,tsurf,
     U                 WESN, HTSNN, SNDZ,
     o                 evap, shflux, runoff, runsrf, smelt,
     o                 hlwup, hlatn, qinfil, bflow, ar1, ar2, rzeq,
     o                 eint,esoi,eveg, esno, ghflux, 
     O                 tmsnow, asnow0, tp1,tp2,tp3,tp4,tp5,tp6,
     o                 a11,a12,a21,a22,a23,a32,a33,q11,q22,q33,
     o                 es,ev,et,qin,srflw,rzflw,bflw,adj,cor,exc1,exc2
     &                )
c****
c****       this subroutine calls chipx, which is the chip subroutine
c**** to be used with the gcm.  this subroutine acts as an intermediary
c**** between chipx and the chip calls made by our old offline forcing
c**** programs.
c**** 
c      implicit none
      integer  ntyps, frstch, memfac, mxchp
      integer  nlay, sfcly, rootly, rechly

      real  zero, one, pie
      real alhe, alhs, alhm, tf, stefan, rgas, shw, shi, rhow, grav
      real epsilon, nosnow

      parameter (ntyps = 10, frstch = 1, memfac = 5, mxchp=3500)
      parameter (nlay = 3)
      parameter (sfcly = 1, rootly = sfcly + 1, rechly = rootly + 1)

      parameter (zero = 0., one = 1., pie = 3.14159265)
      parameter (alhe = 2.4548e6, alhs = 2.8368e6, alhm = alhs-alhe)
      parameter (tf = 273.16)
      parameter (stefan = 5.669e-8)
      parameter (rgas = .286*1003.5)
      parameter (shw = 4200., shi = 2060.)
      parameter (rhow = 1000.)
      parameter (grav = 9.81)
      parameter (epsilon = 18.01/28.97)
      parameter (nosnow = 0.)

      logical*1 bug,bug1

      integer nch, n
      integer ityp(nch)

      real       dtstep,       sfrac,  tsnow(nch),     um(nch),
     &       eturb(nch),  dedea(nch), hsturb(nch), dhsdtc(nch),
     &     tm(nch),ra1(nch),ra2(nch),ra4(nch),ras(nch),
     &     sunang(nch),dhsdea(nch),
     &       em (nch),  par(nch),  pdir(nch),  swnetf(nch), swnets(nch),
     &      hlwdwn(nch),   psur(nch),   zlai(nch),  green(nch),
     &          z2(nch), sqscat(nch),  dedtc(nch)
      real  rsoil1(nch), rsoil2(nch), satcap(nch),    rdc(nch),
     &       u2fac(nch),     tc1(nch), TC2(NCH), TC4(NCH),  
     &         ea1(nch),     EA2(NCH),     EA4(NCH),
     &       catdef(nch),  rzex(nch), rzeq(nch), capac(nch), 
     &         evap(nch), shflux(nch), runoff(nch),
     &       ghtcnt(6,nch), srfex(nch), TSURF(NCH), wesn(3,nch),
     &       htsnn(3,nch), sndz(3,nch)

      real    bflow(nch), runsrf(nch), smelt(nch), 
     &        hlwup(nch), hlatn(nch),
     &          eint(nch), qinfil(nch), 
     &        trainc(nch),trainl(nch),
     &        esoi(nch),  eveg(nch),  esno(nch), ar1(nch), ar2(nch),
     &        ghflux(nch), tmsnow(nch), asnow0(nch)
      real bf1(nch), bf2(nch), bf3(nch), 
     &      vgwmax(nch), cdcr1(nch), cdcr2(nch),
     &      ars1(nch),ars2(nch),ars3(nch),ara1(nch),ara2(nch),ara3(nch),
     &      ara4(nch),arw1(nch),arw2(nch),arw3(nch),arw4(nch),
     &      psis(nch), bee(nch), poros(nch), wpwet(nch),
     &      cond(nch), gnu(nch),
     &      tsa1(NCH), tsa2(NCH),tsb1(NCH), tsb2(NCH),
     &      tp1(nch),tp2(nch),tp3(nch),tp4(nch), tp5(nch), tp6(nch)

      REAL ETURB1(NCH),DEDEA1(NCH),HSTURB1(NCH),DHSDTC1(NCH),
     &     DHSDEA1(NCH),DEDTC1(NCH)
      REAL ETURB2(NCH),DEDEA2(NCH),HSTURB2(NCH),DHSDTC2(NCH),
     &     DHSDEA2(NCH),DEDTC2(NCH)
      REAL ETURB4(NCH),DEDEA4(NCH),HSTURB4(NCH),DHSDTC4(NCH),
     &     DHSDEA4(NCH),DEDTC4(NCH)
      REAL ETURBS(NCH),DEDEAS(NCH),HSTURBS(NCH),DHSDTCS(NCH),
     &     DHSDEAS(NCH),DEDTCS(NCH)

      REAL DEDQA1(MXCHP),DHSDQA1(MXCHP)
      REAL DEDQA2(MXCHP),DHSDQA2(MXCHP)
      REAL DEDQA4(MXCHP),DHSDQA4(MXCHP)
      REAL DEDQAS(MXCHP),DHSDQAS(MXCHP)

      REAL QSAT1(MXCHP),DQS1(MXCHP),ALW1(MXCHP),BLW1(MXCHP),
     &     QSAT2(MXCHP),DQS2(MXCHP),ALW2(MXCHP),BLW2(MXCHP),
     &     QSAT4(MXCHP),DQS4(MXCHP),ALW4(MXCHP),BLW4(MXCHP),
     &     QSATS(MXCHP),DQSS(MXCHP),ALWS(MXCHP),BLWS(MXCHP) 

      real dummy, esat, lwattc, dlwdtc, snwfrc

      real  pardir(mxchp), pardif(mxchp), qsattc(mxchp), dqsdtc(mxchp),
     &          qm(mxchp),    qa1(mxchp),    qa2(mxchp),    qa4(mxchp),  
     &       dedqa(mxchp), dhsdqa(mxchp),
     &     cd1(mxchp), cd2(mxchp), cd4(mxchp), cds(mxchp), alhx(mxchp),
     &      strdg1(mxchp),strdg2(mxchp),strdg3(mxchp),strdg4(mxchp)
      real snwmid(ntyps)
      data snwmid /50.,50.,50.,2.,50.,2.,2.,2.,2.,2./
c
cc jpw 8/4/99 ----------------------------------------------------------
c  for soil moisture assimilation
c ----------------------------------------------------------------------
c
      real a11(nch),a12(nch),a21(nch),a22(nch),a23(nch),a32(nch)
      real a33(nch),q11(nch),q22(nch),q33(nch)
      logical*1 assim
c
      real es(nch),ev(nch),et(nch),qin(nch),srflw(nch),rzflw(nch)
      real bflw(nch),adj(nch),cor(nch),exc1(nch),exc2(nch)
c
cc ---------------------------------------------------------------------
c

c****
c**** -----------------------------------------------------------------
c****
      bug1=.false.
      if (bug) then
         bug=.false.
         bug1=.true.
      endif

      do 10 n = 1, nch
         qm(n) = em(n) * epsilon / psur(n)
         qa1(n) = ea1(n) * epsilon / psur(n)
         qa2(n) = ea2(n) * epsilon / psur(n)
         qa4(n) = ea4(n) * epsilon / psur(n)

         qsat1(n) = qsat(tc1(n),psur(n),alhe)
         dqs1(n) = qsat1(n) * 5418. / ( tc1(n) * tc1(n) )
         dedqa1(n) = dedea1(n) * psur(n) / epsilon
         dhsdqa1(n) = dhsdea1(n) * psur(n) / epsilon

         qsat2(n) = qsat(tc2(n),psur(n),alhe)
         dqs2(n) = qsat2(n) * 5418. / ( tc2(n) * tc2(n) )
         dedqa2(n) = dedea2(n) * psur(n) / epsilon
         dhsdqa2(n) = dhsdea2(n) * psur(n) / epsilon

         qsat4(n) = qsat(tc4(n),psur(n),alhe)
         dqs4(n) = qsat4(n) * 5418. / ( tc4(n) * tc4(n) )
         dedqa4(n) = dedea4(n) * psur(n) / epsilon
         dhsdqa4(n) = dhsdea4(n) * psur(n) / epsilon

         qsats(n) = qsat(tmsnow(n),psur(n),alhs)
         dqss(n) = qsats(n) * 5418. / ( tmsnow(n) * tmsnow(n) )
         dedqas(n) = dedeas(n) * psur(n) / epsilon
         dhsdqas(n) = dhsdeas(n) * psur(n) / epsilon 

         pardir(n) = par(n) * pdir(n)
         pardif(n) = par(n) * ( 1. - pdir(n) )

         cd1(n) = 1. / ( um(n) * ra1(n) )
         cd2(n) = 1. / ( um(n) * ra2(n) )
         cd4(n) = 1. / ( um(n) * ra4(n) )
         cds(n) = 1. / ( um(n) * ras(n) )

         tsnow(n) = tsnow(n) / dtstep
         trainc(n) = trainc(n) / dtstep
         trainl(n) = trainl(n) / dtstep
   10 continue
c****
c**** compute constants for longwave radiation linearization
      do 12 n = 1, nch

         lwattc = 5.67e-8 * tc1(n) * tc1(n) * tc1(n) * tc1(n)
         dlwdtc =  4. * lwattc / tc1(n)
         alw1(n) = lwattc - dlwdtc * tc1(n)
         blw1(n) = dlwdtc

         lwattc = 5.67e-8 * tc2(n) * tc2(n) * tc2(n) * tc2(n)
         dlwdtc =  4. * lwattc / tc2(n)
         alw2(n) = lwattc - dlwdtc * tc2(n)
         blw2(n) = dlwdtc

         lwattc = 5.67e-8 * tc4(n) * tc4(n) * tc4(n) * tc4(n)
         dlwdtc =  4. * lwattc / tc4(n)
         alw4(n) = lwattc - dlwdtc * tc4(n)
         blw4(n) = dlwdtc

         lwattc = 5.67e-8 *tmsnow(n)*tmsnow(n)*tmsnow(n)*tmsnow(n)
         dlwdtc =  4. * lwattc / tmsnow(n)
         alws(n) = lwattc - dlwdtc * tmsnow(n)
         blws(n) = dlwdtc

   12 continue
c****

      if (bug) then 
         WRITE(*,*) 'In chip, before catchment OK'
      endif
c
cc jpw 8/4/99 ----------------------------------------------------------
c  added variables assim, a11, a12, a21, a22, a23, a32, a33, q11, q22
c  and q33
c ----------------------------------------------------------------------
c
      call catchment (
     I               NCH, DTSTEP,SFRAC,ITYP, TRAINC, TRAINL, TSNOW,  UM,
     I               ETURB1, DEDQA1, DEDTC1, HSTURB1,DHSDQA1, DHSDTC1,
     I               ETURB2, DEDQA2, DEDTC2, HSTURB2,DHSDQA2, DHSDTC2,
     I               ETURB4, DEDQA4, DEDTC4, HSTURB4,DHSDQA4, DHSDTC4,
     I               ETURBS, DEDQAS, DEDTCS, HSTURBS,DHSDQAS, DHSDTCS,
     I               TM, QM, CD1, CD2, CD4, CDS, SUNANG, PARDIR, PARDIF,
     I               SWNETF, SWNETS, HLWDWN, PSUR,  ZLAI,   GREEN,  Z2,
     I               SQSCAT, RSOIL1, RSOIL2,   RDC,    U2FAC,
     I               QSAT1, DQS1, ALW1, BLW1,
     I               QSAT2, DQS2, ALW2, BLW2,
     I               QSAT4, DQS4, ALW4, BLW4,
     I               QSATS, DQSS, ALWS, BLWS,
     I               BF1, BF2, BF3, VGWMAX, 
     I               cdcr1, cdcr2, psis, bee, poros, wpwet, cond, gnu,
     I               ars1,ars2,ars3,ara1,ara2,ara3,ara4,
     I               arw1,arw2,arw3,arw4,
     I               tsa1,tsa2,tsb1,tsb2,BUG1,assim,
     U               TC1, TC2, TC4, QA1, QA2, QA4, CAPAC, 
     U               CATDEF, RZEX, srfex, ghtcnt, TSURF,
     U               WESN, HTSNN, SNDZ,
     O               EVAP, SHFLUX, RUNOFF,
     O               EINT,   ESOI,   EVEG,   ESNO,
     O               BFLOW,RUNSRF,SMELT, 
     O               HLWUP,HLATN,QINFIL,AR1,  AR2, RZEQ, 
     O               GHFLUX, TMSNOW, ASNOW0, tp1,tp2,tp3,
     O               tp4,tp5,tp6,
     o               a11,a12,a21,a22,a23,a32,a33,q11,q22,q33,
     o               es,ev,et,qin,srflw,rzflw,bflw,adj,cor,exc1,exc2
     &           )

      if (bug) then
         WRITE(*,*) 'after catchment OK'
      endif

c****
      do 20 n = 1, nch
         ea1(n) = qa1(n) * psur(n) / epsilon
         ea2(n) = qa2(n) * psur(n) / epsilon
         ea4(n) = qa4(n) * psur(n) / epsilon
   20 continue
c****
      return
      end
c****
c**** ------------------------------------------------------------------
c**** //////////////////////////////////////////////////////////////////
c**** ------------------------------------------------------------------
c****
      real function esat(t)
c jpw real function esat(t,alhx)
c**** 
      implicit none
      integer  ntyps, frstch, memfac
      integer  nlay, sfcly, rootly, rechly

      real  zero, one, pie
      real alhe, alhs, alhm, tf, stefan, rgas, shw, shi, rhow, grav
      real epsilon, nosnow

      parameter (ntyps = 10, frstch = 1, memfac = 5)
      parameter (nlay = 3)
      parameter (sfcly = 1, rootly = sfcly + 1, rechly = rootly + 1)

      parameter (zero = 0., one = 1., pie = 3.14159265)
      parameter (alhe = 2.4548e6, alhs = 2.8368e6, alhm = alhs-alhe)
      parameter (tf = 273.16)
      parameter (stefan = 5.669e-8)
      parameter (rgas = .286*1003.5)
      parameter (shw = 4200., shi = 2060.)
      parameter (rhow = 1000.)
      parameter (grav = 9.81)
      parameter (epsilon = 18.01/28.97)
      parameter (nosnow = 0.)

      real  t, alhx
c****
      esat = exp(21.18123 - 5418./t) / epsilon
c****
      return
      end
