
      subroutine snowrt(
     I                  ts,t1,area,pr,snowf,dtsi,eturb,dedtc,dedea,
     I                  dhsdtc,dhsdea,esattc,desdtc,dhlwtc,ea,
     I                  swnet,hlwdwn,hlwtc,hsturb,
     I                  xklhw,xklh,zc,bug,
     U                  wesn,htsnn,sndz,
     O                  tpsn,areasc,pre,fh31w,fh31i,fh31d,hlwout,
     O                  evap0,shflux,areasc0
     &                 )
      parameter (ngm=7)
      dimension fices(3),wesn(3),sndz(3),rhos(3),sdnew(3),htsno(3),
     *htsnn(3),sdold(3),cpsn(3),fhsn(3),tpsn(3),wcont(3),tksn(3),
     *fold(3),t1(3),area(3)
      dimension xklh(ngm),zc(ngm)
      logical*1 bug
c
      real tstart,tfinal,tenergy,tsnow,tfluxes,tground
      common/timing/tcatch,tenergy,tsnow,tfluxes,tground

c      call cpu_time(tstart)

      dts=dtsi
      tpsn(1)=0.0
      tpsn(2)=0.0
      tpsn(3)=0.0
      wcont(1)=0.0
      wcont(2)=0.0
      wcont(3)=0.0
      areasc=0.0
      zfsh=0.0
      pre=0.0
      fh31w=0.0
      fh31i=0.0
      fh31d=0.0
      hlwout=0.0
      evap=0.0
      shflux=0.0
      evap0=0.0
      areasc0=0.
      fresh=0.
      if(wesn(1)+wesn(2)+wesn(3) .eq. 0. .and. snowf.eq. 0.)goto 3000

      iflag0=0
      xsmin=0.013
      lhv=2.4548e6
      lhs=2.8368e6
      fsn=1000.*(lhs-lhv)
      alhx=lhs
      tfrz=273.16
      effpr=pr
      snowdo=wesn(1)+wesn(2)+wesn(3)
      if (snowf .gt. 0) ts=min(ts,tfrz)
      snowd=wesn(1)+wesn(2)+wesn(3)
      if (snowd .ge. xsmin)then
         do i=1,3
            rhos(i)=wesn(i)*1000/sndz(i)
            cpsn(i)=rhos(i)*1.9e6/920
            tksn(i)=3.2217e-12*(rhos(i)*1e3)*(rhos(i)*1e3)
         enddo
         deltg= 2.5e6*(0.05*0.05)/xklhw
         if (deltg .lt. 1200.)then
            xklhw=2.5e6*(0.05*0.05)/1200.
            xklh(2)=xklhw
         endif
         dtm=1200
         DO I=1,3
            dtm2=CPSN(I)*(SNDZ(I)*SNDZ(I))/(TKSN(I)+1E-12)
            DTM=MIN(DTM,dtm2)
         enddo
         if (dtm .lt. 1200)then
            xsmin=snowd
            iflag0=1
         endif
      endif
c**** set flag
      if (snowf .gt. 0. .and. snowd .eq. 0.) then
         iflag=1
      else
         iflag=2
      endif
c**** initialize snow parameters
      if (snowf .gt. 0. .and. snowd .eq. 0) then
         wesn(1)=snowf*dts
         snowd=wesn(1)
         areasc=snowd/xsmin
         areasc=min(areasc,1.0)
         rhos(1)=150.
         sndz(1)=(xsmin*1000./rhos(1))
         totsn=sndz(1)
         tpsns=ts-tfrz
         fices(1)=1
         cpsn(1)=rhos(1)*1.9e6/920
c1231         htsnn(1)=-fsn*wesn(1)+tpsns*cpsn(1)*sndz(1)
         htsnn(1)=-fsn*wesn(1)+tpsns*cpsn(1)*(wesn(1)*1000./rhos(1))
         fhsn(1)=0
         wcont(1)=0
         wesn(2)=0.
         wesn(3)=0.
         htsnn(2)=0.
         htsnn(3)=0.
         sndz(2)=0.
         sndz(3)=0.
         rhos(2)=150.
         rhos(3)=150.
         pre=0.
       endif
c
c**** determine the fractional snow coverage
      snowd=wesn(1)+wesn(2)+wesn(3)
      areasc=snowd/xsmin
      areasc=min(areasc,1.0)
      areasc0=areasc
c
      enin=htsnn(1)+htsnn(2)+htsnn(3)
c**** if the water equivalent snow depth is less than xsmin *****
c******************we use a single layer snow model***************
c
      if (snowd .gt. 0. .and. snowd .le. xsmin)then
c**** lets initialize all of our snow parameters*****
c1/6/98         sndz(1)=(sndz(1)+sndz(2)+sndz(3))
         wesn(1)=wesn(1)+wesn(2)+wesn(3)
         rhos(1)=150.
         sndz(1)=(xsmin*1000.)/rhos(1)
         snowd=wesn(1)
         totsn=sndz(1)
         htsnn(1)=htsnn(1)+htsnn(2)+htsnn(3)
         cpsn(1)=rhos(1)*1.9e6/920
         tksn(1)=3.2217e-12*(rhos(1)*1e3)*(rhos(1)*1e3)
c******* from the layer heat and water contents, calculate old layer
c        temperatures
         tpsn(1)=0
         fices(1)=-htsnn(1)/(fsn*wesn(1))
         if(fsn*wesn(1)+htsnn(1).ge.0)go to 6
         tpsn(1)=(htsnn(1)+wesn(1)*fsn)/(cpsn(1)*
     1           (wesn(1)*1000./rhos(1)))
         fices(1)=1
         go to 9
   6     if(htsnn(1).le.0.)go to 9
         tpsn(1)=htsnn(1)/(wesn(1)*cpsn(1))
         fices(1)=0
         if(tpsn(1) .gt. 0)tpsn(1)=0
   9     continue
c
c******* calculate the heat flow in each layer and determine surface
c        energy fluxes
         factw=area(1)
         facti=area(2)
         factd=area(3)
         denom=(sndz(1)/2-zc(1))
         term1=sqrt(xklhw*tksn(1))
         fh31w=-term1*factw*(tpsn(1)-t1(1))/denom
         term2=sqrt(xklh(2)*tksn(1))
         fh31i=-term2*facti*(tpsn(1)-t1(2))/denom
         fh31d=-term2*factd*(tpsn(1)-t1(3))/denom
c
         fh31wo=fh31w
         fh31io=fh31i
         fh31do=fh31d
c
c***Randy
         fh0=fh31w+fh31i+fh31d
         dfh31w=-term1*factw/denom
         dfh31i=-term2*facti/denom
         dfh31d=-term2*factd/denom
         dfh0=dfh31w+dfh31i+dfh31d

         a11=(cpsn(1)*sndz(1))/dts+dhlwtc+dhsdtc+alhx*dedtc-dfh0
         a12=dhsdea+alhx*dedea
         a21=-desdtc
         a22=1.
         f0=esattc-ea
         q0=swnet+hlwdwn-hlwtc-hsturb-alhx*eturb+fh0

         dea=(q0*a21-a11*f0)/(a12*a21-a11*a22)
         dtc=(q0-a12*dea)/a11
         evap=eturb+dedea*dea+dedtc*dtc
         shflux=hsturb+dhsdea*dea+dhsdtc*dtc
         hlwout=hlwtc+dhlwtc*dtc
         fh1=fh0+dfh0*dtc
         fhsn(1)=swnet+hlwdwn-hlwout-evap*alhx-shflux+fh1
         fh31w=fh31wo+dfh31w*dtc
         fh31i=fh31io+dfh31i*dtc
         fh31d=fh31do+dfh31d*dtc

         if(tpsn(1).le.0. .and. tpsn(1)+dtc.gt.0) then
            dtc=0.-tpsn(1)
            dea=(f0-a21*dtc)/a22
            evap=eturb+dedea*dea+dedtc*dtc
            shflux=hsturb+dhsdea*dea+dhsdtc*dtc
            hlwout=hlwtc+dhlwtc*dtc
            fh1=fh0+dfh0*dtc
            fhsn(1)=swnet+hlwdwn-hlwout-evap*alhx-shflux+fh1
            fh31w=fh31wo+dfh31w*dtc
            fh31i=fh31io+dfh31i*dtc
            fh31d=fh31do+dfh31d*dtc
         endif
c
c * test
c         fh1=fh0
c         fhsn(1)=fbsw+fblw-evp-sen+fh1
c         evap=evp/alhx
c * end of test
c1/5/98         evap=(evap/1000.)*dts
         evap0=evap*alhx
         evap=areasc*(evap/1000.)*dts
c**** take care of evporation/condensation from the snow surface
         if(evap .ge. 0.)then
            if (evap .gt. wesn(1))then
               shflux=((evap-wesn(1))*1000.*alhx)/dts
               evap=wesn(1)
               wesn(1)=0.
            else
               wesn(1)=wesn(1)-evap
            endif
         else
            wesn(1)=wesn(1)-evap
         endif
c1231         htsnn(1)=htsnn(1)+evap*fsn
         pre=0
         smelt=0
         tmelt=0
         flowth=0
         fresh=0
c******* take care of meltwater moving through the pack
         wcont(1)=(1-fices(1))*wesn(1)
         smelt=wcont(1)
         wesn(1)=wesn(1)-smelt
         pre=smelt/dts
c030199         htsnn(1)=htsnn(1)+smelt*fsn
c******* take care of water flow through due to precip
         if(effpr-snowf .gt. 0)pre=pre+areasc*(effpr-snowf)
c******* update snowdepth(water equiv.) for snowf
         fresh=snowf*dts
c******* scale the fluxes by the snow covered fraction
c1231         fh31w=fh31w*areasc
c1231         fh31i=fh31i*areasc
c1231         fh31d=fh31d*areasc
c1231         hlwout=hlwout*areasc
c1231         evap=evap*areasc
c1231         shflux=shflux*areasc
c******* update our new boundaries, heat and water contents
         if (iflag .eq. 2) then
            wesn(1)=wesn(1)+fresh
            cpnew=150.*1.9e6/920
            hnewp=(ts-tfrz)*cpnew-fsn*150./1000.
            zfsh=fresh*1000./150.
            htsnn(1)=htsnn(1)+hnewp*zfsh
         endif
         snowd=wesn(1)
         areasc0=snowd/xsmin
         areasc0=min(areasc0,1.0)
         if (wesn(1) .gt. 0.)then
            sndz(1)=(xsmin*1000.)/rhos(1)
            totsn=sndz(1)
            if(iflag0 .eq. 1)areasc=1.
         else
            tpsn(1)=0.
            sndz(1)=0.
            totsn=0.
            goto 18
         endif
c
c******* from surface-atm flux and snow-ground flux determine the
c        snow layer heat content
c1231         htsnn(1)=htsnn(1)+(fhsn(1))*dts
         htsnn(1)=htsnn(1)+((fhsn(1))*dts)*areasc
c******* from updated heat and water content, determine fraction of ice
c        and temp of snow layers
         tpsn(1)=0
         fices(1)=-htsnn(1)/(fsn*wesn(1))
         if(fsn*wesn(1)+htsnn(1).ge.0.)go to 12
         tpsn(1)=(htsnn(1)+wesn(1)*fsn)/(cpsn(1)*
     1           (wesn(1)*1000./rhos(1)))
         fices(1)=1
         go to 16
   12    if(htsnn(1).le.0)go to 16
         fices(1)=0
   16    tempu=tpsn(1)+tfrz
c******* determine the terrestrial temperature
         if (tpsn(1) .eq. 0. .and. fices(1) .eq. 0.) then
            pre=pre+wesn(1)/dts
            wesn(1)=0
            snowd=0.
            totsn=0.
            areasc=0.
         endif
c******* lets initialize all of our snow parameters*****
   18    sdnew(1)=totsn
         sdnew(2)=totsn/2 + totsn/4
         sdnew(3)=totsn/4
         sndz(1)=sdnew(1)-sdnew(2)
         sndz(2)=sdnew(2)-sdnew(3)
         sndz(3)=sdnew(3)
         wesn(1)=.25*snowd
         wesn(2)=.5*snowd
         wesn(3)=.25*snowd
         do i=1,3
c2/22/99            fices(i)=1
            if (wesn(1) .eq. 0.) fices(i)=0.
            cpsn(i)=rhos(1)*1.9e6/920
c2/22/99            htsnn(i)=-fsn*wesn(i)+tpsn(1)*cpsn(i)*
c2/22/99     1                (wesn(i)*1000./rhos(1))
            htsnn(i)=-fices(1)*fsn*wesn(i)+tpsn(1)*cpsn(i)*
     1                (wesn(i)*1000./rhos(1))
            fhsn(i)=0
            wcont(i)=0
            sdold(i)=sdnew(i)
            htsno(i)=htsnn(i)
         enddo
      goto 2000
      endif
c************* end of single layer snow model************************
c
c**** if the water equivalent snow depth is greater than xsmin*******
c************************we use a three layer snow model*************
c
      if (snowd .gt. xsmin) then
c**** initialize snow parameters
         totsn=sndz(1)+sndz(2)+sndz(3)
         if(totsn .gt. 0.2) then
            sdnew(1)=totsn
            sdnew(2)=totsn-0.05
            sdnew(3)=.66*sdnew(2)
         else
            sdnew(1)=totsn
            sdnew(2)=totsn/2 + totsn/4
            sdnew(3)=totsn/4
         endif
         do i=1,3
            rhos(i)=wesn(i)*1000/sndz(i)
            cpsn(i)=rhos(i)*1.9e6/920
            tksn(i)=3.2217e-12*(rhos(i)*1e3)*(rhos(i)*1e3)
c            htsno(i)=htsnn(i)
         enddo
c******* from layer heat and water content determine fraction of ice and

c        temp of snow layers
         do 20 i=1,3
            tpsn(i)=0
            fices(i)=-htsnn(i)/(fsn*wesn(i))
            if(fsn*wesn(i)+htsnn(i).ge.0)go to 19
            tpsn(i)=(htsnn(i)+wesn(i)*fsn)/(cpsn(i)*sndz(i))
            fices(i)=1
            go to 20
  19        if(htsnn(i).le.0)go to 20
            fices(i)=0
  20     continue
         tempu=tpsn(1)+tfrz
c
c******* calculate the heat flow in each layer and determine surface
c        energy fluxes
         sndist=(sdnew(1)-sdnew(2))/2+sdnew(2)
         denom=sndist-sdnew(1)/2
         term1=sqrt(tksn(2)*tksn(1))
         fh0=-term1*(tpsn(1)-tpsn(2))/denom
c
c***Randy
         dfh0=-term1/denom
         a11=(cpsn(1)*sndz(1))/dts+dhlwtc+dhsdtc+alhx*dedtc-dfh0
         a12=dhsdea+alhx*dedea
         a21=-desdtc
         a22=1.
         f0=esattc-ea
         q0=swnet+hlwdwn-hlwtc-hsturb-alhx*eturb+fh0

         dea=(q0*a21-a11*f0)/(a12*a21-a11*a22)
         dtc=(q0-a12*dea)/a11
         evap=eturb+dedea*dea+dedtc*dtc
         shflux=hsturb+dhsdea*dea+dhsdtc*dtc
         hlwout=hlwtc+dhlwtc*dtc
         fhsn(2)=fh0+dfh0*dtc
         fhsn(1)=swnet+hlwdwn-hlwout-evap*alhx-shflux+fhsn(2)

         if(tpsn(1).le.0. .and. tpsn(1)+dtc.gt.0) then
            dtc=0.-tpsn(1)
            dea=(f0-a21*dtc)/a22
            evap=eturb+dedea*dea+dedtc*dtc
            shflux=hsturb+dhsdea*dea+dhsdtc*dtc
            hlwout=hlwtc+dhlwtc*dtc
            fhsn(2)=fh0+dfh0*dtc
            fhsn(1)=swnet+hlwdwn-hlwout-evap*alhx-shflux+fhsn(2)
         endif
c
c******* calculate the ground-snow energy flux
         denom=(sdnew(3)/2-zc(1))
         fh31w=-sqrt(xklhw*tksn(3))*area(1)*(tpsn(3)-t1(1))/denom
         fh31i=-sqrt(xklh(2)*tksn(3))*area(2)*(tpsn(3)-t1(2))/denom
         fh31d=-sqrt(xklh(2)*tksn(3))*area(3)*(tpsn(3)-t1(3))/denom
         fh1=fh31w+fh31i+fh31d
c******* calculate the heat flux between snow layers
         fhsn(3)=-sqrt(tksn(2)*tksn(3))*(tpsn(2)-tpsn(3))/
     1   (sdnew(1)/2-sdnew(3)/2)
c******* from the energy fluxes calculate layer heat content
         htsnn(3)=htsnn(3)+(fh1-fhsn(3))*dts
         htsnn(2)=htsnn(2)+(fhsn(3)-fhsn(2))*dts
         htsnn(1)=htsnn(1)+(fhsn(1))*dts
         do i=1,3
            htsno(i)=htsnn(i)
         enddo
c******* from layer heat and water content determine fraction of ice and

c        temp of snow layers
         do 140 i=1,3
            tpsn(i)=0
            fices(i)=-htsnn(i)/(fsn*wesn(i))
            if(fsn*wesn(i)+htsnn(i).ge.0)go to 130
            tpsn(i)=(htsnn(i)+wesn(i)*fsn)/(cpsn(i)*sndz(i))
            fices(i)=1
            go to 140
  130       if(htsnn(i).le.0)go to 140
            fices(i)=0
  140    continue
         tempu=tpsn(1)+tfrz
c
         evap0=evap*alhx
         evap=(evap/1000.)*dts
c**** take care of evporation/condensation from the snow surface
         if(evap .ge. 0.)then
            wesn(1)=wesn(1)-evap
            sndz(1)=sndz(1)-evap*(1000./rhos(1))
         else
            wesn(1)=wesn(1)-evap
            sndz(1)=sndz(1)-evap
         endif
c1231         htsnn(1)=htsnn(1)+evap*fsn
         pre=0
         smelt=0
         tmelt=0
         flowth=0
         fresh=0
c******* take care of meltwater moving through the pack
         wcont(1)=(1-fices(1))*wesn(1)
         do 50 i=1,3
            smelt=0.
            if(wcont(i) .gt. .055*sndz(i))then
               smelt=wcont(i) - .055*sndz(i)
               wesn(i)=wesn(i)-smelt
               wcont(i)=.055*sndz(i)
               if (i .eq. 1)sndz(i)=sndz(i)-smelt*(1000./rhos(i))
c030199               htsnn(i)=htsnn(i)+smelt*fsn
c
  35           if(i .eq. 3 ) go to 60
               wesn(i+1)=wesn(i+1)+smelt
               fices(i+1)=-htsnn(i+1)/(fsn*wesn(i+1))
               if(fices(i+1) .gt. 1) fices(i+1)=1
               if(fices(i+1) .lt. 0) fices(i+1)=0
               wcont(i+1)=(1-fices(1+1))*wesn(i+1)
            endif
   50    continue
   60    pre=smelt/dts
c******* take care of water flow through due to precip
         if(effpr-snowf .gt. 0.)then
            wesn(1)=wesn(1)+(effpr-snowf)*dts
            wcont(1)=wcont(1)+(effpr-snowf)*dts
            do 70 i=1,3
               flowth=0.
               if( wcont(i) .gt. .055*sndz(i))then
                  flowth=(wcont(i) - .055*sndz(i))
                  wesn(i)=wesn(i)-flowth
                  wcont(i)=.055*sndz(i)
c
  75              if(i .eq. 3 ) go to 80
                  wesn(i+1)=wesn(i+1)+flowth
                  fices(i+1)=-htsnn(i+1)/(fsn*wesn(i+1))
                  if(fices(i+1) .gt. 1) fices(i+1)=1
                  if(fices(i+1) .lt. 0) fices(i+1)=0
                  wcont(i+1)=(1-fices(1+1))*wesn(i+1)
               endif
   70       continue
   80       pre=pre+flowth/dts
         endif
c******* update water equivalent snowdepth from snowf
         fresh=snowf*dts
c******* update the snow density due to compaction. taken from verseghy
c         international journal of climatology,vol 11, 111-133,1991
         do 95 i=1,3
            rhos(i)=wesn(i)*1000./sndz(i)
            mass=0.5*wesn(1)*1000.
            if(i.ge.2)mass=2.0*mass+0.5*wesn(2)*1000.
            if(i.eq.3)mass=mass+0.5*wesn(2)*1000.+0.5*wesn(3)*1000.
            tplk=min(tpsn(i)+273.16,273.16)
            rhos(i)=rhos(i)+(.5*9.81*1e-7*rhos(i)*mass*exp(14.643-
     1      (4000./tplk)-.02*rhos(i)))*dts
            sndz(i)=wesn(i)*1000./rhos(i)
   95    continue
c******* update our new boundaries
         sdold(1)=sndz(1)+sndz(2)+sndz(3)
         sdold(2)=sndz(2)+sndz(3)
         sdold(3)=sndz(3)
         totsn=sndz(1)+sndz(2)+sndz(3)
         totsn=totsn+(fresh*1000./150.)
         if(totsn .gt. 0.2) then
            sdnew(1)=totsn
            sdnew(2)=totsn-0.05
            sdnew(3)=.66*sdnew(2)
         else
            sdnew(1)=totsn
            sdnew(2)=totsn/2 + totsn/4
            sdnew(3)=totsn/4
         endif
         sndz(1)=sdnew(1)-sdnew(2)
         sndz(2)=sdnew(2)-sdnew(3)
         sndz(3)=sdnew(3)
         snowd=wesn(1)+wesn(2)+wesn(3)
c
c******* since the snow boundary is a moving boundary the heat contents
c        and water equivalents of the snow layers must be updated.
c        From old layer heat and mass contents and old boundaries
c        redistribute new heatand water contents to the new layers.
c******* dtermine the heat capacity, heat content, and water content of
c        fresh snow
         cpnew=150.*1.9e6/920
         hnewp=(ts-tfrz)*cpnew-fsn*150./1000.
         wfrsh=150./1000.
c******* determine heat and water contents of the old layer distribution

         if(sdold(1) .gt. sdold(2))then
            h1p=htsno(1)/(sdold(1)-sdold(2))
            we1p=wesn(1)/(sdold(1)-sdold(2))
         endif
         if(sdold(2) .gt. sdold(3))then
            h2p=htsno(2)/(sdold(2)-sdold(3))
            we2p=wesn(2)/(sdold(2)-sdold(3))
         endif
         if(sdold(3) .gt. 0.)then
            h3p=htsno(3)/sdold(3)
            we3p=wesn(3)/sdold(3)
         endif
c******* calculate the heat contents and water equiv for the new layer 3

         if(sdold(3) .le. sdnew(3))then
            dz3=sdold(3)
         else
            dz3=sdnew(3)
         endif
         if(sdold(2) .le. sdnew(3))then
            dz2=sdold(2)-sdold(3)
         else if(sdold(3) .le. sdnew(3))then
            dz2=sdnew(3)-sdold(3)
         else
            dz2=0
         endif
         if(sdold(1) .le. sdnew(3))then
            dz1=sdold(1)-sdold(2)
         else if(sdold(2) .le. sdnew(3))then
            dz1=sdnew(3)-sdold(2)
         else
            dz1=0
         endif
         if(sdold(1) .le. sdnew(3))then
            zfsh=sdnew(3)-sdold(1)
         else
            zfsh=0
         endif
            htsnn(3)=h3p*dz3+h2p*dz2+h1p*dz1+hnewp*zfsh
            wesn(3)=we3p*dz3+we2p*dz2+we1p*dz1+wfrsh*zfsh
         if(abs(sndz(3)-dz1-dz2-dz3-zfsh).gt.1e-6)write(6,*)
     1   'bad dz layer 3',sndz(2),dz1,dz2,dz3,zfsh
c******* calculate the heat contents and water equiv for the new layer 2

         if(sdold(3) .ge. sdnew(2))then
            dz3=sdnew(2)-sdnew(3)
         else if(sdold(3) .ge. sdnew(3))then
            dz3=sdold(3)-sdnew(3)
         else
            dz3=0
         endif
         if(sdold(2) .le. sdnew(2) .and. sdold(3) .ge. sdnew(3))then
            dz2=sdold(2)-sdold(3)
         else if(sdold(2).le.sdnew(2) .and. sdold(3) .le. sdnew(3))then
            dz2=sdold(2)-sdnew(3)
         else if(sdold(2).ge.sdnew(2) .and. sdold(3) .le. sdnew(2))then
            dz2=sdnew(2)-sdold(3)
            if(sdold(3) .le. sdnew(3))dz2=sdnew(2)-sdnew(3)
         else
            dz2=0
         endif
         if(sdold(1) .le. sdnew(2) .and. sdold(2) .ge. sdnew(3))then
            dz1=sdold(1)-sdold(2)
         else if(sdold(1).le.sdnew(2) .and. sdold(2) .le. sdnew(3))then
            dz1=sdold(1)-sdnew(3)
         else if(sdold(1).ge.sdnew(2) .and. sdold(2) .le. sdnew(2))then
            dz1=sdnew(2)-sdold(2)
            if(sdold(2) .le. sdnew(3))dz1=sdnew(2)-sdnew(3)
         else
            dz1=0
         endif
         if(sdold(1) .le. sdnew(3))then
            zfsh=sdnew(2)-sdnew(3)
         else if (sdold(1) .le. sdnew(2))then
            zfsh=sdnew(2)-sdold(1)
         else
            zfsh=0
         endif
         htsnn(2)=h3p*dz3+h2p*dz2+h1p*dz1+hnewp*zfsh
         wesn(2)=we3p*dz3+we2p*dz2+we1p*dz1+wfrsh*zfsh
         if(abs(sndz(2)-dz1-dz2-dz3-zfsh).gt.1e-6) then
            write(6,*) 'bad dz layer 2',sndz(2),dz1,dz2,dz3,zfsh
            endif
c******* calculate the heat contents and water equiv for the new layer 1

         if(sdold(3) .ge. sdnew(2))then
            dz3=sdold(3)-sdnew(2)
         else
            dz3=0
         endif
         if(sdold(2) .le. sdnew(1) .and. sdold(3) .ge. sdnew(2))then
            dz2=sdold(2)-sdold(3)
         else if(sdold(2) .ge. sdnew(2))then
            dz2=sdold(2)-sdnew(2)
         else
            dz2=0
         endif
         if(sdold(1) .le. sdnew(1) .and. sdold(2) .ge. sdnew(2))then
            dz1=sdold(1)-sdold(2)
         else if(sdold(1).ge.sdnew(2) .and. sdold(2) .le.sdnew(2))then
            dz1=sdold(1)-sdnew(2)
         else
            dz1=0
         endif
         if(sdold(1) .ge. sdnew(2))then
            zfsh=sdnew(1)-sdold(1)
         else
            zfsh=sdnew(1)-sdnew(2)
         endif
         htsnn(1)=h3p*dz3+h2p*dz2+h1p*dz1+hnewp*zfsh
         wesn(1)=we3p*dz3+we2p*dz2+we1p*dz1+wfrsh*zfsh
         if(abs(sndz(1)-dz1-dz2-dz3-zfsh).gt.1e-6)write(6,*)
     1   'bad dz layer 1',sndz(1),dz1,dz2,dz3,zfsh
c******* determine layer density(rhos(i)), specific heat
c        content(cpsn(i)) and thermal conductivity(tksn(i))
         do 100 i=1,3
            rhos(i)=wesn(i)*1000/sndz(i)
            cpsn(i)=rhos(i)*1.9e6/920
            tksn(i)=3.2217e-12*(rhos(i)*1e3)*(rhos(i)*1e3)
            if(rhos(i) .ge. 1000)then
               write(6,*)'rhos too large for timestep'
c              stop 'rhos'
            endif
 100     continue
         totsn=sndz(1)+sndz(2)+sndz(3)
         snowd=wesn(1)+wesn(2)+wesn(3)
c******* from the new layer heat and water contents, calculate layer
c        temperatures
         do 115 i=1,3
            fold(i)=fices(i)
            tpsn(i)=0
            fices(i)=-htsnn(i)/(fsn*wesn(i))
            if(fsn*wesn(i)+htsnn(i).ge.0)go to 112
               tpsn(i)=(htsnn(i)+wesn(i)*fsn)/(cpsn(i)*sndz(i))
               fices(i)=1
               go to 115
 112           if(htsnn(i).le.0)go to 115
               tpsn(i)=htsnn(i)/(wesn(i)*cpsn(i))
               fices(i)=0
               if(tpsn(i) .gt. 0)tpsn(i)=0
 115     continue
         tempu=tpsn(1)+tfrz
         do i=1,3
            if(fices(i) .lt. fold(i))sndz(i)=(fices(i)/fold(i))*sndz(i)
         enddo
      endif
c
c****************** end of three layer snow model ************
 2000 continue
      cpnews=150.*1.9e6/920
      hnews=((ts-tfrz)*cpnews-fsn*150./1000.)*zfsh
      if (iflag .eq. 1)hnews=0.
c030199      hmelt=smelt*fsn
      flxtg=(fh31w+fh31i+fh31d)
      surfflx=(swnet+hlwdwn-hlwout-evap0-shflux)
      enout=htsnn(1)+htsnn(2)+htsnn(3)
c030199 bal=enin/dts+hmelt/dts+hnews/dts+areasc*(surfflx+flxtg)-enout/dts
      balh=enin/dts+hnews/dts+areasc*(surfflx+flxtg)-enout/dts
      balw=snowdo+(snowf*dts)-evap+areasc*(effpr-snowf)*dts-pre*dts-
     1     (wesn(1)+wesn(2)+wesn(3))
c      if (kk .eq. 50000)then
c      cccx=6
c      endif
 3000 continue

c      call cpu_time(tfinal)
c      tsnow=tsnow+(tfinal-tstart)



      return
      end
