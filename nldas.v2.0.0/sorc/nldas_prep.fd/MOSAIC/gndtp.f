      subroutine gndtp0(t1,zbar,thetaf,ht,fh21w,fh21i,fh21d,
     1                dfh21w,dfh21i,dfh21d,tp,xklh,zc,xklhw)
c using a diffusion equation this code generates ground temperatures
c with depth given t1
c            *****************************************
c        input
c        dts     timestep in seconds
c        t1      terrestrial (layer 1) temperature in deg C
c        zbar    mean depth to the water table. 
c        thetaf  mean vadose zone soil moisture factor (0-1) 
c        output,
c        ht      heat content in layers 2-7
c        tp      ground temperatures in layers 2-7 
c        tdeep   the temperature of the "deep"
c        f21     heat flux between layer 2 and the terrestrial layer (1)
c        df21    derivative of f21 with respect to temperature
c        xfice   a total soil column ice factor (0-1)
c             ***********************************
      parameter (ngm=7,ng=ngm+1)  
      dimension ht(ngm),tp(ngm),fh(ng),fice(ngm),zc(ngm),zb(ng),shc(ngm)
      dimension t1(3),dz(ngm)
      dimension xklh(ngm)
      data dz/0.05,0.0988,0.1952,0.3859,0.7626,1.5071,10.0/

      real tstart,tfinal,tenergy,tsnow,tfluxes,tground
      common/timing/tcatch,tenergy,tsnow,tfluxes,tground

c      call cpu_time(tstart)

c initialize parameters
      n=ngm
      fsn=3.34e+8
      shw=4.185e+6
      shi=2.06e+6
      phi=0.45

c calculate the boundaries, based on the layer thicknesses(dz)
      zb(1)=0.
      do l=1,3
        zb(l+1)=zb(l)-dz(l)
        shc(l)=(2.4e+6)*(1-phi)*dz(l)
      enddo
      do l=1,2
        zc(l)=0.5*(zb(l)+zb(l+1))
      enddo

c 
c evuates the temperatures in the soil  layers based on the heat
c values.  
c             ***********************************
c             input:
c             xw - water in soil layers, m
c             ht - heat in soil layers
c             fsn - heat of fusion of water
c             shc - specific heat capacity of soil
c             shi - specific heat capacity of ice
c             shw - specific heat capcity of water
c             snowd - snow depth, equivalent water m
c             output:
c             tp - temperature of layers, c
c             fice - fraction of ice of layers
c             pre - extra precipitation, i.e. snowmelt, m s-1
c             snowd - snow depth after melting, equivalent water m.
c             ***********************************
c determine fraction of ice and temp of soil layers based on layer
c heat and water content
      do 10 l=2,2
           ws=phi*dz(l)
           xw=0.5*ws 
           tp(l)=0.
           if(xw.ge.1.e-12)then
              fice(l)=-ht(l)/(fsn*xw)
           else
              fice(l)=0.
           endif
           if(fsn*xw+ht(l).ge.0.)go to 6
           tp(l)=(ht(l)+xw*fsn)/(shc(l)+xw*shi)
           fice(l)=1.
           go to 10
    6      if(ht(l).le.0.)go to 10
           tp(l)=ht(l)/(shc(l)+xw*shw)
           fice(l)=0.
   10 continue
c
c            
c evaluates:  layer thermal conductivities
c *****************************************
c             from farouki(cold regions sci and tech, 5, 1981,
c             67-75) the values for tk1,tk2,tk3 are as follows:
c             tk2=2.2**(phi(l,ibv)-xw), tk3=.57**xw, and
c             tk1=3**(1-phi(l,ibv) for %sand<.5, and
c             for computation purposes i have fit these eqs.
c             to 2nd order polynomials.
c             ***********************************
c             input:
c             sklh - soil heat conductivities of layers
c             zb - soil layer boundaries, m
c             zc - soil layer centers, m
c             dz - layer thickness, m
c             w - soil water content, m
c             phi - soil porosity, dimensionless
c             q - % sand, silt, clay, peat
c             fice - fraction of ice in layers
c             output:
c             xklh - thermal conductivity, w m-2 k-1
c             ***********************************
c lets get the thermal conductivity for the layers
      do l=2,2
         a1=1-phi
         tk1=1.01692+a1*(0.89865+1.06211*a1)
         xw=phi*(1.-fice(l))
         a2=phi-xw
         tk2=1.00543+a2*(0.723371+.464342*a2)
         tk3=0.998899+xw*(-0.548043+0.120291*xw)
         tksat=tk1*tk2*tk3
c*****
         xwi=1.0
         if (zbar .le. zb(l+1))then
            xwi=thetaf
         elseif (zbar .ge. zb(l+1) .and. zbar .le. zb(l))then
            xd1=zb(l)-zbar
            xd2=zbar-zb(l+1)
            xwi=((xd1*thetaf)+xd2)/(xd1+xd2)
         endif 
c*****  
         xwi=min(xwi,1.)
         xklh(l)=(tksat-tk1)*xwi + tk1
         if(l .eq. 2)then
            xklhw=(tksat-tk1) + tk1
         endif   
      enddo                
c
c evaluates heat flux between layers due to heat diffussion
c             ***********************************
c             input:
c             zb - soil layer boundaries, m
c             zc - soil layer centers, m
c             dz - layer thickness, m
c             fice - fraction of ice in layers
c             tp - temperature of layers, c
c             shw - specific heat of water
c             shi - specific heat of ice
c             fsn - heat of fusion
c             output:
c             fh - heat flux between layers
c             ***********************************
c
c total heat flux is via diffusion along the temperature gradient
      fh(n+1)=0.
      denom=zc(1)-zc(2)
      fh21w=-xklhw*(t1(1)-273.16-tp(2))/denom
      fh21i=-xklh(2)*(t1(2)-273.16-tp(2))/denom
      fh21d=-xklh(2)*(t1(3)-273.16-tp(2))/denom
      dfh21w=-xklhw/denom
      dfh21i=-xklh(2)/denom
      dfh21d=dfh21i
c      
c      call cpu_time(tfinal)
c      tground=tground+(tfinal-tstart)

      return
      end





      subroutine gndtmp(dts,zbar,thetaf,fh21,ht,xfice,tp)
c using a diffusion equation this code generates ground temperatures
c with depth given t1
c            *****************************************
c        input
c        dts     timestep in seconds
c        t1      terrestrial (layer 1) temperature in deg C
c        zbar    mean depth to the water table. 
c        thetaf  mean vadose zone soil moisture factor (0-1) 
c        output,
c        ht      heat content in layers 2-7
c        tp      ground temperatures in layers 2-7 
c        tdeep   the temperature of the "deep"
c        f21     heat flux between layer 2 and the terrestrial layer (1)
c        df21    derivative of f21 with respect to temperature
c        xfice   a total soil column ice factor (0-1)
c             ***********************************
      parameter (ngm=7,ng=ngm+1)  
      dimension ht(ngm),tp(ngm),fh(ng),fice(ngm),zc(ngm),zb(ng),shc(ngm)
      dimension dz(ngm)
      dimension xklh(ngm)
      real tstart,tfinal,tenergy,tsnow,tfluxes,tground
      common/timing/tcatch,tenergy,tsnow,tfluxes,tground
      data dz/0.05,0.0988,0.1952,0.3859,0.7626,1.5071,10.0/

c      call cpu_time(tstart)

c initialize parameters
      n=ngm
      fsn=3.34e+8
      shw=4.185e+6
      shi=2.06e+6
      phi=0.45

c calculate the boundaries, based on the layer thicknesses(dz)
      zb(1)=0.
      do l=1,n
        zb(l+1)=zb(l)-dz(l)
        shc(l)=(2.4e+6)*(1-phi)*dz(l)
      enddo
      do l=1,n
        zc(l)=0.5*(zb(l)+zb(l+1))
      enddo

c 
c evuates the temperatures in the soil  layers based on the heat
c values.  
c             ***********************************
c             input:
c             xw - water in soil layers, m
c             ht - heat in soil layers
c             fsn - heat of fusion of water
c             shc - specific heat capacity of soil
c             shi - specific heat capacity of ice
c             shw - specific heat capcity of water
c             snowd - snow depth, equivalent water m
c             output:
c             tp - temperature of layers, c
c             fice - fraction of ice of layers
c             pre - extra precipitation, i.e. snowmelt, m s-1
c             snowd - snow depth after melting, equivalent water m.
c             ***********************************
c determine fraction of ice and temp of soil layers based on layer
c heat and water content
      do 10 l=2,n
            ws=phi*dz(l)
c            xl=l        
c            xw=(1/(7-xl))*ws
           xw=0.5*ws 
           tp(l)=0.
           if(xw.ge.1.e-12)then
              fice(l)=-ht(l)/(fsn*xw)
           else
              fice(l)=0.
           endif
           if(fsn*xw+ht(l).ge.0.)go to 6
           tp(l)=(ht(l)+xw*fsn)/(shc(l)+xw*shi)
           fice(l)=1.
           go to 10
    6      if(ht(l).le.0.)go to 10
           tp(l)=ht(l)/(shc(l)+xw*shw)
           fice(l)=0.
   10 continue
c
c            
c evaluates:  layer thermal conductivities
c *****************************************
c             from farouki(cold regions sci and tech, 5, 1981,
c             67-75) the values for tk1,tk2,tk3 are as follows:
c             tk2=2.2**(phi(l,ibv)-xw), tk3=.57**xw, and
c             tk1=3**(1-phi(l,ibv) for %sand<.5, and
c             for computation purposes i have fit these eqs.
c             to 2nd order polynomials.
c             ***********************************
c             input:
c             sklh - soil heat conductivities of layers
c             zb - soil layer boundaries, m
c             zc - soil layer centers, m
c             dz - layer thickness, m
c             w - soil water content, m
c             phi - soil porosity, dimensionless
c             q - % sand, silt, clay, peat
c             fice - fraction of ice in layers
c             output:
c             xklh - thermal conductivity, w m-2 k-1
c             ***********************************
c lets get the thermal conductivity for the layers
      do l=2,n
         a1=1-phi
         tk1=1.01692+a1*(0.89865+1.06211*a1)
         xw=phi*(1.-fice(l))
         a2=phi-xw
         tk2=1.00543+a2*(0.723371+.464342*a2)
         tk3=0.998899+xw*(-0.548043+0.120291*xw)
         tksat=tk1*tk2*tk3
c*****
         xwi=1.0
         if (zbar .le. zb(l+1))then
            xwi=thetaf
         elseif (zbar .ge. zb(l+1) .and. zbar .le. zb(l))then
            xd1=zb(l)-zbar
            xd2=zbar-zb(l+1)
            xwi=((xd1*thetaf)+xd2)/(xd1+xd2)
         endif 
c*****  
         xwi=min(xwi,1.)
         xklh(l)=(tksat-tk1)*xwi + tk1   
      enddo                
c
c evaluates heat flux between layers due to heat diffussion
c             ***********************************
c             input:
c             zb - soil layer boundaries, m
c             zc - soil layer centers, m
c             dz - layer thickness, m
c             fice - fraction of ice in layers
c             tp - temperature of layers, c
c             shw - specific heat of water
c             shi - specific heat of ice
c             fsn - heat of fusion
c             output:
c             fh - heat flux between layers
c             ***********************************
c
c total heat flux is via diffusion along the temperature gradient
      fh(n+1)=0.
      fh(2)=fh21
      do l=3,n
C THIS xkth is NEW (ie., Agnes corrected) - it should be fixed in all
c codes I'm using      
         xkth=((zb(l)-zc(l-1))*xklh(l-1)
     &   +(zc(l)-zb(l))*xklh(l))
     &   /(zc(l)-zc(l-1))     
         fh(l)=-xkth*(tp(l-1)-tp(l))/(zc(l-1)-zc(l))
      enddo
c
c update the heat contents in the model layers; ht(l)
c IF THERE'S SNOW THIS WILL HAVE TO BE MODIFIED L=1,N
      do l=2,n
        ht(l)=ht(l)+(fh(l+1)-fh(l))*dts
      enddo
c 
c evuates the temperatures in the soil  layers based on the heat
c values.  
c             ***********************************
c             input:
c             xw - water in soil layers, m
c             ht - heat in soil layers
c             fsn - heat of fusion of water
c             shc - specific heat capacity of soil
c             shi - specific heat capacity of ice
c             shw - specific heat capcity of water
c             snowd - snow depth, equivalent water m
c             output:
c             tp - temperature of layers, c
c             fice - fraction of ice of layers
c             pre - extra precipitation, i.e. snowmelt, m s-1
c             snowd - snow depth after melting, equivalent water m.
c             ***********************************
c determine fraction of ice and temp of soil layers based on layer
c heat and water content
      do 1000 l=2,n
           ws=phi*dz(l)
c            xl=l        
c            xw=(1/(7-xl))*ws
           xw=0.5*ws 
           tp(l)=0.
           if(xw.ge.1.e-12)then
              fice(l)=-ht(l)/(fsn*xw)
           else
              fice(l)=0.
           endif
           if(fsn*xw+ht(l).ge.0.)go to 600
           tp(l)=(ht(l)+xw*fsn)/(shc(l)+xw*shi)
           fice(l)=1.
           go to 1000
 600       if(ht(l).le.0.)go to 1000
           tp(l)=ht(l)/(shc(l)+xw*shw)
           fice(l)=0.
 1000   continue
c   
c determine the value of xfice
      xfice=0.0
      do l=1,n
         if(zbar .ge. zb(l+1))then
            lstart=l
            go to 20
         endif
      enddo   
   20 do l=lstart,n
         xfice=xfice+fice(l)
      enddo
      xfice=xfice/((n+1)-lstart)      
c       
c      call cpu_time(tfinal)
c      tground=tground+(tfinal-tstart)

      Return
      end           
