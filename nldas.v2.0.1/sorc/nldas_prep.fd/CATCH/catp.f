c***********************************************************************
c
      subroutine vegetation(
     i                      ncat,cat_mvfile,
     o                      green1,green2,green3,green4,green5,green6,
     o                      green7,green8,green9,green10,green11,
     o                      green12,lai1,lai2,lai3,lai4,lai5,lai6,lai7,
     o                      lai8,lai9,lai10,lai11,lai12,zol1,zol2,zol3,
     o                      zol4,zol5,zol6,zol7,zol8,zol9,zol10,zol11,
     o                      zol12,albed1,albed2,albed3,albed4,albed5,
     o                      albed6,albed7,albed8,albed9,albed10,albed11,
     o                      albed12,vgd1,vgd2,vgd3,vgd4,vgd5,vgd6,vgd7,
     o                      vgd8,vgd9,vgd10,vgd11,vgd12)
c
c Reads in a climatology of vegetation data.
c
      implicit none
      integer ncat
      real green1(*),green2(*),green3(*),green4(*),green5(*),green6(*)
      real green7(*),green8(*),green9(*),green10(*),green11(*)
      real green12(*),lai1(*),lai2(*),lai3(*),lai4(*),lai5(*),lai6(*)
      real lai7(*),lai8(*),lai9(*),lai10(*),lai11(*),lai12(*),zol1(*)
      real zol2(*),zol3(*),zol4(*),zol5(*),zol6(*),zol7(*),zol8(*)
      real zol9(*),zol10(*),zol11(*),zol12(*),vgd1(*),vgd2(*),vgd3(*)
      real vgd4(*),vgd5(*),vgd6(*),vgd7(*),vgd8(*),vgd9(*),vgd10(*)
      real vgd11(*),vgd12(*),albed1(*),albed2(*),albed3(*),albed4(*)
      real albed5(*),albed6(*),albed7(*),albed8(*),albed9(*),albed10(*)
      real albed11(*),albed12(*)
      character*40 cat_mvfile
c
      integer n,mon
c
c-----------------------------------------------------------------------
c
      open (10,file=cat_mvfile(1:index(cat_mvfile,' ')-1)//'/green.dat',
     &     form='formatted',status='old')
      open (11,file=cat_mvfile(1:index(cat_mvfile,' ')-1)//'/lai.dat',
     &     form='formatted',status='old')
      open (12,file=cat_mvfile(1:index(cat_mvfile,' ')-1)//'/zol.dat',
     &     form='formatted',status='old')
      open (13,file=cat_mvfile(1:index(cat_mvfile,' ')-1)//'/albed.dat',
     &     form='formatted',status='old')
      open (14,file=cat_mvfile(1:index(cat_mvfile,' ')-1)//'/vgd.dat',
     &     form='formatted',status='old')
c
      write (*,*) 'Reading vegetation climatology'
      write (*,*)
      read (10,*) (green1(n), n=1,ncat)
      read (10,*) (green2(n), n=1,ncat)
      read (10,*) (green3(n), n=1,ncat)
      read (10,*) (green4(n), n=1,ncat)
      read (10,*) (green5(n), n=1,ncat)
      read (10,*) (green6(n), n=1,ncat)
      read (10,*) (green7(n), n=1,ncat)
      read (10,*) (green8(n), n=1,ncat)
      read (10,*) (green9(n), n=1,ncat)
      read (10,*) (green10(n), n=1,ncat)
      read (10,*) (green11(n), n=1,ncat)
      read (10,*) (green12(n), n=1,ncat)
      read (11,*) (lai1(n), n=1,ncat)
      read (11,*) (lai2(n), n=1,ncat)
      read (11,*) (lai3(n), n=1,ncat)
      read (11,*) (lai4(n), n=1,ncat)
      read (11,*) (lai5(n), n=1,ncat)
      read (11,*) (lai6(n), n=1,ncat)
      read (11,*) (lai7(n), n=1,ncat)
      read (11,*) (lai8(n), n=1,ncat)
      read (11,*) (lai9(n), n=1,ncat)
      read (11,*) (lai10(n), n=1,ncat)
      read (11,*) (lai11(n), n=1,ncat)
      read (11,*) (lai12(n), n=1,ncat)
      read (12,*) (zol1(n), n=1,ncat)
      read (12,*) (zol2(n), n=1,ncat)
      read (12,*) (zol3(n), n=1,ncat)
      read (12,*) (zol4(n), n=1,ncat)
      read (12,*) (zol5(n), n=1,ncat)
      read (12,*) (zol6(n), n=1,ncat)
      read (12,*) (zol7(n), n=1,ncat)
      read (12,*) (zol8(n), n=1,ncat)
      read (12,*) (zol9(n), n=1,ncat)
      read (12,*) (zol10(n), n=1,ncat)
      read (12,*) (zol11(n), n=1,ncat)
      read (12,*) (zol12(n), n=1,ncat)
      read (13,*) (albed1(n), n=1,ncat)
      read (13,*) (albed2(n), n=1,ncat)
      read (13,*) (albed3(n), n=1,ncat)
      read (13,*) (albed4(n), n=1,ncat)
      read (13,*) (albed5(n), n=1,ncat)
      read (13,*) (albed6(n), n=1,ncat)
      read (13,*) (albed7(n), n=1,ncat)
      read (13,*) (albed8(n), n=1,ncat)
      read (13,*) (albed9(n), n=1,ncat)
      read (13,*) (albed10(n), n=1,ncat)
      read (13,*) (albed11(n), n=1,ncat)
      read (13,*) (albed12(n), n=1,ncat)
      read (14,*) (vgd1(n), n=1,ncat)
      read (14,*) (vgd2(n), n=1,ncat)
      read (14,*) (vgd3(n), n=1,ncat)
      read (14,*) (vgd4(n), n=1,ncat)
      read (14,*) (vgd5(n), n=1,ncat)
      read (14,*) (vgd6(n), n=1,ncat)
      read (14,*) (vgd7(n), n=1,ncat)
      read (14,*) (vgd8(n), n=1,ncat)
      read (14,*) (vgd9(n), n=1,ncat)
      read (14,*) (vgd10(n), n=1,ncat)
      read (14,*) (vgd11(n), n=1,ncat)
      read (14,*) (vgd12(n), n=1,ncat)
c
      do n=0,4
        close (10+n,status='keep')
      end do
c
      return
      end
c
c***********************************************************************
c
      subroutine land(
     i                cat_sfile,cat_vfile,ncat,
     o                vegcls,dpth,poros,bf1,bf2,bf3,ars1,ars2,ars3,ara1,
     o                ara2,ara3,ara4,arw1,arw2,arw3,arw4,tsa1,tsa2,tsb1,
     o                tsb2)
c
c Reads in the vegetation class, soil properties and topographic
c parameters.
c
      implicit none
      integer ncat,vegcls(*)
      real dpth(*),poros(*),bf1(*),bf2(*),bf3(*),ars1(*),ars2(*),ars3(*)
      real ara1(*),ara2(*),ara3(*),ara4(*),arw1(*),arw2(*),arw3(*)
      real arw4(*),tsa1(*),tsa2(*),tsb1(*),tsb2(*)
      character*40 cat_sfile,cat_vfile
c
      integer n
      character*100 filename
c
c-----------------------------------------------------------------------
c
c Vegetation class
c
      filename=cat_vfile(1:index(cat_sfile,' ')-1)//'/veg_mosaic'
      write (*,*) 'Reading vegetation class from: ',filename
      open (10,file=filename,form='formatted',status='old')
      read (10,*) (vegcls(n),n=1,ncat)
      close (10,status='keep')
c
c Soil porosity
c
      filename=cat_sfile(1:index(cat_sfile,' ')-1)//'/poros.sol'
      write (*,*) 'Reading soil porosity from: ',filename
      open (10,file=filename,form='formatted',status='old')
      read (10,*) (poros(n),n=1,ncat)
      close (10,status='keep')
c
c Soil depth
c
      filename=cat_sfile(1:index(cat_sfile,' ')-1)//'/profdpth.sol'
      write (*,*) 'Reading soil depth from: ',filename
      open (10,file=filename,form='formatted',status='old')
      read (10,*) (dpth(n),n=1,ncat)
      close (10,status='keep')
c
c Topological parameters
c
      filename=cat_sfile(1:index(cat_sfile,' ')-1)//'/bf.regpar'
      write (*,*) 'Reading topological parameter bf from: ',filename
      open (10,file=filename,form='formatted',status='old')
      read (10,*) (bf1(n),n=1,ncat)
      read (10,*) (bf2(n),n=1,ncat)
      read (10,*) (bf3(n),n=1,ncat)
      close (10,status='keep')
c
      filename=cat_sfile(1:index(cat_sfile,' ')-1)//'/ar.regpar'
      write (*,*) 'Reading topological parameter ar from: ',filename
      open (10,file=filename,form='formatted',status='old')
      read (10,*) (ars1(n),n=1,ncat)
      read (10,*) (ars2(n),n=1,ncat)
      read (10,*) (ars3(n),n=1,ncat)
      read (10,*) (ara1(n),n=1,ncat)
      read (10,*) (ara2(n),n=1,ncat)
      read (10,*) (ara3(n),n=1,ncat)
      read (10,*) (ara4(n),n=1,ncat)
      read (10,*) (arw1(n),n=1,ncat)
      read (10,*) (arw2(n),n=1,ncat)
      read (10,*) (arw3(n),n=1,ncat)
      read (10,*) (arw4(n),n=1,ncat)
      close (10,status='keep')
c
      filename=cat_sfile(1:index(cat_sfile,' ')-1)//'/ts.regpar'
      write (*,*) 'Reading topological parameter ts from: ',filename
      open (10,file=filename,form='formatted',status='old')
      read (10,*) (tsa1(n),n=1,ncat)
      read (10,*) (tsa2(n),n=1,ncat)
      read (10,*) (tsb1(n),n=1,ncat)
      read (10,*) (tsb2(n),n=1,ncat)
      close (10,status='keep')
c
      return
      end
c
c***********************************************************************
c
      subroutine zoom_setup(
     i                      cat_sfile,vegclslu,
     o                      ncatm,nff,latt,lonn)
c
c Determines the zoom pointer arrays and catchment lattitude/longitude.
c
      implicit none
      integer ncatm,nff(*),vegclslu(*)
      real latt(*),lonn(*)
      character*40 cat_sfile
c
      integer n,ncat,ncont,inzoom,iaux,ldas_mask
      real maxx,maxy,minx,miny
      character*100 filename
c
c-----------------------------------------------------------------------
c
c Open catchment definition file
c
      filename=cat_sfile(1:index(cat_sfile,' ')-1)//'/catchment.def'
      open (10,file=filename,form='formatted',status='old')
      read (10,*) ncat
c
c Open LDAS mask file
c 
      filename=cat_sfile(1:index(cat_sfile,' ')-1)//'/ldas.c.mask'
      open (11,file=filename,form='formatted',status='old')
c
c Initialise arrays
c
cc    do n=1,ncat
cc      nff(n)=0
cc      nrst(n)=0
cc    end do
      inzoom=1
      ncont=0
c
c Read LDAS mask and catchment definition files
c
      do n=1,ncat
        read (10,*) iaux,minx,maxx,miny,maxy,iaux
        read (11,*) ldas_mask
c
c Check vegetation class for ice catchments
c
        if ((vegclslu(n) .gt. 0) .and. (vegclslu(n) .ne. 13)) then
          ncont=ncont+1
c
c Check LDAS mask for if catchment in LDAS domain
c
          if (ldas_mask .eq. 1) then
c
c Pointer Indices
c
*           nrst(inzoom)=ncont
            nff(inzoom)=n
            latt(inzoom)=miny+(maxy-miny)/2.
            lonn(inzoom)=minx+(maxx-minx)/2.
            inzoom = inzoom+1
          end if
        end if
      end do
      close (10,status='keep')
      close (11,status='keep')
      ncatm=inzoom-1
c
      return
      end
c
c***********************************************************************
c
      subroutine zoom(
     i                ncatm,nff,tsa1lu,tsa2lu,tsb1lu,tsb2lu,
     i                arw1lu,arw2lu,arw3lu,arw4lu,ara1lu,ara2lu,
     i                ara3lu,ara4lu,ars1lu,ars2lu,ars3lu,bf1lu,
     i                bf2lu,bf3lu,poroslu,dpthlu,vegclslu,
     o                tsa1,tsa2,tsb1,tsb2,arw1,arw2,arw3,arw4,
     o                ara1,ara2,ara3,ara4,ars1,ars2,ars3,bf1,bf2,
     o                bf3,poros,wpwet,bee,zdep1,zdep2,zdep3,
     o                cdcr1,cdcr2,psis,cond,gnu,vgwmax,vegcls)
c
c Zooms the catchment parameters to the LDAS domain.
c
      implicit none
      integer ncatm,nff(*),vegclslu(*),vegcls(*)
      real tsa1lu(*),tsa1(*),tsa2lu(*),tsa2(*),tsb1lu(*),gnu(*)
      real tsb1(*),tsb2lu(*),tsb2(*),arw1lu(*),arw1(*),arw2lu(*)
      real arw2(*),arw3lu(*),arw3(*),arw4lu(*),arw4(*),ara1lu(*)
      real ara1(*),ara2lu(*),ara2(*),ara3lu(*),ara3(*),ara4lu(*)
      real ara4(*),ars1lu(*),ars1(*),ars2lu(*),ars2(*),ars3lu(*)
      real ars3(*),bf1lu(*),bf1(*),bf2lu(*),bf2(*),bf3(*),bf3lu(*)
      real poros(*),poroslu(*),dpthlu(*),wpwet(*),bee(*),cond(*)
      real zdep1(*),zdep2(*),zdep3(*),cdcr1(*),cdcr2(*),psis(*)
      real vgwmax(*)
c
      integer n
      real zmet,term1,term2
c
c-----------------------------------------------------------------------
c
      do n=1,ncatm
c
c 1. Vegetation class
c
        vegcls(n)=vegclslu(nff(n))
c
c 2. Soil parameters
c
c Three soil depths: zdep2 <-> root zone -> water capacity of the root zone
c                    zdep3 <-> unsaturated zone -> bedrock depth approximation
c
        zdep1(n)=20.
        zdep2(n)=1000.
        zdep3(n)=amax1(1000.,dpthlu(nff(n)))
        if (zdep2(n) .gt. 0.75*zdep3(n)) then
          zdep2(n)=0.75*zdep3(n)
        end if
c
c Soil properties
c
        poros(n)=poroslu(nff(n))
        wpwet(n)=0.148/poros(n)
*
* values used to derive the ar and bf parameters: we could also
* spatially variable parameters, but we might have to recompute
* the ar and bf parameters
*
        psis(n)=-.281
        bee(n)=4.
        cond(n)=2.2E-3
        gnu(n)=3.26
c
c Soil storages
c
        vgwmax(n)=poros(n)*zdep2(n)
        zmet=zdep3(n)/1000.
        term1=-1.+((psis(n)-zmet)/psis(n))**
     &        ((bee(n)-1.)/bee(n))
        term2=psis(n)*bee(n)/(bee(n)-1)
        cdcr1(n)=1000.*poros(n)*(zmet-(-term2*term1))
        cdcr2(n)=(1.-wpwet(n))*poros(n)*
     &                zdep3(n)
c
c 3. Topological parameters
c
        bf1(n)=bf1lu(nff(n))
        bf2(n)=bf2lu(nff(n))
        bf3(n)=bf3lu(nff(n))
        ars1(n)=ars1lu(nff(n))
        ars2(n)=ars2lu(nff(n))
        ars3(n)=ars3lu(nff(n))
        ara1(n)=ara1lu(nff(n))
        ara2(n)=ara2lu(nff(n))
        ara3(n)=ara3lu(nff(n))
        ara4(n)=ara4lu(nff(n))
        arw1(n)=arw1lu(nff(n))
        arw2(n)=arw2lu(nff(n))
        arw3(n)=arw3lu(nff(n))
        arw4(n)=arw4lu(nff(n))
        tsa1(n)=tsa1lu(nff(n))
        tsa2(n)=tsa2lu(nff(n))
        tsb1(n)=tsb1lu(nff(n))
        tsb2(n)=tsb2lu(nff(n))
      end do
c
      return
      end
c
c***********************************************************************
c
      subroutine calc_moist(
     i                      ncat,srfexc,rzexc,catdef,
     i                      cdcr1,cdcr2,wpwet,vgwmax,psis,bee,poros,
     i                      d1,d2,d3,ara1,ara2,ara3,ara4,
     i                      arw1,arw2,arw3,arw4,
     o                      srfmc,rzmc,promc)
c
c Calculates the surface soil moisture content, root zone soil moisture
c content and profile soil moisture content from the catchment defecit,
c root zone excess and surface excess
c
      implicit none
      integer ncat
      real srfexc(*),rzexc(*),catdef(*),cdcr1(*),cdcr2(*)
      real wpwet(*),vgwmax(*),psis(*),bee(*),poros(*)
      real ara1(*),ara2(*),ara3(*),ara4(*),arw1(*)
      real arw2(*),arw3(*),arw4(*),srfmc(*),rzmc(*)
      real promc(*),d1(*),d2(*),d3(*)
c
      integer n
      real catdefx,cdi,alpha,wmin,rzwet,rzsuc,srfwet
c
c ----------------------------------------------------------------------
c
      do n=1,ncat
c
c calculate the minimum soil wetness in the root zone moisture
c distribution at equilibrium based on the catchment defecit
c
        catdefx=min(catdef(n),cdcr1(n))
        wmin=arw4(n)+(1.-arw4(n))*(1.+arw1(n)*catdefx)/
     &       (1.+arw2(n)*catdefx+arw3(n)*catdefx*catdefx)
        wmin=min(1.,max(0.,wmin))
c
c calculate the alpha parameter
c
        if (ara1(n) .ne. ara3(n)) then
          cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
        else
          cdi=0.
        end if
        if (catdefx .ge. cdi) then
          alpha=ara3(n)*catdefx+ara4(n)
        else
          alpha=ara1(n)*catdefx+ara2(n)
        end if
c
c calculate the mean root zone wetness at equilibrium by integrating the
c distribution from wmin to infinity
c
        rzwet=exp(-alpha*(1.-wmin))*(wmin-1.-2./alpha)+wmin+2./alpha
        rzwet=min(1.,max(0.,rzwet))
c
c if the catchment defecit is such that a water table no longer exists,
c ramp the equilibrium mean root zone wetness by a scaling factor
c
        if (catdef(n) .gt. cdcr1(n)) then
          rzwet=wpwet(n)+(rzwet-wpwet(n))*(cdcr2(n)-catdef(n))/
     &          (cdcr2(n)-cdcr1(n))
          rzwet=min(1.,max(0.,rzwet))
        end if
c
c calculate the mean soil wetness in the root zone by adding the root
c zone excess to the equilibrium value
c
        rzwet=rzwet+rzexc(n)/vgwmax(n)
        rzwet=min(1.,max(0.,rzwet))
c
c calculate the root zone matric head
c
        rzsuc=psis(n)*rzwet**(-bee(n))
c
c calculate the surface wetness
c
        srfwet=((rzsuc-(d2(n)-d1(n))/2000.)/psis(n))**(-1./bee(n))+
     &         srfexc(n)/(d1(n)*poros(n))
        if (srfwet .lt. 0.) write(*,*) 'negative moisture srf',
     &                      srfwet*poros(n)
        if (srfwet .gt. 1.) write(*,*) 'porosity exceeded srf',
     &                      srfwet*poros(n)
        srfwet=min(1.,max(0.,srfwet))
c
c calculate the surace volumetric soil moisture
c
        srfmc(n)=srfwet*poros(n)
c
c calculate the rootzone soil moisture
c
        rzmc(n)=rzwet*poros(n)+srfexc(n)/d2(n)
        if (rzmc(n) .lt. 0.) write(*,*) 'negative moisture rz',rzmc(n)
        if (rzmc(n) .gt. poros(n)) write(*,*) 'porosity exceeded rz',
     &                             rzmc(n),catdef(n),rzexc(n),srfexc(n)
c
c calculate the total profile soil moisture
c
        promc(n)=((poros(n)*d3(n)-catdef(n))+rzexc(n)+srfexc(n))/d3(n)
        if (promc(n) .lt. 0.) write(*,*) 'negative moisture pro',
     &                        promc(n),catdef(n),rzexc(n),srfexc(n),
     &                        poros(n),d3(n)
        if (promc(n) .gt. poros(n)) write(*,*) 'porosity exceeded pro',
     &                              promc(n)
      end do
c
      return
      end

