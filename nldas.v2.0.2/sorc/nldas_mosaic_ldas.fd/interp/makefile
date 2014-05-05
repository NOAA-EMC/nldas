iplib:	../iplib(gausslat.o)\
	../iplib(gcdist.o)\
	../iplib(gdsawt.o)\
	../iplib(gdswiz.o)\
	../iplib(gdswiz00.o)\
	../iplib(gdswiz01.o)\
	../iplib(gdswiz03.o)\
	../iplib(gdswiz04.o)\
	../iplib(gdswiz05.o)\
	../iplib(gdswizc9.o)\
	../iplib(gdswizca.o)\
	../iplib(gdswzd.o)\
	../iplib(gdswzd00.o)\
	../iplib(gdswzd01.o)\
	../iplib(gdswzd03.o)\
	../iplib(gdswzd04.o)\
	../iplib(gdswzd05.o)\
	../iplib(gdswzdc9.o)\
	../iplib(gdswzdca.o)\
	../iplib(ijkgds.o)\
	../iplib(ipolates.o)\
	../iplib(ipolatev.o)\
	../iplib(ipxetas.o)\
	../iplib(ipxwafs.o)\
	../iplib(ipxwafs2.o)\
	../iplib(makgds.o)\
	../iplib(movect.o)\
	../iplib(polates0.o)\
	../iplib(polates1.o)\
	../iplib(polates2.o)\
	../iplib(polates3.o)\
	../iplib(polates4.o)\
	../iplib(polates5.o)\
	../iplib(polates6.o)\
	../iplib(polatev0.o)\
	../iplib(polatev1.o)\
	../iplib(polatev2.o)\
	../iplib(polatev3.o)\
	../iplib(polatev4.o)\
	../iplib(polatev6.o)\
	../iplib(polfixs.o)\
	../iplib(polfixv.o)
../iplib(polates0.o):	ijkgds.f
../iplib(polates1.o):	ijkgds.f
../iplib(polates2.o):	ijkgds.f
../iplib(polates3.o):	ijkgds.f
../iplib(polates6.o):	ijkgds.f
../iplib(polatev0.o):	ijkgds.f movect.f
../iplib(polatev1.o):	ijkgds.f movect.f
../iplib(polatev2.o):	ijkgds.f movect.f
../iplib(polatev3.o):	ijkgds.f movect.f
../iplib(polatev6.o):	ijkgds.f movect.f
.f.a:
	f90 -c -O nofastint $<
	ar -ru $@ $*.o
	rm -f $*.o
