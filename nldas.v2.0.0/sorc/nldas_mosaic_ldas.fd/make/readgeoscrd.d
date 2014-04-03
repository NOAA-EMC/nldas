readgeoscrd.o readgeoscrd.d : readgeoscrd.F90
readgeoscrd.o : geosdrv_module.o
readgeoscrd.o : geosopendap_module.o
