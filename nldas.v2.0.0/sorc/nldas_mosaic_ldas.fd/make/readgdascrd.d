readgdascrd.o readgdascrd.d : readgdascrd.F90
readgdascrd.o : gdasdrv_module.o
readgdascrd.o : gdasopendap_module.o
