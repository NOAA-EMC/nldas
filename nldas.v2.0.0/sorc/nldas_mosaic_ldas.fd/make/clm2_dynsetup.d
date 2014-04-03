clm2_dynsetup.o clm2_dynsetup.d : clm2_dynsetup.F90
clm2_dynsetup.o : misc.h
clm2_dynsetup.o : lisdrv_module.o
clm2_dynsetup.o : spmdMod.o
clm2_dynsetup.o : clm2pardef_module.o
clm2_dynsetup.o : clm2_laitable.o
