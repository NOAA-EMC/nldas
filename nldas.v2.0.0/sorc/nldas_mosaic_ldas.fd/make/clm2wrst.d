clm2wrst.o clm2wrst.d : clm2wrst.F90
clm2wrst.o : spmdMod.o
clm2wrst.o : restFileMod.o
clm2wrst.o : clm_varctl.o
clm2wrst.o : lisdrv_module.o
