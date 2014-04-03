initializeMod.o initializeMod.d : initializeMod.F90
initializeMod.o : misc.h
initializeMod.o : spmdMod.o
initializeMod.o : lisdrv_module.o
initializeMod.o : clm_varcon.o
initializeMod.o : clm_varmap.o
initializeMod.o : clm_varctl.o
initializeMod.o : controlMod.o
initializeMod.o : fileutils.o
initializeMod.o : pftcFileMod.o
initializeMod.o : lis_indices_module.o
initializeMod.o : precision.o
