histFileMod.o histFileMod.d : histFileMod.F90
histFileMod.o : misc.h
histFileMod.o : precision.o
histFileMod.o : clmtype.o
histFileMod.o : clm_varpar.o
histFileMod.o : clm_varmap.o
histFileMod.o : shr_const_mod.o
histFileMod.o : fileutils.o
histFileMod.o : clm_varctl.o
histFileMod.o : spmdMod.o
histFileMod.o : clm_varder.o
histFileMod.o : clm_varsur.o
histFileMod.o : time_manager.o
histFileMod.o : infnan.o
