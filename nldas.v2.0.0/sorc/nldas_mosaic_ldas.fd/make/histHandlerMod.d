histHandlerMod.o histHandlerMod.d : histHandlerMod.F90
histHandlerMod.o : misc.h
histHandlerMod.o : precision.o
histHandlerMod.o : clm_varpar.o
histHandlerMod.o : clm_varctl.o
histHandlerMod.o : histFileMod.o
histHandlerMod.o : fileutils.o
histHandlerMod.o : spmdMod.o
histHandlerMod.o : time_manager.o
histHandlerMod.o : shr_const_mod.o
histHandlerMod.o : clm_csmMod.o
histHandlerMod.o : lisdrv_module.o
