inicFileMod.o inicFileMod.d : inicFileMod.F90
inicFileMod.o : misc.h
inicFileMod.o : precision.o
inicFileMod.o : clm_varder.o
inicFileMod.o : clm_varmap.o
inicFileMod.o : clm_varpar.o
inicFileMod.o : clm_varcon.o
inicFileMod.o : fileutils.o
inicFileMod.o : spmdMod.o
inicFileMod.o : mpishorthand.o
inicFileMod.o : RtmMod.o
inicFileMod.o : clm_varctl.o
inicFileMod.o : time_manager.o
inicFileMod.o : lisdrv_module.o
