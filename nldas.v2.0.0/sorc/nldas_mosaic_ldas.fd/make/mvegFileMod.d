mvegFileMod.o mvegFileMod.d : mvegFileMod.F90
mvegFileMod.o : misc.h
mvegFileMod.o : precision.o
mvegFileMod.o : clm_varpar.o
mvegFileMod.o : clm_varmap.o
mvegFileMod.o : fileutils.o
mvegFileMod.o : spmdMod.o
mvegFileMod.o : mpishorthand.o
mvegFileMod.o : time_manager.o
mvegFileMod.o : infnan.o
