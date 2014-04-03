surfFileMod.o surfFileMod.d : surfFileMod.F90
surfFileMod.o : misc.h
surfFileMod.o : precision.o
surfFileMod.o : clm_varpar.o
surfFileMod.o : clm_varctl.o
surfFileMod.o : pft_varcon.o
surfFileMod.o : fileutils.o
surfFileMod.o : spmdMod.o
surfFileMod.o : areaMod.o
