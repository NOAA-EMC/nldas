histUpdate.o histUpdate.d : histUpdate.F90
histUpdate.o : misc.h
histUpdate.o : precision.o
histUpdate.o : clmtype.o
histUpdate.o : clm_varder.o
histUpdate.o : clm_varmap.o
histUpdate.o : clm_varcon.o
histUpdate.o : shr_const_mod.o
histUpdate.o : accumulMod.o
histUpdate.o : pft_varcon.o
histUpdate.o : time_manager.o
histUpdate.o : histFileMod.o
