accumDGVM.o accumDGVM.d : accumDGVM.F90
accumDGVM.o : misc.h
accumDGVM.o : precision.o
accumDGVM.o : clm_varder.o
accumDGVM.o : accumulMod.o
accumDGVM.o : clm_varmap.o
accumDGVM.o : time_manager.o
