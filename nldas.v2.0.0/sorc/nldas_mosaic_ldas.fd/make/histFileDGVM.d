histFileDGVM.o histFileDGVM.d : histFileDGVM.F90
histFileDGVM.o : misc.h
histFileDGVM.o : precision.o
histFileDGVM.o : fileutils.o
histFileDGVM.o : clm_varctl.o
histFileDGVM.o : time_manager.o
histFileDGVM.o : spmdMod.o
histFileDGVM.o : shr_const_mod.o
histFileDGVM.o : lisdrv_module.o
