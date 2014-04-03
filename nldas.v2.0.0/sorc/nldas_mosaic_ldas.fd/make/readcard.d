readcard.o readcard.d : readcard.F90
readcard.o : lisdrv_module.o
readcard.o : time_manager.o
readcard.o : soils_pluginMod.o
readcard.o : lai_pluginMod.o
readcard.o : lis_module.o
