vic_setup.o vic_setup.d : vic_setup.F90
vic_setup.o : misc.h
vic_setup.o : lisdrv_module.o
vic_setup.o : vic_varder.o
vic_setup.o : spmdMod.o
vic_setup.o : opendap_module.o
