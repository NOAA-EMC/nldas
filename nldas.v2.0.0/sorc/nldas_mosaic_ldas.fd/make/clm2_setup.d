clm2_setup.o clm2_setup.d : clm2_setup.F90
clm2_setup.o : lisdrv_module.o
clm2_setup.o : spmdMod.o
clm2_setup.o : time_manager.o
clm2_setup.o : clm_varder.o
clm2_setup.o : clm_varcon.o
clm2_setup.o : clm_varctl.o
