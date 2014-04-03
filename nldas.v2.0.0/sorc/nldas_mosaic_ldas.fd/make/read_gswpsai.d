read_gswpsai.o read_gswpsai.d : read_gswpsai.F90
read_gswpsai.o : lisdrv_module.o
read_gswpsai.o : time_manager.o
read_gswpsai.o : clm_varder.o
