read_gswplai.o read_gswplai.d : read_gswplai.F90
read_gswplai.o : misc.h
read_gswplai.o : lisdrv_module.o
read_gswplai.o : time_manager.o
read_gswplai.o : gswp_module.o
