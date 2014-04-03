mos_gribout.o mos_gribout.d : mos_gribout.F90
mos_gribout.o : lisdrv_module.o
mos_gribout.o : drv_output_mod.o
mos_gribout.o : mos_varder.o
mos_gribout.o : time_manager.o
