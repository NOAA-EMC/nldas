mos_binout.o mos_binout.d : mos_binout.F90
mos_binout.o : lisdrv_module.o
mos_binout.o : drv_output_mod.o
mos_binout.o : mos_varder.o
