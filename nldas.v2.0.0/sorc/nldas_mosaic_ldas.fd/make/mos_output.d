mos_output.o mos_output.d : mos_output.F90
mos_output.o : lisdrv_module.o
mos_output.o : mos_varder.o
mos_output.o : spmdMod.o
