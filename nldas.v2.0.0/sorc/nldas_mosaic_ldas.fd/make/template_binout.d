template_binout.o template_binout.d : template_binout.F90
template_binout.o : lisdrv_module.o
template_binout.o : drv_output_mod.o
template_binout.o : template_varder.o
