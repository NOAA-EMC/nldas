ssib_binout.o ssib_binout.d : ssib_binout.F90
ssib_binout.o : lisdrv_module.o
ssib_binout.o : drv_output_mod.o
ssib_binout.o : ssib_varder.o
