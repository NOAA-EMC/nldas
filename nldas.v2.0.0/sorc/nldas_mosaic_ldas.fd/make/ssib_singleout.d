ssib_singleout.o ssib_singleout.d : ssib_singleout.F90
ssib_singleout.o : lis_module.o
ssib_singleout.o : tile_module.o
ssib_singleout.o : ssib_varder.o
ssib_singleout.o : drv_output_mod.o
