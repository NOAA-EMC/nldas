ssib_main.o ssib_main.d : ssib_main.F90
ssib_main.o : lisdrv_module.o
ssib_main.o : ssib_varder.o
ssib_main.o : tile_spmdMod.o
