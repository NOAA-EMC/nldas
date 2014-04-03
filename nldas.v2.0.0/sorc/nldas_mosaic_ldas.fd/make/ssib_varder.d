ssib_varder.o ssib_varder.d : ssib_varder.F90
ssib_varder.o : misc.h
ssib_varder.o : ssib_module.o
ssib_varder.o : tile_spmdMod.o
ssib_varder.o : ssibpardef_module.o
ssib_varder.o : ssibdrv_module.o
ssib_varder.o : opendap_module.o
