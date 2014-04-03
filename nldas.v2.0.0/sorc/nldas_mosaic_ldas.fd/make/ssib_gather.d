ssib_gather.o ssib_gather.d : ssib_gather.F90
ssib_gather.o : misc.h
ssib_gather.o : tile_spmdMod.o
ssib_gather.o : ssib_varder.o
ssib_gather.o : ssibpardef_module.o
