ssib_singlegather.o ssib_singlegather.d : ssib_singlegather.F90
ssib_singlegather.o : misc.h
ssib_singlegather.o : lisdrv_module.o
ssib_singlegather.o : tile_spmdMod.o
ssib_singlegather.o : ssib_varder.o
ssib_singlegather.o : ssibpardef_module.o
