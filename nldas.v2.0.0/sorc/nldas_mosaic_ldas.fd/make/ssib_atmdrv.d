ssib_atmdrv.o ssib_atmdrv.d : ssib_atmdrv.F90
ssib_atmdrv.o : lisdrv_module.o
ssib_atmdrv.o : spmdMod.o
ssib_atmdrv.o : tile_spmdMod.o
ssib_atmdrv.o : ssib_varder.o
