ssib_coldstart.o ssib_coldstart.d : ssib_coldstart.F90
ssib_coldstart.o : lisdrv_module.o
ssib_coldstart.o : ssib_varder.o
ssib_coldstart.o : time_manager.o
ssib_coldstart.o : tile_spmdMod.o
