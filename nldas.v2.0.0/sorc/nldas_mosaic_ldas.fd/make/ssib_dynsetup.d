ssib_dynsetup.o ssib_dynsetup.d : ssib_dynsetup.F90
ssib_dynsetup.o : lisdrv_module.o
ssib_dynsetup.o : ssib_varder.o
ssib_dynsetup.o : spmdMod.o
ssib_dynsetup.o : ssibpardef_module.o
