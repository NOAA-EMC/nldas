hyssib_setup.o hyssib_setup.d : hyssib_setup.F90
hyssib_setup.o : lisdrv_module.o
hyssib_setup.o : hyssib_varder.o
hyssib_setup.o : spmdMod.o
