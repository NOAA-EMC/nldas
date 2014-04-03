hyssib_dynsetup.o hyssib_dynsetup.d : hyssib_dynsetup.F90
hyssib_dynsetup.o : lisdrv_module.o
hyssib_dynsetup.o : hyssib_varder.o
hyssib_dynsetup.o : spmdMod.o
hyssib_dynsetup.o : hyssibpardef_module.o
