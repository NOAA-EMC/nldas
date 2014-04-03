hyssib_output.o hyssib_output.d : hyssib_output.F90
hyssib_output.o : lisdrv_module.o
hyssib_output.o : hyssib_varder.o
hyssib_output.o : spmdMod.o
