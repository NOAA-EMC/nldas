hyssib_totinit.o hyssib_totinit.d : hyssib_totinit.F90
hyssib_totinit.o : hyssib_varder.o
hyssib_totinit.o : tile_spmdMod.o
hyssib_totinit.o : lisdrv_module.o
