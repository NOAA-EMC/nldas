hyssib_main.o hyssib_main.d : hyssib_main.F90
hyssib_main.o : lisdrv_module.o
hyssib_main.o : hyssib_varder.o
hyssib_main.o : tile_spmdMod.o
