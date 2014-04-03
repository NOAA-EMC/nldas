hyssib_scatter.o hyssib_scatter.d : hyssib_scatter.F90
hyssib_scatter.o : misc.h
hyssib_scatter.o : tile_spmdMod.o
hyssib_scatter.o : hyssib_varder.o
hyssib_scatter.o : hyssibpardef_module.o
