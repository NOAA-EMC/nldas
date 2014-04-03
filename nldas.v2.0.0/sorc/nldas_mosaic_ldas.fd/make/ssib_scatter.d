ssib_scatter.o ssib_scatter.d : ssib_scatter.F90
ssib_scatter.o : misc.h
ssib_scatter.o : tile_spmdMod.o
ssib_scatter.o : ssib_varder.o
ssib_scatter.o : ssibpardef_module.o
