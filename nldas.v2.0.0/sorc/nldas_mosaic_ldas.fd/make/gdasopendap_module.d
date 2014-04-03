gdasopendap_module.o gdasopendap_module.d : gdasopendap_module.F90
gdasopendap_module.o : lisdrv_module.o
gdasopendap_module.o : grid_spmdMod.o
gdasopendap_module.o : tile_spmdMod.o
gdasopendap_module.o : spmdMod.o
gdasopendap_module.o : gdasdrv_module.o
gdasopendap_module.o : opendap_module.o
