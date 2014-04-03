mosrst.o mosrst.d : mosrst.F90
mosrst.o : lisdrv_module.o
mosrst.o : mos_varder.o
mosrst.o : time_manager.o
mosrst.o : tile_spmdMod.o
