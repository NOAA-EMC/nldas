mos_main.o mos_main.d : mos_main.F90
mos_main.o : lisdrv_module.o
mos_main.o : mos_varder.o
mos_main.o : sibalb_module.o
mos_main.o : tile_spmdMod.o
