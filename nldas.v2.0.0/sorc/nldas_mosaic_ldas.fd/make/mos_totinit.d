mos_totinit.o mos_totinit.d : mos_totinit.F90
mos_totinit.o : mos_varder.o
mos_totinit.o : tile_spmdMod.o
mos_totinit.o : lisdrv_module.o
