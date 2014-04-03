mos_atmdrv.o mos_atmdrv.d : mos_atmdrv.F90
mos_atmdrv.o : lisdrv_module.o
mos_atmdrv.o : spmdMod.o
mos_atmdrv.o : tile_spmdMod.o
mos_atmdrv.o : mos_varder.o
