mos_writerst.o mos_writerst.d : mos_writerst.F90
mos_writerst.o : lisdrv_module.o
mos_writerst.o : mos_varder.o
mos_writerst.o : time_manager.o
mos_writerst.o : tile_spmdMod.o
