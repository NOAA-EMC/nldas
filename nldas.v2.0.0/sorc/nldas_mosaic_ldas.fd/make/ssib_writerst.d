ssib_writerst.o ssib_writerst.d : ssib_writerst.F90
ssib_writerst.o : lisdrv_module.o
ssib_writerst.o : ssib_varder.o
ssib_writerst.o : time_manager.o
ssib_writerst.o : tile_spmdMod.o
