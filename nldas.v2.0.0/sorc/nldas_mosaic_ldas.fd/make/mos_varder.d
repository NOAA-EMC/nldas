mos_varder.o mos_varder.d : mos_varder.F90
mos_varder.o : misc.h
mos_varder.o : mos_module.o
mos_varder.o : mosdrv_module.o
mos_varder.o : mospardef_module.o
mos_varder.o : tile_spmdMod.o
mos_varder.o : opendap_module.o
