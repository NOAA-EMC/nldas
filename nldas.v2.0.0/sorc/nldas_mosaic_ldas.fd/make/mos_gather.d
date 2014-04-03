mos_gather.o mos_gather.d : mos_gather.F90
mos_gather.o : misc.h
mos_gather.o : tile_spmdMod.o
mos_gather.o : mos_varder.o
mos_gather.o : mospardef_module.o
