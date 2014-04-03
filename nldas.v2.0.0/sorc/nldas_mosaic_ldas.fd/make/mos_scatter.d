mos_scatter.o mos_scatter.d : mos_scatter.F90
mos_scatter.o : misc.h
mos_scatter.o : tile_spmdMod.o
mos_scatter.o : mos_varder.o
mos_scatter.o : mospardef_module.o
