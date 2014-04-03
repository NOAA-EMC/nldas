clm2_singlegather.o clm2_singlegather.d : clm2_singlegather.F90
clm2_singlegather.o : misc.h
clm2_singlegather.o : lisdrv_module.o
clm2_singlegather.o : clm_varcon.o
clm2_singlegather.o : clm_varpar.o
clm2_singlegather.o : clm_varder.o
clm2_singlegather.o : tile_spmdMod.o
clm2_singlegather.o : clm2pardef_module.o
