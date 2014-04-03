clm_varder.o clm_varder.d : clm_varder.F90
clm_varder.o : misc.h
clm_varder.o : clmtype.o
clm_varder.o : spmdMod.o
clm_varder.o : infnan.o
clm_varder.o : tile_spmdMod.o
clm_varder.o : clm_varmap.o
clm_varder.o : clm_varcon.o
clm_varder.o : clm_varctl.o
clm_varder.o : shr_orb_mod.o
clm_varder.o : initializeMod.o
clm_varder.o : clm2pardef_module.o
