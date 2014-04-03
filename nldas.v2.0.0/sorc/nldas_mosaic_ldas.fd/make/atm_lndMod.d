atm_lndMod.o atm_lndMod.d : atm_lndMod.F90
atm_lndMod.o : misc.h
atm_lndMod.o : precision.o
atm_lndMod.o : shr_const_mod.o
atm_lndMod.o : initializeMod.o
atm_lndMod.o : lnd_atmMod.o
atm_lndMod.o : mpishorthand.o
atm_lndMod.o : time_manager.o
