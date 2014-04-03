lnd_atmMod.o lnd_atmMod.d : lnd_atmMod.F90
lnd_atmMod.o : misc.h
lnd_atmMod.o : precision.o
lnd_atmMod.o : infnan.o
lnd_atmMod.o : clm_varpar.o
lnd_atmMod.o : clm_varmap.o
lnd_atmMod.o : spmdMod.o
lnd_atmMod.o : mpishorthand.o
lnd_atmMod.o : clm_varder.o
lnd_atmMod.o : clm_varcon.o
lnd_atmMod.o : clm_varsur.o
