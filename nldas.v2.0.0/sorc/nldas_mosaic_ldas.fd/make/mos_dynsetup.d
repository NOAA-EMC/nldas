mos_dynsetup.o mos_dynsetup.d : mos_dynsetup.F90
mos_dynsetup.o : lisdrv_module.o
mos_dynsetup.o : mos_varder.o
mos_dynsetup.o : spmdMod.o
mos_dynsetup.o : mospardef_module.o
