noah_dynsetup.o noah_dynsetup.d : noah_dynsetup.F90
noah_dynsetup.o : misc.h
noah_dynsetup.o : lisdrv_module.o
noah_dynsetup.o : noah_varder.o
noah_dynsetup.o : spmdMod.o
noah_dynsetup.o : noahpardef_module.o
