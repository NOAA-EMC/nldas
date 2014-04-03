mos_coldstart.o mos_coldstart.d : mos_coldstart.F90
mos_coldstart.o : lisdrv_module.o
mos_coldstart.o : mos_varder.o
mos_coldstart.o : time_manager.o
mos_coldstart.o : spmdMod.o
