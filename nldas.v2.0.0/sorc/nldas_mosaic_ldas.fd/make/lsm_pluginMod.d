lsm_pluginMod.o lsm_pluginMod.d : lsm_pluginMod.F90
lsm_pluginMod.o : template_varder.o
lsm_pluginMod.o : noah_varder.o
lsm_pluginMod.o : clm_varder.o
lsm_pluginMod.o : vic_varder.o
lsm_pluginMod.o : atmdrvMod.o
lsm_pluginMod.o : mos_varder.o
lsm_pluginMod.o : hyssib_varder.o
lsm_pluginMod.o : ssib_varder.o
