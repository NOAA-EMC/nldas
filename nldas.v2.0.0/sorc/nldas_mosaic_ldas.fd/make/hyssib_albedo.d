hyssib_albedo.o hyssib_albedo.d : hyssib_albedo.F90
hyssib_albedo.o : time_manager.o
hyssib_albedo.o : hyssib_varder.o
hyssib_albedo.o : lisdrv_module.o
hyssib_albedo.o : opendap_module.o
