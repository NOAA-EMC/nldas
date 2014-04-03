getgdas.o getgdas.d : getgdas.F90
getgdas.o : misc.h
getgdas.o : lisdrv_module.o
getgdas.o : baseforcing_module.o
getgdas.o : time_manager.o
getgdas.o : gdasdomain_module.o
getgdas.o : bilinear_interpMod.o
getgdas.o : conserv_interpMod.o
getgdas.o : spmdMod.o
getgdas.o : lis_indices_module.o
getgdas.o : grid_spmdMod.o
