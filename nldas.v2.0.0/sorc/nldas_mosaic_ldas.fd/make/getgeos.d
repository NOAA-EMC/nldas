getgeos.o getgeos.d : getgeos.F90
getgeos.o : misc.h
getgeos.o : lisdrv_module.o
getgeos.o : time_manager.o
getgeos.o : spmdMod.o
getgeos.o : tile_spmdMod.o
getgeos.o : baseforcing_module.o
getgeos.o : geosdomain_module.o
getgeos.o : bilinear_interpMod.o
getgeos.o : conserv_interpMod.o
getgeos.o : lis_indices_module.o
getgeos.o : grid_spmdMod.o
