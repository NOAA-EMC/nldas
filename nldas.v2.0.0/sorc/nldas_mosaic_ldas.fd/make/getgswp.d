getgswp.o getgswp.d : getgswp.F90
getgswp.o : misc.h
getgswp.o : lisdrv_module.o
getgswp.o : time_manager.o
getgswp.o : spmdMod.o
getgswp.o : tile_spmdMod.o
getgswp.o : baseforcing_module.o
getgswp.o : gswpdomain_module.o
