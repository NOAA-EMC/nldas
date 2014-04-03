getecmwf.o getecmwf.d : getecmwf.F90
getecmwf.o : lisdrv_module.o
getecmwf.o : baseforcing_module.o
getecmwf.o : time_manager.o
getecmwf.o : ecmwfdomain_module.o
