getnldas.o getnldas.d : getnldas.F90
getnldas.o : lisdrv_module.o
getnldas.o : baseforcing_module.o
getnldas.o : nldasdomain_module.o
getnldas.o : time_manager.o
