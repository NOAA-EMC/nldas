getnldas2.o getnldas2.d : getnldas2.F90
getnldas2.o : lisdrv_module.o
getnldas2.o : baseforcing_module.o
getnldas2.o : nldas2domain_module.o
getnldas2.o : time_manager.o
