retnldas.o retnldas.d : retnldas.F90
retnldas.o : lisdrv_module.o
retnldas.o : nldasdomain_module.o
retnldas.o : baseforcing_module.o
retnldas.o : bilinear_interpMod.o
retnldas.o : conserv_interpMod.o
