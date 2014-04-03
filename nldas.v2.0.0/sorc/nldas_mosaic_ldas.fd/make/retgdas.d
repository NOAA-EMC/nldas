retgdas.o retgdas.d : retgdas.F90
retgdas.o : misc.h
retgdas.o : lisdrv_module.o
retgdas.o : time_manager.o
retgdas.o : baseforcing_module.o
retgdas.o : gdasdomain_module.o
retgdas.o : spmdMod.o
retgdas.o : opendap_module.o
retgdas.o : gdasopendap_module.o
retgdas.o : lis_openfileMod.o
retgdas.o : lis_indices_module.o
retgdas.o : bilinear_interpMod.o
retgdas.o : conserv_interpMod.o
