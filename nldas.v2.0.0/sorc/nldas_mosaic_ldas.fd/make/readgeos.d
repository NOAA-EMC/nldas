readgeos.o readgeos.d : readgeos.F90
readgeos.o : misc.h
readgeos.o : lisdrv_module.o
readgeos.o : spmdMod.o
readgeos.o : baseforcing_module.o
readgeos.o : bilinear_interpMod.o
readgeos.o : conserv_interpMod.o
readgeos.o : geosdomain_module.o
readgeos.o : geosopendap_module.o
readgeos.o : lis_indices_module.o
readgeos.o : lis_openfileMod.o
