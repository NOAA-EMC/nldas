readgswp.o readgswp.d : readgswp.F90
readgswp.o : misc.h
readgswp.o : lisdrv_module.o
readgswp.o : baseforcing_module.o
readgswp.o : gswp_module.o
readgswp.o : gswpdomain_module.o
