clm_varctl.o clm_varctl.d : clm_varctl.F90
clm_varctl.o : misc.h
clm_varctl.o : precision.o
clm_varctl.o : clm_varpar.o
clm_varctl.o : clm2drv_module.o
