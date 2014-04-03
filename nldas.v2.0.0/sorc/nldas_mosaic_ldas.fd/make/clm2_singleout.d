clm2_singleout.o clm2_singleout.d : clm2_singleout.F90
clm2_singleout.o : lisdrv_module.o
clm2_singleout.o : clm_varcon.o
clm2_singleout.o : clm_varpar.o
clm2_singleout.o : clm_varmap.o
clm2_singleout.o : clm_varctl.o
clm2_singleout.o : drv_output_mod.o
