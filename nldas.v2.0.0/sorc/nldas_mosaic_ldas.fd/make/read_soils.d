read_soils.o read_soils.d : read_soils.F90
read_soils.o : lisdrv_module.o
read_soils.o : vic_varder.o
read_soils.o : lis_openfileMod.o
read_soils.o : lis_indices_module.o
